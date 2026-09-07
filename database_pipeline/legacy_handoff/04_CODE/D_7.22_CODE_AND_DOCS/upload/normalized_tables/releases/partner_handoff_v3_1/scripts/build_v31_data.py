import csv
import hashlib
import io
import json
import re
import shutil
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"C:\github-repos\upload\normalized_tables")
OUT = ROOT / "outputs" / "v3_1_completion"
OUT.mkdir(parents=True, exist_ok=True)
csv.field_size_limit(100_000_000)

INPUTS = {
    "binary": ROOT / "drug_protein_binary_relationships.tsv",
    "protein": ROOT / "protein_database.tsv",
    "molecule": ROOT / "small_molecule_database.tsv",
    "pocket": ROOT / "pocket_instances.tsv",
}

OUTPUTS = {
    "binary": OUT / "drug_protein_binary_relationships_v3_1.tsv",
    "protein": OUT / "protein_database_v3_1.tsv",
    "molecule": OUT / "small_molecule_database_v3_1.tsv",
    "pocket": OUT / "pocket_instances_v3_1.tsv",
}

UNIPROT_CACHE = OUT / "uniprot_v3_1_enrichment_cache.json"
RETRIEVED_AT = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
UNIPROT_FIELDS = [
    "accession", "id", "reviewed", "length", "cc_similarity", "ft_transmem",
    "ft_signal", "cc_subcellular_location", "xref_alphafolddb", "xref_pdb", "keyword",
]


def load_tsv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def clean_cell(value):
    return str(value or "").replace("\r", " ").replace("\n", " ").replace("\t", " ").strip()


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({field: clean_cell(row.get(field, "")) for field in fields})


def split_semicolon(value):
    return [x.strip() for x in str(value or "").split(";") if x.strip()]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_url(url, attempts=4):
    req = urllib.request.Request(url, headers={"User-Agent": "membrane-protein-database-v3.1/1.0"})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)


def fetch_uniprot(accessions):
    if UNIPROT_CACHE.exists():
        cached = json.loads(UNIPROT_CACHE.read_text(encoding="utf-8"))
    else:
        cached = {}
    missing = [x for x in accessions if x not in cached]
    for start in range(0, len(missing), 40):
        batch = missing[start:start + 40]
        query = "(" + " OR ".join(f"accession:{x}" for x in batch) + ")"
        params = urllib.parse.urlencode({
            "query": query,
            "format": "tsv",
            "fields": ",".join(UNIPROT_FIELDS),
            "size": 500,
        })
        text = fetch_url("https://rest.uniprot.org/uniprotkb/search?" + params)
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        returned = set()
        for row in reader:
            acc = row.get("Entry", "")
            if acc:
                cached[acc] = dict(row)
                returned.add(acc)
        for acc in batch:
            if acc not in returned:
                cached[acc] = {"_retrieval_status": "not_returned"}
        UNIPROT_CACHE.write_text(json.dumps(cached, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"UNIPROT {min(start + len(batch), len(missing))}/{len(missing)}", flush=True)
        time.sleep(0.15)
    return cached


def transmem_ranges(raw):
    return re.findall(r"TRANSMEM\s+(\d+)\.\.(\d+)", raw or "")


def signal_ranges(raw):
    return re.findall(r"SIGNAL\s+(\d+)\.\.(\d+)", raw or "")


def strip_evidence(text):
    text = re.sub(r"\s*\{ECO:[^}]+\}", "", text or "")
    text = re.sub(r"^SIMILARITY:\s*", "", text)
    return re.sub(r"\s+", " ", text).strip()


def family_parts(similarity):
    clean = strip_evidence(similarity)
    family = ""
    subfamily = ""
    m = re.search(r"Belongs to the ([^.]+? family)\.", clean, re.I)
    if m:
        family = m.group(1).strip()
    m = re.search(r"([^.]+? subfamily)\.", clean, re.I)
    if m:
        subfamily = m.group(1).strip()
    return clean, family, subfamily


def protein_class(keywords, protein_name, family):
    text = " ".join([keywords or "", protein_name or "", family or ""]).lower()
    if "g protein-coupled receptor" in text:
        return "GPCR", "UniProt keyword/family rule"
    if "ion channel" in text or re.search(r"\bchannel\b", text):
        return "ion_channel", "UniProt keyword/name rule"
    if "transporter" in text or "transport protein" in text or "solute carrier" in text:
        return "transporter", "UniProt keyword/name rule"
    if "receptor" in text and ("kinase" in text or "tyrosine-protein kinase" in text):
        return "receptor_kinase", "UniProt keyword/name rule"
    if "receptor" in text:
        return "other_receptor", "UniProt keyword/name rule"
    if "enzyme" in text or any(x in text for x in ["protease", "peptidase", "hydrolase", "oxidase", "reductase", "transferase", "phosphatase"]):
        return "membrane_associated_enzyme", "UniProt keyword/name rule"
    return "unclassified_by_current_rule", "conservative rule; manual review needed"


def membrane_type(tm_count, location, keywords):
    text = f"{location} {keywords}".lower()
    if tm_count >= 2:
        return "multi_pass_integral_membrane", "UniProt TRANSMEM features"
    if tm_count == 1:
        return "single_pass_integral_membrane", "UniProt TRANSMEM feature"
    if "gpi-anchor" in text or "lipid-anchor" in text or "lipid anchor" in text:
        return "lipid_anchored", "UniProt keyword/location"
    if "membrane" in text:
        return "membrane_associated_no_transmem_feature", "UniProt location/keyword; review recommended"
    return "no_current_uniprot_membrane_feature", "record retrieved; manual scope review recommended"


def build_protein(rows, cache):
    added = [
        "sequence_length", "uniprot_family_annotation", "protein_family", "protein_subfamily",
        "protein_class", "protein_class_rule", "membrane_protein_type", "membrane_type_rule",
        "transmembrane_topology_summary", "transmembrane_annotation_status", "signal_peptide_summary",
        "alphafold_model_id", "experimental_pdb_ids", "uniprot_keywords",
        "annotation_source", "annotation_retrieved_at_utc", "record_qc_status",
    ]
    fields = list(rows[0].keys()) + added
    stats = Counter()
    protein_meta = {}
    for row in rows:
        acc = row["target_uniprot_id"]
        u = cache.get(acc, {})
        if u.get("_retrieval_status") == "not_returned" or not u:
            row.update({x: "" for x in added})
            row["annotation_source"] = "UniProt REST API"
            row["annotation_retrieved_at_utc"] = RETRIEVED_AT
            row["record_qc_status"] = "uniprot_record_not_returned"
            stats["uniprot_not_returned"] += 1
            protein_meta[acc] = {}
            continue
        tm = transmem_ranges(u.get("Transmembrane", ""))
        sig = signal_ranges(u.get("Signal peptide", ""))
        tm_count = len(tm)
        topology = ";".join(f"{a}-{b}" for a, b in tm)
        signal = ";".join(f"{a}-{b}" for a, b in sig)
        fam_ann, fam, subfam = family_parts(u.get("Sequence similarities", ""))
        pclass, pclass_rule = protein_class(u.get("Keywords", ""), row.get("protein_name", ""), fam_ann)
        mtype, mtype_rule = membrane_type(tm_count, u.get("Subcellular location [CC]", ""), u.get("Keywords", ""))
        af_ids = split_semicolon(u.get("AlphaFoldDB", ""))
        pdb_ids = split_semicolon(u.get("PDB", ""))
        row["sequence_length"] = u.get("Length", "")
        row["uniprot_family_annotation"] = fam_ann
        row["protein_family"] = fam
        row["protein_subfamily"] = subfam
        row["protein_class"] = pclass
        row["protein_class_rule"] = pclass_rule
        row["membrane_protein_type"] = mtype
        row["membrane_type_rule"] = mtype_rule
        row["transmembrane_topology_summary"] = topology
        row["transmembrane_annotation_status"] = "feature_present" if tm_count else "record_retrieved_no_transmem_feature"
        row["signal_peptide_summary"] = signal
        row["alphafold_model_id"] = ";".join(af_ids)
        row["experimental_pdb_ids"] = ";".join(pdb_ids)
        row["uniprot_keywords"] = u.get("Keywords", "")
        row["annotation_source"] = "UniProt REST API"
        row["annotation_retrieved_at_utc"] = RETRIEVED_AT
        row["record_qc_status"] = "review_membrane_scope" if mtype == "no_current_uniprot_membrane_feature" else "ok"
        row["transmembrane_count"] = str(tm_count)
        row["pdb_structures"] = str(len(pdb_ids))
        row["reviewed"] = u.get("Reviewed", row.get("reviewed", ""))
        stats[f"membrane_type:{mtype}"] += 1
        stats[f"protein_class:{pclass}"] += 1
        if fam:
            stats["family_filled"] += 1
        if af_ids:
            stats["alphafold_filled"] += 1
        protein_meta[acc] = {
            "alphafold_model_id": af_ids[0] if af_ids else "",
            "uniprot_returned": True,
        }
    return fields, rows, stats, protein_meta


def relationship_semantics(row):
    source = row.get("compound_source", "")
    if source == "bioassay_active_relation":
        return "assay_activity", "assay_result", "assay_dependent_unspecified_directness"
    if source == "structure_ligand_context":
        return "structure_ligand_observation", "structure_instance_context", "not_inferred_as_pharmacological_target"
    if source == "uniprot_binding_site_annotation":
        return "uniprot_binding_site_annotation", "annotated_ligand_site", "curated_site_ligand"
    return "curated_target_assertion", "compound_target_claim", "curated_mechanism_or_target_claim"


def normalize_activity(activity_type, value):
    raw = (activity_type or "").strip()
    relation = "=" if str(value or "").strip() else ""
    for op in ("<=", ">=", "<", ">", "~"):
        if raw.endswith(op):
            raw = raw[:-len(op)].strip()
            relation = op
            break
        if raw.startswith(op):
            raw = raw[len(op):].strip()
            relation = op
            break
    canonical = {"KI": "Ki", "KD": "Kd", "IC50": "IC50", "EC50": "EC50"}.get(raw.upper(), raw)
    return canonical, relation


def build_binary(rows):
    added = [
        "relationship_context_id", "evidence_type", "evidence_scope", "evidence_directness",
        "activity_relation", "activity_value_unit", "record_qc_status",
    ]
    fields = list(rows[0].keys()) + added
    digest_counts = Counter()
    stats = Counter()
    for row in rows:
        original = "\x1f".join(row.get(k, "") for k in rows[0].keys())
        digest = hashlib.sha1(original.encode("utf-8")).hexdigest()[:14].upper()
        digest_counts[digest] += 1
        row["relationship_context_id"] = f"RCTX_{digest}_{digest_counts[digest]:02d}"
        etype, scope, directness = relationship_semantics(row)
        row["evidence_type"] = etype
        row["evidence_scope"] = scope
        row["evidence_directness"] = directness
        normalized, relation = normalize_activity(row.get("activity_type", ""), row.get("activity_value_uM", ""))
        old = row.get("activity_type", "")
        row["activity_type"] = normalized
        row["activity_relation"] = relation
        row["activity_value_unit"] = "uM" if row.get("activity_value_uM", "") else ""
        qc = []
        if old in {"A2", "D2"}:
            qc.append("review_suspected_activity_type")
            stats[f"suspected_activity_type:{old}"] += 1
        if normalized != old:
            stats[f"activity_type_normalized:{old}->{normalized}"] += 1
        row["record_qc_status"] = ";".join(qc) if qc else "ok"
        stats[f"evidence_type:{etype}"] += 1
    return fields, rows, stats


def relationship_summary(binary_rows):
    by_drug = defaultdict(lambda: {"rows": 0, "targets": set(), "sources": set(), "compound_sources": set()})
    for r in binary_rows:
        d = r["drug_id"]
        x = by_drug[d]
        x["rows"] += 1
        x["targets"].add(r["target_uniprot_id"])
        x["sources"].update(split_semicolon(r["source_database"]))
        x["compound_sources"].update(split_semicolon(r["compound_source"]))
    return by_drug


def build_molecules(rows, summaries):
    added = ["record_qc_status", "record_qc_note"]
    fields = list(rows[0].keys()) + added
    stats = Counter()
    for row in rows:
        d = row["drug_id"]
        summary = summaries.get(d)
        if summary:
            row["source_row_count"] = str(summary["rows"])
            row["unique_target_count"] = str(len(summary["targets"]))
            row["source_databases"] = ";".join(sorted(summary["sources"]))
            row["compound_source"] = ";".join(sorted(summary["compound_sources"]))
            row["record_qc_status"] = "ok"
            row["record_qc_note"] = ""
        else:
            previous = f"previous_source_databases={row.get('source_databases','')};previous_source_row_count={row.get('source_row_count','')};previous_unique_target_count={row.get('unique_target_count','')}"
            row["source_row_count"] = "0"
            row["unique_target_count"] = "0"
            row["source_databases"] = ""
            row["compound_source"] = ""
            row["record_qc_status"] = "orphan_no_current_relationship"
            row["record_qc_note"] = previous
            stats["orphan_no_current_relationship"] += 1
    return fields, rows, stats


def build_pockets(rows, protein_meta):
    added = [
        "structure_model_id", "structure_model_source", "structure_model_version",
        "pocket_method", "pocket_method_version", "pocket_qc_status",
    ]
    fields = list(rows[0].keys()) + added
    stats = Counter()
    for row in rows:
        acc = row["target_uniprot_id"]
        model = protein_meta.get(acc, {}).get("alphafold_model_id", "")
        row["structure_model_id"] = model or f"AF-{acc}-F1"
        row["structure_model_source"] = "AlphaFoldDB"
        row["structure_model_version"] = "AlphaFold_DB_v6"
        row["pocket_method"] = "geometry_based_concavity_detection"
        row["pocket_method_version"] = "legacy_pipeline_unversioned"
        row["pocket_qc_status"] = "model_xref_confirmed" if model else "model_id_constructed_review"
        stats[row["pocket_qc_status"]] += 1
    return fields, rows, stats


def main():
    proteins = load_tsv(INPUTS["protein"])
    cache = fetch_uniprot([r["target_uniprot_id"] for r in proteins])
    pfields, proteins, pstats, protein_meta = build_protein(proteins, cache)

    binary = load_tsv(INPUTS["binary"])
    bfields, binary, bstats = build_binary(binary)
    summaries = relationship_summary(binary)

    molecules = load_tsv(INPUTS["molecule"])
    mfields, molecules, mstats = build_molecules(molecules, summaries)

    pockets = load_tsv(INPUTS["pocket"])
    kfields, pockets, kstats = build_pockets(pockets, protein_meta)

    write_tsv(OUTPUTS["protein"], pfields, proteins)
    write_tsv(OUTPUTS["binary"], bfields, binary)
    write_tsv(OUTPUTS["molecule"], mfields, molecules)
    write_tsv(OUTPUTS["pocket"], kfields, pockets)

    report = {
        "generated_at_utc": RETRIEVED_AT,
        "scope": "v3.1 targeted completion preserving four-table design",
        "inputs_unchanged": {k: {"path": str(v), "sha256": sha256(v), "bytes": v.stat().st_size} for k, v in INPUTS.items()},
        "outputs": {k: {"path": str(v), "sha256": sha256(v), "bytes": v.stat().st_size} for k, v in OUTPUTS.items()},
        "row_counts": {"binary": len(binary), "protein": len(proteins), "molecule": len(molecules), "pocket": len(pockets)},
        "column_counts": {"binary": len(bfields), "protein": len(pfields), "molecule": len(mfields), "pocket": len(kfields)},
        "protein_enrichment": dict(pstats),
        "relationship_changes": dict(bstats),
        "molecule_qc": dict(mstats),
        "pocket_metadata": dict(kstats),
        "known_limits": [
            "Protein class is a conservative rule derived from UniProt keywords/name/family and retains its rule column.",
            "Membrane protein type is evidence-based from current UniProt features; records without features are flagged for review.",
            "Pocket pLDDT/PAE and physicochemical descriptors are not included in this first v3.1 pass.",
            "Original assay value/unit cannot be reconstructed where v3 retained only standardized micromolar values.",
        ],
    }
    (OUT / "V3_1_QC_REPORT.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
