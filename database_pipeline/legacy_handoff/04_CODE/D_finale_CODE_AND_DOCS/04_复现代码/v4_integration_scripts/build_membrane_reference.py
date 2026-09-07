from __future__ import annotations

import csv
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
RAW = WORK / "raw"
INTERMEDIATE = WORK / "intermediate"
REPORTS = WORK / "reports"
CURRENT_PROTEINS = (
    ROOT
    / "upload"
    / "normalized_tables"
    / "releases"
    / "partner_handoff_v3_1"
    / "v3_1"
    / "protein_database_v3_1.tsv"
)
RETRIEVAL_DATE = datetime.now(timezone.utc).date().isoformat()
UNIPROT_RELEASE_ENDPOINT = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = [
    "accession",
    "id",
    "reviewed",
    "protein_name",
    "gene_primary",
    "gene_names",
    "organism_id",
    "length",
    "ft_transmem",
    "ft_intramem",
    "ft_signal",
    "cc_subcellular_location",
    "keyword",
    "protein_families",
    "go_c",
    "xref_pdb",
    "xref_alphafolddb",
    "xref_hgnc",
    "xref_geneid",
    "xref_ensembl_full",
]
QUERY = "(proteome:UP000005640) AND (organism_id:9606) AND (reviewed:true)"
UNREVIEWED_QUERY = (
    "(proteome:UP000005640) AND (organism_id:9606) AND (reviewed:false) "
    "AND (keyword:KW-0812)"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_tsv(path: Path, query: str = QUERY) -> dict[str, str | int]:
    params = urllib.parse.urlencode(
        {
            "query": query,
            "format": "tsv",
            "fields": ",".join(FIELDS),
        }
    )
    url = f"{UNIPROT_RELEASE_ENDPOINT}?{params}"
    headers = {"User-Agent": "mempro1-v4-build/1.0"}
    request = urllib.request.Request(url, headers=headers)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                payload = response.read()
            path.write_bytes(payload)
            with path.open(encoding="utf-8", newline="") as handle:
                total = max(0, sum(1 for _ in handle) - 1)
            return {"query": query, "url": url, "total_results": total, "pages": 1}
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("UniProt download failed")


def text_has(value: str, phrases: tuple[str, ...]) -> bool:
    low = value.lower()
    return any(phrase in low for phrase in phrases)


def count_features(value: str, feature: str) -> int:
    return len(re.findall(rf"\b{re.escape(feature)}\b", value))


def classify_membrane(row: dict[str, str]) -> dict[str, str | int]:
    transmem = row["Transmembrane"]
    intramem = row["Intramembrane"]
    location = row["Subcellular location [CC]"]
    keywords = row["Keywords"]
    go_cc = row["Gene Ontology (cellular component)"]
    combined = " ; ".join([location, keywords, go_cc])
    tm_count = count_features(transmem, "TRANSMEM")
    intramem_count = count_features(intramem, "INTRAMEM")
    beta_barrel = text_has(
        " ".join([transmem, keywords, row["Protein families"]]),
        ("beta-barrel", "transmembrane beta strand", "porin"),
    )
    lipid_anchor = text_has(
        combined,
        (
            "gpi-anchor",
            "gpi anchor",
            "lipid-anchor",
            "lipid anchor",
            "prenylated",
            "myristoylated",
        ),
    )
    membrane_location = text_has(
        combined,
        (
            "membrane",
            "sarcolemma",
            "postsynaptic density",
            "synaptic vesicle",
        ),
    )

    if tm_count:
        scope = "integral_membrane"
        if beta_barrel:
            topology = "beta_barrel"
        elif tm_count == 1:
            location_low = location.lower()
            if "single-pass type i membrane protein" in location_low:
                topology = "single_pass_type_i"
            elif "single-pass type ii membrane protein" in location_low:
                topology = "single_pass_type_ii"
            elif "single-pass type iii membrane protein" in location_low:
                topology = "single_pass_type_iii"
            elif "single-pass type iv membrane protein" in location_low:
                topology = "single_pass_type_iv"
            else:
                topology = "single_pass_unresolved_type"
        else:
            topology = "multi_pass"
        status = "accepted_core"
        evidence = "uniprot_transmembrane_feature"
    elif beta_barrel:
        scope = "integral_membrane"
        topology = "beta_barrel_no_explicit_feature"
        status = "review"
        evidence = "uniprot_beta_barrel_annotation"
    elif intramem_count:
        scope = "integral_monotopic"
        topology = "intramembrane"
        status = "accepted_core"
        evidence = "uniprot_intramembrane_feature"
    elif lipid_anchor:
        scope = "lipid_anchored"
        topology = "lipid_anchor"
        status = "accepted_extended"
        evidence = "uniprot_lipid_anchor_annotation"
    elif membrane_location:
        scope = "peripheral_membrane_associated"
        topology = "no_transmembrane_feature"
        status = "accepted_extended"
        evidence = "uniprot_membrane_location_annotation"
    else:
        scope = "non_membrane"
        topology = "none"
        status = "excluded"
        evidence = "no_qualifying_membrane_annotation"

    return {
        "membrane_scope": scope,
        "membrane_topology": topology,
        "transmembrane_count": tm_count,
        "intramembrane_count": intramem_count,
        "membrane_inclusion_status": status,
        "membrane_evidence_basis": evidence,
    }


def classify_function(row: dict[str, str]) -> tuple[str, str, str]:
    text = " ".join(
        [
            row["Protein names"],
            row["Gene Names (primary)"],
            row["Keywords"],
            row["Protein families"],
        ]
    ).lower()
    gene = row["Gene Names (primary)"].upper()
    if "g protein-coupled receptor" in text or "g-protein coupled receptor" in text:
        return "receptor", "GPCR", "UniProt family/keyword rule"
    if text_has(text, ("ion channel", "voltage-gated channel", "ligand-gated ion channel")):
        return "ion_channel", "ion_channel", "UniProt family/keyword rule"
    if gene.startswith(("SLC", "ABC")) or text_has(
        text,
        ("transporter", "solute carrier", "translocase", "antiporter", "symporter"),
    ):
        return "transporter", "transporter", "gene/family/name rule"
    if text_has(text, ("receptor tyrosine kinase", "receptor serine/threonine-protein kinase")):
        return "receptor", "receptor_kinase", "UniProt name rule"
    if "receptor" in text:
        return "receptor", "other_receptor", "UniProt name/keyword rule"
    if text_has(
        text,
        (
            "kinase",
            "phosphatase",
            "protease",
            "peptidase",
            "hydrolase",
            "transferase",
            "oxidoreductase",
            "synthase",
        ),
    ):
        return "enzyme", "membrane_associated_enzyme", "UniProt name/keyword rule"
    if text_has(text, ("cell adhesion", "adhesion protein", "cadherin", "integrin")):
        return "adhesion_recognition", "adhesion_or_recognition", "UniProt name/keyword rule"
    return "other", "unclassified", "no_conservative_function_rule"


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    raw_path = RAW / f"uniprot_human_reviewed_reference_proteome_identity_v2_{RETRIEVAL_DATE}.tsv"
    metadata_path = raw_path.with_suffix(".metadata.json")
    if not raw_path.exists():
        metadata = download_tsv(raw_path)
        metadata.update(
            {
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                "sha256": sha256(raw_path),
                "bytes": raw_path.stat().st_size,
            }
        )
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    unreviewed_path = RAW / f"uniprot_human_unreviewed_transmembrane_candidates_identity_v2_{RETRIEVAL_DATE}.tsv"
    unreviewed_metadata_path = unreviewed_path.with_suffix(".metadata.json")
    if not unreviewed_path.exists():
        metadata = download_tsv(unreviewed_path, UNREVIEWED_QUERY)
        metadata.update(
            {
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                "sha256": sha256(unreviewed_path),
                "bytes": unreviewed_path.stat().st_size,
            }
        )
        unreviewed_metadata_path.write_text(
            json.dumps(metadata, indent=2) + "\n",
            encoding="utf-8",
        )

    with raw_path.open(encoding="utf-8", newline="") as handle:
        uniprot_rows = list(csv.DictReader(handle, delimiter="\t"))
    with unreviewed_path.open(encoding="utf-8", newline="") as handle:
        unreviewed_rows = list(csv.DictReader(handle, delimiter="\t"))
    with CURRENT_PROTEINS.open(encoding="utf-8", newline="") as handle:
        current_rows = list(csv.DictReader(handle, delimiter="\t"))
    current_by_id = {row["target_uniprot_id"]: row for row in current_rows}
    reviewed_gene_symbols = {
        row["Gene Names (primary)"].upper()
        for row in uniprot_rows
        if row["Gene Names (primary)"]
    }

    master = []
    evidence = []
    classification = []
    for row in uniprot_rows:
        membrane = classify_membrane(row)
        if membrane["membrane_inclusion_status"] == "excluded":
            continue
        primary_class, subclass, class_rule = classify_function(row)
        accession = row["Entry"]
        ensembl_gene_ids = ";".join(
            sorted(set(re.findall(r"ENSG\d+", row["Ensembl"])))
        )
        current = current_by_id.get(accession, {})
        master_row = {
            "membrane_protein_id": f"HMP_{accession}",
            "target_uniprot_id": accession,
            "uniprot_entry_name": row["Entry Name"],
            "reviewed": row["Reviewed"],
            "protein_name": row["Protein names"],
            "approved_symbol": row["Gene Names (primary)"],
            "gene_names": row["Gene Names"],
            "hgnc_ids": row["HGNC"].strip(";"),
            "ncbi_gene_ids": row["GeneID"].strip(";"),
            "ensembl_gene_ids": ensembl_gene_ids,
            "organism_id": row["Organism (ID)"],
            "sequence_length": row["Length"],
            **membrane,
            "transmembrane_features": row["Transmembrane"],
            "intramembrane_features": row["Intramembrane"],
            "signal_peptide_features": row["Signal peptide"],
            "subcellular_location": row["Subcellular location [CC]"],
            "go_cellular_component": row["Gene Ontology (cellular component)"],
            "uniprot_keywords": row["Keywords"],
            "protein_families": row["Protein families"],
            "pdb_ids": row["PDB"],
            "alphafolddb_ids": row["AlphaFoldDB"],
            "functional_primary_class": primary_class,
            "functional_subclass": subclass,
            "functional_class_rule": class_rule,
            "present_in_v3_1": "1" if accession in current_by_id else "0",
            "v3_1_unique_drug_count": current.get("unique_drug_count", ""),
            "source": "UniProt",
            "source_query": QUERY,
            "retrieval_date": RETRIEVAL_DATE,
        }
        master.append(master_row)
        evidence.append(
            {
                "target_uniprot_id": accession,
                "source": "UniProt",
                "source_query": QUERY,
                "evidence_type": membrane["membrane_evidence_basis"],
                "membrane_scope": membrane["membrane_scope"],
                "membrane_topology": membrane["membrane_topology"],
                "transmembrane_features": row["Transmembrane"],
                "intramembrane_features": row["Intramembrane"],
                "subcellular_location": row["Subcellular location [CC]"],
                "go_cellular_component": row["Gene Ontology (cellular component)"],
                "retrieval_date": RETRIEVAL_DATE,
            }
        )
        classification.append(
            {
                "target_uniprot_id": accession,
                "primary_class": primary_class,
                "subclass": subclass,
                "classification_rule": class_rule,
                "protein_family_annotation": row["Protein families"],
                "membrane_scope": membrane["membrane_scope"],
                "membrane_topology": membrane["membrane_topology"],
                "subcellular_location": row["Subcellular location [CC]"],
            }
        )

    master_fields = list(master[0])
    write_tsv(INTERMEDIATE / "human_membrane_protein_master.tsv", master, master_fields)
    write_tsv(
        INTERMEDIATE / "membrane_annotation_evidence.tsv",
        evidence,
        list(evidence[0]),
    )
    write_tsv(
        INTERMEDIATE / "protein_classification.tsv",
        classification,
        list(classification[0]),
    )
    review_rows = [
        row
        for row in master
        if row["membrane_inclusion_status"] == "review"
        or row["target_uniprot_id"] in current_by_id
        and row["membrane_scope"] == "peripheral_membrane_associated"
    ]
    write_tsv(
        INTERMEDIATE / "membrane_candidate_review.tsv",
        review_rows,
        master_fields,
    )

    master_ids = {row["target_uniprot_id"] for row in master}
    missing_current = [
        {
            "target_uniprot_id": accession,
            "approved_symbol": row["approved_symbol"],
            "protein_name": row["protein_name"],
            "v3_1_membrane_protein_type": row["membrane_protein_type"],
            "v3_1_record_qc_status": row["record_qc_status"],
            "v3_1_unique_drug_count": row["unique_drug_count"],
            "comparison_status": "not_in_reviewed_uniprot_membrane_reference",
        }
        for accession, row in current_by_id.items()
        if accession not in master_ids
    ]
    write_tsv(
        INTERMEDIATE / "v3_1_proteins_not_in_reference.tsv",
        missing_current,
        list(missing_current[0]) if missing_current else ["target_uniprot_id"],
    )
    unreviewed_candidates = []
    for row in unreviewed_rows:
        gene = row["Gene Names (primary)"].upper()
        membrane = classify_membrane(row)
        current = current_by_id.get(row["Entry"], {})
        if gene and gene in reviewed_gene_symbols and not current:
            canonical_resolution = "reviewed_gene_symbol_exists"
        elif current:
            canonical_resolution = "present_in_v3_1_requires_resolution"
        elif gene:
            canonical_resolution = "no_reviewed_gene_symbol_match"
        else:
            canonical_resolution = "missing_primary_gene_symbol"
        unreviewed_candidates.append(
            {
                "target_uniprot_id": row["Entry"],
                "uniprot_entry_name": row["Entry Name"],
                "approved_symbol": row["Gene Names (primary)"],
                "protein_name": row["Protein names"],
                "sequence_length": row["Length"],
                **membrane,
                "protein_families": row["Protein families"],
                "subcellular_location": row["Subcellular location [CC]"],
                "present_in_v3_1": "1" if current else "0",
                "v3_1_unique_drug_count": current.get("unique_drug_count", ""),
                "canonical_resolution_status": canonical_resolution,
                "source_query": UNREVIEWED_QUERY,
                "retrieval_date": RETRIEVAL_DATE,
            }
        )
    write_tsv(
        INTERMEDIATE / "unreviewed_transmembrane_candidate_review.tsv",
        unreviewed_candidates,
        list(unreviewed_candidates[0]),
    )

    scope_counts = Counter(str(row["membrane_scope"]) for row in master)
    topology_counts = Counter(str(row["membrane_topology"]) for row in master)
    class_counts = Counter(str(row["functional_subclass"]) for row in master)
    summary = {
        "status": "passed",
        "retrieval_date": RETRIEVAL_DATE,
        "uniprot_reviewed_human_reference_rows": len(uniprot_rows),
        "membrane_reference_rows": len(master),
        "present_in_v3_1": sum(int(row["present_in_v3_1"]) for row in master),
        "v3_1_total_rows": len(current_rows),
        "v3_1_not_in_reviewed_membrane_reference": len(missing_current),
        "unreviewed_transmembrane_candidate_rows": len(unreviewed_candidates),
        "unreviewed_candidates_present_in_v3_1": sum(
            int(row["present_in_v3_1"]) for row in unreviewed_candidates
        ),
        "unreviewed_candidates_without_reviewed_gene_symbol_match": sum(
            row["canonical_resolution_status"] == "no_reviewed_gene_symbol_match"
            for row in unreviewed_candidates
        ),
        "membrane_scope_counts": dict(sorted(scope_counts.items())),
        "membrane_topology_counts": dict(sorted(topology_counts.items())),
        "functional_subclass_counts": dict(class_counts.most_common()),
        "raw_uniprot_sha256": sha256(raw_path),
    }
    (REPORTS / "MEMBRANE_REFERENCE_BUILD_REPORT.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
