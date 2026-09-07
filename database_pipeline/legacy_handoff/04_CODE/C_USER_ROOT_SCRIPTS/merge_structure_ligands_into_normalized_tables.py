import csv
import hashlib
import json
import math
import statistics
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")
BINARY = BASE / "drug_protein_binary_relationships.tsv"
PROTEIN = BASE / "protein_database.tsv"
SMALL = BASE / "small_molecule_database.tsv"
REVIEW = BASE / "protein_structure_ligand_context_review.tsv"
MANIFEST = BASE / "MANIFEST.json"
REPORT = BASE / "structure_context_merge_report.json"
PUBCHEM_PROPS_CACHE = BASE / "xlsx_build" / "pubchem_structure_context_properties_cache.json"

csv.field_size_limit(100_000_000)

PUBCHEM_PROPERTIES = [
    "MolecularFormula",
    "MolecularWeight",
    "CanonicalSMILES",
    "IsomericSMILES",
    "InChIKey",
    "IUPACName",
    "XLogP",
    "TPSA",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
    "Complexity",
    "FormalCharge",
]


def load_tsv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def split_values(value):
    return [part.strip() for part in str(value or "").split(";") if part.strip()]


def join_values(values):
    unique = sorted({str(value).strip() for value in values if str(value).strip()})
    return ";".join(unique)


def append_semicolon(existing, additions):
    values = split_values(existing)
    values.extend(additions)
    return join_values(values)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def update_manifest(paths):
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    files = manifest.get("files", [])
    by_name = {Path(item["path"]).name: item for item in files}
    for path in paths:
        item = by_name.get(path.name)
        if item is None:
            item = {"path": f"upload/normalized_tables/{path.name}"}
            files.append(item)
        item["sha256"] = sha256(path)
        item["bytes"] = path.stat().st_size
    manifest["files"] = files
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def stable_unmapped_id(row):
    code = row.get("pdb_ligand_code", "").strip().upper()
    if code:
        return f"UNMAPPED_PDB_LIGAND:{code}"
    chembl_id = row.get("chembl_id", "").strip().upper()
    if chembl_id:
        return f"UNMAPPED_CHEMBL:{chembl_id}"
    return f"UNMAPPED_STRUCTURE_LIGAND:{row.get('structure_ligand_id', '').strip()}"


def ligand_display_name(row):
    name = row.get("ligand_name_raw", "").strip()
    if name:
        return name
    label = row.get("raw_ligand_label", "").strip()
    if label and label != "()":
        return label
    code = row.get("pdb_ligand_code", "").strip().upper()
    if code:
        return f"PDB ligand {code}"
    return row.get("structure_ligand_id", "").strip()


def parse_activity_uM(value, unit):
    value = str(value or "").strip()
    unit = str(unit or "").strip().lower().replace("μ", "u")
    if not value:
        return ""
    try:
        numeric = float(value)
    except ValueError:
        return ""
    if not math.isfinite(numeric):
        return ""
    if unit in {"um", "µm", "microm", "micromolar"}:
        converted = numeric
    elif unit in {"nm", "nanom", "nanomolar"}:
        converted = numeric / 1000.0
    elif unit in {"mm", "millim", "millimolar"}:
        converted = numeric * 1000.0
    elif unit in {"m"}:
        converted = numeric * 1_000_000.0
    else:
        return ""
    return f"{converted:.6g}"


def site_value(row, drug_id):
    prefix = f"CID{drug_id}|" if str(drug_id).isdigit() else f"{drug_id}|"
    return prefix + row.get("raw_entry", "")


def mechanism(row):
    parts = [
        "structure_ligand_context",
        f"category={row.get('ligand_category', '')}",
        f"cid_status={row.get('cid_mapping_status', '')}",
        f"cid_method={row.get('cid_mapping_method', '')}",
        f"keep_small={row.get('keep_for_small_molecule_binding_site', '')}",
        f"keep_context={row.get('keep_for_binding_site_context', '')}",
        f"context_id={row.get('structure_ligand_id', '')}",
    ]
    if row.get("exclusion_reason"):
        parts.append(f"exclusion={row['exclusion_reason']}")
    if row.get("affinity_value"):
        parts.append(f"affinity_original={row.get('affinity_value', '')}{row.get('affinity_unit', '')}")
    return "; ".join(part for part in parts if part)


def make_binary_row(row, binary_fields, protein_by_uid, small_by_id):
    uid = row["target_uniprot_id"]
    protein = protein_by_uid.get(uid, {})
    cid = row.get("ligand_pubchem_cid", "").strip()
    drug_id = cid if cid.isdigit() else stable_unmapped_id(row)
    small = small_by_id.get(drug_id, {})
    source = row.get("structure_source", "")

    out = {field: "" for field in binary_fields}
    out["target_uniprot_id"] = uid
    out["approved_symbol"] = row.get("approved_symbol") or protein.get("approved_symbol", "")
    out["ncbi_gene_id"] = protein.get("ncbi_gene_id", "")
    out["protein_name"] = protein.get("protein_name", "")
    out["drug_id"] = drug_id
    out["compound_id_type"] = "PubChem CID" if cid.isdigit() else "unmapped_structure_ligand"
    out["drug_name"] = small.get("drug_name") or ligand_display_name(row)
    out["drug_name_type"] = small.get("drug_name_type") or row.get("ligand_category", "")
    out["compound_source"] = "structure_context"
    out["source_database"] = source
    out["assay_or_mechanism"] = mechanism(row)
    out["activity_type"] = row.get("affinity_type", "")
    out["activity_value_uM"] = parse_activity_uM(row.get("affinity_value", ""), row.get("affinity_unit", ""))
    out["clinical_or_approval_status"] = ""

    out["has_pdb_biolip_matched"] = "1" if source == "BioLiP" else "0"
    out["has_pdb_scpdb_matched"] = "1" if source == "sc-PDB" else "0"
    out["has_pdbbind_matched"] = "1" if source == "PDBbind" else "0"
    if source == "BioLiP":
        out["pdb_biolip_matched_sites"] = site_value(row, drug_id)
    if source == "sc-PDB":
        out["pdb_scpdb_matched_sites"] = site_value(row, drug_id)
    if source == "PDBbind":
        out["pdbbind_matched_sites"] = site_value(row, drug_id)

    out["has_stitch_matched"] = "0"
    out["has_exp_binding_site"] = "0"
    out["target_group_id"] = protein.get("target_group_id", "")
    return out


def fetch_pubchem_properties(cids):
    if PUBCHEM_PROPS_CACHE.exists():
        cache = json.loads(PUBCHEM_PROPS_CACHE.read_text(encoding="utf-8"))
    else:
        cache = {}

    missing = [cid for cid in cids if cid.isdigit() and cid not in cache]
    props = ",".join(PUBCHEM_PROPERTIES)
    for index in range(0, len(missing), 100):
        batch = missing[index:index + 100]
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            + ",".join(batch)
            + f"/property/{props}/JSON"
        )
        ok = False
        for attempt in range(3):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "codex-structure-context-merge/1.0"})
                with urllib.request.urlopen(req, timeout=45) as response:
                    data = json.loads(response.read().decode("utf-8", errors="ignore"))
                for item in data.get("PropertyTable", {}).get("Properties", []):
                    cid = str(item.get("CID", ""))
                    if cid:
                        cache[cid] = item
                ok = True
                break
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    break
            except Exception:
                pass
            time.sleep(1.5 * (attempt + 1))
        if not ok:
            for cid in batch:
                cache.setdefault(cid, {})
        if index % 500 == 0:
            PUBCHEM_PROPS_CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    PUBCHEM_PROPS_CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    return cache


def aggregate_binary(rows):
    agg = defaultdict(lambda: {
        "sources": set(),
        "compound_sources": set(),
        "targets": set(),
        "activity_types": set(),
        "clinical": set(),
        "activity_values": [],
        "structure_context": False,
        "row_count": 0,
    })
    for row in rows:
        drug_id = row.get("drug_id", "").strip()
        if not drug_id:
            continue
        item = agg[drug_id]
        item["row_count"] += 1
        if row.get("source_database"):
            item["sources"].add(row["source_database"])
        if row.get("compound_source"):
            for part in split_values(row["compound_source"]):
                item["compound_sources"].add(part)
        if row.get("target_uniprot_id"):
            item["targets"].add(row["target_uniprot_id"])
        if row.get("activity_type"):
            item["activity_types"].add(row["activity_type"])
        if row.get("clinical_or_approval_status"):
            for part in split_values(row["clinical_or_approval_status"]):
                item["clinical"].add(part)
        try:
            value = float(str(row.get("activity_value_uM", "")).strip())
            if math.isfinite(value):
                item["activity_values"].append(value)
        except ValueError:
            pass
        if row.get("compound_source") == "structure_context":
            item["structure_context"] = True
    return agg


def apply_small_summaries(small_rows, small_fields, binary_rows, review_rows):
    small_by_id = {row["drug_id"]: row for row in small_rows}
    review_by_drug = defaultdict(list)
    for row in review_rows:
        if row.get("keep_for_binding_site_context") != "1":
            continue
        drug_id = row.get("ligand_pubchem_cid", "").strip()
        if not drug_id.isdigit():
            drug_id = stable_unmapped_id(row)
        review_by_drug[drug_id].append(row)

    agg = aggregate_binary(binary_rows)
    new_drug_ids = [drug_id for drug_id in sorted(agg) if drug_id not in small_by_id]
    pubchem_props = fetch_pubchem_properties([drug_id for drug_id in new_drug_ids if drug_id.isdigit()])

    output = list(small_rows)
    for drug_id in new_drug_ids:
        examples = review_by_drug.get(drug_id, [])
        example = examples[0] if examples else {}
        props = pubchem_props.get(drug_id, {}) if drug_id.isdigit() else {}
        row = {field: "" for field in small_fields}
        row["drug_id"] = drug_id
        row["compound_id_type"] = "PubChem CID" if drug_id.isdigit() else "unmapped_structure_ligand"
        row["drug_name"] = ligand_display_name(example) if example else drug_id
        row["drug_name_type"] = example.get("ligand_category", "structure_ligand_context") if example else "structure_ligand_context"
        row["compound_source"] = "structure_context"
        row["cid_mapping_method"] = example.get("cid_mapping_method", "mapped" if drug_id.isdigit() else "unmapped")
        row["evidence_levels"] = "structure_ligand_context"
        row["molecular_formula"] = str(props.get("MolecularFormula", ""))
        row["molecular_weight"] = str(props.get("MolecularWeight", ""))
        row["canonical_smiles"] = str(props.get("CanonicalSMILES", ""))
        row["isomeric_smiles"] = str(props.get("IsomericSMILES", ""))
        row["inchikey"] = str(props.get("InChIKey", ""))
        row["iupac_name"] = str(props.get("IUPACName", ""))
        row["xlogp"] = str(props.get("XLogP", ""))
        row["tpsa"] = str(props.get("TPSA", ""))
        row["hbond_donor_count"] = str(props.get("HBondDonorCount", ""))
        row["hbond_acceptor_count"] = str(props.get("HBondAcceptorCount", ""))
        row["rotatable_bond_count"] = str(props.get("RotatableBondCount", ""))
        row["complexity"] = str(props.get("Complexity", ""))
        row["formal_charge"] = str(props.get("FormalCharge", ""))
        output.append(row)
        small_by_id[drug_id] = row

    for row in output:
        drug_id = row.get("drug_id", "")
        item = agg.get(drug_id)
        if not item:
            continue
        row["source_databases"] = join_values(item["sources"]) or row.get("source_databases", "")
        row["source_row_count"] = str(item["row_count"])
        row["unique_target_count"] = str(len(item["targets"]))
        row["activity_types"] = join_values(item["activity_types"]) or row.get("activity_types", "")
        row["clinical_statuses"] = join_values(item["clinical"]) or row.get("clinical_statuses", "")
        if item["activity_values"]:
            values = sorted(item["activity_values"])
            row["activity_value_uM_min"] = f"{values[0]:.6g}"
            row["activity_value_uM_median"] = f"{statistics.median(values):.6g}"
        row["compound_source"] = append_semicolon(row.get("compound_source", ""), item["compound_sources"])
        if item["structure_context"]:
            row["evidence_levels"] = append_semicolon(row.get("evidence_levels", ""), ["structure_ligand_context"])
    return output, new_drug_ids


def apply_protein_summaries(protein_rows, binary_rows):
    by_target = defaultdict(lambda: {
        "drug_ids": set(),
        "sources": set(),
        "evidence": set(),
        "has_match": False,
    })
    for row in binary_rows:
        uid = row.get("target_uniprot_id", "")
        if not uid:
            continue
        item = by_target[uid]
        if row.get("drug_id"):
            item["drug_ids"].add(row["drug_id"])
        if row.get("source_database"):
            item["sources"].add(row["source_database"])
        if row.get("compound_source") == "structure_context":
            item["evidence"].add("structure_ligand_context")
        for flag in [
            "has_pdb_biolip_matched",
            "has_pdb_scpdb_matched",
            "has_pdbbind_matched",
            "has_stitch_matched",
            "has_exp_binding_site",
        ]:
            if str(row.get(flag, "")).strip() == "1":
                item["has_match"] = True
    for row in protein_rows:
        item = by_target.get(row["target_uniprot_id"])
        if not item:
            row["unique_drug_count"] = "0"
            row["unique_drug_cid_sample"] = ""
            continue
        row["unique_drug_count"] = str(len(item["drug_ids"]))
        row["unique_drug_cid_sample"] = ";".join(sorted(item["drug_ids"])[:50])
        row["source_databases"] = join_values(item["sources"])
        row["evidence_levels"] = append_semicolon(row.get("evidence_levels", ""), item["evidence"])
        row["has_matched_evidence"] = "1" if item["has_match"] else "0"
    return protein_rows


def main():
    binary_rows = load_tsv(BINARY)
    binary_fields = list(binary_rows[0].keys())
    protein_rows = load_tsv(PROTEIN)
    protein_fields = list(protein_rows[0].keys())
    small_rows = load_tsv(SMALL)
    small_fields = list(small_rows[0].keys())
    review_rows = load_tsv(REVIEW)

    protein_by_uid = {row["target_uniprot_id"]: row for row in protein_rows}
    small_by_id = {row["drug_id"]: row for row in small_rows}

    base_binary_rows = [
        row for row in binary_rows
        if not (
            row.get("compound_source") == "structure_context"
            and "structure_ligand_context" in row.get("assay_or_mechanism", "")
        )
    ]
    merge_review_rows = [
        row for row in review_rows
        if row.get("keep_for_binding_site_context") == "1"
    ]
    structure_binary_rows = [
        make_binary_row(row, binary_fields, protein_by_uid, small_by_id)
        for row in merge_review_rows
    ]
    final_binary_rows = base_binary_rows + structure_binary_rows

    final_small_rows, new_drug_ids = apply_small_summaries(
        small_rows, small_fields, final_binary_rows, merge_review_rows
    )
    final_protein_rows = apply_protein_summaries(protein_rows, final_binary_rows)

    write_tsv(BINARY, binary_fields, final_binary_rows)
    write_tsv(SMALL, small_fields, final_small_rows)
    write_tsv(PROTEIN, protein_fields, final_protein_rows)
    update_manifest([BINARY, SMALL, PROTEIN])

    report = {
        "binary_rows_before": len(binary_rows),
        "binary_existing_structure_context_rows_removed": len(binary_rows) - len(base_binary_rows),
        "structure_context_rows_added": len(structure_binary_rows),
        "structure_context_rows_added_mapped_pubchem": sum(1 for row in structure_binary_rows if row["compound_id_type"] == "PubChem CID"),
        "structure_context_rows_added_unmapped": sum(1 for row in structure_binary_rows if row["compound_id_type"] == "unmapped_structure_ligand"),
        "binary_rows_after": len(final_binary_rows),
        "small_molecule_rows_before": len(small_rows),
        "small_molecule_rows_added": len(new_drug_ids),
        "small_molecule_rows_after": len(final_small_rows),
        "new_pubchem_cid_rows": sum(1 for drug_id in new_drug_ids if drug_id.isdigit()),
        "new_unmapped_structure_ligand_rows": sum(1 for drug_id in new_drug_ids if not drug_id.isdigit()),
        "new_drug_id_sample": new_drug_ids[:50],
        "protein_rows_updated": len(final_protein_rows),
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    update_manifest([REPORT, PUBCHEM_PROPS_CACHE])
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
