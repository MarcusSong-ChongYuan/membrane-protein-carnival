from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22")
RELEASE = ROOT / "membrane_master_v5_working" / "releases" / "release_v5_4_pgd_binding"
WORK = ROOT / "small_molecule_v1_working"
REPORTS = WORK / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)

INDEX = RELEASE / "small_molecule_index_v4_1.tsv"
EVIDENCE = RELEASE / "binding_evidence_master_v4_1.tsv"


def split_values(value: str) -> list[str]:
    if not value:
        return []
    normalized = value.replace("|", ";").replace(",", ";")
    return [item.strip() for item in normalized.split(";") if item.strip()]


def main() -> None:
    csv.field_size_limit(100_000_000)
    counts = Counter()
    type_counts = Counter()
    source_counts = Counter()
    scope_counts = Counter()
    qc_counts = Counter()
    inchikey_groups: dict[str, int] = Counter()
    pubchem_ids: set[str] = set()
    chembl_ids: set[str] = set()
    names: set[str] = set()
    compound_ids: set[str] = set()
    structure_by_type = Counter()
    missing_structure_by_type = Counter()
    mapped_both = 0

    with INDEX.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            counts["index_rows"] += 1
            compound_id = row.get("compound_id", "").strip()
            compound_ids.add(compound_id)
            compound_type = row.get("compound_id_type", "").strip() or "missing"
            type_counts[compound_type] += 1
            scope_counts[row.get("small_molecule_scope_status", "").strip() or "missing"] += 1
            qc_counts[row.get("record_qc_status", "").strip() or "missing"] += 1
            for source in split_values(row.get("source_databases", "")):
                source_counts[source] += 1
            smiles = row.get("canonical_smiles", "").strip()
            inchikey = row.get("inchikey", "").strip()
            if smiles:
                counts["rows_with_smiles"] += 1
                structure_by_type[compound_type] += 1
            else:
                missing_structure_by_type[compound_type] += 1
            if inchikey:
                counts["rows_with_inchikey"] += 1
                inchikey_groups[inchikey] += 1
            if row.get("compound_name", "").strip():
                counts["rows_with_name"] += 1
                names.add(row["compound_name"].strip().casefold())
            row_pubchem = set(split_values(row.get("pubchem_cids", "")))
            row_chembl = set(split_values(row.get("chembl_ids", "")))
            if compound_type == "PubChem CID" and compound_id:
                row_pubchem.add(compound_id.replace("PUBCHEM:", ""))
            if compound_type == "ChEMBL ID" and compound_id:
                row_chembl.add(compound_id.replace("CHEMBL:", ""))
            pubchem_ids.update(row_pubchem)
            chembl_ids.update(row_chembl)
            if row_pubchem and row_chembl:
                mapped_both += 1

    evidence_compounds: set[str] = set()
    evidence_default_compounds: set[str] = set()
    evidence_tiers = Counter()
    with EVIDENCE.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            counts["evidence_rows"] += 1
            compound_id = row.get("compound_id", "").strip()
            evidence_compounds.add(compound_id)
            evidence_tiers[row.get("evidence_tier", "").strip() or "missing"] += 1
            if row.get("default_release_inclusion", "").strip() == "1":
                evidence_default_compounds.add(compound_id)

    duplicate_inchikey_groups = {
        key: size for key, size in inchikey_groups.items() if size > 1
    }
    report = {
        "input_index": str(INDEX),
        "input_evidence": str(EVIDENCE),
        "index_rows": counts["index_rows"],
        "distinct_compound_ids": len(compound_ids),
        "rows_with_name": counts["rows_with_name"],
        "distinct_casefolded_names": len(names),
        "rows_with_smiles": counts["rows_with_smiles"],
        "rows_with_inchikey": counts["rows_with_inchikey"],
        "structure_coverage": (
            counts["rows_with_smiles"] / counts["index_rows"]
            if counts["index_rows"]
            else 0
        ),
        "pubchem_id_count": len(pubchem_ids),
        "chembl_id_count": len(chembl_ids),
        "rows_with_both_pubchem_and_chembl": mapped_both,
        "compound_id_type_counts": dict(type_counts),
        "structure_by_type": dict(structure_by_type),
        "missing_structure_by_type": dict(missing_structure_by_type),
        "source_database_counts": dict(source_counts),
        "scope_status_counts": dict(scope_counts),
        "record_qc_counts": dict(qc_counts),
        "inchikey_group_count": len(inchikey_groups),
        "duplicate_inchikey_group_count": len(duplicate_inchikey_groups),
        "rows_in_duplicate_inchikey_groups": sum(duplicate_inchikey_groups.values()),
        "evidence_rows": counts["evidence_rows"],
        "evidence_distinct_compounds": len(evidence_compounds),
        "default_evidence_distinct_compounds": len(evidence_default_compounds),
        "index_missing_evidence_compound_ids": len(evidence_compounds - compound_ids),
        "evidence_tier_counts": dict(evidence_tiers),
    }
    output = REPORTS / "SMALL_MOLECULE_V1_INPUT_INVENTORY.json"
    output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
