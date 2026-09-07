from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path


ROOT = Path(r"D:\7.22")
INDEX = (
    ROOT
    / "membrane_master_v5_working"
    / "releases"
    / "release_v5_4_pgd_binding"
    / "small_molecule_index_v4_1.tsv"
)
CHEMREPS = (
    ROOT
    / "small_molecule_v1_working"
    / "raw"
    / "chembl37"
    / "chembl_37_chemreps_full.txt.gz"
)
METADATA = (
    ROOT
    / "small_molecule_v1_working"
    / "raw"
    / "chembl37"
    / "chembl37_selected_metadata.jsonl"
)
OUT = (
    ROOT
    / "small_molecule_v1_working"
    / "raw"
    / "chembl37"
    / "chembl37_selected_chemreps.tsv"
)
REPORT = (
    ROOT
    / "small_molecule_v1_working"
    / "reports"
    / "CHEMBL37_STRUCTURE_EXTRACTION_REPORT.json"
)
CHEMBL_PATTERN = re.compile(r"CHEMBL\d+", re.IGNORECASE)


def extract_ids(value: str) -> set[str]:
    return {match.upper() for match in CHEMBL_PATTERN.findall(value or "")}


def main() -> None:
    wanted: set[str] = set()
    with INDEX.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            wanted.update(extract_ids(row.get("compound_id", "")))
            wanted.update(extract_ids(row.get("chembl_ids", "")))
    if METADATA.exists():
        with METADATA.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    wanted.update(extract_ids(json.loads(line).get("parent_chembl_id", "")))

    found: set[str] = set()
    with gzip.open(CHEMREPS, "rt", encoding="utf-8", newline="") as source, OUT.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        fields = [
            "chembl_id",
            "canonical_smiles",
            "standard_inchi",
            "standard_inchi_key",
        ]
        writer = csv.DictWriter(
            destination,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            chembl_id = row["chembl_id"].upper()
            if chembl_id not in wanted:
                continue
            writer.writerow({field: row.get(field, "") for field in fields})
            found.add(chembl_id)

    report = {
        "wanted_chembl_ids": len(wanted),
        "found_chembl_structures": len(found),
        "missing_chembl_ids": len(wanted - found),
        "missing_examples": sorted(wanted - found)[:100],
        "source_file": str(CHEMREPS),
        "source_version": "ChEMBL 37",
        "source_sha256": "ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7",
        "output_file": str(OUT),
    }
    REPORT.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
