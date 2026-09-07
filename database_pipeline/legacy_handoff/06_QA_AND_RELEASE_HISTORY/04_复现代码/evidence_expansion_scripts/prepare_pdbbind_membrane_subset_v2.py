#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import re
import tarfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
SOURCE_DIR = ROOT / "raw" / "manual_downloads" / "pdbbind_2020r1"
INDEX_FILE = (
    SOURCE_DIR
    / "extracted_index"
    / "index"
    / "INDEX_general_PL.2020R1.lst"
)
ARCHIVE = SOURCE_DIR / "P-L.tar.gz"
PROTEIN_MASTER = Path(
    r"D:\7.22\small_molecule_v1_working\releases\release_small_molecule_v1_0"
    r"\human_membrane_audit_master_v5_5.tsv"
)
SUBSET_DIR = ROOT / "raw" / "pdbbind_2020r1" / "membrane_subset"
MANIFEST = ROOT / "intermediate" / "pdbbind_2020r1_membrane_subset_manifest.tsv"
REPORT = ROOT / "reports" / "PDBBIND_2020R1_MEMBRANE_SUBSET_REPORT.json"
LOG = ROOT / "logs" / "pdbbind_2020r1_subset.log"


def log(message: str) -> None:
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def load_membrane_pdb_map() -> dict[str, dict[str, set[str]]]:
    mapping: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"uniprot_ids": set(), "symbols": set()}
    )
    with PROTEIN_MASTER.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("membrane_class_v52") not in {"A", "B", "C"}:
                continue
            for pdb_id in re.split(r"[;,\s]+", row.get("pdb_ids", "")):
                pdb_id = pdb_id.strip().lower()
                if len(pdb_id) != 4:
                    continue
                mapping[pdb_id]["uniprot_ids"].add(row["target_uniprot_id"].strip())
                symbol = row.get("approved_symbol", "").strip()
                if symbol:
                    mapping[pdb_id]["symbols"].add(symbol)
    return mapping


def load_index() -> dict[str, dict[str, str]]:
    records = {}
    with INDEX_FILE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            pdb_id = parts[0].lower()
            ligand_match = re.search(r"\(([^()]*)\)\s*(.*)$", line)
            records[pdb_id] = {
                "pdb_id": pdb_id,
                "resolution": parts[1],
                "release_year": parts[2],
                "binding_data_raw": parts[3],
                "ligand_annotation": ligand_match.group(1).strip()
                if ligand_match
                else "",
                "index_comment": ligand_match.group(2).strip()
                if ligand_match
                else "",
                "index_line_raw": line.rstrip("\n"),
            }
    return records


def archive_pdb_id(member_name: str) -> str:
    parts = Path(member_name).parts
    if len(parts) >= 3 and parts[0] == "P-L":
        candidate = parts[2].lower()
        if len(candidate) == 4:
            return candidate
    return ""


def main() -> None:
    SUBSET_DIR.mkdir(parents=True, exist_ok=True)
    pdb_map = load_membrane_pdb_map()
    index = load_index()
    selected = set(pdb_map) & set(index)
    log(
        f"{datetime.now(timezone.utc).isoformat()} "
        f"selected_pdb_ids={len(selected)}"
    )

    extracted_files = 0
    extracted_bytes = 0
    extracted_pdb_ids = set()
    with tarfile.open(ARCHIVE, "r:gz") as archive:
        for member_number, member in enumerate(archive, 1):
            pdb_id = archive_pdb_id(member.name)
            if pdb_id not in selected:
                continue
            archive.extract(member, path=SUBSET_DIR, filter="data")
            extracted_pdb_ids.add(pdb_id)
            if member.isfile():
                extracted_files += 1
                extracted_bytes += member.size
            if extracted_files and extracted_files % 5000 == 0:
                log(
                    f"{datetime.now(timezone.utc).isoformat()} "
                    f"files={extracted_files} pdb_ids={len(extracted_pdb_ids)}"
                )

    fields = [
        "pdb_id",
        "target_uniprot_ids",
        "approved_symbols",
        "resolution",
        "release_year",
        "binding_data_raw",
        "ligand_annotation",
        "index_comment",
        "structure_directory",
        "index_line_raw",
    ]
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for pdb_id in sorted(selected):
            record = index[pdb_id]
            protein = pdb_map[pdb_id]
            writer.writerow(
                {
                    **record,
                    "target_uniprot_ids": ";".join(sorted(protein["uniprot_ids"])),
                    "approved_symbols": ";".join(sorted(protein["symbols"])),
                    "structure_directory": str(
                        next(SUBSET_DIR.glob(f"P-L/*/{pdb_id}"), "")
                    ),
                }
            )

    report = {
        "status": "passed"
        if extracted_pdb_ids == selected
        else "incomplete",
        "pdbbind_protein_ligand_index_records": len(index),
        "membrane_master_pdb_ids": len(pdb_map),
        "selected_intersecting_pdb_ids": len(selected),
        "selected_membrane_proteins": len(
            {
                uniprot
                for pdb_id in selected
                for uniprot in pdb_map[pdb_id]["uniprot_ids"]
            }
        ),
        "extracted_pdb_ids": len(extracted_pdb_ids),
        "extracted_files": extracted_files,
        "extracted_uncompressed_bytes": extracted_bytes,
        "missing_selected_pdb_ids": sorted(selected - extracted_pdb_ids),
        "source_archive": str(ARCHIVE),
        "manifest": str(MANIFEST),
        "subset_directory": str(SUBSET_DIR),
    }
    with REPORT.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    log(f"{datetime.now(timezone.utc).isoformat()} complete {json.dumps(report)}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
