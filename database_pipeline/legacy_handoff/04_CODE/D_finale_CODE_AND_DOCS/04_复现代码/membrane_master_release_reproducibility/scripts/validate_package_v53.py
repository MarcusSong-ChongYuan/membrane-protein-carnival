from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import sys
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
RELEASE = ROOT / "releases" / "release_v5_3_sequences"
PREVIOUS = ROOT / "releases" / "release_v5_2_evidence_tiered"
OUTPUT = Path(r"C:\Users\Administrator\outputs\membrane_v53_sequences_20260724")
DESKTOP = Path(r"C:\Users\Administrator\Desktop\semifinal")
WORKBOOK_NAME = "human_membrane_protein_master_v5_3_with_sequences.xlsx"
WORKBOOK_SOURCE = OUTPUT / WORKBOOK_NAME
AUDIT = RELEASE / "human_membrane_audit_master_v5_3.tsv"
EXPECTED_SHEETS = [
    "Summary",
    "E1 Experimental",
    "E2 Supported",
    "E1+E2 Default",
    "E3 Candidates",
    "E0 Excluded",
    "Sequence & Evidence Rules",
]
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_file_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return sum(1 for _ in handle) - 1


def verify_previous_release() -> dict[str, object]:
    manifest = PREVIOUS / "SHA256SUMS_v5_2.txt"
    checks: list[dict[str, object]] = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        path = PREVIOUS / relative.strip()
        actual = sha256(path) if path.exists() else None
        checks.append(
            {
                "path": relative.strip(),
                "exists": path.exists(),
                "matches": actual == expected,
            }
        )
    return {
        "manifest": manifest.name,
        "checked_files": len(checks),
        "all_match": all(item["exists"] and item["matches"] for item in checks),
        "failures": [item for item in checks if not (item["exists"] and item["matches"])],
    }


def parse_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current = ""
    with path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                match = re.match(r">[^|]+\|([^|]+)\|", line)
                if not match:
                    raise ValueError(f"Invalid FASTA header in {path.name}: {line}")
                current = match.group(1)
                sequences[current] = []
            elif current:
                sequences[current].append(line)
            else:
                raise ValueError(f"Sequence before header in {path.name}")
    return {accession: "".join(parts) for accession, parts in sequences.items()}


def workbook_sheet_targets(
    archive: zipfile.ZipFile,
) -> tuple[list[str], dict[str, str]]:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    rel_targets = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")
    }
    names: list[str] = []
    targets: dict[str, str] = {}
    for sheet in workbook.findall(f".//{{{MAIN_NS}}}sheet"):
        name = sheet.attrib["name"]
        relation_id = sheet.attrib[f"{{{REL_NS}}}id"]
        target = rel_targets[relation_id].replace("\\", "/")
        if target.startswith("/"):
            target = target.lstrip("/")
        elif not target.startswith("xl/"):
            target = "xl/" + target
        names.append(name)
        targets[name] = target
    return names, targets


def cell_text(cell: ET.Element) -> str:
    cell_type = cell.attrib.get("t", "")
    if cell_type == "inlineStr":
        return "".join(
            node.text or "" for node in cell.findall(f".//{{{MAIN_NS}}}t")
        )
    value = cell.find(f"{{{MAIN_NS}}}v")
    return value.text if value is not None and value.text is not None else ""


def worksheet_row(
    archive: zipfile.ZipFile, target: str, row_number: int
) -> dict[str, str]:
    root = ET.fromstring(archive.read(target))
    row = root.find(f".//{{{MAIN_NS}}}row[@r='{row_number}']")
    if row is None:
        raise ValueError(f"Row {row_number} not found in {target}")
    return {
        re.match(r"[A-Z]+", cell.attrib["r"]).group(0): cell_text(cell)
        for cell in row.findall(f"{{{MAIN_NS}}}c")
    }


def inspect_workbook(path: Path, audit_rows: list[dict[str, str]]) -> dict[str, object]:
    titin = next(row for row in audit_rows if row["target_uniprot_id"] == "Q8WZ42")
    evidence_level = titin["evidence_level_v52"]
    sheet_name = {
        "E1": "E1 Experimental",
        "E2": "E2 Supported",
        "E3": "E3 Candidates",
        "E0": "E0 Excluded",
    }[evidence_level]
    subset_file = {
        "E1": RELEASE / "human_membrane_experimental_core_E1_v5_3.tsv",
        "E2": RELEASE / "human_membrane_strong_supported_E2_v5_3.tsv",
        "E3": RELEASE / "human_membrane_prediction_candidates_E3_v5_3.tsv",
        "E0": RELEASE / "human_membrane_excluded_E0_v5_3.tsv",
    }[evidence_level]
    subset = read_tsv(subset_file)
    titin_row_number = next(
        index + 2
        for index, row in enumerate(subset)
        if row["target_uniprot_id"] == "Q8WZ42"
    )

    with zipfile.ZipFile(path) as archive:
        sheets, targets = workbook_sheet_targets(archive)
        summary_row = worksheet_row(archive, targets["Summary"], 10)
        titin_excel = worksheet_row(
            archive, targets[sheet_name], titin_row_number
        )
    combined = titin_excel.get("L", "") + titin_excel.get("M", "")
    return {
        "sheet_names": sheets,
        "summary_excel_split_records": summary_row.get("E", ""),
        "titin_sheet": sheet_name,
        "titin_row": titin_row_number,
        "titin_accession": titin_excel.get("B", ""),
        "titin_sequence_length_cell": titin_excel.get("H", ""),
        "titin_part_1_length": len(titin_excel.get("L", "")),
        "titin_part_2_length": len(titin_excel.get("M", "")),
        "titin_combined_length": len(combined),
        "titin_combined_sha256": sequence_sha256(combined),
        "titin_matches_full_sequence": combined == titin["canonical_sequence"],
    }


def main() -> int:
    DESKTOP.mkdir(parents=True, exist_ok=True)
    release_workbook = RELEASE / WORKBOOK_NAME
    shutil.copy2(WORKBOOK_SOURCE, release_workbook)
    shutil.copy2(
        ROOT / "v53_sequence_working" / "uniprot_sequence_retrieval_report_v53.json",
        RELEASE / "reports" / "UNIPROT_SEQUENCE_RETRIEVAL_REPORT_v5_3.json",
    )
    shutil.copy2(
        OUTPUT / "build_v53_workbook.mjs",
        RELEASE / "reproducibility" / "scripts" / "build_v53_workbook.mjs",
    )
    shutil.copy2(
        Path(__file__),
        RELEASE / "reproducibility" / "scripts" / Path(__file__).name,
    )

    rows = read_tsv(AUDIT)
    ids = [row["target_uniprot_id"] for row in rows]
    levels = Counter(row["evidence_level_v52"] for row in rows)
    default_rows = [
        row for row in rows if row["evidence_level_v52"] in {"E1", "E2"}
    ]
    invalid_characters = [
        row["target_uniprot_id"]
        for row in rows
        if not re.fullmatch(r"[A-Z]+", row["canonical_sequence"])
    ]
    bad_sha = [
        row["target_uniprot_id"]
        for row in rows
        if sequence_sha256(row["canonical_sequence"]) != row["sequence_sha256_v53"]
    ]
    bad_length = [
        row["target_uniprot_id"]
        for row in rows
        if len(row["canonical_sequence"]) != int(row["sequence_length_v53"])
    ]
    prior_length_mismatch = [
        row["target_uniprot_id"]
        for row in rows
        if row["sequence_length"] != row["sequence_length_v53"]
    ]

    subset_files = {
        "E1": RELEASE / "human_membrane_experimental_core_E1_v5_3.tsv",
        "E2": RELEASE / "human_membrane_strong_supported_E2_v5_3.tsv",
        "E1_E2": RELEASE / "human_membrane_default_E1_E2_v5_3.tsv",
        "E3": RELEASE / "human_membrane_prediction_candidates_E3_v5_3.tsv",
        "E0": RELEASE / "human_membrane_excluded_E0_v5_3.tsv",
        "A_default": RELEASE / "human_membrane_class_A_E1_E2_v5_3.tsv",
        "B_default": RELEASE / "human_membrane_class_B_E1_E2_v5_3.tsv",
        "C_default": RELEASE / "human_membrane_class_C_E1_E2_v5_3.tsv",
    }
    subset_counts = {name: count_file_rows(path) for name, path in subset_files.items()}
    expected_subset_counts = {
        "E1": 3451,
        "E2": 4453,
        "E1_E2": 7904,
        "E3": 2652,
        "E0": 441,
        "A_default": 5490,
        "B_default": 440,
        "C_default": 1974,
    }

    all_fasta = parse_fasta(RELEASE / "human_membrane_all_audit_v5_3.fasta")
    default_fasta = parse_fasta(
        RELEASE / "human_membrane_default_E1_E2_v5_3.fasta"
    )
    master_by_id = {row["target_uniprot_id"]: row for row in rows}
    all_fasta_mismatches = [
        accession
        for accession, sequence in all_fasta.items()
        if sequence != master_by_id[accession]["canonical_sequence"]
    ]
    default_ids = {row["target_uniprot_id"] for row in default_rows}
    default_fasta_mismatches = [
        accession
        for accession, sequence in default_fasta.items()
        if accession not in default_ids
        or sequence != master_by_id[accession]["canonical_sequence"]
    ]

    workbook = inspect_workbook(release_workbook, rows)
    previous_release = verify_previous_release()
    checks = {
        "audit_rows_10997": len(rows) == 10_997,
        "unique_accessions_10997": len(ids) == len(set(ids)) == 10_997,
        "evidence_counts_unchanged": dict(levels)
        == {"E0": 441, "E1": 3451, "E2": 4453, "E3": 2652},
        "all_sequences_present": all(bool(row["canonical_sequence"]) for row in rows),
        "all_sequences_valid_uppercase": not invalid_characters,
        "all_sequence_sha256_match": not bad_sha,
        "all_sequence_lengths_match_content": not bad_length,
        "all_sequence_lengths_match_v52": not prior_length_mismatch,
        "all_accessions_resolved_primary": all(
            row["sequence_status_v53"] == "retrieved_primary"
            and row["sequence_resolved_accession_v53"] == row["target_uniprot_id"]
            for row in rows
        ),
        "all_sequences_from_uniprot_2026_02": all(
            row["uniprot_release_v53"] == "2026_02" for row in rows
        ),
        "subset_counts_match": subset_counts == expected_subset_counts,
        "all_audit_fasta_count_10997": len(all_fasta) == 10_997,
        "all_audit_fasta_exact": not all_fasta_mismatches
        and set(all_fasta) == set(ids),
        "default_fasta_count_7904": len(default_fasta) == 7_904,
        "default_fasta_exact": not default_fasta_mismatches
        and set(default_fasta) == default_ids,
        "workbook_valid_zip": zipfile.is_zipfile(release_workbook),
        "workbook_sheet_names_match": workbook["sheet_names"] == EXPECTED_SHEETS,
        "workbook_summary_split_count_1": workbook[
            "summary_excel_split_records"
        ]
        in {"1", "1.0"},
        "workbook_titin_parts_32000_2350": workbook["titin_part_1_length"] == 32_000
        and workbook["titin_part_2_length"] == 2_350,
        "workbook_titin_full_sequence_exact": bool(
            workbook["titin_matches_full_sequence"]
        ),
        "previous_v5_2_release_unchanged": bool(previous_release["all_match"]),
    }

    report = {
        "release": "v5.3-sequence-enhanced",
        "validation_date": date.today().isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": {
            "audit_rows": len(rows),
            "unique_accessions": len(set(ids)),
            "evidence_levels": dict(sorted(levels.items())),
            "default_E1_E2_rows": len(default_rows),
            "total_amino_acids": sum(
                int(row["sequence_length_v53"]) for row in rows
            ),
            "longest_sequence_length": max(
                int(row["sequence_length_v53"]) for row in rows
            ),
            "subset_files": subset_counts,
            "all_fasta_records": len(all_fasta),
            "default_fasta_records": len(default_fasta),
        },
        "exceptions": {
            "invalid_sequence_characters": invalid_characters,
            "bad_sequence_sha256": bad_sha,
            "bad_sequence_length": bad_length,
            "prior_length_mismatch": prior_length_mismatch,
            "all_fasta_mismatches": all_fasta_mismatches,
            "default_fasta_mismatches": default_fasta_mismatches,
        },
        "workbook": {
            "name": release_workbook.name,
            "bytes": release_workbook.stat().st_size,
            "sha256": sha256(release_workbook),
            **workbook,
        },
        "previous_release_integrity": previous_release,
    }

    report_path = RELEASE / "reports" / "V53_VALIDATION_REPORT.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    manifest_files = sorted(
        path
        for path in RELEASE.rglob("*")
        if path.is_file()
        and path.name not in {"SHA256SUMS_v5_3.txt", "RELEASE_MANIFEST_v5_3.json"}
    )
    manifest = {
        "release": "v5.3-sequence-enhanced",
        "release_date": date.today().isoformat(),
        "validation_status": report["status"],
        "file_count": len(manifest_files),
        "files": [
            {
                "path": path.relative_to(RELEASE).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in manifest_files
        ],
    }
    (RELEASE / "RELEASE_MANIFEST_v5_3.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    checksum_files = sorted(path for path in RELEASE.rglob("*") if path.is_file())
    checksum_lines = [
        f"{sha256(path)}  {path.relative_to(RELEASE).as_posix()}"
        for path in checksum_files
        if path.name != "SHA256SUMS_v5_3.txt"
    ]
    (RELEASE / "SHA256SUMS_v5_3.txt").write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8"
    )

    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
