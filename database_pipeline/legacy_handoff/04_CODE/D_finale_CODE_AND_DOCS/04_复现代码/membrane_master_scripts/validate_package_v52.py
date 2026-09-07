from __future__ import annotations

import csv
import hashlib
import json
import shutil
import sys
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
RELEASE = ROOT / "releases" / "release_v5_2_evidence_tiered"
PREVIOUS = ROOT / "releases" / "release_v5_1"
OUTPUT = Path(r"C:\Users\Administrator\outputs\membrane_v52_20260724")
DESKTOP = Path(r"C:\Users\Administrator\Desktop\semifinal")
WORKBOOK_NAME = "human_membrane_protein_master_v5_2_evidence_tiered.xlsx"
WORKBOOK_SOURCE = OUTPUT / WORKBOOK_NAME
AUDIT = RELEASE / "human_membrane_audit_master_v5_2.tsv"
EXPECTED_SHEETS = [
    "Summary",
    "E1 Experimental",
    "E2 Supported",
    "E1+E2 Default",
    "E3 Candidates",
    "E0 Excluded",
    "Evidence Rules",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def workbook_sheet_names(path: Path) -> list[str]:
    namespace = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        tree = ET.fromstring(archive.read("xl/workbook.xml"))
    return [sheet.attrib["name"] for sheet in tree.findall("m:sheets/m:sheet", namespace)]


def verify_previous_release() -> dict[str, object]:
    manifest = PREVIOUS / "SHA256SUMS_v5_1.txt"
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


def count_file_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return sum(1 for _ in handle) - 1


def main() -> int:
    RELEASE.mkdir(parents=True, exist_ok=True)
    (RELEASE / "reports").mkdir(parents=True, exist_ok=True)
    (RELEASE / "reproducibility" / "scripts").mkdir(parents=True, exist_ok=True)
    DESKTOP.mkdir(parents=True, exist_ok=True)

    release_workbook = RELEASE / WORKBOOK_NAME
    shutil.copy2(WORKBOOK_SOURCE, release_workbook)

    rows = read_tsv(AUDIT)
    ids = [row["membrane_protein_id"] for row in rows]
    levels = Counter(row["evidence_level_v52"] for row in rows)
    default_rows = [row for row in rows if row["website_default_v52"] == "1"]
    default_classes = Counter(row["membrane_class_v52"] for row in default_rows)
    direct_flags = Counter(
        (row["evidence_level_v52"], row["direct_experimental_flag_v52"]) for row in rows
    )
    class_by_level = {
        level: dict(
            sorted(
                Counter(
                    row["membrane_class_v52"]
                    for row in rows
                    if row["evidence_level_v52"] == level
                ).items()
            )
        )
        for level in ("E0", "E1", "E2", "E3")
    }

    expected_level_counts = {"E0": 441, "E1": 3451, "E2": 4453, "E3": 2652}
    expected_default_classes = {"A": 5490, "B": 440, "C": 1974}
    subset_files = {
        "E1": RELEASE / "human_membrane_experimental_core_E1_v5_2.tsv",
        "E2": RELEASE / "human_membrane_strong_supported_E2_v5_2.tsv",
        "E1_E2": RELEASE / "human_membrane_default_E1_E2_v5_2.tsv",
        "E3": RELEASE / "human_membrane_prediction_candidates_E3_v5_2.tsv",
        "E0": RELEASE / "human_membrane_excluded_E0_v5_2.tsv",
        "A_default": RELEASE / "human_membrane_class_A_E1_E2_v5_2.tsv",
        "B_default": RELEASE / "human_membrane_class_B_E1_E2_v5_2.tsv",
        "C_default": RELEASE / "human_membrane_class_C_E1_E2_v5_2.tsv",
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

    sheets = workbook_sheet_names(release_workbook)
    previous_release = verify_previous_release()
    checks = {
        "audit_rows_10997": len(rows) == 10997,
        "unique_membrane_protein_ids": len(ids) == len(set(ids)),
        "evidence_levels_partition_audit": sum(levels.values()) == len(rows),
        "evidence_level_counts_match": dict(levels) == expected_level_counts,
        "default_rows_7904": len(default_rows) == 7904,
        "default_class_counts_match": dict(default_classes) == expected_default_classes,
        "default_flag_only_E1_E2": all(
            (row["website_default_v52"] == "1")
            == (row["evidence_level_v52"] in {"E1", "E2"})
            for row in rows
        ),
        "E1_all_direct_experimental": all(
            row["direct_experimental_flag_v52"] == "1"
            for row in rows
            if row["evidence_level_v52"] == "E1"
        ),
        "non_E1_not_direct_experimental": all(
            row["direct_experimental_flag_v52"] == "0"
            for row in rows
            if row["evidence_level_v52"] != "E1"
        ),
        "included_levels_have_ABC_class": all(
            row["membrane_class_v52"] in {"A", "B", "C"}
            for row in rows
            if row["evidence_level_v52"] in {"E1", "E2", "E3"}
        ),
        "E0_has_unknown_class": all(
            row["membrane_class_v52"] == "unknown"
            for row in rows
            if row["evidence_level_v52"] == "E0"
        ),
        "subset_file_counts_match": subset_counts == expected_subset_counts,
        "workbook_is_valid_zip": zipfile.is_zipfile(release_workbook),
        "workbook_sheet_names_match": sheets == EXPECTED_SHEETS,
        "previous_v5_1_release_unchanged": bool(previous_release["all_match"]),
    }

    required_files = [
        "README_v5_2.md",
        "EVIDENCE_POLICY_v5_2.md",
        "DATA_DICTIONARY_v5_2.md",
        "SOURCE_REGISTRY_v5_2.tsv",
        WORKBOOK_NAME,
        *[path.name for path in subset_files.values()],
        AUDIT.name,
    ]
    checks["required_release_files_present"] = all(
        (RELEASE / relative).exists() for relative in required_files
    )

    report = {
        "release": "v5.2-evidence-tiered",
        "validation_date": date.today().isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": {
            "audit_rows": len(rows),
            "evidence_levels": dict(sorted(levels.items())),
            "default_rows": len(default_rows),
            "default_classes": dict(sorted(default_classes.items())),
            "class_by_evidence_level": class_by_level,
            "direct_experimental_flags": {
                f"{level}:{flag}": count
                for (level, flag), count in sorted(direct_flags.items())
            },
            "subset_files": subset_counts,
        },
        "workbook": {
            "name": release_workbook.name,
            "sha256": sha256(release_workbook),
            "sheets": sheets,
        },
        "previous_release_integrity": previous_release,
    }

    report_path = RELEASE / "reports" / "V52_VALIDATION_REPORT.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    shutil.copy2(Path(__file__), RELEASE / "reproducibility" / "scripts" / Path(__file__).name)

    manifest_files = sorted(
        path
        for path in RELEASE.rglob("*")
        if path.is_file()
        and path.name not in {"SHA256SUMS_v5_2.txt", "RELEASE_MANIFEST_v5_2.json"}
    )
    manifest = {
        "release": "v5.2-evidence-tiered",
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
    (RELEASE / "RELEASE_MANIFEST_v5_2.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    checksum_files = sorted(path for path in RELEASE.rglob("*") if path.is_file())
    checksum_lines = [
        f"{sha256(path)}  {path.relative_to(RELEASE).as_posix()}"
        for path in checksum_files
        if path.name != "SHA256SUMS_v5_2.txt"
    ]
    (RELEASE / "SHA256SUMS_v5_2.txt").write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8"
    )

    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
