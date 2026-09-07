from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
V5 = ROOT / "releases" / "release_v5"
RELEASE = ROOT / "releases" / "release_v5_1"
REPORT = RELEASE / "reports" / "V51_VALIDATION_REPORT.json"
MANIFEST = RELEASE / "RELEASE_MANIFEST_v5_1.json"
CHECKSUMS = RELEASE / "SHA256SUMS_v5_1.txt"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check(condition: bool, name: str, detail: object, checks: list[dict]) -> None:
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


def validate_v5_immutable(checks: list[dict]) -> None:
    checksum_file = V5 / "SHA256SUMS_v5.txt"
    mismatches: list[str] = []
    missing: list[str] = []
    for line in checksum_file.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        path = V5 / relative
        if not path.exists():
            missing.append(relative)
        elif sha256(path) != expected:
            mismatches.append(relative)
    check(
        not missing and not mismatches,
        "immutable_v5_release_checksums",
        {"missing": missing, "mismatched": mismatches},
        checks,
    )


def workbook_sheet_names(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as archive:
        archive.testzip()
        xml = archive.read("xl/workbook.xml")
    root = ET.fromstring(xml)
    namespace = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    return [
        element.attrib["name"]
        for element in root.findall("m:sheets/m:sheet", namespace)
    ]


def main() -> None:
    checks: list[dict] = []
    required = [
        "human_membrane_protein_audit_master_v5_1.tsv",
        "human_membrane_protein_master_v5_1.tsv",
        "human_integral_membrane_view_v5_1.tsv",
        "human_lipid_anchored_membrane_view_v5_1.tsv",
        "human_peripheral_membrane_view_v5_1.tsv",
        "human_membrane_candidate_view_v5_1.tsv",
        "human_membrane_excluded_audit_v5_1.tsv",
        "human_membrane_protein_master_v5_1.xlsx",
        "review/r1109_final_review_v51.tsv",
        "normalized/pdbtm_target_chain_topology_v51.tsv",
        "README_v5_1.md",
        "METHODS_v5_1.md",
        "INCLUSION_POLICY_v5_1.md",
        "DATA_DICTIONARY_v5_1.md",
        "REVIEW_REPORT_v5_1.md",
        "SOURCE_REGISTRY_v5_1.tsv",
    ]
    missing = [name for name in required if not (RELEASE / name).exists()]
    check(not missing, "required_files_present", {"missing": missing}, checks)

    audit = read_tsv(RELEASE / "human_membrane_protein_audit_master_v5_1.tsv")
    included = read_tsv(RELEASE / "human_membrane_protein_master_v5_1.tsv")
    class_a = read_tsv(RELEASE / "human_integral_membrane_view_v5_1.tsv")
    class_b = read_tsv(RELEASE / "human_lipid_anchored_membrane_view_v5_1.tsv")
    class_c = read_tsv(RELEASE / "human_peripheral_membrane_view_v5_1.tsv")
    candidates = read_tsv(RELEASE / "human_membrane_candidate_view_v5_1.tsv")
    excluded = read_tsv(RELEASE / "human_membrane_excluded_audit_v5_1.tsv")
    review = read_tsv(RELEASE / "review" / "r1109_final_review_v51.tsv")

    def ids(rows: list[dict[str, str]]) -> list[str]:
        return [row["target_uniprot_id"] for row in rows]

    check(len(audit) == 10997, "audit_master_row_count", len(audit), checks)
    check(len(set(ids(audit))) == len(audit), "audit_master_unique_accessions", len(set(ids(audit))), checks)
    check(len(included) == 10076, "published_master_row_count", len(included), checks)
    check(len(set(ids(included))) == len(included), "published_master_unique_accessions", len(set(ids(included))), checks)
    check(
        (len(class_a), len(class_b), len(class_c)) == (5505, 440, 4131),
        "published_class_counts",
        {"A": len(class_a), "B": len(class_b), "C": len(class_c)},
        checks,
    )
    check(
        set(ids(included)) == set(ids(class_a)) | set(ids(class_b)) | set(ids(class_c)),
        "class_views_partition_published_master",
        True,
        checks,
    )
    check(
        not (set(ids(class_a)) & set(ids(class_b)))
        and not (set(ids(class_a)) & set(ids(class_c)))
        and not (set(ids(class_b)) & set(ids(class_c))),
        "class_views_are_mutually_exclusive",
        True,
        checks,
    )
    check(len(review) == 1109, "former_R_review_row_count", len(review), checks)
    check(len(set(ids(review))) == 1109, "former_R_unique_accessions", len(set(ids(review))), checks)
    check(
        (len(candidates), len(excluded)) == (480, 441),
        "candidate_and_excluded_counts",
        {"candidate": len(candidates), "excluded": len(excluded)},
        checks,
    )
    included_r = [row for row in review if row["inclusion_status_v51"] == "included"]
    check(len(included_r) == 188, "former_R_included_count", len(included_r), checks)
    check(
        len(included_r) + len(candidates) + len(excluded) == 1109,
        "former_R_outcomes_cover_all_rows",
        len(included_r) + len(candidates) + len(excluded),
        checks,
    )
    check(
        not (set(ids(candidates)) & set(ids(excluded))),
        "candidate_excluded_disjoint",
        True,
        checks,
    )
    check(
        all(row["final_membrane_class_v51"] in {"A", "B", "C"} for row in included),
        "published_rows_have_final_ABC_class",
        True,
        checks,
    )
    check(
        all(
            row["final_membrane_class_v51"] == "unknown"
            and row["evidence_status_v51"] == "uncertain"
            for row in candidates
        ),
        "candidate_rows_are_unknown_uncertain",
        True,
        checks,
    )
    check(
        all(
            row["final_membrane_class_v51"] == "unknown"
            and row["evidence_status_v51"] == "excluded"
            for row in excluded
        ),
        "excluded_rows_are_unknown_excluded",
        True,
        checks,
    )
    check(
        sum(row["pdbtm_target_chain_tm_flag"] == "1" for row in review) == 0,
        "no_former_R_target_is_a_PDBTM_TM_chain",
        0,
        checks,
    )
    check(
        all(row["decision_code_v51"] and row["decision_basis_v51"] for row in review),
        "all_former_R_rows_have_auditable_decisions",
        True,
        checks,
    )

    workbook = RELEASE / "human_membrane_protein_master_v5_1.xlsx"
    expected_sheets = [
        "Summary",
        "Included",
        "Former R Audit",
        "Candidates",
        "Excluded",
        "Data Dictionary",
    ]
    try:
        sheets = workbook_sheet_names(workbook)
        workbook_ok = sheets == expected_sheets
        workbook_detail = {"sheets": sheets, "bytes": workbook.stat().st_size}
    except Exception as exc:  # pragma: no cover - release guard
        workbook_ok = False
        workbook_detail = {"error": repr(exc)}
    check(workbook_ok, "xlsx_structure_and_sheets", workbook_detail, checks)
    validate_v5_immutable(checks)

    status = "PASS" if all(item["passed"] for item in checks) else "FAIL"
    report = {
        "release": "v5.1",
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "checks": checks,
    }
    REPORT.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )
    if status != "PASS":
        failed = [item for item in checks if not item["passed"]]
        raise RuntimeError(f"v5.1 validation failed: {failed}")

    destination = RELEASE / "reproducibility" / "scripts" / Path(__file__).name
    shutil.copy2(Path(__file__), destination)

    excluded_from_manifest = {MANIFEST.name, CHECKSUMS.name}
    files = sorted(
        path
        for path in RELEASE.rglob("*")
        if path.is_file() and path.name not in excluded_from_manifest
    )
    entries = [
        {
            "path": path.relative_to(RELEASE).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in files
    ]
    MANIFEST.write_text(
        json.dumps(
            {
                "release": "v5.1",
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "file_count_excluding_manifests": len(entries),
                "files": entries,
            },
            indent=2,
            ensure_ascii=False,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    CHECKSUMS.write_text(
        "".join(f'{entry["sha256"]}  {entry["path"]}\n' for entry in entries),
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
    print(f"Wrote {MANIFEST}")
    print(f"Wrote {CHECKSUMS}")


if __name__ == "__main__":
    main()
