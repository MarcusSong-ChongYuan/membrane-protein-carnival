#!/usr/bin/env python3
"""Repair V6.0 legacy blockers, revalidate, and freeze formal V6.1."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
import sqlite3

import finalize_v61_release as v
from resume_v61_after_canonical import load_slim_master


def gzip_writer(path, fields):
    handle = gzip.open(path, "wt", encoding="utf-8", newline="")
    writer = csv.DictWriter(
        handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    return handle, writer


def rekey_duplicate_evidence_ids() -> int:
    path = v.CANDIDATE / "binding_evidence_master_v6_1.tsv.gz"
    temp = path.with_suffix(path.suffix + ".repair")
    audit_path = (
        v.CANDIDATE / "duplicate_evidence_id_rekey_audit_v1_0.tsv.gz"
    )
    audit_temp = audit_path.with_suffix(audit_path.suffix + ".tmp")
    if temp.exists():
        temp.unlink()
    if audit_temp.exists():
        audit_temp.unlink()
    seen: set[str] = set()
    rekeyed = 0
    audit_fields = [
        "original_evidence_id",
        "replacement_evidence_id",
        "target_uniprot_id",
        "compound_internal_id",
        "source_database",
        "source_record_id",
        "repair_action",
    ]
    audit_handle, audit_writer = gzip_writer(audit_temp, audit_fields)
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or [])
        output_handle, writer = gzip_writer(temp, fields)
        try:
            for row in reader:
                evidence_id = row["evidence_id"]
                if evidence_id in seen:
                    rekeyed += 1
                    payload = "\x1f".join(
                        [
                            evidence_id,
                            row.get("target_uniprot_id", ""),
                            row.get("compound_internal_id", ""),
                            row.get("source_database", ""),
                            row.get("source_record_id", ""),
                            str(rekeyed),
                        ]
                    )
                    new_id = "BE61D-" + hashlib.sha256(
                        payload.encode("utf-8")
                    ).hexdigest()[:24].upper()
                    while new_id in seen:
                        payload += "\x1f"
                        new_id = "BE61D-" + hashlib.sha256(
                            payload.encode("utf-8")
                        ).hexdigest()[:24].upper()
                    row["evidence_id"] = new_id
                    row["duplicate_status_v60"] = (
                        "legacy_duplicate_evidence_id_rekeyed_v61"
                    )
                    row["evidence_action_v61"] = (
                        "baseline_preserved_duplicate_id_rekeyed"
                    )
                    audit_writer.writerow(
                        {
                            "original_evidence_id": evidence_id,
                            "replacement_evidence_id": new_id,
                            "target_uniprot_id": row.get(
                                "target_uniprot_id", ""
                            ),
                            "compound_internal_id": row.get(
                                "compound_internal_id", ""
                            ),
                            "source_database": row.get(
                                "source_database", ""
                            ),
                            "source_record_id": row.get(
                                "source_record_id", ""
                            ),
                            "repair_action": (
                                "deterministic_rekey_preserve_evidence"
                            ),
                        }
                    )
                seen.add(row["evidence_id"])
                writer.writerow(row)
        finally:
            output_handle.close()
            audit_handle.close()
    temp.replace(path)
    audit_temp.replace(audit_path)
    return rekeyed


def split_unmapped_negative(
    valid_compounds: set[str],
) -> tuple[int, int]:
    path = v.CANDIDATE / "negative_binding_evidence_v1_1.tsv.gz"
    mapped_temp = path.with_suffix(path.suffix + ".mapped")
    unmapped_path = (
        v.CANDIDATE / "negative_binding_evidence_unmapped_review_v1_1.tsv.gz"
    )
    unmapped_temp = unmapped_path.with_suffix(unmapped_path.suffix + ".tmp")
    for candidate in (mapped_temp, unmapped_temp):
        if candidate.exists():
            candidate.unlink()
    mapped = 0
    unmapped = 0
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or [])
        review_fields = list(fields)
        if "unmapped_reason_v61" not in review_fields:
            review_fields.append("unmapped_reason_v61")
        mapped_handle, mapped_writer = gzip_writer(mapped_temp, fields)
        review_handle, review_writer = gzip_writer(
            unmapped_temp, review_fields
        )
        try:
            for row in reader:
                compound_id = row.get("compound_internal_id", "")
                if compound_id and compound_id in valid_compounds:
                    mapped_writer.writerow(row)
                    mapped += 1
                else:
                    row["unmapped_reason_v61"] = (
                        "missing_compound_internal_id"
                        if not compound_id
                        else "compound_internal_id_not_in_v6_1_master"
                    )
                    row["release_action_v61"] = (
                        "retained_unmapped_negative_review"
                    )
                    review_writer.writerow(row)
                    unmapped += 1
                if (mapped + unmapped) % 1000000 == 0:
                    v.update_progress(
                        "split_unmapped_negative_evidence",
                        rows_processed=mapped + unmapped,
                        mapped_rows=mapped,
                        unmapped_review_rows=unmapped,
                    )
        finally:
            mapped_handle.close()
            review_handle.close()
    mapped_temp.replace(path)
    unmapped_temp.replace(unmapped_path)
    return mapped, unmapped


def main() -> None:
    report_path = v.CANDIDATE / "V61_VALIDATION_REPORT.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    all_master = load_slim_master()
    valid_compounds = set(all_master)

    v.update_progress("rekey_duplicate_evidence_ids")
    rekeyed = rekey_duplicate_evidence_ids()
    v.update_progress(
        "split_unmapped_negative_evidence",
        duplicate_evidence_ids_rekeyed=rekeyed,
    )
    mapped_negative, unmapped_negative = split_unmapped_negative(
        valid_compounds
    )

    evidence_stats = report["evidence_stats"]
    evidence_stats["counts"]["baseline_negative"] = mapped_negative
    evidence_stats["counts"][
        "baseline_negative_unmapped_review"
    ] = unmapped_negative
    evidence_stats["counts"][
        "duplicate_evidence_ids_rekeyed"
    ] = rekeyed

    con = sqlite3.connect(v.DB)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    v.update_progress("validate_referential_integrity_after_repairs")
    validation = v.duplicate_and_fk_checks(
        con, all_master, report["reconciliation"]
    )
    validation.update(
        {
            key: value
            for key, value in report.items()
            if key
            not in {
                "blocking_errors",
                "blocking_error_total",
                "status",
                "evidence_stats",
                "generated_utc",
            }
        }
    )
    validation.update(
        {
            "status": (
                "PASS" if validation["blocking_error_total"] == 0 else "FAIL"
            ),
            "generated_utc": v.now(),
            "evidence_stats": evidence_stats,
            "duplicate_evidence_ids_rekeyed": rekeyed,
            "mapped_negative_evidence_rows": mapped_negative,
            "unmapped_negative_review_rows": unmapped_negative,
        }
    )
    v.write_json(report_path, validation)
    v.write_release_docs(
        validation["identity_counts"],
        validation["reconciliation"],
        validation["canonicalization"],
        validation["master_stats"],
        evidence_stats,
        validation["output_stats"],
        validation,
    )
    v.manifest_and_hashes()
    if validation["blocking_error_total"]:
        v.update_progress(
            "validation_failed",
            status="failed",
            validation_report=str(report_path),
            blocking_errors=validation["blocking_errors"],
        )
        raise RuntimeError(
            f"V6.1 not frozen; blocking errors: {validation['blocking_errors']}"
        )

    shutil.copytree(v.CANDIDATE, v.FINAL)
    freeze = {
        "status": "frozen",
        "release_version": v.RELEASE_VERSION,
        "release_path": str(v.FINAL),
        "frozen_utc": v.now(),
        "validation_report_sha256": v.sha256_file(
            v.FINAL / "V61_VALIDATION_REPORT.json"
        ),
        "release_manifest_sha256": v.sha256_file(
            v.FINAL / "RELEASE_MANIFEST_v6_1.tsv"
        ),
        "resumed_from_checkpoint": True,
    }
    v.write_json(v.FINAL / "FREEZE_v6_1.json", freeze)
    v.update_progress(
        "complete",
        status="complete",
        release_path=str(v.FINAL),
        validation_status="PASS",
        freeze_manifest=str(v.FINAL / "FREEZE_v6_1.json"),
    )
    con.close()
    print(json.dumps(freeze, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
