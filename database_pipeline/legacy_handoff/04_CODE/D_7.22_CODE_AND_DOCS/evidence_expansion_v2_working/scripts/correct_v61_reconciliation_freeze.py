#!/usr/bin/env python3
"""Correct review-row reconciliation, revalidate, and replace pre-handoff freeze."""

from __future__ import annotations

import json
import shutil
import sqlite3

import finalize_v61_release as v
from resume_v61_after_canonical import load_slim_master


SUPERSEDED = (
    v.ROOT
    / "releases"
    / "superseded_pre_handoff_release_mempro_v6_1_20260729"
)


def main() -> None:
    report_path = v.CANDIDATE / "V61_VALIDATION_REPORT.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    con = sqlite3.connect(v.DB)
    all_master = load_slim_master()

    unique_dispositions = sum(
        report["reconciliation"]["disposition_counts"].values()
    )
    review_total = v.count_rows(v.REVIEW0)
    duplicate_source_ids = review_total - unique_dispositions
    if duplicate_source_ids < 0:
        raise RuntimeError("Disposition count exceeds V6.0 review row count")
    report["reconciliation"].update(
        {
            "review_rows_total": review_total,
            "disposition_unique_source_evidence_ids": unique_dispositions,
            "duplicate_source_evidence_id_rows": duplicate_source_ids,
        }
    )

    v.update_progress(
        "correct_review_reconciliation_and_revalidate",
        review_rows_total=review_total,
        disposition_unique_source_evidence_ids=unique_dispositions,
        duplicate_source_evidence_id_rows=duplicate_source_ids,
    )
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
                "generated_utc",
                "review_rows_accounted",
                "review_rows_total",
            }
        }
    )
    validation.update(
        {
            "status": (
                "PASS" if validation["blocking_error_total"] == 0 else "FAIL"
            ),
            "generated_utc": v.now(),
            "pre_handoff_reconciliation_correction": True,
        }
    )
    v.write_json(report_path, validation)
    v.write_release_docs(
        validation["identity_counts"],
        validation["reconciliation"],
        validation["canonicalization"],
        validation["master_stats"],
        validation["evidence_stats"],
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

    if SUPERSEDED.exists():
        raise FileExistsError(f"Superseded audit path already exists: {SUPERSEDED}")
    if v.FINAL.exists():
        shutil.move(str(v.FINAL), str(SUPERSEDED))
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
        "pre_handoff_reconciliation_correction": True,
        "superseded_pre_handoff_snapshot": str(SUPERSEDED),
    }
    v.write_json(v.FINAL / "FREEZE_v6_1.json", freeze)
    v.update_progress(
        "complete",
        status="complete",
        release_path=str(v.FINAL),
        validation_status="PASS",
        freeze_manifest=str(v.FINAL / "FREEZE_v6_1.json"),
        duplicate_source_evidence_id_rows=duplicate_source_ids,
    )
    con.close()
    print(json.dumps(freeze, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
