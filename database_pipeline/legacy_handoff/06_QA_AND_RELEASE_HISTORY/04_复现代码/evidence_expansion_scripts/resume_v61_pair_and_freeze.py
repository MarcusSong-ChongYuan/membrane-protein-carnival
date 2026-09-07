#!/usr/bin/env python3
"""Resume V6.1 from repaired evidence through pair outputs and freeze."""

from __future__ import annotations

import json
import shutil
import sqlite3

import finalize_v61_release as v
from resume_v61_after_canonical import (
    load_slim_master,
    reconstruct_audit_counts,
)


def main() -> None:
    v.require_inputs()
    checkpoint = json.loads(v.PROGRESS.read_text(encoding="utf-8"))
    evidence_stats = checkpoint.get("evidence_stats")
    demoted_count = int(checkpoint.get("legacy_rows_demoted", 0))
    if not evidence_stats or demoted_count <= 0:
        raise RuntimeError("Repaired evidence checkpoint metadata is unavailable")

    con = sqlite3.connect(v.DB)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA temp_store=FILE")
    master_by_id, _, _, max_compound, max_form = v.load_baseline_identity()
    all_master = load_slim_master()
    identity_counts, reconciliation, canonical = reconstruct_audit_counts(
        con, max_compound, max_form
    )
    master_stats = {
        "master_rows": len(all_master),
        "new_master_rows": len(all_master) - len(master_by_id),
        "new_form_rows": canonical["new_exact_forms"],
    }

    v.update_progress(
        "clean_invalid_pair_foreign_keys",
        legacy_rows_demoted=demoted_count,
    )
    con.executescript(
        """
        DROP TABLE IF EXISTS valid_compound_v61;
        CREATE TABLE valid_compound_v61 (
          compound_internal_id TEXT PRIMARY KEY
        );
        """
    )
    con.executemany(
        "INSERT INTO valid_compound_v61 VALUES (?)",
        ((compound_id,) for compound_id in all_master),
    )
    invalid_pair_rows = con.execute(
        """
        SELECT COUNT(*) FROM pair_raw
        WHERE compound NOT IN (
          SELECT compound_internal_id FROM valid_compound_v61
        )
        """
    ).fetchone()[0]
    invalid_quant_rows = con.execute(
        """
        SELECT COUNT(*) FROM quant_raw
        WHERE compound NOT IN (
          SELECT compound_internal_id FROM valid_compound_v61
        )
        """
    ).fetchone()[0]
    con.executescript(
        """
        DELETE FROM pair_raw
        WHERE compound NOT IN (
          SELECT compound_internal_id FROM valid_compound_v61
        );
        DELETE FROM quant_raw
        WHERE compound NOT IN (
          SELECT compound_internal_id FROM valid_compound_v61
        );
        DROP INDEX IF EXISTS idx_pair_raw_key;
        DROP INDEX IF EXISTS idx_quant_raw_key;
        DROP INDEX IF EXISTS idx_pair_agg_key;
        DROP TABLE IF EXISTS pair_agg;
        DROP TABLE IF EXISTS target_agg;
        """
    )
    con.commit()
    for path in v.CANDIDATE.glob("protein_compound_summary_v6_1.tsv.gz*"):
        path.unlink()

    v.update_progress(
        "build_pair_site_protein_outputs",
        invalid_pair_rows_removed=invalid_pair_rows,
        invalid_quant_rows_removed=invalid_quant_rows,
    )
    output_stats = v.build_pair_site_protein_outputs(con, all_master)
    v.update_progress("validate_referential_integrity")
    validation = v.duplicate_and_fk_checks(con, all_master, reconciliation)
    validation.update(
        {
            "status": (
                "PASS" if validation["blocking_error_total"] == 0 else "FAIL"
            ),
            "release_version": v.RELEASE_VERSION,
            "release_date": v.RELEASE_DATE,
            "generated_utc": v.now(),
            "identity_counts": identity_counts,
            "reconciliation": reconciliation,
            "canonicalization": canonical,
            "master_stats": master_stats,
            "evidence_stats": evidence_stats,
            "output_stats": output_stats,
            "legacy_default_evidence_demoted_for_missing_compound_fk": (
                demoted_count
            ),
            "invalid_pair_rows_removed": invalid_pair_rows,
            "invalid_quant_rows_removed": invalid_quant_rows,
            "resumed_from_checkpoint": True,
        }
    )
    v.write_json(v.CANDIDATE / "V61_VALIDATION_REPORT.json", validation)
    v.write_release_docs(
        identity_counts,
        reconciliation,
        canonical,
        master_stats,
        evidence_stats,
        output_stats,
        validation,
    )
    v.manifest_and_hashes()
    if validation["blocking_error_total"]:
        v.update_progress(
            "validation_failed",
            status="failed",
            validation_report=str(
                v.CANDIDATE / "V61_VALIDATION_REPORT.json"
            ),
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
