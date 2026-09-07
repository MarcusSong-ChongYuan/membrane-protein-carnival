#!/usr/bin/env python3
"""Resume V6.1 after evidence rebuild, repairing legacy blank compound FKs."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
import sqlite3
from collections import Counter

import finalize_v61_release as v
from resume_v61_after_canonical import (
    load_slim_master,
    reconstruct_audit_counts,
)


def sanitize_legacy_blank_compounds(
    all_master: dict[str, dict[str, str]],
) -> int:
    evidence_path = v.CANDIDATE / "binding_evidence_master_v6_1.tsv.gz"
    review_path = v.CANDIDATE / "binding_evidence_review_queue_v6_1.tsv.gz"
    demoted: list[dict[str, str]] = []
    stale_atomic_temp = evidence_path.with_suffix(evidence_path.suffix + ".tmp")
    if stale_atomic_temp.exists():
        stale_atomic_temp.unlink()
    evidence_temp = evidence_path.with_suffix(evidence_path.suffix + ".repair")
    if evidence_temp.exists():
        evidence_temp.unlink()
    with gzip.open(
        evidence_path, "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        evidence_fields = list(reader.fieldnames or [])
        with gzip.open(
            evidence_temp, "wt", encoding="utf-8", newline=""
        ) as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=evidence_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                compound_id = row.get("compound_internal_id", "")
                if compound_id and compound_id in all_master:
                    writer.writerow(row)
                    continue
                demoted.append(row)
    evidence_temp.replace(evidence_path)

    review_temp = review_path.with_suffix(review_path.suffix + ".repair")
    if review_temp.exists():
        review_temp.unlink()
    with gzip.open(
        review_path, "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        review_fields = list(reader.fieldnames or [])
        with gzip.open(
            review_temp, "wt", encoding="utf-8", newline=""
        ) as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=review_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                writer.writerow(row)
            for row in demoted:
                output = {field: row.get(field, "") for field in review_fields}
                output.update(
                    {
                        "source_evidence_id": (
                            row.get("source_evidence_id_v61")
                            or row.get("evidence_id", "")
                        ),
                        "target_uniprot_id": row.get("target_uniprot_id", ""),
                        "approved_symbol": row.get("approved_symbol", ""),
                        "source_database": row.get("source_database", ""),
                        "source_record_id": row.get("source_record_id", ""),
                        "compound_source_id": row.get("compound_id", ""),
                        "compound_name": row.get("compound_name", ""),
                        "identity_resolution_v61": "baseline_missing_compound_fk",
                        "resolved_compound_internal_id_v61": "",
                        "resolved_compound_form_id_v61": "",
                        "release_disposition_v61": "review",
                        "release_reason_v61": (
                            "legacy V6.0 default evidence lacked a valid "
                            "compound_internal_id"
                        ),
                        "release_version_v61": v.RELEASE_VERSION,
                        "release_version": v.RELEASE_VERSION,
                    }
                )
                writer.writerow(output)
    review_temp.replace(review_path)
    return len(demoted)


def reconstruct_evidence_stats(
    con: sqlite3.Connection,
    demoted_count: int,
) -> dict:
    action_counts = {
        action: count
        for action, count in con.execute(
            "SELECT action,COUNT(*) FROM disposition GROUP BY action"
        )
    }
    duplicate_count = v.count_rows(
        v.CANDIDATE / "duplicate_evidence_audit_v1_0.tsv.gz"
    )
    counts = {
        "baseline_positive": v.count_rows(v.EVIDENCE0) - demoted_count,
        "baseline_positive_demoted_missing_compound_fk": demoted_count,
        "baseline_negative": v.count_rows(v.NEGATIVE0),
        "new_positive": max(
            0, action_counts.get("positive", 0) - duplicate_count
        ),
        "new_negative": action_counts.get("negative", 0),
        "retained_review": action_counts.get("review", 0) + demoted_count,
        "retained_excluded": action_counts.get("excluded", 0),
        "duplicate": duplicate_count,
        "retained_duplicate": duplicate_count,
    }
    type_to_source = {
        "PUBCHEM_CID": "PubChem BioAssay",
        "PDB_CCD": "PDBe",
        "PDBBIND_COMPLEX": "PDBbind",
        "BRENDA_NAME": "BRENDA",
    }
    source_counts = {
        f"{type_to_source.get(source_type, source_type)}:{action}": count
        for source_type, action, count in con.execute(
            """
            SELECT source_type,action,COUNT(*)
            FROM disposition GROUP BY source_type,action
            """
        )
    }
    return {
        "counts": counts,
        "source_counts": source_counts,
        "manual_validation_sample_rows": v.count_rows(
            v.CANDIDATE / "manual_validation_sample_pending_v6_1.tsv"
        ),
    }


def main() -> None:
    v.require_inputs()
    if not v.DB.exists():
        raise FileNotFoundError(f"Evidence checkpoint is missing: {v.DB}")
    con = sqlite3.connect(v.DB)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA temp_store=FILE")
    master_by_id, _, form_by_key, max_compound, max_form = (
        v.load_baseline_identity()
    )
    all_master = load_slim_master()
    identity_counts, reconciliation, canonical = reconstruct_audit_counts(
        con, max_compound, max_form
    )
    master_stats = {
        "master_rows": len(all_master),
        "new_master_rows": len(all_master) - len(master_by_id),
        "new_form_rows": canonical["new_exact_forms"],
    }

    v.update_progress("repair_legacy_blank_compound_fk")
    demoted_count = sanitize_legacy_blank_compounds(all_master)
    con.execute("DELETE FROM pair_raw WHERE compound='' OR compound IS NULL")
    con.execute("DELETE FROM quant_raw WHERE compound='' OR compound IS NULL")
    con.commit()
    for path in v.CANDIDATE.glob("protein_compound_summary_v6_1.tsv.gz*"):
        path.unlink()

    evidence_stats = reconstruct_evidence_stats(con, demoted_count)
    v.update_progress(
        "build_pair_site_protein_outputs",
        legacy_rows_demoted=demoted_count,
        evidence_stats=evidence_stats,
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
