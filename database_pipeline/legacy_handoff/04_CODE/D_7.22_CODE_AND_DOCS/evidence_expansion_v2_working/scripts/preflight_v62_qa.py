#!/usr/bin/env python3
"""Stratified preflight audit required before the V6.2 freeze step."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import random
from collections import Counter
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RUN = ROOT / "runs" / "v62_completion_20260729"
QA = RUN / "qa"
EXPR = RUN / "staging" / "expression_location_v2"
SUBUNIT = RUN / "staging" / "subunit_assembly_v2"
NEG = RUN / "staging" / "negative_identity"
GATE = QA / "V62_PREFLIGHT_QA_GATE.json"
SAMPLE = QA / "V62_STRATIFIED_VALIDATION_SAMPLE.tsv"


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def reservoir(path: Path, size: int, seed: int) -> list[dict]:
    rng = random.Random(seed)
    sample = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for index, row in enumerate(
            csv.DictReader(handle, delimiter="\t"), start=1
        ):
            if len(sample) < size:
                sample.append(row)
            else:
                position = rng.randint(1, index)
                if position <= size:
                    sample[position - 1] = row
    return sample


def main() -> None:
    expression_report = json.loads(
        (QA / "EXPRESSION_LOCATION_V2_VALIDATION.json").read_text(
            encoding="utf-8"
        )
    )
    subunit_report = json.loads(
        (QA / "SUBUNIT_ASSEMBLY_V2_VALIDATION.json").read_text(
            encoding="utf-8"
        )
    )
    negative_report = json.loads(
        (QA / "NEGATIVE_EVIDENCE_RESOLUTION_V62_VALIDATION.json").read_text(
            encoding="utf-8"
        )
    )
    blocking = Counter()
    if expression_report.get("summary_rows") != 10_997:
        blocking["expression_summary_row_count"] += 1
    if subunit_report.get("outputs", {}).get("summary_rows") != 10_997:
        blocking["subunit_summary_row_count"] += 1
    negative_counts = negative_report.get("negative_table_counts", {})
    source_rows = negative_counts.get("source_rows", 0)
    accounted = sum(
        negative_counts.get(key, 0)
        for key in [
            "mapped_release_rows",
            "unmapped_review_rows",
            "duplicates_removed",
        ]
    )
    if source_rows != 10_246_948 or accounted != source_rows:
        blocking["negative_row_reconciliation"] += abs(source_rows - accounted) or 1

    unsafe_heteromer_target_copy = 0
    direct_stoichiometry_overclaim = 0
    with gzip.open(
        SUBUNIT / "protein_biological_assembly_v2.tsv.gz",
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                int(row.get("distinct_protein_entity_count") or 0) > 1
                and row.get("target_copy_count")
            ):
                unsafe_heteromer_target_copy += 1
            if (
                row.get("assembly_direct_experimental_confirmation_flag")
                == "1"
                and row.get("target_state_inference_status")
                != "direct_experimental_stoichiometry"
            ):
                direct_stoichiometry_overclaim += 1
    if unsafe_heteromer_target_copy:
        blocking["heteromer_target_copy_overclaim"] += unsafe_heteromer_target_copy
    if direct_stoichiometry_overclaim:
        blocking["direct_stoichiometry_overclaim"] += direct_stoichiometry_overclaim

    unsafe_identity = 0
    crosswalk_rows = 0
    with gzip.open(
        NEG / "negative_cid_identity_crosswalk_v62.tsv.gz",
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            crosswalk_rows += 1
            resolution = row["identity_resolution"]
            mapped = bool(row["compound_internal_id"])
            safe = resolution in {
                "exact_full_inchikey_unique_form",
                "canonical_parent_full_inchikey_unique",
                "new_unique_canonical_parent",
            }
            if mapped != safe:
                unsafe_identity += 1
            if (
                resolution == "new_unique_canonical_parent"
                and row["exact_inchikey"]
                != row["canonical_parent_inchikey"]
            ):
                unsafe_identity += 1
    if crosswalk_rows != 1_622_591:
        blocking["cid_crosswalk_row_count"] += abs(
            crosswalk_rows - 1_622_591
        )
    if unsafe_identity:
        blocking["unsafe_identity_resolution"] += unsafe_identity

    expression_sample = reservoir(
        EXPR / "expression_location_summary_v2.tsv.gz", 200, 6201
    )
    subunit_sample = reservoir(
        SUBUNIT / "protein_subunit_summary_v2.tsv.gz", 150, 6202
    )
    negative_sample = reservoir(
        NEG / "negative_cid_identity_crosswalk_v62.tsv.gz", 150, 6203
    )
    sample_fields = [
        "sample_id",
        "module",
        "primary_id",
        "secondary_id",
        "status_or_grade",
        "review_flags",
        "automated_audit_result",
        "audit_rule",
    ]
    with SAMPLE.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=sample_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        index = 0
        for row in expression_sample:
            index += 1
            writer.writerow(
                {
                    "sample_id": f"V62-QA-{index:04d}",
                    "module": "expression_location",
                    "primary_id": row["target_uniprot_id"],
                    "secondary_id": row["hpa_ensembl_gene_ids_v62"],
                    "status_or_grade": row["hpa_mapping_status_v62"],
                    "review_flags": row[
                        "expression_location_conflict_reason_v62"
                    ],
                    "automated_audit_result": (
                        "pass"
                        if row["target_uniprot_id"]
                        else "fail"
                    ),
                    "audit_rule": "protein_key_present_and_summary_typed",
                }
            )
        for row in subunit_sample:
            index += 1
            writer.writerow(
                {
                    "sample_id": f"V62-QA-{index:04d}",
                    "module": "subunit_assembly",
                    "primary_id": row["target_uniprot_id"],
                    "secondary_id": row["biological_assembly_ids_v62"],
                    "status_or_grade": row["subunit_evidence_grade_v62"],
                    "review_flags": row["subunit_conflict_reason_v62"],
                    "automated_audit_result": (
                        "pass"
                        if row["subunit_evidence_grade_v62"].startswith("SA")
                        else "fail"
                    ),
                    "audit_rule": "evidence_grade_and_conflict_status_explicit",
                }
            )
        for row in negative_sample:
            index += 1
            resolution = row["identity_resolution"]
            mapped = bool(row["compound_internal_id"])
            safe = resolution in {
                "exact_full_inchikey_unique_form",
                "canonical_parent_full_inchikey_unique",
                "new_unique_canonical_parent",
            }
            writer.writerow(
                {
                    "sample_id": f"V62-QA-{index:04d}",
                    "module": "negative_identity",
                    "primary_id": row["cid"],
                    "secondary_id": row["compound_internal_id"],
                    "status_or_grade": resolution,
                    "review_flags": row["review_reason"],
                    "automated_audit_result": (
                        "pass" if mapped == safe else "fail"
                    ),
                    "audit_rule": "only_safe_full_structure_resolution_is_mapped",
                }
            )

    gate = {
        "status": "PASS" if not blocking else "FAIL",
        "completed_utc": now(),
        "blocking_errors": dict(blocking),
        "blocking_error_total": sum(blocking.values()),
        "checks": {
            "negative_source_rows": source_rows,
            "negative_accounted_rows": accounted,
            "cid_crosswalk_rows": crosswalk_rows,
            "heteromer_target_copy_overclaims": unsafe_heteromer_target_copy,
            "direct_stoichiometry_overclaims": direct_stoichiometry_overclaim,
            "unsafe_identity_resolutions": unsafe_identity,
            "stratified_sample_rows": 500,
        },
        "sample_file": str(SAMPLE),
        "review_scope": "agent-assisted stratified rule audit; unresolved scientific conflicts remain explicitly flagged",
    }
    GATE.write_text(
        json.dumps(gate, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(gate, ensure_ascii=False, indent=2))
    if gate["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
