#!/usr/bin/env python3
"""Validate V6.1 compound/form/evidence preview referential integrity."""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
RUN = ROOT / "runs" / "incremental_v61_20260727"
DELTA = RUN / "merge_candidate" / "release_delta_preview"
REPORT = RUN / "qa" / "V61_RELEASE_DELTA_PREVIEW_REPORT.json"
OUT = RUN / "qa" / "V61_RELEASE_DELTA_PREVIEW_VALIDATION.json"


def write_json(payload: dict) -> None:
    temporary = OUT.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    temporary.replace(OUT)


def main() -> None:
    preview = json.loads(REPORT.read_text(encoding="utf-8"))
    counts: Counter[str] = Counter()
    blockers: list[str] = []

    baseline_compounds: set[str] = set()
    with (
        RELEASE / "small_molecule_master_v1_1.tsv"
    ).open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            baseline_compounds.add(row["compound_internal_id"])

    proteins: set[str] = set()
    with (
        RELEASE / "human_membrane_protein_master_v6_0.tsv"
    ).open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            proteins.add(row["target_uniprot_id"])

    preview_compounds: set[str] = set()
    preview_keys: set[str] = set()
    with gzip.open(
        DELTA / "small_molecule_master_delta_preview_v61.tsv.gz",
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["compound_delta_rows"] += 1
            compound_id = row["provisional_compound_internal_id"]
            key = row["standard_inchikey"]
            if not compound_id or compound_id in preview_compounds:
                counts["duplicate_or_blank_preview_compound_id"] += 1
            if not key or key in preview_keys:
                counts["duplicate_or_blank_preview_inchikey"] += 1
            if compound_id in baseline_compounds:
                counts["preview_compound_id_collides_baseline"] += 1
            preview_compounds.add(compound_id)
            preview_keys.add(key)

    form_ids: set[str] = set()
    with gzip.open(
        DELTA / "compound_form_hierarchy_delta_preview_v61.tsv.gz",
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["form_delta_rows"] += 1
            form_id = row["provisional_form_id"]
            parent = row["provisional_compound_internal_id"]
            if not form_id or form_id in form_ids:
                counts["duplicate_or_blank_preview_form_id"] += 1
            if parent not in preview_compounds:
                counts["form_missing_preview_parent"] += 1
            form_ids.add(form_id)

    evidence_ids: set[str] = set()
    action_counts: Counter[str] = Counter()
    allowed_compounds = baseline_compounds | preview_compounds
    with gzip.open(
        DELTA / "binding_evidence_resolution_delta_preview_v61.tsv.gz",
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["evidence_delta_rows"] += 1
            evidence_id = row["source_evidence_id"]
            compound = row["resolved_compound_internal_id_preview"]
            action = row["proposed_release_action"]
            default = row["proposed_default_inclusion"]
            if not evidence_id or evidence_id in evidence_ids:
                counts["duplicate_or_blank_preview_evidence_id"] += 1
            evidence_ids.add(evidence_id)
            if row["target_uniprot_id"] not in proteins:
                counts["evidence_missing_protein"] += 1
            if compound and compound not in allowed_compounds:
                counts["evidence_missing_compound"] += 1
            if default == "1":
                counts["preview_default_rows"] += 1
                if action != "positive_release_candidate" or not compound:
                    counts["invalid_preview_default_action"] += 1
            action_counts[action] += 1

    if counts["compound_delta_rows"] != preview["candidate_compound_rows"]:
        counts["compound_report_difference"] = (
            counts["compound_delta_rows"] - preview["candidate_compound_rows"]
        )
    if counts["evidence_delta_rows"] != preview[
        "evidence_resolution_preview_rows"
    ]:
        counts["evidence_report_difference"] = (
            counts["evidence_delta_rows"]
            - preview["evidence_resolution_preview_rows"]
        )
    if dict(action_counts) != preview["evidence_action_counts"]:
        counts["evidence_action_count_mismatch"] = 1

    blocking_fields = [
        "duplicate_or_blank_preview_compound_id",
        "duplicate_or_blank_preview_inchikey",
        "preview_compound_id_collides_baseline",
        "duplicate_or_blank_preview_form_id",
        "form_missing_preview_parent",
        "duplicate_or_blank_preview_evidence_id",
        "evidence_missing_protein",
        "evidence_missing_compound",
        "invalid_preview_default_action",
        "compound_report_difference",
        "evidence_report_difference",
        "evidence_action_count_mismatch",
    ]
    for field in blocking_fields:
        if counts[field]:
            blockers.append(f"{field}: {counts[field]}")

    payload = {
        "status": "passed_preview" if not blockers else "failed_preview",
        "scope": "current completed sources only; provisional identifiers",
        "counts": dict(counts),
        "action_counts": dict(action_counts),
        "blockers": blockers,
        "release_ready": False,
        "waiting_for": ["PubChem P2", "PubChem P3", "BRENDA"],
    }
    write_json(payload)
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if blockers:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
