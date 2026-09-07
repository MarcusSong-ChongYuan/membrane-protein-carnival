#!/usr/bin/env python3
"""Validate the current-input V6.1 premerge and the immutable V6.0 baseline."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import time
from collections import Counter
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
RUN = ROOT / "runs" / "incremental_v61_20260727"
PRE = RUN / "merge_candidate" / "premerge_existing_inputs"
BUILD_REPORT = RUN / "qa" / "V61_PREMERGE_EXISTING_INPUTS_REPORT.json"
OUT = RUN / "qa" / "V61_PREMERGE_VALIDATION_REPORT.json"
FROZEN = RUN / "manifests" / "V60_FROZEN_INPUT_MANIFEST.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def main() -> None:
    deadline = time.time() + 3 * 60 * 60
    while time.time() < deadline:
        if BUILD_REPORT.exists():
            try:
                build = json.loads(BUILD_REPORT.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                build = {}
            if build.get("status") == "complete_premerge_not_release":
                break
        time.sleep(15)
    else:
        raise TimeoutError("Premerge report was not completed within three hours")

    blockers: list[str] = []
    counts: Counter[str] = Counter()
    candidate_keys: set[str] = set()
    provisional_ids: set[str] = set()
    candidate_path = PRE / "global_identity_clusters_current_inputs_v61.tsv.gz"
    with gzip.open(
        candidate_path, "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["candidate_rows"] += 1
            key = row["inchikey"]
            provisional = row["provisional_cluster_id"]
            if not key or key in candidate_keys:
                counts["duplicate_or_blank_candidate_key"] += 1
            candidate_keys.add(key)
            if row["candidate_action"] == "new_candidate_pending_full_refresh":
                if not provisional or provisional in provisional_ids:
                    counts["duplicate_or_blank_provisional_id"] += 1
                provisional_ids.add(provisional)
            if row["structure_parse_status"] != "ok":
                counts["candidate_structure_manual_review"] += 1
            if row["full_key_matches_source"] != "1":
                counts["candidate_computed_key_mismatch"] += 1

    form_keys: set[str] = set()
    form_path = PRE / "parent_form_stereo_candidates_current_inputs_v61.tsv.gz"
    with gzip.open(form_path, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["parent_form_rows"] += 1
            key = row["exact_inchikey"]
            if not key or key in form_keys:
                counts["duplicate_or_blank_parent_form_key"] += 1
            form_keys.add(key)
            if not row["required_release_action"]:
                counts["blank_parent_form_action"] += 1

    if candidate_keys != form_keys:
        counts["candidate_parent_key_symmetric_difference"] = len(
            candidate_keys.symmetric_difference(form_keys)
        )

    lookup_keys: set[tuple[str, str]] = set()
    lookup_path = PRE / "source_identity_resolution_lookup_current_inputs_v61.tsv.gz"
    with gzip.open(lookup_path, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["source_lookup_rows"] += 1
            key = (row["source_object_type"], row["source_object_id"])
            if not key[0] or not key[1] or key in lookup_keys:
                counts["duplicate_or_blank_source_lookup_key"] += 1
            lookup_keys.add(key)
            if (
                row["provisional_cluster_id"]
                and row["provisional_cluster_id"] not in provisional_ids
            ):
                counts["lookup_missing_provisional_cluster"] += 1

    review_total = sum(build["review_queue_premerge_actions"].values())
    counts["review_actions_total"] = review_total
    if review_total != build["review_queue_rows_scanned"]:
        counts["review_action_reconciliation_difference"] = (
            review_total - build["review_queue_rows_scanned"]
        )

    frozen = json.loads(FROZEN.read_text(encoding="utf-8"))
    for item in frozen["files"]:
        path = RELEASE / item["file"]
        counts["v60_files_checked"] += 1
        if (
            not path.exists()
            or path.stat().st_size != item["size_bytes"]
            or sha256(path) != item["sha256"]
        ):
            counts["v60_hash_mismatches"] += 1

    blocking_fields = [
        "duplicate_or_blank_candidate_key",
        "duplicate_or_blank_provisional_id",
        "duplicate_or_blank_parent_form_key",
        "blank_parent_form_action",
        "candidate_parent_key_symmetric_difference",
        "duplicate_or_blank_source_lookup_key",
        "lookup_missing_provisional_cluster",
        "review_action_reconciliation_difference",
        "v60_hash_mismatches",
    ]
    for field in blocking_fields:
        if counts[field]:
            blockers.append(f"{field}: {counts[field]}")

    report = {
        "status": "passed_premerge" if not blockers else "failed_premerge",
        "scope": "current completed inputs only; not a V6.1 release",
        "counts": dict(counts),
        "blockers": blockers,
        "nonblocking_manual_review_counts": {
            "candidate_structure_manual_review": counts[
                "candidate_structure_manual_review"
            ],
            "candidate_computed_key_mismatch": counts[
                "candidate_computed_key_mismatch"
            ],
        },
        "release_ready": False,
        "release_waiting_for": ["PubChem P2", "PubChem P3", "BRENDA"],
    }
    write_json(OUT, report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if blockers:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
