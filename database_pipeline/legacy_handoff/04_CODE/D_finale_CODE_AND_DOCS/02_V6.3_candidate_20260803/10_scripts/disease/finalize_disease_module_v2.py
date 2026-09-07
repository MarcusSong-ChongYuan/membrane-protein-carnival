#!/usr/bin/env python3
"""Add release-disposition fields and perform stricter disease-module QA."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


SAFE = {"direct_mondo_id", "unique_official_exact_mapping"}


def read(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write(path: Path, fields: list[str], rows) -> int:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        count = 0
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--qa-dir", type=Path, required=True)
    args = p.parse_args()
    source_evidence = list(read(args.output_dir / "protein_disease_source_evidence_v2.tsv"))
    source_status_by_relation = defaultdict(list)
    for row in source_evidence:
        source_status_by_relation[row["disease_relation_id_v2"]].append(row["canonicalization_status"])
    original = list(read(args.output_dir / "protein_disease_relation_v2.tsv"))
    final_rows = []
    for row in original:
        statuses = source_status_by_relation[row["disease_relation_id_v2"]]
        safe_count = sum(status in SAFE for status in statuses)
        review_count = len(statuses) - safe_count
        row = dict(row)
        row.update({
            "canonical_mapping_statuses": ";".join(sorted(set(statuses))),
            "safe_source_evidence_count": safe_count,
            "review_source_evidence_count": review_count,
            "mapping_review_flag": int(review_count > 0),
            "default_release_inclusion_v63": int(row.get("default_relation_inclusion") == "1" and safe_count > 0),
            "release_disposition_v63": "default" if row.get("default_relation_inclusion") == "1" and safe_count > 0 else "review_or_nondefault",
        })
        final_rows.append(row)
    fields = list(final_rows[0].keys())
    final_path = args.output_dir / "protein_disease_relation_v2_1.tsv"
    write(final_path, fields, final_rows)
    keys = [(row["target_uniprot_id"], row["canonical_disease_id"]) for row in final_rows]
    checks = {
        "canonical_relation_keys_unique": len(keys) == len(set(keys)),
        "source_evidence_rows_preserved": len(source_evidence) == 10264,
        "every_relation_has_source_evidence": all(row["disease_relation_id_v2"] in source_status_by_relation for row in final_rows),
        "unsafe_only_relations_excluded_from_default": all(not (int(row["safe_source_evidence_count"]) == 0 and int(row["default_release_inclusion_v63"]) == 1) for row in final_rows),
        "review_evidence_flagged": all(not (int(row["review_source_evidence_count"]) > 0 and int(row["mapping_review_flag"]) == 0) for row in final_rows),
    }
    validation = {
        "module": "disease_identity_hierarchy_classification_anatomy",
        "module_version": "0.2.0",
        "candidate_release": "V6.3",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": {
            "canonical_relations": len(final_rows),
            "default_relations": sum(int(row["default_release_inclusion_v63"]) for row in final_rows),
            "relations_with_mapping_review_evidence": sum(int(row["mapping_review_flag"]) for row in final_rows),
            "release_disposition": dict(Counter(row["release_disposition_v63"] for row in final_rows)),
        },
        "authoritative_relation_table": final_path.name,
        "supersedes_relation_table": "protein_disease_relation_v2.tsv",
    }
    validation_path = args.qa_dir / "DISEASE_MODULE_V0_2_VALIDATION.json"
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest = {
        "generated_at_utc": validation["generated_at_utc"],
        "files": [
            {"relative_path": final_path.name, "bytes": final_path.stat().st_size, "sha256": sha(final_path)},
            {"relative_path": validation_path.name, "bytes": validation_path.stat().st_size, "sha256": sha(validation_path)},
        ],
    }
    (args.qa_dir / "DISEASE_MODULE_V0_2_MANIFEST.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
