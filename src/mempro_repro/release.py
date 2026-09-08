from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path


TABLES = {
    "proteins": "01_core_v72/01_release_tables/protein_master_v72.tsv.gz",
    "pairs": "01_core_v72/01_release_tables/protein_compound_pair_v72.tsv.gz",
    "positive_evidence": "01_core_v72/01_release_tables/positive_interaction_evidence_v72.tsv.gz",
    "compounds": "01_core_v72/01_release_tables/compound_master_v72.tsv.gz",
    "disease_relations": "01_core_v72/01_release_tables/protein_disease_relation_v72.tsv.gz",
    "binding_sites": "01_core_v72/01_release_tables/interaction_site_v72.tsv.gz",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(data_root: Path) -> list[dict[str, str]]:
    manifest = data_root / "05_metadata" / "FORMAL_MANIFEST_SHA256.tsv"
    with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify_manifest(data_root: Path, sample_size: int | None = None) -> dict:
    rows = read_manifest(data_root)
    if sample_size is not None and sample_size < len(rows):
        # Deterministic coverage: beginning, middle and end of sorted manifest.
        indices = sorted({round(i * (len(rows) - 1) / (sample_size - 1)) for i in range(sample_size)}) if sample_size > 1 else [0]
        rows = [rows[index] for index in indices]

    failures: list[dict[str, str]] = []
    for row in rows:
        path = data_root / row["relative_path"]
        if not path.exists():
            failures.append({"relative_path": row["relative_path"], "reason": "missing"})
        elif path.stat().st_size != int(row["bytes"]):
            failures.append({"relative_path": row["relative_path"], "reason": "byte_count_mismatch"})
        elif sha256(path).lower() != row["sha256"].lower():
            failures.append({"relative_path": row["relative_path"], "reason": "sha256_mismatch"})
    return {"checked_files": len(rows), "failures": failures, "status": "PASS" if not failures else "FAIL"}


def tsv_row_count(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        # Some biological names contain embedded newlines.  Count parsed TSV
        # records, never physical text lines, so the result has the same unit
        # as a release-table row and the FORMAL QA report.
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def build_release_summary(data_root: Path) -> dict:
    counts = {label: tsv_row_count(data_root / relative_path) for label, relative_path in TABLES.items()}
    qa_path = data_root / "07_QA" / "FORMAL_VALIDATION_REPORT.json"
    qa = json.loads(qa_path.read_text(encoding="utf-8"))
    return {
        "release": "MemPro FORMAL",
        "counts_from_frozen_tables": counts,
        "release_qa_status": qa.get("status"),
        "release_qa_fixed_counts": qa.get("fixed_counts", {}),
        "interpretation": {
            "alphaFold": "Predicted models are not experimental structures.",
            "negative_evidence": "Archived separately; retained for audit and conflict analysis.",
            "evidence_lineage": "Contributing database count is not an independent-experiment count.",
        },
    }
