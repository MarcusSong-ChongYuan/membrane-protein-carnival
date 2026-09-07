#!/usr/bin/env python3
"""Initialize an isolated, append-only incremental expansion run.

The published V5.5/V5.4/V4.2/V1.0 inputs remain read-only.  This script only
creates a run directory and records the already-frozen baseline plus current
source/download state so every later transformation is reproducible.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import pathlib
import shutil
import sys
from collections import Counter


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN_ID = os.environ.get("MEMPRO_INCREMENTAL_RUN_ID", "incremental_20260727")
RUN_ROOT = ROOT / "runs" / RUN_ID
BASELINE_CONFIG = ROOT / "config" / "baseline_paths.json"
BASELINE_MANIFEST = ROOT / "reports" / "baseline_file_manifest.tsv"
PUBCHEM_CHECKPOINT = (
    ROOT / "raw" / "pubchem_202607" / "pubchem_protein_concise_v2.jsonl"
)


def sha256(path: pathlib.Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def read_baseline_manifest() -> list[dict[str, str]]:
    with BASELINE_MANIFEST.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pubchem_state() -> dict[str, object]:
    state: dict[str, object] = {
        "checkpoint_path": str(PUBCHEM_CHECKPOINT),
        "checkpoint_exists": PUBCHEM_CHECKPOINT.exists(),
        "processed_targets": 0,
        "status_counts": {},
        "row_count": 0,
        "quantitative_row_count": 0,
        "last_retrieved_at_utc": None,
    }
    if not PUBCHEM_CHECKPOINT.exists():
        return state
    statuses: Counter[str] = Counter()
    last: dict[str, object] | None = None
    with PUBCHEM_CHECKPOINT.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            item = json.loads(line)
            statuses[str(item.get("status") or "unknown")] += 1
            state["processed_targets"] = int(state["processed_targets"]) + 1
            state["row_count"] = int(state["row_count"]) + int(
                item.get("row_count") or 0
            )
            state["quantitative_row_count"] = int(
                state["quantitative_row_count"]
            ) + int(item.get("quantitative_row_count") or 0)
            last = item
    state["status_counts"] = dict(sorted(statuses.items()))
    if last:
        state["last_retrieved_at_utc"] = last.get("retrieved_at_utc")
        state["last_target_uniprot_id"] = last.get("target_uniprot_id")
    return state


def source_file(path: pathlib.Path, include_hash: bool = True) -> dict[str, object]:
    result: dict[str, object] = {
        "path": str(path),
        "exists": path.exists(),
    }
    if path.exists():
        stat = path.stat()
        result.update(
            {
                "size_bytes": stat.st_size,
                "mtime_utc": dt.datetime.fromtimestamp(
                    stat.st_mtime, tz=dt.timezone.utc
                ).isoformat(),
            }
        )
        if include_hash:
            result["sha256"] = sha256(path)
    return result


def atomic_json(path: pathlib.Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    temporary.replace(path)


def main() -> int:
    if not BASELINE_CONFIG.exists() or not BASELINE_MANIFEST.exists():
        raise FileNotFoundError("Frozen baseline configuration/manifest is missing")

    directories = [
        RUN_ROOT / "staging" / "brenda",
        RUN_ROOT / "staging" / "pdbe",
        RUN_ROOT / "staging" / "pdbbind",
        RUN_ROOT / "staging" / "pubchem",
        RUN_ROOT / "staging" / "gpcrdb",
        RUN_ROOT / "merge_candidate",
        RUN_ROOT / "qa",
        RUN_ROOT / "logs",
        RUN_ROOT / "manifests",
    ]
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)

    baseline_paths = json.loads(BASELINE_CONFIG.read_text(encoding="utf-8"))
    baseline_rows = read_baseline_manifest()
    baseline_validation = []
    for row in baseline_rows:
        path = pathlib.Path(row["absolute_path"])
        baseline_validation.append(
            {
                "input_id": row["input_id"],
                "path": str(path),
                "exists": path.exists(),
                "size_matches": (
                    path.exists() and path.stat().st_size == int(row["size_bytes"])
                ),
                "expected_sha256": row["sha256"],
            }
        )

    sources = {
        "brenda_2026_1": source_file(
            ROOT
            / "raw"
            / "manual_downloads"
            / "brenda_2026_1"
            / "brenda_2026_1.json.tar.gz"
        ),
        "pdbbind_index": source_file(
            ROOT
            / "raw"
            / "manual_downloads"
            / "pdbbind_2020r1"
            / "index.tar.gz"
        ),
        "pdbbind_protein_ligand": source_file(
            ROOT
            / "raw"
            / "manual_downloads"
            / "pdbbind_2020r1"
            / "P-L.tar.gz",
            include_hash=False,
        ),
        "pdbe_fetch_report": source_file(
            ROOT / "reports" / "PDBE_LIGAND_SITES_FETCH_REPORT.json"
        ),
        "gpcrdb_qa_report": source_file(
            ROOT / "reports" / "GPCRDB_NORMALIZATION_V2_QA.json"
        ),
    }

    manifest = {
        "run_id": RUN_ID,
        "run_type": "incremental_supplement_not_rebuild",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "workspace_root": str(ROOT),
        "run_root": str(RUN_ROOT),
        "baseline_policy": "read_only_no_overwrite",
        "baseline_paths": baseline_paths,
        "baseline_validation": baseline_validation,
        "source_inputs": sources,
        "pubchem_state_at_run_start": pubchem_state(),
        "disk_free_bytes_at_run_start": shutil.disk_usage(ROOT.drive + "\\").free,
        "concurrency_policy": {
            "maximum_heavy_io_workers": 2,
            "minimum_free_disk_bytes": 25 * 1024**3,
            "source_staging_isolated": True,
            "shared_output_writes_during_source_normalization": False,
        },
    }
    atomic_json(RUN_ROOT / "RUN_MANIFEST.json", manifest)

    readme = RUN_ROOT / "README.md"
    if not readme.exists():
        readme.write_text(
            "\n".join(
                [
                    f"# {RUN_ID}",
                    "",
                    "Incremental supplement to the frozen MemProDB baseline.",
                    "Published V5.5/V5.4/V4.2/V1.0 inputs are read-only.",
                    "Source-specific staging outputs are merged only after QA.",
                    "",
                ]
            ),
            encoding="utf-8",
        )

    print(
        json.dumps(
            {
                "status": "initialized",
                "run_id": RUN_ID,
                "run_root": str(RUN_ROOT),
                "baseline_inputs": len(baseline_validation),
                "baseline_files_present": sum(
                    1 for item in baseline_validation if item["exists"]
                ),
                "baseline_sizes_match": sum(
                    1 for item in baseline_validation if item["size_matches"]
                ),
                "pubchem_processed_targets": manifest[
                    "pubchem_state_at_run_start"
                ]["processed_targets"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
