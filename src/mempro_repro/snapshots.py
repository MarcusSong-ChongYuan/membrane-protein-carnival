"""Auditable inventory of frozen upstream inputs.

The release package is a *release snapshot*, not an assertion that every
historical upstream API can still recreate 2026-era inputs.  This module makes
that distinction machine-readable: a missing historical input is a documented
limitation, whereas a missing or corrupted input declared as PRESENT is an
engineering error.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Iterable


REGISTRY_COLUMNS = {
    "source_id",
    "module",
    "source_database",
    "source_version",
    "release_or_retrieval_date",
    "record_scope",
    "historical_snapshot_status",
    "legacy_parser_or_workflow",
}
SNAPSHOT_COLUMNS = {
    "snapshot_id",
    "source_id",
    "snapshot_relpath",
    "bytes",
    "sha256",
    "acquisition_mode",
    "status",
    "note",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _missing_columns(rows: list[dict[str, str]], expected: set[str]) -> list[str]:
    if not rows:
        return sorted(expected)
    return sorted(expected.difference(rows[0]))


def resolve_snapshot_root(data_root: Path, override: str | None = None) -> Path:
    raw = override or os.environ.get("MEMPRO_RAW_SNAPSHOT_ROOT")
    return Path(raw).expanduser().resolve() if raw else (data_root / "06_scripts" / "legacy_handoff").resolve()


def default_registry_path(repo_root: Path) -> Path:
    return repo_root / "config" / "sources" / "source_registry_v1.tsv"


def default_manifest_path(repo_root: Path) -> Path:
    return repo_root / "config" / "sources" / "snapshot_manifest_v1.tsv"


def _duplicate_values(rows: Iterable[dict[str, str]], column: str) -> list[str]:
    values: set[str] = set()
    duplicates: set[str] = set()
    for row in rows:
        value = row.get(column, "")
        if value in values:
            duplicates.add(value)
        values.add(value)
    return sorted(duplicates)


def audit_snapshots(
    repo_root: Path,
    data_root: Path,
    snapshot_root_override: str | None = None,
    registry_path: Path | None = None,
    manifest_path: Path | None = None,
) -> dict:
    """Audit source coverage and hashes without downloading any live data."""
    registry_path = registry_path or default_registry_path(repo_root)
    manifest_path = manifest_path or default_manifest_path(repo_root)
    snapshot_root = resolve_snapshot_root(data_root, snapshot_root_override)
    failures: list[dict[str, str]] = []

    if not registry_path.exists():
        return {"status": "FAIL", "failures": [{"reason": "registry_missing", "path": str(registry_path)}]}
    if not manifest_path.exists():
        return {"status": "FAIL", "failures": [{"reason": "manifest_missing", "path": str(manifest_path)}]}

    registry = read_tsv(registry_path)
    snapshots = read_tsv(manifest_path)
    missing_registry_columns = _missing_columns(registry, REGISTRY_COLUMNS)
    missing_snapshot_columns = _missing_columns(snapshots, SNAPSHOT_COLUMNS)
    if missing_registry_columns:
        failures.append({"reason": "registry_columns_missing", "columns": ",".join(missing_registry_columns)})
    if missing_snapshot_columns:
        failures.append({"reason": "manifest_columns_missing", "columns": ",".join(missing_snapshot_columns)})

    registry_ids = {row["source_id"] for row in registry if row.get("source_id")}
    duplicates = _duplicate_values(registry, "source_id")
    if duplicates:
        failures.append({"reason": "duplicate_source_id", "source_ids": ",".join(duplicates)})
    snapshot_duplicates = _duplicate_values(snapshots, "snapshot_id")
    if snapshot_duplicates:
        failures.append({"reason": "duplicate_snapshot_id", "snapshot_ids": ",".join(snapshot_duplicates)})

    snapshot_by_source: dict[str, list[dict[str, str]]] = {}
    snapshot_results: list[dict[str, str]] = []
    for snapshot in snapshots:
        source_id = snapshot.get("source_id", "")
        if source_id not in registry_ids:
            failures.append({"reason": "snapshot_source_not_registered", "source_id": source_id})
            continue
        snapshot_by_source.setdefault(source_id, []).append(snapshot)
        expected_path = snapshot_root / snapshot["snapshot_relpath"]
        result = {
            "snapshot_id": snapshot["snapshot_id"],
            "source_id": source_id,
            "relative_path": snapshot["snapshot_relpath"],
            "expected_bytes": snapshot["bytes"],
            "expected_sha256": snapshot["sha256"],
            "observed_bytes": "",
            "observed_sha256": "",
            "verification": "",
        }
        if snapshot.get("status") != "PRESENT":
            failures.append({"reason": "manifest_snapshot_not_present", "snapshot_id": snapshot.get("snapshot_id", "")})
            result["verification"] = "INVALID_MANIFEST_STATUS"
        elif not expected_path.exists():
            failures.append({"reason": "declared_present_snapshot_missing", "snapshot_id": snapshot["snapshot_id"]})
            result["verification"] = "MISSING"
        else:
            result["observed_bytes"] = str(expected_path.stat().st_size)
            result["observed_sha256"] = sha256(expected_path)
            if result["observed_bytes"] != snapshot["bytes"]:
                failures.append({"reason": "snapshot_byte_count_mismatch", "snapshot_id": snapshot["snapshot_id"]})
                result["verification"] = "BYTE_COUNT_MISMATCH"
            elif result["observed_sha256"].lower() != snapshot["sha256"].lower():
                failures.append({"reason": "snapshot_sha256_mismatch", "snapshot_id": snapshot["snapshot_id"]})
                result["verification"] = "SHA256_MISMATCH"
            else:
                result["verification"] = "VERIFIED"
        snapshot_results.append(result)

    source_results: list[dict[str, str]] = []
    unavailable = 0
    for row in registry:
        source_id = row["source_id"]
        declared = row["historical_snapshot_status"]
        source_snapshots = snapshot_by_source.get(source_id, [])
        source_snapshot_results = [x for x in snapshot_results if x["source_id"] == source_id]
        if declared == "SNAPSHOT_PRESENT":
            state = "SNAPSHOT_VERIFIED" if source_snapshot_results and all(x["verification"] == "VERIFIED" for x in source_snapshot_results) else "SNAPSHOT_ERROR"
            rebuild = "INPUT_AVAILABLE_FOR_SCHEMA_SPECIFIC_REPLAY"
        elif declared == "DERIVED_FROM_FROZEN_FORMAL_TABLES":
            state = "DERIVED_FROM_FORMAL"
            rebuild = "REPRODUCIBLE_FROM_FORMAL_WORKFLOWS"
        else:
            state = "HISTORICAL_SNAPSHOT_UNAVAILABLE"
            rebuild = "EXACT_REBUILD_BLOCKED_PENDING_ORIGINAL_SNAPSHOT"
            unavailable += 1
        source_results.append({
            "source_id": source_id,
            "module": row["module"],
            "source_database": row["source_database"],
            "source_version": row["source_version"],
            "declared_snapshot_status": declared,
            "audit_status": state,
            "rebuild_status": rebuild,
            "snapshot_count": str(len(source_snapshots)),
            "legacy_parser_or_workflow": row["legacy_parser_or_workflow"],
        })

    verified = sum(item["verification"] == "VERIFIED" for item in snapshot_results)
    status = "FAIL" if failures else ("PASS_WITH_LIMITATIONS" if unavailable else "PASS")
    return {
        "status": status,
        "registry_path": str(registry_path),
        "snapshot_manifest_path": str(manifest_path),
        "snapshot_root": str(snapshot_root),
        "source_counts": {
            "registered": len(registry),
            "snapshot_available": sum(item["historical_snapshot_status"] == "SNAPSHOT_PRESENT" for item in registry),
            "derived_from_formal": sum(item["historical_snapshot_status"] == "DERIVED_FROM_FROZEN_FORMAL_TABLES" for item in registry),
            "historical_snapshot_unavailable": unavailable,
        },
        "snapshot_counts": {"declared": len(snapshots), "verified": verified},
        "failures": failures,
        "sources": source_results,
        "snapshots": snapshot_results,
        "interpretation": (
            "Unavailable historical snapshots are documented provenance limitations, not a reason to replace them with current live API data. "
            "Only inputs declared PRESENT are required to pass hash validation."
        ),
    }


def write_snapshot_audit(output_root: Path, payload: dict) -> tuple[Path, Path]:
    output_dir = output_root / "reproducibility"
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "raw_snapshot_audit.json"
    tsv_path = output_dir / "raw_snapshot_audit.tsv"
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    rows = payload.get("sources", [])
    fieldnames = list(rows[0]) if rows else ["source_id", "audit_status"]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return json_path, tsv_path
