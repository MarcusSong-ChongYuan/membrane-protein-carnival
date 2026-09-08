"""Recover and checksum legacy source artifacts without changing FORMAL data."""

from __future__ import annotations

import csv
import hashlib
import shutil
from datetime import datetime, timezone
from pathlib import Path


PLAN_FIELDNAMES = (
    "snapshot_id", "source_id", "source_database", "source_relpath", "target_relpath",
    "artifact_type", "release_equivalence", "note",
)
PLAN_COLUMNS = set(PLAN_FIELDNAMES)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_plan(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    missing = PLAN_COLUMNS.difference(rows[0] if rows else {})
    if missing:
        raise ValueError(f"Recovery plan lacks columns: {', '.join(sorted(missing))}")
    if len({row['snapshot_id'] for row in rows}) != len(rows):
        raise ValueError("Recovery plan contains duplicate snapshot_id values")
    return rows


def default_plan_path(repo_root: Path) -> Path:
    return repo_root / "config" / "sources" / "legacy_context_recovery_v1.tsv"


def recover_legacy_context(
    repo_root: Path,
    legacy_source_root: Path,
    destination: Path,
    plan_path: Path | None = None,
) -> dict:
    """Copy available legacy inputs into a new, checksum-indexed context bundle.

    This intentionally refuses a non-empty destination.  The bundle is an
    archival recovery layer and is never used to overwrite a frozen V7.2 table.
    """
    plan_path = plan_path or default_plan_path(repo_root)
    plan = read_plan(plan_path)
    destination = destination.resolve()
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(f"Destination must be empty: {destination}")
    destination.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, str]] = []
    for row in plan:
        source = legacy_source_root / row["source_relpath"]
        target = destination / row["target_relpath"]
        record = dict(row)
        record.update({"bytes": "", "sha256": "", "recovery_status": ""})
        if not source.is_file():
            record["recovery_status"] = "SOURCE_NOT_FOUND"
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            record["bytes"] = str(target.stat().st_size)
            record["sha256"] = sha256(target)
            record["recovery_status"] = "RECOVERED_AND_VERIFIED"
        records.append(record)

    manifest = destination / "LEGACY_CONTEXT_MANIFEST.tsv"
    fieldnames = list(PLAN_FIELDNAMES) + ["bytes", "sha256", "recovery_status"]
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    recovered = sum(row["recovery_status"] == "RECOVERED_AND_VERIFIED" for row in records)
    unavailable = len(records) - recovered
    (destination / "README.md").write_text(
        "# MemPro legacy source-context recovery\n\n"
        "This package preserves raw exports and API caches recovered from an earlier MemPro pipeline. "
        "It is intentionally separate from the frozen V7.2 FORMAL release. Every record is tagged "
        "NOT_EQUIVALENT_TO_V72 or NOT_IN_V72_SOURCE_REGISTRY unless exact V7.2 equivalence has been independently demonstrated.\n\n"
        f"Created: {datetime.now(timezone.utc).isoformat()}\n\n"
        f"Recovered: {recovered}/{len(records)} artifacts. Source not found: {unavailable}.\n\n"
        "Use LEGACY_CONTEXT_MANIFEST.tsv to validate every retained artifact. Never overwrite frozen FORMAL tables with this package.\n",
        encoding="utf-8",
    )
    return {
        "status": "PASS" if not unavailable else "PASS_WITH_LIMITATIONS",
        "legacy_source_root": str(legacy_source_root),
        "destination": str(destination),
        "plan_path": str(plan_path),
        "planned": len(records),
        "recovered": recovered,
        "source_not_found": unavailable,
        "manifest": str(manifest),
        "records": records,
    }


def verify_legacy_context(destination: Path) -> dict:
    """Verify every file in a recovered context bundle against its manifest."""
    destination = destination.resolve()
    manifest = destination / "LEGACY_CONTEXT_MANIFEST.tsv"
    if not manifest.exists():
        return {"status": "FAIL", "failures": [{"reason": "manifest_missing", "path": str(manifest)}]}
    with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = PLAN_COLUMNS | {"bytes", "sha256", "recovery_status"}
    missing_columns = expected.difference(rows[0] if rows else {})
    if missing_columns:
        return {"status": "FAIL", "failures": [{"reason": "manifest_columns_missing", "columns": ",".join(sorted(missing_columns))}]}

    failures: list[dict[str, str]] = []
    checked = 0
    for row in rows:
        if row["recovery_status"] != "RECOVERED_AND_VERIFIED":
            continue
        checked += 1
        target = destination / row["target_relpath"]
        if not target.exists():
            failures.append({"snapshot_id": row["snapshot_id"], "reason": "missing"})
        elif target.stat().st_size != int(row["bytes"]):
            failures.append({"snapshot_id": row["snapshot_id"], "reason": "byte_count_mismatch"})
        elif sha256(target).lower() != row["sha256"].lower():
            failures.append({"snapshot_id": row["snapshot_id"], "reason": "sha256_mismatch"})
    return {
        "status": "PASS" if not failures else "FAIL",
        "destination": str(destination),
        "checked": checked,
        "failures": failures,
        "interpretation": "A passing context bundle preserves recovered legacy artifacts only; it does not establish V7.2 source equivalence.",
    }
