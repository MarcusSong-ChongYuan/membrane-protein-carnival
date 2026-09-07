#!/usr/bin/env python3
"""Inventory all unresolved negative-evidence CIDs and reuse V6.1 structures."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import os
import re
import sqlite3
from collections import Counter
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_1_20260729"
RUN = ROOT / "runs" / "v62_completion_20260729"
STAGING = RUN / "staging" / "negative_identity"
QA = RUN / "qa"
DB = STAGING / "negative_cid_identity_v62.sqlite"
PROGRESS = QA / "NEGATIVE_CID_INVENTORY_V62_PROGRESS.json"
NEGATIVE = RELEASE / "negative_binding_evidence_unmapped_review_v1_1.tsv.gz"
PROPERTY_FILES = [
    ROOT
    / "runs"
    / "incremental_v61_20260727"
    / "staging"
    / "pubchem_identity"
    / f"pubchem_properties_{tag}_v61.tsv.gz"
    for tag in ("p1", "p2", "p3")
]
CID_RE = re.compile(r"(?:CID[: ]?)(\d+)", re.IGNORECASE)


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(payload: dict) -> None:
    tmp = PROGRESS.with_suffix(PROGRESS.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    os.replace(tmp, PROGRESS)


def counts(conn: sqlite3.Connection) -> dict[str, int]:
    result = {
        status: count
        for status, count in conn.execute(
            "SELECT fetch_status,COUNT(*) FROM cid_inventory GROUP BY fetch_status"
        )
    }
    result["total"] = sum(result.values())
    return result


def initialize(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        PRAGMA temp_store=MEMORY;
        CREATE TABLE IF NOT EXISTS cid_inventory (
            cid INTEGER PRIMARY KEY,
            evidence_rows INTEGER NOT NULL DEFAULT 0,
            fetch_status TEXT NOT NULL DEFAULT 'pending',
            attempts INTEGER NOT NULL DEFAULT 0,
            molecular_formula TEXT,
            molecular_weight TEXT,
            smiles TEXT,
            connectivity_smiles TEXT,
            inchi TEXT,
            inchikey TEXT,
            charge INTEGER,
            http_status INTEGER,
            error_message TEXT,
            structure_source TEXT,
            updated_utc TEXT
        );
        """
    )
    conn.commit()


def flush_counter(conn: sqlite3.Connection, counter: Counter[int]) -> None:
    conn.executemany(
        """
        INSERT INTO cid_inventory(cid,evidence_rows)
        VALUES (?,?)
        ON CONFLICT(cid) DO UPDATE
        SET evidence_rows=evidence_rows+excluded.evidence_rows
        """,
        counter.items(),
    )
    conn.commit()
    counter.clear()


def inventory_negative(conn: sqlite3.Connection, started: str) -> int:
    existing = conn.execute(
        "SELECT COALESCE(SUM(evidence_rows),0) FROM cid_inventory"
    ).fetchone()[0]
    if existing:
        return int(existing)
    counter: Counter[int] = Counter()
    scanned = 0
    missing_cid = 0
    with gzip.open(NEGATIVE, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            scanned += 1
            value = (row.get("pubchem_cid") or "").strip()
            if not value.isdigit():
                match = CID_RE.search(
                    (row.get("source_record_id") or "")
                    + " "
                    + (row.get("compound_source_id") or "")
                )
                value = match.group(1) if match else ""
            if value.isdigit():
                counter[int(value)] += 1
            else:
                missing_cid += 1
            if scanned % 250_000 == 0:
                flush_counter(conn, counter)
                write_json(
                    {
                        "status": "running",
                        "stage": "inventory_negative_cids",
                        "started_utc": started,
                        "updated_utc": now(),
                        "evidence_rows_scanned": scanned,
                        "rows_without_cid": missing_cid,
                        "counts": counts(conn),
                    }
                )
    if counter:
        flush_counter(conn, counter)
    write_json(
        {
            "status": "running",
            "stage": "reuse_v61_property_cache",
            "started_utc": started,
            "updated_utc": now(),
            "evidence_rows_scanned": scanned,
            "rows_without_cid": missing_cid,
            "counts": counts(conn),
        }
    )
    return scanned


def seed_properties(conn: sqlite3.Connection) -> dict[str, int]:
    reused: dict[str, int] = {}
    statement = """
        UPDATE cid_inventory SET
            fetch_status=?,
            attempts=?,
            molecular_formula=?,
            molecular_weight=?,
            smiles=?,
            connectivity_smiles=?,
            inchi=?,
            inchikey=?,
            charge=?,
            http_status=?,
            error_message=?,
            structure_source=?,
            updated_utc=?
        WHERE cid=?
    """
    for path in PROPERTY_FILES:
        if not path.exists():
            reused[path.name] = 0
            continue
        changed_before = conn.total_changes
        batch = []
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                status = row.get("status", "")
                if status != "resolved":
                    continue
                batch.append(
                    (
                        "resolved_cache",
                        int(row.get("attempts") or 0),
                        row.get("molecular_formula", ""),
                        row.get("molecular_weight", ""),
                        row.get("smiles", ""),
                        row.get("connectivity_smiles", ""),
                        row.get("inchi", ""),
                        row.get("inchikey", "").upper(),
                        row.get("charge") or None,
                        int(row.get("http_status") or 200),
                        row.get("error_message", ""),
                        path.name,
                        row.get("updated_utc", ""),
                        int(row["cid"]),
                    )
                )
                if len(batch) >= 20_000:
                    conn.executemany(statement, batch)
                    conn.commit()
                    batch.clear()
        if batch:
            conn.executemany(statement, batch)
            conn.commit()
        reused[path.name] = conn.total_changes - changed_before
    return reused


def main() -> None:
    STAGING.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    started = now()
    conn = sqlite3.connect(DB)
    initialize(conn)
    scanned = inventory_negative(conn, started)
    reused = seed_properties(conn)
    state = {
        "status": "complete",
        "stage": "inventory_and_cache_reuse_complete",
        "started_utc": started,
        "completed_utc": now(),
        "negative_review_file": str(NEGATIVE),
        "negative_evidence_rows_scanned": scanned,
        "unique_cid_expected_from_prior_audit": 1_622_591,
        "counts": counts(conn),
        "reused_property_rows_by_file": reused,
        "database": str(DB),
    }
    state["unique_cid_count_matches_expected"] = (
        state["counts"]["total"] == state["unique_cid_expected_from_prior_audit"]
    )
    write_json(state)
    conn.close()
    print(json.dumps(state, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
