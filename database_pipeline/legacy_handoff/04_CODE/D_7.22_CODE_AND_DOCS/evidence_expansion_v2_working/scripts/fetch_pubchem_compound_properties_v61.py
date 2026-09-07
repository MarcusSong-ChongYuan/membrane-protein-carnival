#!/usr/bin/env python3
"""Checkpointed PubChem PUG REST structure retrieval for V6.1."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import random
import re
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


DEFAULT_RUN = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\incremental_v61_20260727"
)
PROPERTIES = (
    "CanonicalSMILES,IsomericSMILES,InChI,InChIKey,"
    "MolecularFormula,MolecularWeight,Charge"
)
API_PREFIX = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def api_request(cids: list[int], timeout: int) -> tuple[int, dict]:
    cid_string = ",".join(str(cid) for cid in cids)
    url = f"{API_PREFIX}/compound/cid/{cid_string}/property/{PROPERTIES}/JSON"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "MemPro-V6.1-identity-resolution/1.0",
            "Accept": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status, json.load(response)
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        try:
            payload = json.loads(body)
        except json.JSONDecodeError:
            payload = {"error": body[:1000]}
        return exc.code, payload


def initialize_database(
    conn: sqlite3.Connection, inventory: Path, priority: str
) -> int:
    conn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        CREATE TABLE IF NOT EXISTS compound (
            cid INTEGER PRIMARY KEY,
            priority TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'pending',
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
            updated_utc TEXT
        );
        """
    )
    added = 0
    with gzip.open(inventory, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            if row["priority"] != priority:
                continue
            rows.append((int(row["cid"]), priority))
            if len(rows) >= 20_000:
                before = conn.total_changes
                conn.executemany(
                    "INSERT OR IGNORE INTO compound(cid,priority) VALUES (?,?)", rows
                )
                added += conn.total_changes - before
                rows.clear()
        if rows:
            before = conn.total_changes
            conn.executemany(
                "INSERT OR IGNORE INTO compound(cid,priority) VALUES (?,?)", rows
            )
            added += conn.total_changes - before
    conn.commit()
    return added


def counts(conn: sqlite3.Connection) -> dict[str, int]:
    result = {
        status: count
        for status, count in conn.execute(
            "SELECT status,COUNT(*) FROM compound GROUP BY status"
        )
    }
    result["total"] = sum(result.values())
    return result


def mark_error(
    conn: sqlite3.Connection,
    cids: list[int],
    http_status: int,
    message: str,
    terminal: bool,
) -> None:
    status = "not_found" if terminal else "retry"
    timestamp = now()
    conn.executemany(
        """
        UPDATE compound
        SET status=?,attempts=attempts+1,http_status=?,error_message=?,updated_utc=?
        WHERE cid=?
        """,
        [(status, http_status, message[:1000], timestamp, cid) for cid in cids],
    )


def save_batch(conn: sqlite3.Connection, cids: list[int], payload: dict) -> None:
    properties = payload.get("PropertyTable", {}).get("Properties", [])
    found: set[int] = set()
    timestamp = now()
    for item in properties:
        cid = int(item["CID"])
        found.add(cid)
        conn.execute(
            """
            UPDATE compound SET
                status='resolved',
                attempts=attempts+1,
                molecular_formula=?,
                molecular_weight=?,
                smiles=?,
                connectivity_smiles=?,
                inchi=?,
                inchikey=?,
                charge=?,
                http_status=200,
                error_message=NULL,
                updated_utc=?
            WHERE cid=?
            """,
            (
                item.get("MolecularFormula", ""),
                str(item.get("MolecularWeight", "")),
                item.get("SMILES", item.get("IsomericSMILES", "")),
                item.get(
                    "ConnectivitySMILES", item.get("CanonicalSMILES", "")
                ),
                item.get("InChI", ""),
                item.get("InChIKey", ""),
                item.get("Charge"),
                timestamp,
                cid,
            ),
        )
    missing = [cid for cid in cids if cid not in found]
    if missing:
        mark_error(conn, missing, 200, "CID absent from successful response", True)


def export_results(conn: sqlite3.Connection, output: Path) -> int:
    header = [
        "cid",
        "priority",
        "status",
        "attempts",
        "molecular_formula",
        "molecular_weight",
        "smiles",
        "connectivity_smiles",
        "inchi",
        "inchikey",
        "charge",
        "http_status",
        "error_message",
        "updated_utc",
    ]
    count = 0
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in conn.execute(
            "SELECT " + ",".join(header) + " FROM compound ORDER BY cid"
        ):
            writer.writerow(row)
            count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--priority", default="P1_ACTIVE_BE2")
    parser.add_argument(
        "--job-tag",
        default="p1",
        help="Short unique suffix used for checkpoint, progress, and output files.",
    )
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--requests-per-second", type=float, default=3.0)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--max-attempts", type=int, default=6)
    args = parser.parse_args()

    run = args.run.resolve()
    inventory = (
        run
        / "staging"
        / "pubchem_identity"
        / "pubchem_cid_priority_inventory_v61.tsv.gz"
    )
    output_dir = run / "staging" / "pubchem_identity"
    qa_dir = run / "qa"
    output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir.mkdir(parents=True, exist_ok=True)
    if not inventory.exists():
        raise FileNotFoundError(
            f"Priority inventory is not ready: {inventory}. Run initializer first."
        )

    job_tag = re.sub(r"[^a-z0-9_-]+", "_", args.job_tag.casefold()).strip("_")
    if not job_tag:
        raise ValueError("--job-tag must contain at least one letter or digit")
    db_path = output_dir / f"pubchem_property_fetch_{job_tag}_v61.sqlite"
    progress_path = qa_dir / f"PUBCHEM_IDENTITY_FETCH_{job_tag.upper()}_V61_PROGRESS.json"
    output_path = output_dir / f"pubchem_properties_{job_tag}_v61.tsv.gz"
    conn = sqlite3.connect(db_path)
    added = initialize_database(conn, inventory, args.priority)
    started = now()
    min_interval = 1.0 / max(args.requests_per_second, 0.1)
    batches = 0

    write_json(
        progress_path,
        {
            "status": "running",
            "started_utc": started,
            "updated_utc": now(),
            "priority": args.priority,
            "new_cids_added": added,
            "counts": counts(conn),
            "endpoint": API_PREFIX,
            "requests_per_second": args.requests_per_second,
        },
    )

    while True:
        cids = [
            row[0]
            for row in conn.execute(
                """
                SELECT cid FROM compound
                WHERE status IN ('pending','retry') AND attempts<?
                ORDER BY attempts,cid LIMIT ?
                """,
                (args.max_attempts, args.batch_size),
            )
        ]
        if not cids:
            break

        request_started = time.monotonic()
        try:
            status, payload = api_request(cids, args.timeout)
        except Exception as exc:
            status, payload = 0, {"error": f"{type(exc).__name__}: {exc}"}

        if status == 200 and "PropertyTable" in payload:
            save_batch(conn, cids, payload)
        elif status in (400, 404) and len(cids) > 1:
            # Requeue the batch as smaller requests so one invalid CID cannot hide
            # valid compounds. A one-CID 400/404 becomes terminal not_found.
            conn.executemany(
                """
                UPDATE compound
                SET attempts=attempts+1,status='retry',http_status=?,error_message=?,
                    updated_utc=?
                WHERE cid=?
                """,
                [
                    (
                        status,
                        "batch split after PubChem 400/404",
                        now(),
                        cid,
                    )
                    for cid in cids
                ],
            )
            args.batch_size = max(1, min(args.batch_size // 2, 25))
        else:
            terminal = status in (400, 404) and len(cids) == 1
            message = json.dumps(payload, ensure_ascii=False)[:1000]
            mark_error(conn, cids, status, message, terminal)

        conn.commit()
        batches += 1
        if batches % 10 == 0:
            state = counts(conn)
            write_json(
                progress_path,
                {
                    "status": "running",
                    "started_utc": started,
                    "updated_utc": now(),
                    "priority": args.priority,
                    "batches_completed": batches,
                    "counts": state,
                    "endpoint": API_PREFIX,
                    "requests_per_second": args.requests_per_second,
                },
            )

        if status in (429, 500, 502, 503, 504, 0):
            time.sleep(min(60.0, 2.0 + random.random() * 3.0))
        elapsed = time.monotonic() - request_started
        if elapsed < min_interval:
            time.sleep(min_interval - elapsed)

    # Anything that exhausted the retry budget gets a distinct terminal state.
    conn.execute(
        """
        UPDATE compound SET status='failed_after_retries',updated_utc=?
        WHERE status='retry' AND attempts>=?
        """,
        (now(), args.max_attempts),
    )
    conn.commit()
    exported = export_results(conn, output_path)
    final_counts = counts(conn)
    terminal = final_counts.get("pending", 0) == 0 and final_counts.get("retry", 0) == 0
    write_json(
        progress_path,
        {
            "status": "complete" if terminal else "incomplete",
            "started_utc": started,
            "completed_utc": now(),
            "priority": args.priority,
            "batches_completed": batches,
            "counts": final_counts,
            "exported_rows": exported,
            "output": str(output_path),
            "database": str(db_path),
            "endpoint": API_PREFIX,
            "requests_per_second": args.requests_per_second,
        },
    )
    conn.close()
    print(json.dumps(json.loads(progress_path.read_text(encoding="utf-8")), indent=2))


if __name__ == "__main__":
    main()
