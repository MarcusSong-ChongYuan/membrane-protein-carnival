#!/usr/bin/env python3
"""Checkpointed PubChem structure retrieval for unresolved negative CIDs."""

from __future__ import annotations

import datetime as dt
import gzip
import csv
import json
import random
import sqlite3
import time
import urllib.error
import urllib.request
from pathlib import Path


RUN = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\v62_completion_20260729"
)
DB = RUN / "staging" / "negative_identity" / "negative_cid_identity_v62.sqlite"
OUTPUT = (
    RUN
    / "staging"
    / "negative_identity"
    / "negative_cid_properties_v62.tsv.gz"
)
PROGRESS = RUN / "qa" / "NEGATIVE_CID_FETCH_V62_PROGRESS.json"
API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPERTIES = (
    "CanonicalSMILES,IsomericSMILES,InChI,InChIKey,"
    "MolecularFormula,MolecularWeight,Charge"
)


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(payload: dict) -> None:
    tmp = PROGRESS.with_suffix(PROGRESS.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(PROGRESS)


def counts(conn: sqlite3.Connection) -> dict[str, int]:
    result = {
        status: count
        for status, count in conn.execute(
            "SELECT fetch_status,COUNT(*) FROM cid_inventory GROUP BY fetch_status"
        )
    }
    result["total"] = sum(result.values())
    return result


def request(cids: list[int]) -> tuple[int, dict]:
    url = (
        f"{API}/compound/cid/{','.join(map(str, cids))}"
        f"/property/{PROPERTIES}/JSON"
    )
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "MemPro-V6.2-negative-identity/1.0",
            "Accept": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=90) as response:
            return response.status, json.load(response)
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        try:
            return exc.code, json.loads(body)
        except json.JSONDecodeError:
            return exc.code, {"error": body[:1000]}


def mark_error(
    conn: sqlite3.Connection,
    cids: list[int],
    status: int,
    message: str,
    terminal: bool,
) -> None:
    conn.executemany(
        """
        UPDATE cid_inventory
        SET fetch_status=?,attempts=attempts+1,http_status=?,
            error_message=?,updated_utc=?
        WHERE cid=?
        """,
        [
            (
                "not_found" if terminal else "retry",
                status,
                message[:1000],
                now(),
                cid,
            )
            for cid in cids
        ],
    )


def save(conn: sqlite3.Connection, cids: list[int], payload: dict) -> None:
    found = set()
    for item in payload.get("PropertyTable", {}).get("Properties", []):
        cid = int(item["CID"])
        found.add(cid)
        conn.execute(
            """
            UPDATE cid_inventory SET
                fetch_status='resolved_pubchem',
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
                structure_source='PubChem PUG REST v62',
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
                item.get("InChIKey", "").upper(),
                item.get("Charge"),
                now(),
                cid,
            ),
        )
    missing = [cid for cid in cids if cid not in found]
    if missing:
        mark_error(
            conn,
            missing,
            200,
            "CID absent from successful response",
            True,
        )


def export(conn: sqlite3.Connection) -> int:
    columns = [
        "cid",
        "evidence_rows",
        "fetch_status",
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
        "structure_source",
        "updated_utc",
    ]
    count = 0
    with gzip.open(OUTPUT, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in conn.execute(
            "SELECT " + ",".join(columns) + " FROM cid_inventory ORDER BY cid"
        ):
            writer.writerow(row)
            count += 1
    return count


def main() -> None:
    if not DB.exists():
        raise FileNotFoundError(
            f"CID inventory is not ready: {DB}. Run prepare script first."
        )
    conn = sqlite3.connect(DB)
    started = now()
    batch_size = 100
    batches = 0
    write_json(
        {
            "status": "running",
            "started_utc": started,
            "updated_utc": now(),
            "counts": counts(conn),
            "requests_per_second": 3.0,
        }
    )
    while True:
        cids = [
            row[0]
            for row in conn.execute(
                """
                SELECT cid FROM cid_inventory
                WHERE fetch_status IN ('pending','retry') AND attempts<6
                ORDER BY attempts,cid LIMIT ?
                """,
                (batch_size,),
            )
        ]
        if not cids:
            break
        started_request = time.monotonic()
        try:
            status, payload = request(cids)
        except Exception as exc:
            status, payload = 0, {"error": f"{type(exc).__name__}: {exc}"}
        if status == 200 and "PropertyTable" in payload:
            save(conn, cids, payload)
        elif status in (400, 404) and len(cids) > 1:
            conn.executemany(
                """
                UPDATE cid_inventory
                SET fetch_status='retry',attempts=attempts+1,http_status=?,
                    error_message='batch split after PubChem 400/404',
                    updated_utc=?
                WHERE cid=?
                """,
                [(status, now(), cid) for cid in cids],
            )
            batch_size = max(1, min(batch_size // 2, 25))
        else:
            mark_error(
                conn,
                cids,
                status,
                json.dumps(payload, ensure_ascii=False),
                status in (400, 404) and len(cids) == 1,
            )
        conn.commit()
        batches += 1
        if batches % 10 == 0:
            write_json(
                {
                    "status": "running",
                    "started_utc": started,
                    "updated_utc": now(),
                    "batches_completed": batches,
                    "counts": counts(conn),
                    "requests_per_second": 3.0,
                }
            )
        if status in (0, 429, 500, 502, 503, 504):
            time.sleep(min(60, 2 + random.random() * 3))
        elapsed = time.monotonic() - started_request
        if elapsed < 1 / 3:
            time.sleep(1 / 3 - elapsed)

    conn.execute(
        """
        UPDATE cid_inventory
        SET fetch_status='failed_after_retries',updated_utc=?
        WHERE fetch_status='retry' AND attempts>=6
        """,
        (now(),),
    )
    conn.commit()
    exported = export(conn)
    final_counts = counts(conn)
    complete = (
        final_counts.get("pending", 0) == 0
        and final_counts.get("retry", 0) == 0
    )
    state = {
        "status": "complete" if complete else "incomplete",
        "started_utc": started,
        "completed_utc": now(),
        "batches_completed": batches,
        "counts": final_counts,
        "exported_rows": exported,
        "output": str(OUTPUT),
        "database": str(DB),
    }
    write_json(state)
    conn.close()
    print(json.dumps(state, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
