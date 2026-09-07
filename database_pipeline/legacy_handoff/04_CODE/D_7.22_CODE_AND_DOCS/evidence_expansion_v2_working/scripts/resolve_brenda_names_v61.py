#!/usr/bin/env python3
"""Resolve BRENDA compound names through structure, never by name alone."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RUN = ROOT / "runs" / "incremental_v61_20260727"
RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
QUEUE = (
    ROOT
    / "runs"
    / "incremental_20260727"
    / "staging"
    / "brenda"
    / "brenda_compound_identity_queue_v3.tsv.gz"
)
OUTDIR = RUN / "staging" / "brenda_identity"
PROGRESS = RUN / "qa" / "BRENDA_IDENTITY_V61_PROGRESS.json"
PROPERTIES = (
    "CanonicalSMILES,IsomericSMILES,InChI,InChIKey,"
    "MolecularFormula,MolecularWeight,Charge"
)
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def request_properties(
    namespace: str, value: str, timeout: int = 60
) -> tuple[int, list[dict], str]:
    url = f"{PUG}/compound/{namespace}/property/{PROPERTIES}/JSON"
    body = urllib.parse.urlencode({namespace: value}).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=body,
        method="POST",
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "application/json",
            "User-Agent": "MemPro-V6.1-BRENDA-identity/1.0",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.load(response)
            return (
                response.status,
                payload.get("PropertyTable", {}).get("Properties", []),
                "",
            )
    except urllib.error.HTTPError as exc:
        try:
            body_text = exc.read().decode("utf-8", errors="replace")
        except Exception as body_exc:
            body_text = (
                f"HTTP error body unavailable: "
                f"{type(body_exc).__name__}:{body_exc}"
            )
        return exc.code, [], body_text[:1000]
    except Exception as exc:
        return 0, [], f"{type(exc).__name__}:{exc}"


def load_master(path: Path):
    full: dict[str, list[tuple[str, str]]] = defaultdict(list)
    connectivity: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    by_id: dict[str, tuple[str, str]] = {}
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            compound_id = row["compound_internal_id"]
            name = row["preferred_name"]
            key = row["standard_inchikey"].strip().upper()
            conn_key = row["connectivity_key"].strip().upper()
            by_id[compound_id] = (key, name)
            if key:
                full[key].append((compound_id, name))
            if conn_key:
                connectivity[conn_key].append((compound_id, name, key))
    return full, connectivity, by_id


def split_semicolon(value: str) -> list[str]:
    return [part.strip() for part in value.split(";")]


def all_cids(value: str) -> list[str]:
    return sorted(
        {
            item.strip()
            for group in split_semicolon(value)
            for item in group.replace(",", "|").split("|")
            if item.strip().isdigit()
        },
        key=int,
    )


def initialize(conn: sqlite3.Connection) -> int:
    conn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        CREATE TABLE IF NOT EXISTS name_resolution (
            normalized_name TEXT PRIMARY KEY,
            priority_rank INTEGER,
            preferred_name TEXT,
            evidence_row_count INTEGER,
            target_count INTEGER,
            generic_flag INTEGER,
            original_queue_status TEXT,
            candidate_compound_internal_ids TEXT,
            candidate_pubchem_cids TEXT,
            candidate_standard_inchikeys TEXT,
            resolution_status TEXT NOT NULL DEFAULT 'pending',
            pubchem_result_cids TEXT,
            pubchem_result_inchikeys TEXT,
            matched_compound_internal_ids TEXT,
            structure_smiles TEXT,
            attempts INTEGER NOT NULL DEFAULT 0,
            http_status INTEGER,
            message TEXT,
            updated_utc TEXT
        );
        CREATE TABLE IF NOT EXISTS cid_cache (
            cid INTEGER PRIMARY KEY,
            fetch_status TEXT,
            attempts INTEGER NOT NULL DEFAULT 0,
            inchikey TEXT,
            smiles TEXT,
            inchi TEXT,
            http_status INTEGER,
            message TEXT
        );
        """
    )
    added = 0
    with gzip.open(QUEUE, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            before = conn.total_changes
            conn.execute(
                """
                INSERT OR IGNORE INTO name_resolution(
                    normalized_name,priority_rank,preferred_name,evidence_row_count,
                    target_count,generic_flag,original_queue_status,
                    candidate_compound_internal_ids,candidate_pubchem_cids,
                    candidate_standard_inchikeys
                ) VALUES (?,?,?,?,?,?,?,?,?,?)
                """,
                (
                    row["compound_name_normalized"],
                    int(row["priority_rank"]),
                    row["preferred_source_name"],
                    int(row["evidence_row_count"]),
                    int(row["target_count"]),
                    int(row["generic_or_common"]),
                    row["identity_queue_status"],
                    row["candidate_compound_internal_ids"],
                    row["candidate_pubchem_cids"],
                    row["candidate_standard_inchikeys"],
                ),
            )
            added += conn.total_changes - before
    conn.execute(
        """
        UPDATE name_resolution
        SET resolution_status='excluded_generic_or_common',updated_utc=?
        WHERE generic_flag=1 AND resolution_status='pending'
        """,
        (now(),),
    )
    conn.commit()
    return added


def state_counts(conn: sqlite3.Connection) -> dict[str, int]:
    result = dict(
        conn.execute(
            "SELECT resolution_status,COUNT(*) FROM name_resolution GROUP BY 1"
        ).fetchall()
    )
    result["total"] = sum(result.values())
    return result


def fetch_candidate_cids(conn: sqlite3.Connection) -> None:
    cid_values: set[str] = set()
    for (value,) in conn.execute(
        """
        SELECT candidate_pubchem_cids FROM name_resolution
        WHERE candidate_pubchem_cids IS NOT NULL AND candidate_pubchem_cids<>''
        """
    ):
        cid_values.update(all_cids(value))
    conn.executemany(
        "INSERT OR IGNORE INTO cid_cache(cid,fetch_status) VALUES (?,'pending')",
        [(int(value),) for value in cid_values],
    )
    conn.commit()
    while True:
        batch = [
            str(row[0])
            for row in conn.execute(
                """
                SELECT cid FROM cid_cache
                WHERE fetch_status='pending' AND attempts<5
                ORDER BY attempts,cid LIMIT 100
                """
            )
        ]
        if not batch:
            break
        status, properties, message = request_properties("cid", ",".join(batch))
        found = set()
        if status == 200:
            for item in properties:
                cid = int(item["CID"])
                found.add(str(cid))
                conn.execute(
                    """
                    UPDATE cid_cache SET fetch_status='resolved',inchikey=?,smiles=?,
                        inchi=?,attempts=attempts+1,http_status=200,message=NULL
                    WHERE cid=?
                    """,
                    (
                        item.get("InChIKey", ""),
                        item.get("SMILES", item.get("IsomericSMILES", "")),
                        item.get("InChI", ""),
                        cid,
                    ),
                )
            missing = [int(value) for value in batch if value not in found]
            conn.executemany(
                """
                UPDATE cid_cache SET fetch_status='not_found',http_status=200,
                    attempts=attempts+1,message='CID absent from response'
                WHERE cid=?
                """,
                [(value,) for value in missing],
            )
        elif status in (400, 404):
            conn.executemany(
                """
                UPDATE cid_cache SET fetch_status='not_found',
                    attempts=attempts+1,http_status=?,message=?
                WHERE cid=?
                """,
                [(status, message, int(value)) for value in batch],
            )
        else:
            conn.executemany(
                """
                UPDATE cid_cache SET attempts=attempts+1,http_status=?,message=?
                WHERE cid=?
                """,
                [(status, message, int(value)) for value in batch],
            )
        conn.commit()
        time.sleep(0.35)
    conn.execute(
        """
        UPDATE cid_cache SET fetch_status='failed_after_retries'
        WHERE fetch_status='pending' AND attempts>=5
        """
    )
    conn.commit()


def verify_existing_candidates(
    conn: sqlite3.Connection, by_id: dict[str, tuple[str, str]]
) -> None:
    query = """
        SELECT normalized_name,candidate_compound_internal_ids,
               candidate_pubchem_cids,candidate_standard_inchikeys
        FROM name_resolution
        WHERE resolution_status='pending'
          AND candidate_compound_internal_ids<>''
    """
    for name, ids_text, cid_text, keys_text in conn.execute(query).fetchall():
        ids = split_semicolon(ids_text)
        cid_groups = split_semicolon(cid_text)
        key_groups = split_semicolon(keys_text)
        exact_ids = set()
        connectivity_ids = set()
        result_cids = set()
        result_keys = set()
        for index, compound_id in enumerate(ids):
            if not compound_id:
                continue
            master_key = by_id.get(compound_id, ("", ""))[0].upper()
            stated_key = (
                key_groups[index].upper()
                if index < len(key_groups) and key_groups[index]
                else master_key
            )
            cids = (
                [
                    value
                    for value in cid_groups[index].replace(",", "|").split("|")
                    if value.isdigit()
                ]
                if index < len(cid_groups)
                else []
            )
            for cid in cids:
                cached = conn.execute(
                    "SELECT inchikey FROM cid_cache WHERE cid=?", (int(cid),)
                ).fetchone()
                if not cached or not cached[0]:
                    continue
                fetched_key = cached[0].upper()
                result_cids.add(cid)
                result_keys.add(fetched_key)
                if fetched_key == stated_key or fetched_key == master_key:
                    exact_ids.add(compound_id)
                elif master_key and fetched_key[:14] == master_key[:14]:
                    connectivity_ids.add(compound_id)
        if len(exact_ids) == 1:
            status = "verified_existing_exact_structure"
            matched = exact_ids
        elif len(exact_ids) > 1:
            status = "ambiguous_multiple_verified_structures"
            matched = exact_ids
        elif len(connectivity_ids) == 1:
            status = "verified_connectivity_form_or_stereo_review"
            matched = connectivity_ids
        elif len(connectivity_ids) > 1:
            status = "ambiguous_connectivity_candidates"
            matched = connectivity_ids
        else:
            status = "candidate_structure_not_verified"
            matched = set()
        conn.execute(
            """
            UPDATE name_resolution SET resolution_status=?,pubchem_result_cids=?,
                pubchem_result_inchikeys=?,matched_compound_internal_ids=?,
                attempts=attempts+1,updated_utc=? WHERE normalized_name=?
            """,
            (
                status,
                ";".join(sorted(result_cids, key=int)),
                ";".join(sorted(result_keys)),
                ";".join(sorted(matched)),
                now(),
                name,
            ),
        )
    conn.commit()


def resolve_remaining_names(
    conn: sqlite3.Connection,
    full_map: dict[str, list[tuple[str, str]]],
    connectivity_map: dict[str, list[tuple[str, str, str]]],
) -> None:
    processed = 0
    while True:
        item = conn.execute(
            """
            SELECT normalized_name,preferred_name,attempts FROM name_resolution
            WHERE resolution_status='pending'
              AND attempts<5
            ORDER BY attempts,evidence_row_count DESC,priority_rank LIMIT 1
            """
        ).fetchone()
        if not item:
            break
        normalized_name, preferred_name, attempts = item
        status, properties, message = request_properties("name", preferred_name)
        if status in (0, 429, 500, 502, 503, 504):
            next_status = (
                "failed_after_retries" if attempts + 1 >= 5 else "pending"
            )
            conn.execute(
                """
                UPDATE name_resolution SET resolution_status=?,
                    attempts=attempts+1,http_status=?,message=?,updated_utc=?
                WHERE normalized_name=?
                """,
                (next_status, status, message, now(), normalized_name),
            )
            conn.commit()
            time.sleep(3)
            continue
        keys = {
            value.get("InChIKey", "").upper()
            for value in properties
            if value.get("InChIKey")
        }
        cids = {
            str(value["CID"])
            for value in properties
            if value.get("CID") is not None
        }
        exact_ids = {
            compound_id
            for key in keys
            for compound_id, _name in full_map.get(key, [])
        }
        connection_ids = {
            compound_id
            for key in keys
            for compound_id, _name, _full_key in connectivity_map.get(
                key[:14], []
            )
        }
        if status in (400, 404) or not properties:
            resolution = "unresolved_no_pubchem_name_match"
            matches = set()
        elif len(keys) > 1:
            resolution = "ambiguous_pubchem_name_multiple_structures"
            matches = exact_ids
        elif len(exact_ids) == 1:
            resolution = "resolved_unique_name_to_existing_exact_structure"
            matches = exact_ids
        elif len(exact_ids) > 1:
            resolution = "ambiguous_existing_exact_duplicates"
            matches = exact_ids
        elif len(connection_ids) == 1:
            resolution = "unique_name_connectivity_form_or_stereo_review"
            matches = connection_ids
        elif len(connection_ids) > 1:
            resolution = "ambiguous_name_connectivity_candidates"
            matches = connection_ids
        elif len(keys) == 1:
            resolution = "new_structure_candidate_from_unique_name"
            matches = set()
        else:
            resolution = "unresolved_pubchem_response"
            matches = set()
        smiles = ";".join(
            sorted(
                {
                    value.get("SMILES", value.get("IsomericSMILES", ""))
                    for value in properties
                    if value.get("SMILES") or value.get("IsomericSMILES")
                }
            )
        )
        conn.execute(
            """
            UPDATE name_resolution SET resolution_status=?,pubchem_result_cids=?,
                pubchem_result_inchikeys=?,matched_compound_internal_ids=?,
                structure_smiles=?,attempts=attempts+1,http_status=?,message=?,
                updated_utc=? WHERE normalized_name=?
            """,
            (
                resolution,
                ";".join(sorted(cids, key=int)),
                ";".join(sorted(keys)),
                ";".join(sorted(matches)),
                smiles,
                status,
                message,
                now(),
                normalized_name,
            ),
        )
        conn.commit()
        processed += 1
        if processed % 100 == 0:
            write_json(
                PROGRESS,
                {
                    "status": "running",
                    "updated_utc": now(),
                    "names_resolved_this_run": processed,
                    "counts": state_counts(conn),
                    "endpoint": PUG,
                },
            )
        time.sleep(0.35)


def export(conn: sqlite3.Connection, path: Path) -> int:
    columns = [
        "normalized_name",
        "priority_rank",
        "preferred_name",
        "evidence_row_count",
        "target_count",
        "generic_flag",
        "original_queue_status",
        "candidate_compound_internal_ids",
        "candidate_pubchem_cids",
        "candidate_standard_inchikeys",
        "resolution_status",
        "pubchem_result_cids",
        "pubchem_result_inchikeys",
        "matched_compound_internal_ids",
        "structure_smiles",
        "attempts",
        "http_status",
        "message",
        "updated_utc",
    ]
    count = 0
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in conn.execute(
            "SELECT " + ",".join(columns) + " FROM name_resolution ORDER BY priority_rank"
        ):
            writer.writerow(row)
            count += 1
    return count


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    started = now()
    db_path = OUTDIR / "brenda_name_resolution_v61.sqlite"
    output = OUTDIR / "brenda_name_identity_v61.tsv.gz"
    conn = sqlite3.connect(db_path)
    added = initialize(conn)
    write_json(
        PROGRESS,
        {
            "status": "running",
            "started_utc": started,
            "updated_utc": now(),
            "new_names_added": added,
            "counts": state_counts(conn),
            "endpoint": PUG,
        },
    )
    full_map, connectivity_map, by_id = load_master(
        RELEASE / "small_molecule_master_v1_1.tsv"
    )
    fetch_candidate_cids(conn)
    verify_existing_candidates(conn, by_id)
    resolve_remaining_names(conn, full_map, connectivity_map)
    exported = export(conn, output)
    final_counts = state_counts(conn)
    evidence_by_status = dict(
        conn.execute(
            """
            SELECT resolution_status,SUM(evidence_row_count)
            FROM name_resolution GROUP BY resolution_status
            """
        ).fetchall()
    )
    conn.close()
    summary = {
        "status": "complete",
        "started_utc": started,
        "completed_utc": now(),
        "name_counts": final_counts,
        "evidence_rows_by_resolution": evidence_by_status,
        "exported_names": exported,
        "output": str(output),
        "database": str(db_path),
        "endpoint": PUG,
    }
    write_json(PROGRESS, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
