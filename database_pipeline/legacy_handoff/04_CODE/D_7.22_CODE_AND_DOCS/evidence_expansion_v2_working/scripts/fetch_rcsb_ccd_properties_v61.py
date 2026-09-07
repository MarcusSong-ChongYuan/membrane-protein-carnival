#!/usr/bin/env python3
"""Resolve PDBe ligand component identities with the official wwPDB/RCSB CCD."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import sqlite3
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
DEFAULT_RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
DEFAULT_RUN = ROOT / "runs" / "incremental_v61_20260727"
GRAPHQL = "https://data.rcsb.org/graphql"


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def load_master(path: Path):
    full: dict[str, list[tuple[str, str]]] = defaultdict(list)
    connectivity: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["standard_inchikey"].strip().upper()
            conn = row["connectivity_key"].strip().upper()
            item = (row["compound_internal_id"], row["preferred_name"])
            if key:
                full[key].append(item)
            if conn:
                connectivity[conn].append((*item, key))
    return full, connectivity


def load_pdbe_review_annotations(review_path: Path) -> dict[str, dict]:
    result: dict[str, dict] = {}
    with gzip.open(review_path, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source_database"] != "PDBe":
                continue
            component = row["compound_source_id"].split(":", 1)[-1].upper()
            item = result.setdefault(
                component,
                {
                    "preferred_name": row["compound_name"],
                    "row_count": 0,
                    "targets": set(),
                    "artifact_rows": 0,
                },
            )
            item["row_count"] += 1
            if row["target_uniprot_id"]:
                item["targets"].add(row["target_uniprot_id"])
            reason = (
                row.get("original_exclusion_reason", "")
                + ";"
                + row.get("review_reason", "")
            ).casefold()
            if any(
                term in reason
                for term in (
                    "solvent",
                    "buffer",
                    "crystal artifact",
                    "non-small molecule",
                    "incomplete",
                )
            ):
                item["artifact_rows"] += 1
    return result


def query_batch(component_ids: list[str], timeout: int) -> tuple[int, dict]:
    quoted = ",".join(json.dumps(value) for value in component_ids)
    query = (
        "query { chem_comps(comp_ids: ["
        + quoted
        + """]) {
          chem_comp {
            id name formula formula_weight type pdbx_formal_charge
          }
          rcsb_chem_comp_descriptor {
            comp_id InChI InChIKey SMILES SMILES_stereo
          }
        } }"""
    )
    body = json.dumps({"query": query}).encode("utf-8")
    request = urllib.request.Request(
        GRAPHQL,
        data=body,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "MemPro-V6.1-CCD-resolution/1.0",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status, json.load(response)
    except urllib.error.HTTPError as exc:
        text = exc.read().decode("utf-8", errors="replace")
        try:
            payload = json.loads(text)
        except json.JSONDecodeError:
            payload = {"error": text[:1000]}
        return exc.code, payload


def initialize_db(conn: sqlite3.Connection, inventory: Path) -> int:
    conn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        CREATE TABLE IF NOT EXISTS component (
            component_id TEXT PRIMARY KEY,
            preferred_name TEXT,
            evidence_row_count INTEGER,
            target_count INTEGER,
            artifact_rows INTEGER,
            fetch_status TEXT NOT NULL DEFAULT 'pending',
            attempts INTEGER NOT NULL DEFAULT 0,
            ccd_name TEXT,
            formula TEXT,
            formula_weight REAL,
            component_type TEXT,
            formal_charge INTEGER,
            smiles TEXT,
            smiles_stereo TEXT,
            inchi TEXT,
            inchikey TEXT,
            http_status INTEGER,
            error_message TEXT,
            updated_utc TEXT
        );
        """
    )
    annotations = load_pdbe_review_annotations(
        inventory.parents[2]
        / ".."
        / ".."
        / "releases"
        / "release_mempro_v6_0_20260727"
        / "binding_evidence_review_queue_v6_0.tsv.gz"
    )
    added = 0
    with gzip.open(inventory, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            component = row["component_id"].upper()
            annotation = annotations.get(component, {})
            before = conn.total_changes
            conn.execute(
                """
                INSERT OR IGNORE INTO component(
                    component_id,preferred_name,evidence_row_count,target_count,artifact_rows
                ) VALUES (?,?,?,?,?)
                """,
                (
                    component,
                    row["preferred_name"],
                    int(row["row_count"]),
                    int(row["target_count"]),
                    int(annotation.get("artifact_rows", 0)),
                ),
            )
            added += conn.total_changes - before
    conn.commit()
    return added


def db_counts(conn: sqlite3.Connection) -> dict[str, int]:
    values = dict(
        conn.execute(
            "SELECT fetch_status,COUNT(*) FROM component GROUP BY fetch_status"
        ).fetchall()
    )
    values["total"] = sum(values.values())
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--requests-per-second", type=float, default=1.5)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--max-attempts", type=int, default=5)
    args = parser.parse_args()

    output_dir = args.run / "staging" / "pdbe_identity"
    qa_dir = args.run / "qa"
    output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir.mkdir(parents=True, exist_ok=True)
    inventory = output_dir / "pdbe_component_priority_inventory_v61.tsv.gz"
    db_path = output_dir / "rcsb_ccd_fetch_v61.sqlite"
    progress = qa_dir / "PDBE_CCD_IDENTITY_V61_PROGRESS.json"
    output = output_dir / "pdbe_ccd_identity_v61.tsv.gz"
    started = now()

    conn = sqlite3.connect(db_path)
    added = initialize_db(conn, inventory)
    write_json(
        progress,
        {
            "status": "running",
            "started_utc": started,
            "updated_utc": now(),
            "new_components_added": added,
            "counts": db_counts(conn),
            "endpoint": GRAPHQL,
        },
    )

    batches = 0
    min_interval = 1.0 / max(args.requests_per_second, 0.1)
    while True:
        component_ids = [
            row[0]
            for row in conn.execute(
                """
                SELECT component_id FROM component
                WHERE fetch_status IN ('pending','retry') AND attempts<?
                ORDER BY attempts,component_id LIMIT ?
                """,
                (args.max_attempts, args.batch_size),
            )
        ]
        if not component_ids:
            break
        started_request = time.monotonic()
        try:
            status, payload = query_batch(component_ids, args.timeout)
        except Exception as exc:
            status, payload = 0, {"error": f"{type(exc).__name__}:{exc}"}

        timestamp = now()
        if status == 200 and payload.get("data", {}).get("chem_comps") is not None:
            found: set[str] = set()
            for item in payload["data"]["chem_comps"]:
                if not item:
                    continue
                chem = item.get("chem_comp") or {}
                descriptor = item.get("rcsb_chem_comp_descriptor") or {}
                component = (
                    descriptor.get("comp_id") or chem.get("id") or ""
                ).upper()
                if not component:
                    continue
                found.add(component)
                conn.execute(
                    """
                    UPDATE component SET fetch_status='resolved',attempts=attempts+1,
                        ccd_name=?,formula=?,formula_weight=?,component_type=?,
                        formal_charge=?,smiles=?,smiles_stereo=?,inchi=?,inchikey=?,
                        http_status=200,error_message=NULL,updated_utc=?
                    WHERE component_id=?
                    """,
                    (
                        chem.get("name", ""),
                        chem.get("formula", ""),
                        chem.get("formula_weight"),
                        chem.get("type", ""),
                        chem.get("pdbx_formal_charge"),
                        descriptor.get("SMILES", ""),
                        descriptor.get("SMILES_stereo", ""),
                        descriptor.get("InChI", ""),
                        descriptor.get("InChIKey", ""),
                        timestamp,
                        component,
                    ),
                )
            missing = [value for value in component_ids if value not in found]
            conn.executemany(
                """
                UPDATE component SET fetch_status='not_found',attempts=attempts+1,
                    http_status=200,error_message='CCD component absent from response',
                    updated_utc=? WHERE component_id=?
                """,
                [(timestamp, value) for value in missing],
            )
        else:
            message = json.dumps(payload, ensure_ascii=False)[:1000]
            conn.executemany(
                """
                UPDATE component SET fetch_status='retry',attempts=attempts+1,
                    http_status=?,error_message=?,updated_utc=? WHERE component_id=?
                """,
                [
                    (status, message, timestamp, value)
                    for value in component_ids
                ],
            )
        conn.commit()
        batches += 1
        if batches % 5 == 0:
            write_json(
                progress,
                {
                    "status": "running",
                    "started_utc": started,
                    "updated_utc": now(),
                    "batches_completed": batches,
                    "counts": db_counts(conn),
                    "endpoint": GRAPHQL,
                },
            )
        if status in (0, 429, 500, 502, 503, 504):
            time.sleep(3)
        elapsed = time.monotonic() - started_request
        if elapsed < min_interval:
            time.sleep(min_interval - elapsed)

    conn.execute(
        """
        UPDATE component SET fetch_status='failed_after_retries',updated_utc=?
        WHERE fetch_status='retry' AND attempts>=?
        """,
        (now(), args.max_attempts),
    )
    conn.commit()

    full_map, connectivity_map = load_master(
        args.release / "small_molecule_master_v1_1.tsv"
    )
    header = [
        "component_id",
        "preferred_name",
        "ccd_name",
        "evidence_row_count",
        "target_count",
        "artifact_rows",
        "fetch_status",
        "attempts",
        "formula",
        "formula_weight",
        "component_type",
        "formal_charge",
        "smiles",
        "smiles_stereo",
        "inchi",
        "inchikey",
        "connectivity_key",
        "identity_resolution",
        "matched_compound_internal_ids",
        "matched_compound_names",
        "release_eligibility",
        "http_status",
        "error_message",
        "updated_utc",
    ]
    resolution_counts: dict[str, int] = defaultdict(int)
    eligibility_counts: dict[str, int] = defaultdict(int)
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        query = "SELECT * FROM component ORDER BY component_id"
        columns = [item[0] for item in conn.execute(query).description]
        for values in conn.execute(query):
            item = dict(zip(columns, values))
            key = (item["inchikey"] or "").upper()
            conn_key = key[:14] if key else ""
            exact = full_map.get(key, [])
            connected = connectivity_map.get(conn_key, [])
            if item["fetch_status"] != "resolved" or not key:
                resolution = "unresolved_ccd"
                matches = []
            elif exact:
                resolution = "exact_full_inchikey"
                matches = exact
            elif connected:
                resolution = "connectivity_only_form_or_stereo"
                matches = [(value[0], value[1]) for value in connected]
            else:
                resolution = "new_structure_candidate"
                matches = []
            component_type = (item["component_type"] or "").casefold()
            if item["artifact_rows"] and item["artifact_rows"] == item["evidence_row_count"]:
                eligibility = "exclude_curated_artifact"
            elif component_type and "non-polymer" not in component_type:
                eligibility = "exclude_or_review_non_small_molecule"
            elif not key:
                eligibility = "manual_review_no_identity"
            else:
                eligibility = "eligible_small_molecule"
            row = {
                **item,
                "connectivity_key": conn_key,
                "identity_resolution": resolution,
                "matched_compound_internal_ids": ";".join(
                    sorted({value[0] for value in matches})
                ),
                "matched_compound_names": ";".join(
                    sorted({value[1] for value in matches})
                ),
                "release_eligibility": eligibility,
            }
            writer.writerow({field: row.get(field, "") for field in header})
            resolution_counts[resolution] += 1
            eligibility_counts[eligibility] += 1

    final_counts = db_counts(conn)
    conn.close()
    summary = {
        "status": "complete",
        "started_utc": started,
        "completed_utc": now(),
        "batches_completed": batches,
        "fetch_counts": final_counts,
        "identity_resolution_counts": resolution_counts,
        "eligibility_counts": eligibility_counts,
        "output": str(output),
        "database": str(db_path),
        "endpoint": GRAPHQL,
    }
    write_json(progress, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
