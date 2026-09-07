#!/usr/bin/env python3
"""Build the V6.1 identity-resolution inventory without modifying V6.0."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import re
import sqlite3
from pathlib import Path


DEFAULT_RELEASE = Path(
    r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_0_20260727"
)
DEFAULT_RUN = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\incremental_v61_20260727"
)
PUBCHEM_RE = re.compile(r"(?:PUBCHEM:|CID[: ]?)(\d+)", re.IGNORECASE)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def open_tsv(path: Path):
    return gzip.open(path, "rt", encoding="utf-8-sig", newline="")


def extract_pubchem_cid(row: dict[str, str]) -> str:
    for field in (
        row.get("compound_source_id", ""),
        row.get("source_record_id", ""),
        row.get("source_evidence_id", ""),
    ):
        match = PUBCHEM_RE.search(field)
        if match:
            return match.group(1)
    return ""


def normalize_name(value: str) -> str:
    return " ".join(value.casefold().split())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--progress-every", type=int, default=100_000)
    args = parser.parse_args()

    release = args.release.resolve()
    run = args.run.resolve()
    dirs = {
        "manifests": run / "manifests",
        "staging": run / "staging",
        "pubchem": run / "staging" / "pubchem_identity",
        "pdbe": run / "staging" / "pdbe_identity",
        "pdbbind": run / "staging" / "pdbbind_identity",
        "brenda": run / "staging" / "brenda_identity",
        "qa": run / "qa",
        "logs": run / "logs",
        "merge": run / "merge_candidate",
    }
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)

    review_path = release / "binding_evidence_review_queue_v6_0.tsv.gz"
    if not review_path.exists():
        raise FileNotFoundError(review_path)

    started = dt.datetime.now(dt.timezone.utc).isoformat()
    progress_path = dirs["qa"] / "V61_IDENTITY_INVENTORY_PROGRESS.json"
    db_path = dirs["staging"] / "v61_identity_inventory.sqlite"
    if db_path.exists():
        db_path.unlink()

    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        PRAGMA temp_store=MEMORY;
        CREATE TABLE pubchem_cid (
            cid INTEGER PRIMARY KEY,
            row_count INTEGER NOT NULL DEFAULT 0,
            active_rows INTEGER NOT NULL DEFAULT 0,
            be2_rows INTEGER NOT NULL DEFAULT 0,
            active_be2_rows INTEGER NOT NULL DEFAULT 0,
            inconclusive_rows INTEGER NOT NULL DEFAULT 0
        );
        CREATE TABLE pubchem_pair (
            cid INTEGER NOT NULL,
            target_uniprot_id TEXT NOT NULL,
            PRIMARY KEY (cid, target_uniprot_id)
        ) WITHOUT ROWID;
        CREATE TABLE pdbe_component (
            component_id TEXT PRIMARY KEY,
            preferred_name TEXT,
            row_count INTEGER NOT NULL DEFAULT 0,
            target_count INTEGER NOT NULL DEFAULT 0
        );
        CREATE TABLE pdbe_pair (
            component_id TEXT NOT NULL,
            target_uniprot_id TEXT NOT NULL,
            PRIMARY KEY (component_id, target_uniprot_id)
        ) WITHOUT ROWID;
        CREATE TABLE pdbbind_complex (
            complex_id TEXT PRIMARY KEY,
            preferred_name TEXT,
            row_count INTEGER NOT NULL DEFAULT 0,
            target_count INTEGER NOT NULL DEFAULT 0
        );
        CREATE TABLE pdbbind_pair (
            complex_id TEXT NOT NULL,
            target_uniprot_id TEXT NOT NULL,
            PRIMARY KEY (complex_id, target_uniprot_id)
        ) WITHOUT ROWID;
        CREATE TABLE brenda_name (
            normalized_name TEXT PRIMARY KEY,
            preferred_name TEXT,
            row_count INTEGER NOT NULL DEFAULT 0,
            active_rows INTEGER NOT NULL DEFAULT 0,
            target_count INTEGER NOT NULL DEFAULT 0
        );
        CREATE TABLE brenda_pair (
            normalized_name TEXT NOT NULL,
            target_uniprot_id TEXT NOT NULL,
            PRIMARY KEY (normalized_name, target_uniprot_id)
        ) WITHOUT ROWID;
        """
    )

    source_counts: dict[str, int] = {}
    total = 0
    with open_tsv(review_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            source = row.get("source_database", "").strip()
            source_counts[source] = source_counts.get(source, 0) + 1
            source_lower = source.casefold()
            target = row.get("target_uniprot_id", "").strip()
            outcome = row.get("activity_outcome", "").strip().casefold()
            tier = row.get("evidence_tier", "").strip().upper()
            name = row.get("compound_name", "").strip()

            if source_lower.startswith("pubchem"):
                cid = extract_pubchem_cid(row)
                if cid:
                    is_active = int(outcome == "active")
                    is_be2 = int(tier == "BE2")
                    is_inconclusive = int(outcome == "inconclusive")
                    conn.execute(
                        """
                        INSERT INTO pubchem_cid
                            (cid,row_count,active_rows,be2_rows,active_be2_rows,inconclusive_rows)
                        VALUES (?,?,?,?,?,?)
                        ON CONFLICT(cid) DO UPDATE SET
                            row_count=row_count+1,
                            active_rows=active_rows+excluded.active_rows,
                            be2_rows=be2_rows+excluded.be2_rows,
                            active_be2_rows=active_be2_rows+excluded.active_be2_rows,
                            inconclusive_rows=inconclusive_rows+excluded.inconclusive_rows
                        """,
                        (
                            int(cid),
                            1,
                            is_active,
                            is_be2,
                            is_active * is_be2,
                            is_inconclusive,
                        ),
                    )
                    if target:
                        conn.execute(
                            "INSERT OR IGNORE INTO pubchem_pair VALUES (?,?)",
                            (int(cid), target),
                        )
            elif source_lower == "pdbe":
                component = row.get("compound_source_id", "").strip()
                if component.upper().startswith("PDBCCD:"):
                    component = component.split(":", 1)[1].upper()
                if component:
                    conn.execute(
                        """
                        INSERT INTO pdbe_component(component_id,preferred_name,row_count)
                        VALUES (?,?,1)
                        ON CONFLICT(component_id) DO UPDATE SET
                            row_count=row_count+1,
                            preferred_name=CASE WHEN preferred_name='' THEN excluded.preferred_name ELSE preferred_name END
                        """,
                        (component, name),
                    )
                    if target:
                        conn.execute(
                            "INSERT OR IGNORE INTO pdbe_pair VALUES (?,?)",
                            (component, target),
                        )
            elif source_lower == "pdbbind":
                complex_id = (
                    row.get("source_record_id", "").split(":", 1)[0].strip().lower()
                    or row.get("pdb_ids", "").split(";", 1)[0].strip().lower()
                )
                if complex_id:
                    conn.execute(
                        """
                        INSERT INTO pdbbind_complex(complex_id,preferred_name,row_count)
                        VALUES (?,?,1)
                        ON CONFLICT(complex_id) DO UPDATE SET
                            row_count=row_count+1,
                            preferred_name=CASE WHEN preferred_name='' THEN excluded.preferred_name ELSE preferred_name END
                        """,
                        (complex_id, name),
                    )
                    if target:
                        conn.execute(
                            "INSERT OR IGNORE INTO pdbbind_pair VALUES (?,?)",
                            (complex_id, target),
                        )
            elif source_lower == "brenda":
                normalized = normalize_name(name)
                if normalized:
                    is_active = int(outcome == "active")
                    conn.execute(
                        """
                        INSERT INTO brenda_name(normalized_name,preferred_name,row_count,active_rows)
                        VALUES (?,?,1,?)
                        ON CONFLICT(normalized_name) DO UPDATE SET
                            row_count=row_count+1,
                            active_rows=active_rows+excluded.active_rows,
                            preferred_name=CASE WHEN preferred_name='' THEN excluded.preferred_name ELSE preferred_name END
                        """,
                        (normalized, name, is_active),
                    )
                    if target:
                        conn.execute(
                            "INSERT OR IGNORE INTO brenda_pair VALUES (?,?)",
                            (normalized, target),
                        )

            if total % args.progress_every == 0:
                conn.commit()
                write_json(
                    progress_path,
                    {
                        "status": "running",
                        "started_utc": started,
                        "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                        "review_rows_scanned": total,
                        "source_rows": source_counts,
                    },
                )

    conn.commit()
    conn.execute(
        """
        UPDATE pdbe_component SET target_count=(
            SELECT COUNT(*) FROM pdbe_pair p WHERE p.component_id=pdbe_component.component_id
        )
        """
    )
    conn.execute(
        """
        UPDATE pdbbind_complex SET target_count=(
            SELECT COUNT(*) FROM pdbbind_pair p WHERE p.complex_id=pdbbind_complex.complex_id
        )
        """
    )
    conn.execute(
        """
        UPDATE brenda_name SET target_count=(
            SELECT COUNT(*) FROM brenda_pair p WHERE p.normalized_name=brenda_name.normalized_name
        )
        """
    )
    conn.commit()

    def export_query(path: Path, header: list[str], query: str) -> int:
        count = 0
        with gzip.open(path, "wt", encoding="utf-8", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for record in conn.execute(query):
                writer.writerow(record)
                count += 1
        return count

    pubchem_count = export_query(
        dirs["pubchem"] / "pubchem_cid_priority_inventory_v61.tsv.gz",
        [
            "priority",
            "cid",
            "row_count",
            "active_rows",
            "be2_rows",
            "active_be2_rows",
            "inconclusive_rows",
            "target_count",
        ],
        """
        SELECT
            CASE
                WHEN active_be2_rows>0 THEN 'P1_ACTIVE_BE2'
                WHEN active_rows>0 THEN 'P2_ACTIVE_OTHER'
                ELSE 'P3_NONACTIVE_OR_UNSPECIFIED'
            END AS priority,
            c.cid,c.row_count,c.active_rows,c.be2_rows,c.active_be2_rows,c.inconclusive_rows,
            (SELECT COUNT(*) FROM pubchem_pair p WHERE p.cid=c.cid) AS target_count
        FROM pubchem_cid c
        ORDER BY
            CASE
                WHEN active_be2_rows>0 THEN 1
                WHEN active_rows>0 THEN 2
                ELSE 3
            END,
            active_be2_rows DESC, active_rows DESC, row_count DESC, c.cid
        """,
    )
    pdbe_count = export_query(
        dirs["pdbe"] / "pdbe_component_priority_inventory_v61.tsv.gz",
        ["component_id", "preferred_name", "row_count", "target_count"],
        "SELECT component_id,preferred_name,row_count,target_count FROM pdbe_component ORDER BY row_count DESC,component_id",
    )
    pdbbind_count = export_query(
        dirs["pdbbind"] / "pdbbind_complex_priority_inventory_v61.tsv.gz",
        ["complex_id", "preferred_name", "row_count", "target_count"],
        "SELECT complex_id,preferred_name,row_count,target_count FROM pdbbind_complex ORDER BY row_count DESC,complex_id",
    )
    brenda_count = export_query(
        dirs["brenda"] / "brenda_name_priority_inventory_v61.tsv.gz",
        [
            "normalized_name",
            "preferred_name",
            "row_count",
            "active_rows",
            "target_count",
        ],
        "SELECT normalized_name,preferred_name,row_count,active_rows,target_count FROM brenda_name ORDER BY row_count DESC,normalized_name",
    )

    priority_counts = dict(
        conn.execute(
            """
            SELECT
                CASE
                    WHEN active_be2_rows>0 THEN 'P1_ACTIVE_BE2'
                    WHEN active_rows>0 THEN 'P2_ACTIVE_OTHER'
                    ELSE 'P3_NONACTIVE_OR_UNSPECIFIED'
                END,
                COUNT(*)
            FROM pubchem_cid GROUP BY 1
            """
        ).fetchall()
    )
    conn.close()

    frozen_files = []
    for path in sorted(release.iterdir()):
        if path.is_file():
            frozen_files.append(
                {
                    "file": path.name,
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    write_json(
        dirs["manifests"] / "V60_FROZEN_INPUT_MANIFEST.json",
        {
            "release_directory": str(release),
            "captured_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "files": frozen_files,
        },
    )

    summary = {
        "status": "complete",
        "started_utc": started,
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "review_rows_scanned": total,
        "source_rows": source_counts,
        "unique_objects": {
            "pubchem_cids": pubchem_count,
            "pdbe_components": pdbe_count,
            "pdbbind_complexes": pdbbind_count,
            "brenda_names": brenda_count,
        },
        "pubchem_priority_counts": priority_counts,
        "v60_files_hashed": len(frozen_files),
        "database": str(db_path),
    }
    write_json(progress_path, summary)
    write_json(dirs["manifests"] / "V61_INITIALIZATION_SUMMARY.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
