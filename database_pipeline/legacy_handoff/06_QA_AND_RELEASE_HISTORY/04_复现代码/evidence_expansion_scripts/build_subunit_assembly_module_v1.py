#!/usr/bin/env python3
"""Build a checkpointed membrane-protein subunit/assembly add-on.

The frozen MemPro release is read-only. This script writes an independent module
from UniProt SUBUNIT comments and PDBe biological assemblies. It intentionally
uses one network request at a time and sleeps between batches so it can run
beside the BRENDA identity resolver without competing aggressively for network
or CPU.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import re
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
MASTER = RELEASE / "human_membrane_protein_master_v6_0.tsv"
RUN = ROOT / "runs" / "subunit_assembly_v1_20260728"
STAGING = RUN / "staging"
QA = RUN / "qa"
LOGS = RUN / "logs"
DB_PATH = STAGING / "subunit_assembly_cache_v1.sqlite"
PROGRESS = QA / "SUBUNIT_ASSEMBLY_V1_PROGRESS.json"
SUMMARY_OUT = STAGING / "membrane_protein_subunit_summary_v1.tsv.gz"
ASSEMBLY_OUT = STAGING / "membrane_protein_biological_assembly_v1.tsv.gz"
UNIPROT_OUT = STAGING / "uniprot_subunit_annotations_v1.tsv.gz"

UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"
PDBE_ENDPOINT = "https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/assembly"
USER_AGENT = "MemPro-subunit-assembly/1.0"
UNIPROT_BATCH = 50
PDBE_BATCH = 50
NETWORK_SLEEP_SECONDS = 0.9
MAX_ATTEMPTS = 5


COUNT_WORDS = {
    "monomer": 1,
    "dimer": 2,
    "trimer": 3,
    "tetramer": 4,
    "pentamer": 5,
    "hexamer": 6,
    "heptamer": 7,
    "octamer": 8,
    "nonamer": 9,
    "decamer": 10,
    "undecamer": 11,
    "dodecamer": 12,
}
COUNT_NAMES = {value: key for key, value in COUNT_WORDS.items()}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def request(
    url: str,
    *,
    body: bytes | None = None,
    content_type: str | None = None,
    timeout: int = 120,
) -> tuple[int, bytes, str]:
    headers = {"Accept": "application/json, text/tab-separated-values", "User-Agent": USER_AGENT}
    if content_type:
        headers["Content-Type"] = content_type
    req = urllib.request.Request(
        url, data=body, method="POST" if body is not None else "GET", headers=headers
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            return response.status, response.read(), ""
    except urllib.error.HTTPError as exc:
        message = exc.read().decode("utf-8", errors="replace")
        return exc.code, b"", message[:2000]
    except Exception as exc:
        return 0, b"", f"{type(exc).__name__}:{exc}"


def split_pdb_ids(value: str) -> list[str]:
    return sorted(
        {
            token.lower()
            for token in re.split(r"[,;|\s]+", value.strip())
            if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", token)
        }
    )


def load_master() -> tuple[list[dict], dict[str, list[str]], list[str]]:
    rows: list[dict] = []
    protein_pdb: dict[str, list[str]] = {}
    all_pdb: set[str] = set()
    with MASTER.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["target_uniprot_id"].strip()
            pdb_ids = split_pdb_ids(row.get("pdb_ids", ""))
            rows.append(
                {
                    "membrane_protein_id": row["membrane_protein_id"],
                    "target_uniprot_id": accession,
                    "approved_symbol": row.get("approved_symbol", ""),
                    "protein_name": row.get("protein_name", ""),
                    "website_default_v52": row.get("website_default_v52", ""),
                    "membrane_class_v52": row.get("membrane_class_v52", ""),
                }
            )
            protein_pdb[accession] = pdb_ids
            all_pdb.update(pdb_ids)
    return rows, protein_pdb, sorted(all_pdb)


def initialize(
    conn: sqlite3.Connection, proteins: list[dict], pdb_ids: list[str]
) -> None:
    conn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        CREATE TABLE IF NOT EXISTS uniprot_subunit (
            accession TEXT PRIMARY KEY,
            fetch_status TEXT NOT NULL DEFAULT 'pending',
            subunit_comment TEXT,
            attempts INTEGER NOT NULL DEFAULT 0,
            http_status INTEGER,
            message TEXT,
            updated_utc TEXT
        );
        CREATE TABLE IF NOT EXISTS pdbe_assembly (
            pdb_id TEXT PRIMARY KEY,
            fetch_status TEXT NOT NULL DEFAULT 'pending',
            assembly_json TEXT,
            attempts INTEGER NOT NULL DEFAULT 0,
            http_status INTEGER,
            message TEXT,
            updated_utc TEXT
        );
        """
    )
    conn.executemany(
        "INSERT OR IGNORE INTO uniprot_subunit(accession) VALUES (?)",
        [(row["target_uniprot_id"],) for row in proteins],
    )
    conn.executemany(
        "INSERT OR IGNORE INTO pdbe_assembly(pdb_id) VALUES (?)",
        [(pdb_id,) for pdb_id in pdb_ids],
    )
    # A stopped run may leave transient network failures with consumed attempts.
    # Only pending transient failures are reset; resolved rows remain immutable.
    conn.execute(
        """
        UPDATE uniprot_subunit
        SET attempts=0,http_status=NULL,message=NULL
        WHERE fetch_status='pending'
          AND (http_status IS NULL OR http_status IN (0,429,500,502,503,504))
        """
    )
    conn.execute(
        """
        UPDATE pdbe_assembly
        SET attempts=0,http_status=NULL,message=NULL
        WHERE fetch_status='pending'
          AND (http_status IS NULL OR http_status IN (0,429,500,502,503,504))
        """
    )
    conn.commit()


def table_counts(conn: sqlite3.Connection, table: str) -> dict[str, int]:
    counts = dict(
        conn.execute(
            f"SELECT fetch_status,COUNT(*) FROM {table} GROUP BY fetch_status"
        ).fetchall()
    )
    counts["total"] = sum(counts.values())
    return counts


def progress_payload(
    conn: sqlite3.Connection, stage: str, started: str, extra: dict | None = None
) -> dict:
    payload = {
        "status": "running",
        "stage": stage,
        "started_utc": started,
        "updated_utc": now(),
        "network_mode": "single_thread_rate_limited",
        "network_sleep_seconds": NETWORK_SLEEP_SECONDS,
        "uniprot": table_counts(conn, "uniprot_subunit"),
        "pdbe": table_counts(conn, "pdbe_assembly"),
        "outputs": {
            "protein_summary": str(SUMMARY_OUT),
            "assembly_detail": str(ASSEMBLY_OUT),
            "uniprot_annotations": str(UNIPROT_OUT),
        },
    }
    if extra:
        payload.update(extra)
    return payload


def fetch_uniprot(conn: sqlite3.Connection, started: str) -> None:
    while True:
        accessions = [
            row[0]
            for row in conn.execute(
                """
                SELECT accession FROM uniprot_subunit
                WHERE fetch_status='pending' AND attempts<?
                ORDER BY attempts,accession LIMIT ?
                """,
                (MAX_ATTEMPTS, UNIPROT_BATCH),
            )
        ]
        if not accessions:
            break
        query = "(" + " OR ".join(f"accession:{item}" for item in accessions) + ")"
        params = urllib.parse.urlencode(
            {
                "query": query,
                "fields": "accession,cc_subunit",
                "format": "tsv",
                "size": str(UNIPROT_BATCH),
            }
        )
        status, raw, message = request(f"{UNIPROT_ENDPOINT}?{params}")
        if status == 200:
            text = raw.decode("utf-8-sig", errors="replace")
            parsed: dict[str, str] = {}
            reader = csv.DictReader(text.splitlines(), delimiter="\t")
            for row in reader:
                accession = (row.get("Entry") or "").strip()
                if accession:
                    parsed[accession] = (row.get("Subunit structure") or "").strip()
            for accession in accessions:
                comment = parsed.get(accession, "")
                conn.execute(
                    """
                    UPDATE uniprot_subunit SET fetch_status=?,subunit_comment=?,
                        attempts=attempts+1,http_status=200,message=NULL,updated_utc=?
                    WHERE accession=?
                    """,
                    (
                        "resolved_with_comment" if comment else "resolved_no_comment",
                        comment,
                        now(),
                        accession,
                    ),
                )
        elif status in (400, 404):
            conn.executemany(
                """
                UPDATE uniprot_subunit SET fetch_status='not_found',
                    attempts=attempts+1,http_status=?,message=?,updated_utc=?
                WHERE accession=?
                """,
                [(status, message, now(), accession) for accession in accessions],
            )
        else:
            conn.executemany(
                """
                UPDATE uniprot_subunit SET attempts=attempts+1,http_status=?,
                    message=?,updated_utc=? WHERE accession=?
                """,
                [(status, message, now(), accession) for accession in accessions],
            )
        conn.commit()
        write_json(PROGRESS, progress_payload(conn, "fetch_uniprot", started))
        time.sleep(3 if status in (0, 429, 500, 502, 503, 504) else NETWORK_SLEEP_SECONDS)
    conn.execute(
        """
        UPDATE uniprot_subunit SET fetch_status='failed_after_retries',updated_utc=?
        WHERE fetch_status='pending' AND attempts>=?
        """,
        (now(), MAX_ATTEMPTS),
    )
    conn.commit()


def fetch_pdbe(conn: sqlite3.Connection, started: str) -> None:
    while True:
        pdb_ids = [
            row[0]
            for row in conn.execute(
                """
                SELECT pdb_id FROM pdbe_assembly
                WHERE fetch_status='pending' AND attempts<?
                ORDER BY attempts,pdb_id LIMIT ?
                """,
                (MAX_ATTEMPTS, PDBE_BATCH),
            )
        ]
        if not pdb_ids:
            break
        body = json.dumps(",".join(pdb_ids)).encode("utf-8")
        status, raw, message = request(
            PDBE_ENDPOINT, body=body, content_type="application/json"
        )
        if status == 200:
            try:
                payload = json.loads(raw.decode("utf-8"))
            except Exception as exc:
                payload = {}
                status = 0
                message = f"invalid_json:{type(exc).__name__}:{exc}"
            if status == 200:
                for pdb_id in pdb_ids:
                    assemblies = payload.get(pdb_id, [])
                    conn.execute(
                        """
                        UPDATE pdbe_assembly SET fetch_status=?,assembly_json=?,
                            attempts=attempts+1,http_status=200,message=NULL,updated_utc=?
                        WHERE pdb_id=?
                        """,
                        (
                            "resolved_with_assembly"
                            if assemblies
                            else "resolved_no_assembly",
                            json.dumps(assemblies, ensure_ascii=False, separators=(",", ":")),
                            now(),
                            pdb_id,
                        ),
                    )
        if status in (400, 404, 422):
            conn.executemany(
                """
                UPDATE pdbe_assembly SET fetch_status='not_found',
                    attempts=attempts+1,http_status=?,message=?,updated_utc=?
                WHERE pdb_id=?
                """,
                [(status, message, now(), pdb_id) for pdb_id in pdb_ids],
            )
        elif status != 200:
            conn.executemany(
                """
                UPDATE pdbe_assembly SET attempts=attempts+1,http_status=?,
                    message=?,updated_utc=? WHERE pdb_id=?
                """,
                [(status, message, now(), pdb_id) for pdb_id in pdb_ids],
            )
        conn.commit()
        write_json(PROGRESS, progress_payload(conn, "fetch_pdbe", started))
        time.sleep(3 if status in (0, 429, 500, 502, 503, 504) else NETWORK_SLEEP_SECONDS)
    conn.execute(
        """
        UPDATE pdbe_assembly SET fetch_status='failed_after_retries',updated_utc=?
        WHERE fetch_status='pending' AND attempts>=?
        """,
        (now(), MAX_ATTEMPTS),
    )
    conn.commit()


def parse_uniprot_states(comment: str) -> tuple[list[str], list[int]]:
    text = comment.lower()
    states: set[str] = set()
    counts: set[int] = set()
    if re.search(r"\bhomo-\s*and/or\s*heterodimer", text):
        states.update({"homodimer", "heterodimer"})
        counts.add(2)
    if re.search(r"\bmonomer(?:ic)?\b", text):
        states.add("monomer")
        counts.add(1)
    for stem, count in COUNT_WORDS.items():
        if stem == "monomer":
            continue
        if re.search(rf"\bhomo[- ]?{re.escape(stem)}(?:ic)?\b", text):
            states.add("homo" + stem)
            counts.add(count)
        if re.search(rf"\bhetero[- ]?{re.escape(stem)}(?:ic)?\b", text):
            states.add("hetero" + stem)
            counts.add(count)
        if re.search(rf"(?<!homo)(?<!hetero)\b{re.escape(stem)}(?:ic)?\b", text):
            states.add(stem + "_unspecified")
            counts.add(count)
    if re.search(r"\bhomo[- ]?oligomer(?:ic)?\b", text):
        states.add("homooligomer_count_unknown")
    if re.search(r"\bhetero[- ]?oligomer(?:ic)?\b", text):
        states.add("heterooligomer_count_unknown")
    if not states and re.search(r"\boligomer(?:ic)?\b", text):
        states.add("oligomer_count_unknown")
    return sorted(states), sorted(counts)


def pdbe_state(assembly: dict) -> tuple[str, int, str, int]:
    proteins = [
        entity
        for entity in assembly.get("entities", [])
        if str(entity.get("molecule_type", "")).lower().startswith("polypeptide")
    ]
    total = sum(int(entity.get("number_of_copies") or 0) for entity in proteins)
    distinct = len(proteins)
    if total <= 0:
        state = "no_polypeptide_entity"
    elif total == 1:
        state = "monomer"
    else:
        stem = COUNT_NAMES.get(total, f"{total}-mer")
        prefix = "homo" if distinct == 1 else "hetero"
        state = prefix + stem
    composition = "; ".join(
        f"{' / '.join(entity.get('molecule_name') or ['unnamed'])} x{int(entity.get('number_of_copies') or 0)}"
        for entity in proteins
    )
    return state, total, composition, distinct


def bool_text(value: bool) -> str:
    return "1" if value else "0"


def export_uniprot(conn: sqlite3.Connection) -> int:
    columns = [
        "target_uniprot_id",
        "fetch_status",
        "subunit_states",
        "subunit_count_values",
        "direct_experimental_flag",
        "pubmed_ids",
        "eco_codes",
        "subunit_comment",
        "retrieval_utc",
    ]
    count = 0
    with gzip.open(UNIPROT_OUT, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for accession, status, comment, updated in conn.execute(
            """
            SELECT accession,fetch_status,COALESCE(subunit_comment,''),updated_utc
            FROM uniprot_subunit ORDER BY accession
            """
        ):
            states, counts = parse_uniprot_states(comment)
            eco = sorted(set(re.findall(r"ECO:\d{7}", comment)))
            pubmed = sorted(set(re.findall(r"PubMed:(\d+)", comment)), key=int)
            writer.writerow(
                {
                    "target_uniprot_id": accession,
                    "fetch_status": status,
                    "subunit_states": ";".join(states),
                    "subunit_count_values": ";".join(map(str, counts)),
                    "direct_experimental_flag": bool_text("ECO:0000269" in eco),
                    "pubmed_ids": ";".join(pubmed),
                    "eco_codes": ";".join(eco),
                    "subunit_comment": comment,
                    "retrieval_utc": updated or "",
                }
            )
            count += 1
    return count


def export_outputs(
    conn: sqlite3.Connection,
    proteins: list[dict],
    protein_pdb: dict[str, list[str]],
) -> tuple[int, int]:
    uniprot_cache = {
        accession: (status, comment or "")
        for accession, status, comment in conn.execute(
            "SELECT accession,fetch_status,subunit_comment FROM uniprot_subunit"
        )
    }
    pdbe_cache = {
        pdb_id: (status, json.loads(payload) if payload else [])
        for pdb_id, status, payload in conn.execute(
            "SELECT pdb_id,fetch_status,assembly_json FROM pdbe_assembly"
        )
    }
    detail_columns = [
        "target_uniprot_id",
        "pdb_id",
        "biological_assembly_id",
        "assembly_details",
        "assembly_assignment_class",
        "structure_experimental_flag",
        "assembly_direct_experimental_confirmation_flag",
        "oligomeric_state",
        "protein_subunit_count",
        "distinct_protein_entity_count",
        "subunit_composition",
        "assembly_composition",
        "polymeric_count",
        "molecular_weight_kda",
        "source_database",
        "source_endpoint",
    ]
    summary_columns = [
        "membrane_protein_id",
        "target_uniprot_id",
        "approved_symbol",
        "protein_name",
        "membrane_class_v52",
        "website_default_v52",
        "subunit_annotation_status",
        "subunit_state_summary",
        "subunit_count_values",
        "subunit_count_min",
        "subunit_count_max",
        "subunit_composition_summary",
        "biological_assembly_ids",
        "pdbe_pdb_count",
        "pdbe_assembly_count",
        "uniprot_subunit_comment",
        "uniprot_direct_experimental_flag",
        "uniprot_pubmed_ids",
        "pdbe_experimental_structure_flag",
        "pdbe_author_defined_assembly_flag",
        "evidence_sources",
        "subunit_state_conflict_flag",
        "subunit_state_conflict_reason",
        "module_version",
    ]

    details_by_accession: dict[str, list[dict]] = defaultdict(list)
    detail_count = 0
    with gzip.open(ASSEMBLY_OUT, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=detail_columns, lineterminator="\n")
        writer.writeheader()
        for protein in proteins:
            accession = protein["target_uniprot_id"]
            for pdb_id in protein_pdb.get(accession, []):
                fetch_status, assemblies = pdbe_cache.get(pdb_id, ("missing_cache", []))
                for assembly in assemblies:
                    state, count, composition, distinct = pdbe_state(assembly)
                    details = str(assembly.get("details") or "")
                    assignment = (
                        "author_defined"
                        if "author" in details.lower()
                        else "software_defined"
                        if any(word in details.lower() for word in ("software", "predicted"))
                        else "unspecified"
                    )
                    row = {
                        "target_uniprot_id": accession,
                        "pdb_id": pdb_id,
                        "biological_assembly_id": str(assembly.get("assembly_id") or ""),
                        "assembly_details": details,
                        "assembly_assignment_class": assignment,
                        "structure_experimental_flag": "1",
                        # A PDB biological assembly is based on an experimental
                        # structure, but its physiological stoichiometry is not
                        # automatically experimentally confirmed.
                        "assembly_direct_experimental_confirmation_flag": "0",
                        "oligomeric_state": state,
                        "protein_subunit_count": count or "",
                        "distinct_protein_entity_count": distinct,
                        "subunit_composition": composition,
                        "assembly_composition": str(assembly.get("assembly_composition") or ""),
                        "polymeric_count": str(assembly.get("polymeric_count") or ""),
                        "molecular_weight_kda": str(assembly.get("molecular_weight") or ""),
                        "source_database": "PDBe",
                        "source_endpoint": PDBE_ENDPOINT,
                    }
                    writer.writerow(row)
                    details_by_accession[accession].append(row)
                    detail_count += 1

    summary_count = 0
    with gzip.open(SUMMARY_OUT, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=summary_columns, lineterminator="\n")
        writer.writeheader()
        for protein in proteins:
            accession = protein["target_uniprot_id"]
            uni_status, comment = uniprot_cache.get(accession, ("missing_cache", ""))
            uni_states, uni_counts = parse_uniprot_states(comment)
            uni_pubmed = sorted(set(re.findall(r"PubMed:(\d+)", comment)), key=int)
            uni_direct = "ECO:0000269" in comment
            pdbe_rows = details_by_accession.get(accession, [])
            pdbe_states = sorted(
                {
                    row["oligomeric_state"]
                    for row in pdbe_rows
                    if row["oligomeric_state"] != "no_polypeptide_entity"
                }
            )
            pdbe_counts = sorted(
                {
                    int(row["protein_subunit_count"])
                    for row in pdbe_rows
                    if str(row["protein_subunit_count"]).isdigit()
                    and int(row["protein_subunit_count"]) > 0
                }
            )
            all_states = sorted(set(uni_states) | set(pdbe_states))
            all_counts = sorted(set(uni_counts) | set(pdbe_counts))
            author_rows = [
                row for row in pdbe_rows if row["assembly_assignment_class"] == "author_defined"
            ]
            author_counts = {
                int(row["protein_subunit_count"])
                for row in author_rows
                if str(row["protein_subunit_count"]).isdigit()
                and int(row["protein_subunit_count"]) > 0
            }
            conflict_reasons: list[str] = []
            if uni_counts and author_counts and not (set(uni_counts) & author_counts):
                conflict_reasons.append("UniProt_vs_PDBe_author_defined_count_disagreement")
            if len(set(uni_counts)) > 1:
                conflict_reasons.append("multiple_UniProt_annotated_states")
            if len(author_counts) > 1:
                conflict_reasons.append("multiple_PDBe_author_defined_states")
            if not conflict_reasons and len(all_counts) > 1:
                conflict_reasons.append("multiple_context_dependent_states")
            evidence_sources = []
            if comment:
                evidence_sources.append("UniProtKB")
            if pdbe_rows:
                evidence_sources.append("PDBe")
            if not comment and not pdbe_rows:
                annotation_status = "no_subunit_evidence_found"
            elif conflict_reasons:
                annotation_status = "multiple_or_conflicting_states"
            elif uni_direct:
                annotation_status = "direct_experimental_annotation"
            elif author_rows:
                annotation_status = "experimental_structure_author_assembly"
            else:
                annotation_status = "supported_annotation"
            compositions = sorted(
                {row["subunit_composition"] for row in pdbe_rows if row["subunit_composition"]}
            )
            assembly_ids = sorted(
                {
                    f"{row['pdb_id']}:{row['biological_assembly_id']}"
                    for row in pdbe_rows
                    if row["biological_assembly_id"]
                }
            )
            writer.writerow(
                {
                    "membrane_protein_id": protein["membrane_protein_id"],
                    "target_uniprot_id": accession,
                    "approved_symbol": protein["approved_symbol"],
                    "protein_name": protein["protein_name"],
                    "membrane_class_v52": protein["membrane_class_v52"],
                    "website_default_v52": protein["website_default_v52"],
                    "subunit_annotation_status": annotation_status,
                    "subunit_state_summary": ";".join(all_states),
                    "subunit_count_values": ";".join(map(str, all_counts)),
                    "subunit_count_min": min(all_counts) if all_counts else "",
                    "subunit_count_max": max(all_counts) if all_counts else "",
                    "subunit_composition_summary": " || ".join(compositions),
                    "biological_assembly_ids": ";".join(assembly_ids),
                    "pdbe_pdb_count": len(set(row["pdb_id"] for row in pdbe_rows)),
                    "pdbe_assembly_count": len(pdbe_rows),
                    "uniprot_subunit_comment": comment,
                    "uniprot_direct_experimental_flag": bool_text(uni_direct),
                    "uniprot_pubmed_ids": ";".join(uni_pubmed),
                    "pdbe_experimental_structure_flag": bool_text(bool(pdbe_rows)),
                    "pdbe_author_defined_assembly_flag": bool_text(bool(author_rows)),
                    "evidence_sources": ";".join(evidence_sources),
                    "subunit_state_conflict_flag": bool_text(bool(conflict_reasons)),
                    "subunit_state_conflict_reason": ";".join(conflict_reasons),
                    "module_version": "subunit_assembly_v1_20260728",
                }
            )
            summary_count += 1
    return summary_count, detail_count


def main() -> None:
    for path in (STAGING, QA, LOGS):
        path.mkdir(parents=True, exist_ok=True)
    started = now()
    proteins, protein_pdb, pdb_ids = load_master()
    conn = sqlite3.connect(DB_PATH)
    initialize(conn, proteins, pdb_ids)
    write_json(
        PROGRESS,
        progress_payload(
            conn,
            "initialized",
            started,
            {
                "protein_rows": len(proteins),
                "proteins_with_pdb": sum(bool(value) for value in protein_pdb.values()),
                "unique_pdb_ids": len(pdb_ids),
            },
        ),
    )
    fetch_uniprot(conn, started)
    fetch_pdbe(conn, started)
    write_json(PROGRESS, progress_payload(conn, "exporting", started))
    uniprot_rows = export_uniprot(conn)
    summary_rows, assembly_rows = export_outputs(conn, proteins, protein_pdb)
    final = progress_payload(
        conn,
        "complete",
        started,
        {
            "status": "complete",
            "completed_utc": now(),
            "protein_rows": len(proteins),
            "proteins_with_pdb": sum(bool(value) for value in protein_pdb.values()),
            "unique_pdb_ids": len(pdb_ids),
            "exported_uniprot_rows": uniprot_rows,
            "exported_protein_summary_rows": summary_rows,
            "exported_assembly_rows": assembly_rows,
        },
    )
    conn.close()
    write_json(PROGRESS, final)
    print(json.dumps(final, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
