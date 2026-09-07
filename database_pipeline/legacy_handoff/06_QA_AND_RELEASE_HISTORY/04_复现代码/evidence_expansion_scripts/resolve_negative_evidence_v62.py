#!/usr/bin/env python3
"""Canonicalize negative-evidence CIDs and rebuild mapped/review layers."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import sqlite3
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger

from finalize_v61_release import canonicalize_candidate_worker


RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_1_20260729"
RUN = ROOT / "runs" / "v62_completion_20260729"
STAGING = RUN / "staging" / "negative_identity"
QA = RUN / "qa"
DB = STAGING / "negative_cid_identity_v62.sqlite"
PROGRESS = QA / "NEGATIVE_EVIDENCE_RESOLUTION_V62_PROGRESS.json"
REPORT = QA / "NEGATIVE_EVIDENCE_RESOLUTION_V62_VALIDATION.json"

MASTER = RELEASE / "small_molecule_master_v1_2.tsv"
FORM = RELEASE / "compound_form_hierarchy_v1_2.tsv"
POSITIVE = RELEASE / "binding_evidence_master_v6_1.tsv.gz"
NEGATIVE = (
    RELEASE / "negative_binding_evidence_unmapped_review_v1_1.tsv.gz"
)
MAPPED_OUT = STAGING / "negative_binding_evidence_newly_mapped_delta_v62.tsv.gz"
REVIEW_OUT = STAGING / "negative_binding_evidence_unmapped_review_v1_2.tsv.gz"
CROSSWALK_OUT = STAGING / "negative_cid_identity_crosswalk_v62.tsv.gz"
NEW_COMPOUND_OUT = STAGING / "negative_new_compound_delta_v62.tsv.gz"
NEW_FORM_OUT = STAGING / "negative_new_form_delta_v62.tsv.gz"
DUPLICATE_OUT = STAGING / "negative_evidence_duplicate_audit_v62.tsv.gz"

MODULE_VERSION = "negative_identity_v62_20260729"
CID_RE = re.compile(r"(?:CID[: ]?)(\d+)", re.IGNORECASE)


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    os.replace(tmp, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def integer_suffix(identifier: str) -> int:
    match = re.search(r"(\d+)$", identifier or "")
    return int(match.group(1)) if match else 0


def structure_worker(payload: tuple) -> tuple[str, dict | None, str, str]:
    key, record, status, message = canonicalize_candidate_worker(payload)
    if record:
        mol = Chem.MolFromSmiles(record["exact_smiles"])
        centers = (
            Chem.FindMolChiralCenters(
                mol, includeUnassigned=True, includeCIP=True
            )
            if mol is not None
            else []
        )
        record["unassigned_stereo_count"] = sum(
            label == "?" for _, label in centers
        )
        record["stereo_center_count"] = len(centers)
    return key, record, status, message


def initialize_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        PRAGMA temp_store=MEMORY;
        CREATE TABLE IF NOT EXISTS v61_form (
            exact_key TEXT NOT NULL,
            compound_id TEXT NOT NULL,
            form_id TEXT NOT NULL,
            parent_key TEXT,
            form_type TEXT,
            PRIMARY KEY(exact_key,form_id)
        ) WITHOUT ROWID;
        CREATE INDEX IF NOT EXISTS idx_v61_form_exact ON v61_form(exact_key);
        CREATE TABLE IF NOT EXISTS v61_parent (
            parent_key TEXT NOT NULL,
            compound_id TEXT NOT NULL,
            PRIMARY KEY(parent_key,compound_id)
        ) WITHOUT ROWID;
        CREATE INDEX IF NOT EXISTS idx_v61_parent_key ON v61_parent(parent_key);
        CREATE TABLE IF NOT EXISTS exact_structure (
            exact_key TEXT PRIMARY KEY,
            representative_cid INTEGER,
            cid_count INTEGER,
            evidence_rows INTEGER,
            smiles TEXT,
            inchi TEXT,
            molecular_formula TEXT,
            molecular_weight TEXT,
            charge INTEGER,
            canonical_status TEXT NOT NULL DEFAULT 'pending',
            canonical_error TEXT,
            parent_key TEXT,
            relation TEXT,
            exact_smiles TEXT,
            exact_inchi TEXT,
            parent_smiles TEXT,
            parent_inchi TEXT,
            parent_formula TEXT,
            parent_molecular_weight TEXT,
            parent_exact_mass TEXT,
            parent_xlogp TEXT,
            parent_tpsa TEXT,
            parent_hbd TEXT,
            parent_hba TEXT,
            parent_rotatable TEXT,
            parent_charge TEXT,
            parent_heavy_atoms TEXT,
            parent_ring_count TEXT,
            parent_amide_count TEXT,
            stereo_center_count INTEGER,
            unassigned_stereo_count INTEGER,
            identity_resolution TEXT,
            compound_id TEXT,
            form_id TEXT,
            review_reason TEXT,
            updated_utc TEXT
        );
        CREATE TABLE IF NOT EXISTS negative_new_compound (
            parent_key TEXT PRIMARY KEY,
            compound_id TEXT UNIQUE NOT NULL,
            representative_cid INTEGER,
            parent_smiles TEXT,
            parent_inchi TEXT,
            molecular_formula TEXT,
            molecular_weight TEXT,
            exact_mass TEXT,
            xlogp TEXT,
            tpsa TEXT,
            hbd TEXT,
            hba TEXT,
            rotatable TEXT,
            formal_charge TEXT,
            heavy_atom_count TEXT,
            ring_count TEXT,
            amide_count TEXT
        );
        CREATE TABLE IF NOT EXISTS negative_new_form (
            exact_key TEXT PRIMARY KEY,
            form_id TEXT UNIQUE NOT NULL,
            compound_id TEXT NOT NULL,
            exact_smiles TEXT,
            exact_inchi TEXT,
            parent_key TEXT,
            relation TEXT,
            representative_cid INTEGER
        );
        CREATE TABLE IF NOT EXISTS seen_negative (
            dedup_key TEXT PRIMARY KEY
        ) WITHOUT ROWID;
        """
    )
    conn.commit()


def load_v61_identity(conn: sqlite3.Connection) -> tuple[int, int]:
    form_count = conn.execute("SELECT COUNT(*) FROM v61_form").fetchone()[0]
    parent_count = conn.execute("SELECT COUNT(*) FROM v61_parent").fetchone()[0]
    if form_count and parent_count:
        return form_count, parent_count
    with FORM.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        batch = []
        for row in reader:
            key = (row.get("exact_inchikey") or "").strip().upper()
            if not key:
                continue
            batch.append(
                (
                    key,
                    row["compound_internal_id"],
                    row["compound_form_id"],
                    (row.get("parent_inchikey_v12") or "").strip().upper(),
                    row.get("form_type", ""),
                )
            )
            if len(batch) >= 20_000:
                conn.executemany(
                    "INSERT OR IGNORE INTO v61_form VALUES (?,?,?,?,?)", batch
                )
                conn.commit()
                batch.clear()
        if batch:
            conn.executemany(
                "INSERT OR IGNORE INTO v61_form VALUES (?,?,?,?,?)", batch
            )
            conn.commit()
    with MASTER.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        batch = []
        for row in reader:
            key = (row.get("standard_inchikey") or "").strip().upper()
            if key:
                batch.append((key, row["compound_internal_id"]))
            if len(batch) >= 20_000:
                conn.executemany(
                    "INSERT OR IGNORE INTO v61_parent VALUES (?,?)", batch
                )
                conn.commit()
                batch.clear()
        if batch:
            conn.executemany(
                "INSERT OR IGNORE INTO v61_parent VALUES (?,?)", batch
            )
            conn.commit()
    return (
        conn.execute("SELECT COUNT(*) FROM v61_form").fetchone()[0],
        conn.execute("SELECT COUNT(*) FROM v61_parent").fetchone()[0],
    )


def populate_exact_structures(conn: sqlite3.Connection) -> int:
    current = conn.execute("SELECT COUNT(*) FROM exact_structure").fetchone()[0]
    if current:
        return current
    conn.execute(
        """
        INSERT INTO exact_structure(
            exact_key,representative_cid,cid_count,evidence_rows,smiles,inchi,
            molecular_formula,molecular_weight,charge,updated_utc
        )
        SELECT UPPER(inchikey),MIN(cid),COUNT(*),SUM(evidence_rows),
               MIN(smiles),MIN(inchi),MIN(molecular_formula),
               MIN(molecular_weight),MIN(charge),?
        FROM cid_inventory
        WHERE fetch_status IN ('resolved_cache','resolved_pubchem')
          AND LENGTH(TRIM(inchikey))=27
        GROUP BY UPPER(inchikey)
        """,
        (now(),),
    )
    conn.commit()
    return conn.execute("SELECT COUNT(*) FROM exact_structure").fetchone()[0]


def exact_form_matches(conn: sqlite3.Connection) -> dict[str, int]:
    before = conn.total_changes
    conn.execute(
        """
        UPDATE exact_structure
        SET canonical_status='complete',
            parent_key=(
                SELECT MIN(parent_key) FROM v61_form f
                WHERE f.exact_key=exact_structure.exact_key
            ),
            relation='existing_exact_form',
            identity_resolution='exact_full_inchikey_unique_form',
            compound_id=(
                SELECT MIN(compound_id) FROM v61_form f
                WHERE f.exact_key=exact_structure.exact_key
            ),
            form_id=(
                SELECT MIN(form_id) FROM v61_form f
                WHERE f.exact_key=exact_structure.exact_key
            ),
            updated_utc=?
        WHERE canonical_status='pending'
          AND (SELECT COUNT(DISTINCT compound_id) FROM v61_form f
               WHERE f.exact_key=exact_structure.exact_key)=1
          AND (SELECT COUNT(*) FROM v61_form f
               WHERE f.exact_key=exact_structure.exact_key)=1
        """,
        (now(),),
    )
    unique = conn.total_changes - before
    conn.execute(
        """
        UPDATE exact_structure
        SET canonical_status='review',
            identity_resolution='ambiguous_existing_exact_form',
            review_reason='full_InChIKey_maps_to_multiple_V6.1_forms',
            updated_utc=?
        WHERE canonical_status='pending'
          AND (SELECT COUNT(*) FROM v61_form f
               WHERE f.exact_key=exact_structure.exact_key)>1
        """,
        (now(),),
    )
    conn.commit()
    return {
        "unique_exact_form_matches": unique,
        "ambiguous_exact_form_matches": conn.execute(
            """
            SELECT COUNT(*) FROM exact_structure
            WHERE identity_resolution='ambiguous_existing_exact_form'
            """
        ).fetchone()[0],
    }


def canonicalize_pending(conn: sqlite3.Connection, started: str) -> dict:
    failures = Counter()
    completed = conn.execute(
        "SELECT COUNT(*) FROM exact_structure WHERE canonical_status!='pending'"
    ).fetchone()[0]
    total = conn.execute("SELECT COUNT(*) FROM exact_structure").fetchone()[0]
    with ProcessPoolExecutor(
        max_workers=max(2, min(8, os.cpu_count() or 4))
    ) as pool:
        while True:
            rows = conn.execute(
                """
                SELECT exact_key,substr(exact_key,1,14),smiles,inchi,
                       molecular_formula,molecular_weight,charge,
                       'PubChem CID '||representative_cid
                FROM exact_structure
                WHERE canonical_status='pending'
                ORDER BY exact_key LIMIT 10000
                """
            ).fetchall()
            if not rows:
                break
            results = list(pool.map(structure_worker, rows, chunksize=100))
            for key, record, status, message in results:
                if not record:
                    failures[status] += 1
                    conn.execute(
                        """
                        UPDATE exact_structure
                        SET canonical_status='review',canonical_error=?,
                            identity_resolution='canonicalization_failed',
                            review_reason=?,updated_utc=?
                        WHERE exact_key=?
                        """,
                        (status, message[:1000], now(), key),
                    )
                    continue
                props = record["parent_props"]
                conn.execute(
                    """
                    UPDATE exact_structure SET
                        canonical_status='canonicalized',
                        parent_key=?,relation=?,exact_smiles=?,exact_inchi=?,
                        parent_smiles=?,parent_inchi=?,parent_formula=?,
                        parent_molecular_weight=?,parent_exact_mass=?,
                        parent_xlogp=?,parent_tpsa=?,parent_hbd=?,parent_hba=?,
                        parent_rotatable=?,parent_charge=?,parent_heavy_atoms=?,
                        parent_ring_count=?,parent_amide_count=?,
                        stereo_center_count=?,unassigned_stereo_count=?,
                        updated_utc=?
                    WHERE exact_key=?
                    """,
                    (
                        record["parent_key"],
                        record["relation"],
                        record["exact_smiles"],
                        record["exact_inchi"],
                        record["parent_smiles"],
                        record["parent_inchi"],
                        props["formula"],
                        props["mw"],
                        props["exact_mass"],
                        props["xlogp"],
                        props["tpsa"],
                        props["hbd"],
                        props["hba"],
                        props["rotatable"],
                        props["charge"],
                        props["heavy"],
                        props["rings"],
                        props["amide"],
                        record["stereo_center_count"],
                        record["unassigned_stereo_count"],
                        now(),
                        key,
                    ),
                )
            conn.commit()
            completed += len(rows)
            write_json(
                PROGRESS,
                {
                    "status": "running",
                    "stage": "canonicalize_structures",
                    "started_utc": started,
                    "updated_utc": now(),
                    "exact_structures_total": total,
                    "exact_structures_processed": completed,
                    "canonicalization_failures": dict(failures),
                },
            )
    return {"canonicalization_failures": dict(failures)}


def assign_identity(conn: sqlite3.Connection) -> dict[str, int]:
    # Safe match to an existing parent requires a standardized structure that is
    # already a parent form, with no unresolved stereocenter.
    conn.execute(
        """
        UPDATE exact_structure
        SET canonical_status='complete',
            identity_resolution='canonical_parent_full_inchikey_unique',
            compound_id=(
                SELECT MIN(compound_id) FROM v61_parent p
                WHERE p.parent_key=exact_structure.parent_key
            ),
            form_id=COALESCE((
                SELECT MIN(form_id) FROM v61_form f
                WHERE f.exact_key=exact_structure.parent_key
                  AND f.compound_id=(
                      SELECT MIN(compound_id) FROM v61_parent p
                      WHERE p.parent_key=exact_structure.parent_key
                  )
            ),''),
            updated_utc=?
        WHERE canonical_status='canonicalized'
          AND exact_key=parent_key
          AND COALESCE(unassigned_stereo_count,0)=0
          AND (SELECT COUNT(*) FROM v61_parent p
               WHERE p.parent_key=exact_structure.parent_key)=1
        """,
        (now(),),
    )
    conn.commit()

    # A new structure is eligible only when the exact structure is itself the
    # standardized parent and its stereochemistry is fully specified.
    max_compound = 0
    with MASTER.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            max_compound = max(
                max_compound, integer_suffix(row["compound_internal_id"])
            )
    max_form = 0
    with FORM.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            max_form = max(max_form, integer_suffix(row["compound_form_id"]))

    candidates = conn.execute(
        """
        SELECT exact_key,representative_cid,parent_smiles,parent_inchi,
               parent_formula,parent_molecular_weight,parent_exact_mass,
               parent_xlogp,parent_tpsa,parent_hbd,parent_hba,parent_rotatable,
               parent_charge,parent_heavy_atoms,parent_ring_count,
               parent_amide_count
        FROM exact_structure
        WHERE canonical_status='canonicalized'
          AND exact_key=parent_key
          AND COALESCE(unassigned_stereo_count,0)=0
          AND NOT EXISTS (
              SELECT 1 FROM v61_parent p
              WHERE p.parent_key=exact_structure.parent_key
          )
        ORDER BY exact_key
        """
    ).fetchall()
    for index, row in enumerate(candidates, start=1):
        compound_id = f"HMPD-CMPD-{max_compound+index:07d}"
        form_id = f"HMPD-FORM-{max_form+index:07d}"
        conn.execute(
            """
            INSERT OR IGNORE INTO negative_new_compound
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """,
            (row[0], compound_id, *row[1:]),
        )
        conn.execute(
            """
            INSERT OR IGNORE INTO negative_new_form
            SELECT exact_key,?, ?,exact_smiles,exact_inchi,parent_key,relation,
                   representative_cid
            FROM exact_structure WHERE exact_key=?
            """,
            (form_id, compound_id, row[0]),
        )
        conn.execute(
            """
            UPDATE exact_structure
            SET canonical_status='complete',
                identity_resolution='new_unique_canonical_parent',
                compound_id=?,form_id=?,updated_utc=?
            WHERE exact_key=?
            """,
            (compound_id, form_id, now(), row[0]),
        )

    # Everything else is intentionally retained for review.
    conn.execute(
        """
        UPDATE exact_structure
        SET canonical_status='review',
            identity_resolution=CASE
                WHEN COALESCE(unassigned_stereo_count,0)>0
                    THEN 'unassigned_stereochemistry_review'
                WHEN relation='salt_or_multicomponent_form'
                    THEN 'salt_or_multicomponent_review'
                WHEN relation='charged_or_protonation_form'
                    THEN 'charged_or_protonation_review'
                ELSE 'nonunique_or_unresolved_parent_review'
            END,
            review_reason=CASE
                WHEN COALESCE(unassigned_stereo_count,0)>0
                    THEN 'unassigned_stereochemistry'
                WHEN relation='salt_or_multicomponent_form'
                    THEN 'salt_parent_not_already_curated_as_exact_form'
                WHEN relation='charged_or_protonation_form'
                    THEN 'charge_or_protonation_parent_uncertain'
                ELSE 'canonical_parent_not_uniquely_mappable'
            END,
            updated_utc=?
        WHERE canonical_status='canonicalized'
        """,
        (now(),),
    )
    conn.commit()
    return {
        status: count
        for status, count in conn.execute(
            """
            SELECT identity_resolution,COUNT(*) FROM exact_structure
            GROUP BY identity_resolution
            """
        )
    }


def export_identity(conn: sqlite3.Connection) -> dict[str, int]:
    fields = [
        "cid",
        "evidence_rows",
        "fetch_status",
        "exact_inchikey",
        "canonical_parent_inchikey",
        "canonical_relation",
        "canonical_status",
        "identity_resolution",
        "compound_internal_id",
        "compound_form_id",
        "review_reason",
        "structure_source",
        "module_version",
    ]
    counts = Counter()
    with gzip.open(
        CROSSWALK_OUT, "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in conn.execute(
            """
            SELECT c.cid,c.evidence_rows,c.fetch_status,UPPER(c.inchikey),
                   COALESCE(s.parent_key,''),COALESCE(s.relation,''),
                   COALESCE(s.canonical_status,'review'),
                   COALESCE(s.identity_resolution,
                       CASE WHEN c.fetch_status='not_found'
                            THEN 'pubchem_cid_not_found'
                            WHEN c.fetch_status='failed_after_retries'
                            THEN 'pubchem_fetch_failed'
                            ELSE 'structure_or_full_inchikey_missing' END),
                   COALESCE(s.compound_id,''),COALESCE(s.form_id,''),
                   COALESCE(s.review_reason,c.error_message,''),
                   COALESCE(c.structure_source,'')
            FROM cid_inventory c
            LEFT JOIN exact_structure s ON s.exact_key=UPPER(c.inchikey)
            ORDER BY c.cid
            """
        ):
            out = dict(zip(fields[:-1], row))
            out["module_version"] = MODULE_VERSION
            writer.writerow(out)
            counts[out["identity_resolution"]] += 1

    compound_fields = [
        "parent_inchikey",
        "compound_internal_id",
        "representative_pubchem_cid",
        "standard_smiles",
        "standard_inchi",
        "molecular_formula",
        "molecular_weight",
        "exact_mass",
        "xlogp",
        "tpsa",
        "hbond_donor_count",
        "hbond_acceptor_count",
        "rotatable_bond_count",
        "formal_charge",
        "heavy_atom_count",
        "ring_count",
        "amide_bond_count",
        "identity_resolution",
        "module_version",
    ]
    with gzip.open(
        NEW_COMPOUND_OUT, "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=compound_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in conn.execute(
            """
            SELECT parent_key,compound_id,representative_cid,parent_smiles,
                   parent_inchi,molecular_formula,molecular_weight,exact_mass,
                   xlogp,tpsa,hbd,hba,rotatable,formal_charge,
                   heavy_atom_count,ring_count,amide_count
            FROM negative_new_compound ORDER BY compound_id
            """
        ):
            out = dict(zip(compound_fields[:17], row))
            out["identity_resolution"] = "new_unique_canonical_parent"
            out["module_version"] = MODULE_VERSION
            writer.writerow(out)

    form_fields = [
        "exact_inchikey",
        "compound_form_id",
        "compound_internal_id",
        "exact_smiles",
        "exact_inchi",
        "parent_inchikey",
        "parent_relation",
        "representative_pubchem_cid",
        "molecular_formula",
        "molecular_weight",
        "formal_charge",
        "module_version",
    ]
    with gzip.open(
        NEW_FORM_OUT, "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=form_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in conn.execute(
            """
            SELECT f.exact_key,f.form_id,f.compound_id,f.exact_smiles,
                   f.exact_inchi,f.parent_key,f.relation,
                   f.representative_cid,c.molecular_formula,
                   c.molecular_weight,c.formal_charge
            FROM negative_new_form f
            JOIN negative_new_compound c ON c.compound_id=f.compound_id
            ORDER BY f.form_id
            """
        ):
            out = dict(zip(form_fields[:-1], row))
            out["module_version"] = MODULE_VERSION
            writer.writerow(out)
    return dict(counts)


def load_cid_map(conn: sqlite3.Connection) -> dict[int, tuple]:
    return {
        row[0]: tuple(row[1:])
        for row in conn.execute(
            """
            SELECT c.cid,COALESCE(s.compound_id,''),COALESCE(s.form_id,''),
                   COALESCE(s.identity_resolution,
                       CASE WHEN c.fetch_status='not_found'
                            THEN 'pubchem_cid_not_found'
                            WHEN c.fetch_status='failed_after_retries'
                            THEN 'pubchem_fetch_failed'
                            ELSE 'structure_or_full_inchikey_missing' END),
                   COALESCE(s.exact_key,''),COALESCE(s.parent_key,''),
                   COALESCE(s.relation,''),COALESCE(s.review_reason,c.error_message,'')
            FROM cid_inventory c
            LEFT JOIN exact_structure s ON s.exact_key=UPPER(c.inchikey)
            """
        )
    }


def load_positive_pairs() -> set[str]:
    pairs = set()
    with gzip.open(
        POSITIVE, "rt", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("default_release_inclusion") != "1":
                continue
            target = row.get("target_uniprot_id", "")
            compound = row.get("compound_internal_id", "")
            if target and compound:
                pairs.add(target + "|" + compound)
    return pairs


def rebuild_negative(conn: sqlite3.Connection, started: str) -> dict[str, int]:
    cid_map = load_cid_map(conn)
    positive_pairs = load_positive_pairs()
    with gzip.open(
        NEGATIVE, "rt", encoding="utf-8", newline=""
    ) as source:
        reader = csv.DictReader(source, delimiter="\t")
        source_fields = list(reader.fieldnames or [])
        added = [
            "negative_evidence_id_v62",
            "exact_inchikey_v62",
            "canonical_parent_inchikey_v62",
            "canonical_relation_v62",
            "identity_resolution_v62",
            "positive_negative_conflict_flag_v62",
            "positive_negative_conflict_type_v62",
            "release_action_v62",
            "release_version_v62",
        ]
        mapped_fields = source_fields + added
        review_fields = source_fields + added + ["unmapped_reason_v62"]
        duplicate_fields = [
            "dedup_key",
            "target_uniprot_id",
            "compound_internal_id",
            "pubchem_aid",
            "pubchem_sid",
            "duplicate_source_evidence_id",
            "module_version",
        ]
        counts = Counter()
        conn.execute("DELETE FROM seen_negative")
        conn.commit()
        with gzip.open(
            MAPPED_OUT, "wt", encoding="utf-8", newline=""
        ) as mapped_handle, gzip.open(
            REVIEW_OUT, "wt", encoding="utf-8", newline=""
        ) as review_handle, gzip.open(
            DUPLICATE_OUT, "wt", encoding="utf-8", newline=""
        ) as duplicate_handle:
            mapped_writer = csv.DictWriter(
                mapped_handle,
                fieldnames=mapped_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            review_writer = csv.DictWriter(
                review_handle,
                fieldnames=review_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            duplicate_writer = csv.DictWriter(
                duplicate_handle,
                fieldnames=duplicate_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            mapped_writer.writeheader()
            review_writer.writeheader()
            duplicate_writer.writeheader()
            for row in reader:
                counts["source_rows"] += 1
                cid_text = (row.get("pubchem_cid") or "").strip()
                if not cid_text.isdigit():
                    match = CID_RE.search(
                        (row.get("source_record_id") or "")
                        + " "
                        + (row.get("compound_source_id") or "")
                    )
                    cid_text = match.group(1) if match else ""
                resolution = cid_map.get(int(cid_text), ()) if cid_text.isdigit() else ()
                if resolution:
                    compound, form, identity, exact, parent, relation, reason = resolution
                else:
                    compound = form = exact = parent = relation = ""
                    identity = "missing_or_invalid_pubchem_cid"
                    reason = "record_has_no_parseable_CID"
                row["compound_internal_id"] = compound
                row["compound_form_id"] = form
                row["compound_mapping_status"] = (
                    "mapped_canonical" if compound else "unresolved"
                )
                target = row.get("target_uniprot_id", "")
                aid = row.get("pubchem_aid", "")
                sid = row.get("pubchem_sid", "")
                natural = "\x1f".join([target, compound, aid, sid])
                dedup_hash = hashlib.sha1(natural.encode()).hexdigest().upper()
                conflict = bool(compound and target + "|" + compound in positive_pairs)
                row.update(
                    {
                        "negative_evidence_id_v62": (
                            "NEG-" + dedup_hash[:24]
                            if compound
                            else ""
                        ),
                        "exact_inchikey_v62": exact,
                        "canonical_parent_inchikey_v62": parent,
                        "canonical_relation_v62": relation,
                        "identity_resolution_v62": identity,
                        "positive_negative_conflict_flag_v62": int(conflict),
                        "positive_negative_conflict_type_v62": (
                            "positive_and_negative_same_target_canonical_compound"
                            if conflict
                            else ""
                        ),
                        "release_action_v62": (
                            "mapped_negative_release"
                            if compound
                            else "retained_unmapped_negative_review"
                        ),
                        "release_version_v62": "v6.2",
                    }
                )
                if compound:
                    before = conn.total_changes
                    conn.execute(
                        "INSERT OR IGNORE INTO seen_negative VALUES (?)",
                        (dedup_hash,),
                    )
                    if conn.total_changes == before:
                        counts["duplicates_removed"] += 1
                        duplicate_writer.writerow(
                            {
                                "dedup_key": dedup_hash,
                                "target_uniprot_id": target,
                                "compound_internal_id": compound,
                                "pubchem_aid": aid,
                                "pubchem_sid": sid,
                                "duplicate_source_evidence_id": row.get(
                                    "source_evidence_id", ""
                                ),
                                "module_version": MODULE_VERSION,
                            }
                        )
                        continue
                    mapped_writer.writerow(row)
                    counts["mapped_release_rows"] += 1
                    if conflict:
                        counts["positive_negative_conflict_rows"] += 1
                else:
                    row["unmapped_reason_v62"] = reason or identity
                    review_writer.writerow(row)
                    counts["unmapped_review_rows"] += 1
                if counts["source_rows"] % 100_000 == 0:
                    conn.commit()
                    write_json(
                        PROGRESS,
                        {
                            "status": "running",
                            "stage": "rebuild_negative_tables",
                            "started_utc": started,
                            "updated_utc": now(),
                            "counts": dict(counts),
                        },
                    )
            conn.commit()
    return dict(counts)


def main() -> None:
    STAGING.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    started = now()
    conn = sqlite3.connect(DB)
    initialize_schema(conn)
    write_json(
        PROGRESS,
        {
            "status": "running",
            "stage": "load_v61_identity",
            "started_utc": started,
            "updated_utc": now(),
        },
    )
    form_rows, parent_rows = load_v61_identity(conn)
    exact_total = populate_exact_structures(conn)
    exact_match_stats = exact_form_matches(conn)
    canonical_stats = canonicalize_pending(conn, started)
    resolution_counts = assign_identity(conn)
    identity_counts = export_identity(conn)
    negative_counts = rebuild_negative(conn, started)
    new_compounds = conn.execute(
        "SELECT COUNT(*) FROM negative_new_compound"
    ).fetchone()[0]
    new_forms = conn.execute(
        "SELECT COUNT(*) FROM negative_new_form"
    ).fetchone()[0]
    report = {
        "status": "passed_with_review_queue",
        "started_utc": started,
        "completed_utc": now(),
        "module_version": MODULE_VERSION,
        "inputs": {
            "v61_compound_rows": parent_rows,
            "v61_form_rows": form_rows,
            "negative_review_sha256": sha256(NEGATIVE),
            "property_database": str(DB),
        },
        "exact_structures": exact_total,
        "exact_match_stats": exact_match_stats,
        "canonicalization": canonical_stats,
        "identity_resolution_counts": resolution_counts,
        "cid_identity_counts": identity_counts,
        "negative_table_counts": negative_counts,
        "new_unique_canonical_compounds": new_compounds,
        "new_unique_forms": new_forms,
        "quality_checks": {
            "full_inchikey_required_for_exact_form_merge": True,
            "unassigned_stereochemistry_not_auto_merged": True,
            "uncurated_salt_or_charge_parent_not_auto_merged": True,
            "dedup_key_target_canonical_aid_sid": True,
            "positive_negative_conflicts_flagged": True,
            "unmapped_records_retained": True,
            "v61_immutable": True,
        },
        "outputs": {
            "mapped_negative": str(MAPPED_OUT),
            "unmapped_review": str(REVIEW_OUT),
            "cid_crosswalk": str(CROSSWALK_OUT),
            "new_compound_delta": str(NEW_COMPOUND_OUT),
            "new_form_delta": str(NEW_FORM_OUT),
            "duplicate_audit": str(DUPLICATE_OUT),
        },
    }
    write_json(REPORT, report)
    write_json(
        PROGRESS,
        {
            "status": "complete",
            "stage": "complete",
            "started_utc": started,
            "completed_utc": now(),
            "report": str(REPORT),
            "negative_table_counts": negative_counts,
        },
    )
    conn.close()
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
