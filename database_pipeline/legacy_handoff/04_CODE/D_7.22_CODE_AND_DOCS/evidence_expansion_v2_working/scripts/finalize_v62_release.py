#!/usr/bin/env python3
"""Integrate the approved V6.2 modules, validate, hash, and freeze the release."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from zoneinfo import ZoneInfo


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
V61 = ROOT / "releases" / "release_mempro_v6_1_20260729"
RUN = ROOT / "runs" / "v62_completion_20260729"
EXPR = RUN / "staging" / "expression_location_v2"
SUBUNIT = RUN / "staging" / "subunit_assembly_v2"
NEG = RUN / "staging" / "negative_identity"
QA = RUN / "qa"
CANDIDATE = RUN / "release_candidate_v6_2"
DATE = dt.datetime.now(ZoneInfo("Asia/Shanghai")).date().isoformat()
DATE_TAG = DATE.replace("-", "")
FINAL = ROOT / "releases" / f"release_mempro_v6_2_{DATE_TAG}"
DB = RUN / "staging" / "v62_release_build.sqlite"
PROGRESS = QA / "V62_FORMAL_BUILD_PROGRESS.json"
REPORT = QA / "V62_VALIDATION_REPORT.json"


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


def stable_negative_id(target: str, compound: str, aid: str, sid: str) -> str:
    natural = "\x1f".join([target, compound, aid, sid])
    return "NEG-" + hashlib.sha1(natural.encode()).hexdigest()[:24].upper()


def initialize_db(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        PRAGMA temp_store=MEMORY;
        CREATE TABLE IF NOT EXISTS seen_negative (
            dedup_key TEXT PRIMARY KEY
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS positive_pair (
            target TEXT NOT NULL,
            compound TEXT NOT NULL,
            PRIMARY KEY(target,compound)
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS negative_pair (
            target TEXT NOT NULL,
            compound TEXT NOT NULL,
            evidence_count INTEGER NOT NULL DEFAULT 0,
            conflict_flag INTEGER NOT NULL DEFAULT 0,
            PRIMARY KEY(target,compound)
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS negative_form_ref (
            form_id TEXT PRIMARY KEY
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS compound_key (
            compound_id TEXT PRIMARY KEY
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS form_key (
            form_id TEXT PRIMARY KEY
        ) WITHOUT ROWID;
        """
    )
    conn.commit()


def load_positive_pairs(conn: sqlite3.Connection) -> int:
    count = conn.execute("SELECT COUNT(*) FROM positive_pair").fetchone()[0]
    if count:
        return count
    source = V61 / "binding_evidence_master_v6_1.tsv.gz"
    batch = []
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("default_release_inclusion") != "1":
                continue
            target = row.get("target_uniprot_id", "")
            compound = row.get("compound_internal_id", "")
            if target and compound:
                batch.append((target, compound))
            if len(batch) >= 50_000:
                conn.executemany(
                    "INSERT OR IGNORE INTO positive_pair VALUES (?,?)", batch
                )
                conn.commit()
                batch.clear()
    if batch:
        conn.executemany(
            "INSERT OR IGNORE INTO positive_pair VALUES (?,?)", batch
        )
        conn.commit()
    return conn.execute("SELECT COUNT(*) FROM positive_pair").fetchone()[0]


def build_combined_negative(conn: sqlite3.Connection, started: str) -> dict:
    baseline = V61 / "negative_binding_evidence_v1_1.tsv.gz"
    delta = NEG / "negative_binding_evidence_newly_mapped_delta_v62.tsv.gz"
    output = CANDIDATE / "negative_binding_evidence_v1_2.tsv.gz"
    duplicate = CANDIDATE / "negative_evidence_duplicate_audit_v6_2.tsv.gz"
    with gzip.open(delta, "rt", encoding="utf-8", newline="") as handle:
        delta_reader = csv.DictReader(handle, delimiter="\t")
        fields = list(delta_reader.fieldnames or [])
    duplicate_fields = [
        "dedup_key",
        "source_layer",
        "target_uniprot_id",
        "compound_internal_id",
        "pubchem_aid",
        "pubchem_sid",
        "source_evidence_id",
        "release_version",
    ]
    conn.execute("DELETE FROM seen_negative")
    conn.execute("DELETE FROM negative_pair")
    conn.execute("DELETE FROM negative_form_ref")
    conn.commit()
    counts = Counter()
    form_cache = {}
    negative_identity_db = sqlite3.connect(
        NEG / "negative_cid_identity_v62.sqlite"
    )
    negative_identity_db.execute(
        "CREATE INDEX IF NOT EXISTS idx_v61_form_formid ON v61_form(form_id)"
    )
    with gzip.open(
        output, "wt", encoding="utf-8", newline=""
    ) as output_handle, gzip.open(
        duplicate, "wt", encoding="utf-8", newline=""
    ) as duplicate_handle:
        writer = csv.DictWriter(
            output_handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        dup_writer = csv.DictWriter(
            duplicate_handle,
            fieldnames=duplicate_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        dup_writer.writeheader()
        for source_layer, path in [
            ("v61_mapped_baseline", baseline),
            ("v62_newly_mapped_delta", delta),
        ]:
            with gzip.open(
                path, "rt", encoding="utf-8", newline=""
            ) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    counts["source_rows"] += 1
                    target = row.get("target_uniprot_id", "")
                    compound = row.get("compound_internal_id", "")
                    form = row.get("compound_form_id", "")
                    aid = row.get("pubchem_aid", "")
                    sid = row.get("pubchem_sid", "")
                    if not target or not compound:
                        counts["unexpected_unmapped_in_mapped_input"] += 1
                        continue
                    dedup = hashlib.sha1(
                        "\x1f".join([target, compound, aid, sid]).encode()
                    ).hexdigest().upper()
                    before = conn.total_changes
                    conn.execute(
                        "INSERT OR IGNORE INTO seen_negative VALUES (?)",
                        (dedup,),
                    )
                    if conn.total_changes == before:
                        counts["duplicates_removed"] += 1
                        dup_writer.writerow(
                            {
                                "dedup_key": dedup,
                                "source_layer": source_layer,
                                "target_uniprot_id": target,
                                "compound_internal_id": compound,
                                "pubchem_aid": aid,
                                "pubchem_sid": sid,
                                "source_evidence_id": row.get(
                                    "source_evidence_id", ""
                                ),
                                "release_version": "v6.2",
                            }
                        )
                        continue
                    if source_layer == "v62_newly_mapped_delta":
                        conflict = (
                            row.get(
                                "positive_negative_conflict_flag_v62", "0"
                            )
                            == "1"
                        )
                    else:
                        conflict = bool(
                            conn.execute(
                                """
                                SELECT 1 FROM positive_pair
                                WHERE target=? AND compound=?
                                """,
                                (target, compound),
                            ).fetchone()
                        )
                    if source_layer == "v61_mapped_baseline":
                        if form not in form_cache:
                            form_cache[form] = negative_identity_db.execute(
                                """
                                SELECT exact_key,parent_key FROM v61_form
                                WHERE form_id=? LIMIT 1
                                """,
                                (form,),
                            ).fetchone() or ("", "")
                        exact, parent = form_cache[form]
                        row.update(
                            {
                                "negative_evidence_id_v62": stable_negative_id(
                                    target, compound, aid, sid
                                ),
                                "exact_inchikey_v62": exact,
                                "canonical_parent_inchikey_v62": parent,
                                "canonical_relation_v62": "baseline_curated_form",
                                "identity_resolution_v62": "baseline_mapped_preserved",
                                "release_action_v62": "mapped_negative_release",
                                "release_version_v62": "v6.2",
                            }
                        )
                    row["positive_negative_conflict_flag_v62"] = int(conflict)
                    row["positive_negative_conflict_type_v62"] = (
                        "positive_and_negative_same_target_canonical_compound"
                        if conflict
                        else ""
                    )
                    writer.writerow({field: row.get(field, "") for field in fields})
                    counts["official_negative_rows"] += 1
                    if conflict:
                        counts["positive_negative_conflict_rows"] += 1
                    conn.execute(
                        """
                        INSERT INTO negative_pair
                        (target,compound,evidence_count,conflict_flag)
                        VALUES (?,?,1,?)
                        ON CONFLICT(target,compound) DO UPDATE SET
                            evidence_count=evidence_count+1,
                            conflict_flag=MAX(conflict_flag,excluded.conflict_flag)
                        """,
                        (target, compound, int(conflict)),
                    )
                    if form:
                        conn.execute(
                            "INSERT OR IGNORE INTO negative_form_ref VALUES (?)",
                            (form,),
                        )
                    if counts["source_rows"] % 100_000 == 0:
                        conn.commit()
                        write_json(
                            PROGRESS,
                            {
                                "status": "running",
                                "stage": "combine_negative_evidence",
                                "started_utc": started,
                                "updated_utc": now(),
                                "counts": dict(counts),
                            },
                        )
        conn.commit()
    negative_identity_db.close()
    return dict(counts)


def negative_compound_counts(conn: sqlite3.Connection) -> dict[str, tuple]:
    return {
        compound: (evidence, targets, conflicts)
        for compound, evidence, targets, conflicts in conn.execute(
            """
            SELECT compound,SUM(evidence_count),COUNT(*),SUM(conflict_flag)
            FROM negative_pair GROUP BY compound
            """
        )
    }


def build_compound_master(conn: sqlite3.Connection) -> tuple[int, int]:
    baseline = V61 / "small_molecule_master_v1_2.tsv"
    delta = NEG / "negative_new_compound_delta_v62.tsv.gz"
    output = CANDIDATE / "small_molecule_master_v1_3.tsv"
    counts = negative_compound_counts(conn)
    added_fields = [
        "negative_evidence_count_v13",
        "negative_target_count_v13",
        "positive_negative_conflict_pair_count_v13",
        "negative_identity_refresh_v13",
        "release_version_v13",
        "release_date_v13",
    ]
    conn.execute("DELETE FROM compound_key")
    baseline_count = 0
    new_count = 0
    with baseline.open(
        "rt", encoding="utf-8-sig", newline=""
    ) as source, output.open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(source, delimiter="\t")
        baseline_fields = list(reader.fieldnames or [])
        fields = baseline_fields + added_fields
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in reader:
            compound = row["compound_internal_id"]
            evidence, targets, conflicts = counts.get(compound, (0, 0, 0))
            row.update(
                {
                    "negative_evidence_count_v13": evidence,
                    "negative_target_count_v13": targets,
                    "positive_negative_conflict_pair_count_v13": conflicts,
                    "negative_identity_refresh_v13": "baseline_preserved",
                    "release_version_v13": "small-molecule-v1.3",
                    "release_date_v13": DATE,
                }
            )
            writer.writerow(row)
            conn.execute(
                "INSERT OR IGNORE INTO compound_key VALUES (?)", (compound,)
            )
            baseline_count += 1
        with gzip.open(
            delta, "rt", encoding="utf-8", newline=""
        ) as delta_handle:
            for item in csv.DictReader(delta_handle, delimiter="\t"):
                compound = item["compound_internal_id"]
                evidence, targets, conflicts = counts.get(compound, (0, 0, 0))
                row = {field: "" for field in fields}
                row.update(
                    {
                        "compound_internal_id": compound,
                        "preferred_name": f"PubChem CID {item['representative_pubchem_cid']}",
                        "preferred_name_source": "PubChem",
                        "compound_scope_status": "negative_evidence_only",
                        "compound_scope_class": "structure_resolved_small_molecule",
                        "identity_confidence": "high",
                        "standard_smiles": item["standard_smiles"],
                        "standard_inchi": item["standard_inchi"],
                        "standard_inchikey": item["parent_inchikey"],
                        "identity_key": "IK:" + item["parent_inchikey"],
                        "connectivity_key": item["parent_inchikey"][:14],
                        "molecular_formula": item["molecular_formula"],
                        "molecular_weight": item["molecular_weight"],
                        "exact_mass": item["exact_mass"],
                        "xlogp": item["xlogp"],
                        "tpsa": item["tpsa"],
                        "hbond_donor_count": item["hbond_donor_count"],
                        "hbond_acceptor_count": item["hbond_acceptor_count"],
                        "rotatable_bond_count": item["rotatable_bond_count"],
                        "formal_charge": item["formal_charge"],
                        "heavy_atom_count": item["heavy_atom_count"],
                        "ring_count": item["ring_count"],
                        "amide_bond_count": item["amide_bond_count"],
                        "computed_structural_class": "structure_resolved_small_molecule",
                        "molecule_types": "negative_assay_compound",
                        "form_count": 1,
                        "source_record_count": evidence,
                        "source_database_count": 1,
                        "source_databases": "PubChem BioAssay",
                        "pubchem_cids": item["representative_pubchem_cid"],
                        "record_qc_status": "ok",
                        "qc_notes": "new unique canonical parent from mapped negative evidence",
                        "protein_target_count": 0,
                        "binding_evidence_count": 0,
                        "release_version": "small-molecule-v1.3",
                        "release_date": DATE,
                        "standardization_policy": "V6.2 full-InChIKey; RDKit cleanup/fragment/charge parent",
                        "identity_refresh_v12": "not_in_v1.2",
                        "release_version_v12": "not_in_v1.2",
                        "negative_evidence_count_v13": evidence,
                        "negative_target_count_v13": targets,
                        "positive_negative_conflict_pair_count_v13": conflicts,
                        "negative_identity_refresh_v13": "new_unique_canonical_parent",
                        "release_version_v13": "small-molecule-v1.3",
                        "release_date_v13": DATE,
                    }
                )
                writer.writerow(row)
                conn.execute(
                    "INSERT OR IGNORE INTO compound_key VALUES (?)",
                    (compound,),
                )
                new_count += 1
    conn.commit()
    return baseline_count, new_count


def build_form_hierarchy(conn: sqlite3.Connection) -> tuple[int, int]:
    baseline = V61 / "compound_form_hierarchy_v1_2.tsv"
    delta = NEG / "negative_new_form_delta_v62.tsv.gz"
    output = CANDIDATE / "compound_form_hierarchy_v1_3.tsv"
    added = [
        "negative_identity_refresh_v13",
        "release_version_v13",
        "release_date_v13",
    ]
    conn.execute("DELETE FROM form_key")
    baseline_count = 0
    new_count = 0
    with baseline.open(
        "rt", encoding="utf-8-sig", newline=""
    ) as source, output.open(
        "wt", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or []) + added
        writer = csv.DictWriter(
            target, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in reader:
            row.update(
                {
                    "negative_identity_refresh_v13": "baseline_preserved",
                    "release_version_v13": "compound-form-v1.3",
                    "release_date_v13": DATE,
                }
            )
            writer.writerow(row)
            conn.execute(
                "INSERT OR IGNORE INTO form_key VALUES (?)",
                (row["compound_form_id"],),
            )
            baseline_count += 1
        with gzip.open(
            delta, "rt", encoding="utf-8", newline=""
        ) as delta_handle:
            for item in csv.DictReader(delta_handle, delimiter="\t"):
                row = {field: "" for field in fields}
                row.update(
                    {
                        "compound_form_id": item["compound_form_id"],
                        "compound_internal_id": item["compound_internal_id"],
                        "form_type": "parent_form",
                        "exact_smiles": item["exact_smiles"],
                        "exact_inchi": item["exact_inchi"],
                        "exact_inchikey": item["exact_inchikey"],
                        "exact_identity_key": "IK:" + item["exact_inchikey"],
                        "molecular_formula": item["molecular_formula"],
                        "molecular_weight": item["molecular_weight"],
                        "formal_charge": item["formal_charge"],
                        "source_record_count": 1,
                        "source_databases": "PubChem BioAssay",
                        "source_compound_ids": item[
                            "representative_pubchem_cid"
                        ],
                        "record_qc_status": "ok",
                        "release_version_v11": "not_in_v1.1",
                        "parent_inchikey_v12": item["parent_inchikey"],
                        "parent_relation_v12": item["parent_relation"],
                        "standardization_status_v12": "not_in_v1.2",
                        "release_version_v12": "not_in_v1.2",
                        "negative_identity_refresh_v13": "new_unique_form",
                        "release_version_v13": "compound-form-v1.3",
                        "release_date_v13": DATE,
                    }
                )
                writer.writerow(row)
                conn.execute(
                    "INSERT OR IGNORE INTO form_key VALUES (?)",
                    (item["compound_form_id"],),
                )
                new_count += 1
    conn.commit()
    return baseline_count, new_count


def load_gzip_index(path: Path, key: str) -> tuple[dict, list[str]]:
    rows = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        for row in reader:
            rows[row[key]] = row
    return rows, fields


def build_protein_master() -> int:
    expression, expression_fields = load_gzip_index(
        EXPR / "expression_location_summary_v2.tsv.gz",
        "target_uniprot_id",
    )
    subunit, subunit_fields = load_gzip_index(
        SUBUNIT / "protein_subunit_summary_v2.tsv.gz",
        "target_uniprot_id",
    )
    source = V61 / "human_membrane_protein_master_v6_1.tsv"
    output = CANDIDATE / "human_membrane_protein_master_v6_2.tsv"
    count = 0
    with source.open(
        "rt", encoding="utf-8-sig", newline=""
    ) as input_handle, output.open(
        "wt", encoding="utf-8", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        added_expression = expression_fields[1:]
        added_subunit = subunit_fields[4:]
        fields = list(reader.fieldnames or []) + added_expression + added_subunit
        writer = csv.DictWriter(
            output_handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            target = row["target_uniprot_id"]
            row.update(
                {field: expression[target][field] for field in added_expression}
            )
            row.update(
                {field: subunit[target][field] for field in added_subunit}
            )
            writer.writerow(row)
            count += 1
    return count


def rewrite_with_release(
    source: Path, output: Path, release_field: str = "release_version_v62"
) -> int:
    opener = gzip.open if source.suffix == ".gz" else open
    out_opener = gzip.open if output.suffix == ".gz" else open
    count = 0
    with opener(
        source, "rt", encoding="utf-8-sig", newline=""
    ) as input_handle, out_opener(
        output, "wt", encoding="utf-8", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if release_field not in fields:
            fields.append(release_field)
        writer = csv.DictWriter(
            output_handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            row[release_field] = "v6.2"
            writer.writerow(row)
            count += 1
    return count


def export_negative_pair_summary(conn: sqlite3.Connection) -> int:
    output = CANDIDATE / "protein_compound_negative_summary_v6_2.tsv.gz"
    fields = [
        "target_uniprot_id",
        "compound_internal_id",
        "negative_evidence_count",
        "positive_negative_conflict_flag",
        "release_version",
    ]
    count = 0
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for target, compound, evidence, conflict in conn.execute(
            """
            SELECT target,compound,evidence_count,conflict_flag
            FROM negative_pair ORDER BY target,compound
            """
        ):
            writer.writerow(
                {
                    "target_uniprot_id": target,
                    "compound_internal_id": compound,
                    "negative_evidence_count": evidence,
                    "positive_negative_conflict_flag": conflict,
                    "release_version": "v6.2",
                }
            )
            count += 1
    return count


def copy_modules() -> None:
    files = [
        EXPR / "protein_tissue_expression_v2.tsv.gz",
        EXPR / "protein_cell_type_expression_v2.tsv.gz",
        EXPR / "protein_cell_cluster_expression_v2.tsv.gz",
        EXPR / "protein_tissue_cell_ihc_expression_v2.tsv.gz",
        EXPR / "protein_subcellular_localization_v2.tsv.gz",
        EXPR / "anatomy_cell_location_ontology_mapping_v2.tsv",
        EXPR / "expression_location_summary_v2.tsv.gz",
        EXPR / "expression_location_review_queue_v2.tsv.gz",
        SUBUNIT / "uniprot_subunit_annotations_v2.tsv.gz",
        SUBUNIT / "protein_biological_assembly_v2.tsv.gz",
        SUBUNIT / "biological_assembly_components_v2.tsv.gz",
        SUBUNIT / "protein_subunit_summary_v2.tsv.gz",
        SUBUNIT / "subunit_conflict_review_v2.tsv.gz",
        NEG / "negative_cid_identity_crosswalk_v62.tsv.gz",
        NEG / "negative_evidence_duplicate_audit_v62.tsv.gz",
        QA / "V62_STRATIFIED_VALIDATION_SAMPLE.tsv",
        QA / "V62_PREFLIGHT_QA_GATE.json",
    ]
    for path in files:
        shutil.copy2(path, CANDIDATE / path.name)
    shutil.copy2(
        NEG / "negative_binding_evidence_unmapped_review_v1_2.tsv.gz",
        CANDIDATE / "negative_binding_evidence_unmapped_review_v1_2.tsv.gz",
    )


def validate_release(conn: sqlite3.Connection, stats: dict) -> dict:
    blocking = Counter()
    compound_rows = conn.execute("SELECT COUNT(*) FROM compound_key").fetchone()[0]
    form_rows = conn.execute("SELECT COUNT(*) FROM form_key").fetchone()[0]
    invalid_negative_compounds = conn.execute(
        """
        SELECT COUNT(*) FROM (
            SELECT DISTINCT compound FROM negative_pair n
            WHERE NOT EXISTS (
                SELECT 1 FROM compound_key c WHERE c.compound_id=n.compound
            )
        )
        """
    ).fetchone()[0]
    invalid_negative_forms = conn.execute(
        """
        SELECT COUNT(*) FROM negative_form_ref n
        WHERE NOT EXISTS (
            SELECT 1 FROM form_key f WHERE f.form_id=n.form_id
        )
        """
    ).fetchone()[0]
    if stats["protein_rows"] != 10_997:
        blocking["protein_row_count"] += 1
    if invalid_negative_compounds:
        blocking["negative_compound_fk"] += invalid_negative_compounds
    if invalid_negative_forms:
        blocking["negative_form_fk"] += invalid_negative_forms
    expression_report = json.loads(
        (QA / "EXPRESSION_LOCATION_V2_VALIDATION.json").read_text(
            encoding="utf-8"
        )
    )
    subunit_report = json.loads(
        (QA / "SUBUNIT_ASSEMBLY_V2_VALIDATION.json").read_text(
            encoding="utf-8"
        )
    )
    negative_report = json.loads(
        (QA / "NEGATIVE_EVIDENCE_RESOLUTION_V62_VALIDATION.json").read_text(
            encoding="utf-8"
        )
    )
    preflight_gate = json.loads(
        (QA / "V62_PREFLIGHT_QA_GATE.json").read_text(encoding="utf-8")
    )
    for name, report in [
        ("expression", expression_report),
        ("subunit", subunit_report),
        ("negative", negative_report),
    ]:
        if not str(report.get("status", "")).startswith("passed"):
            blocking[f"{name}_module_status"] += 1
    if preflight_gate.get("status") != "PASS":
        blocking["preflight_qa_gate"] += 1
    return {
        "status": "PASS" if not blocking else "FAIL",
        "blocking_errors": dict(blocking),
        "blocking_error_total": sum(blocking.values()),
        "protein_rows": stats["protein_rows"],
        "compound_rows": compound_rows,
        "form_rows": form_rows,
        "official_negative_rows": stats["negative"]["official_negative_rows"],
        "negative_pair_rows": stats["negative_pair_rows"],
        "invalid_negative_compound_foreign_keys": invalid_negative_compounds,
        "invalid_negative_form_foreign_keys": invalid_negative_forms,
        "expression_module_status": expression_report["status"],
        "subunit_module_status": subunit_report["status"],
        "negative_module_status": negative_report["status"],
        "preflight_qa_status": preflight_gate["status"],
        "generated_utc": now(),
    }


def write_release_docs(validation: dict, stats: dict) -> None:
    (CANDIDATE / "README_v6_2.md").write_text(
        """# MemPro V6.2

V6.2 is an additive release over frozen V6.1. It adds full selected HPA 25.1
normal-tissue, cell-type, cell-cluster, IHC/MS and subcellular-location
matrices; a conservative UniProt/PDBe subunit and biological-assembly module;
and structure-resolved canonical mapping of the former unmapped PubChem
negative-evidence layer.

Detailed expression and assembly records remain normalized add-on tables.
Only one-row-per-protein summaries are merged into the protein master.
""",
        encoding="utf-8",
    )
    (CANDIDATE / "METHODS_v6_2.md").write_text(
        """# Methods

- V6.1 was treated as immutable.
- HPA RNA, IHC, MS and localization modalities retain their original units;
  missing values and measured zero are distinct.
- Automatic chemical identity merging requires a unique full InChIKey.
- Unassigned stereochemistry and uncurated salt/charge-parent relationships
  remain in review.
- Negative evidence is deduplicated by target UniProt, canonical compound,
  PubChem AID and SID, then compared with positive evidence at target-compound
  level.
- PDB asymmetric-unit chain counts are not used as physiological stoichiometry.
  Heteromer target copy number remains unknown without component-level mapping.
""",
        encoding="utf-8",
    )
    release_info = {
        "release_version": "v6.2",
        "release_date": DATE,
        "baseline_release": "v6.1",
        "status": "frozen",
        "stats": stats,
        "validation": validation,
        "hpa_version": "25.1",
        "small_molecule_master_version": "v1.3",
        "compound_form_version": "v1.3",
    }
    write_json(CANDIDATE / "RELEASE_INFO_v6_2.json", release_info)
    write_json(CANDIDATE / "V62_VALIDATION_REPORT.json", validation)


def manifest_and_hashes() -> None:
    files = sorted(
        path
        for path in CANDIDATE.iterdir()
        if path.is_file()
        and path.name not in {"SHA256SUMS_v6_2.txt", "FREEZE_v6_2.json"}
    )
    with (CANDIDATE / "RELEASE_MANIFEST_v6_2.tsv").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file_name", "bytes", "sha256"])
        for path in files:
            writer.writerow([path.name, path.stat().st_size, sha256(path)])
    files = sorted(
        path
        for path in CANDIDATE.iterdir()
        if path.is_file() and path.name != "FREEZE_v6_2.json"
    )
    with (CANDIDATE / "SHA256SUMS_v6_2.txt").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        for path in files:
            handle.write(f"{sha256(path)}  {path.name}\n")


def main() -> None:
    if CANDIDATE.exists() or FINAL.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing V6.2 candidate/release: {CANDIDATE} / {FINAL}"
        )
    CANDIDATE.mkdir(parents=True)
    started = now()
    write_json(
        PROGRESS,
        {
            "status": "running",
            "stage": "initialize",
            "started_utc": started,
            "updated_utc": now(),
        },
    )
    conn = sqlite3.connect(DB)
    initialize_db(conn)
    positive_pairs = load_positive_pairs(conn)
    negative_stats = build_combined_negative(conn, started)
    baseline_compounds, new_compounds = build_compound_master(conn)
    baseline_forms, new_forms = build_form_hierarchy(conn)
    protein_rows = build_protein_master()
    positive_rows = rewrite_with_release(
        V61 / "binding_evidence_master_v6_1.tsv.gz",
        CANDIDATE / "binding_evidence_master_v6_2.tsv.gz",
    )
    disease_rows = rewrite_with_release(
        V61 / "protein_gene_disease_relations_v6_1.tsv",
        CANDIDATE / "protein_gene_disease_relations_v6_2.tsv",
    )
    pair_rows = rewrite_with_release(
        V61 / "protein_compound_summary_v6_1.tsv.gz",
        CANDIDATE / "protein_compound_summary_v6_2.tsv.gz",
    )
    site_rows = rewrite_with_release(
        V61 / "binding_site_instances_v6_1.tsv.gz",
        CANDIDATE / "binding_site_instances_v6_2.tsv.gz",
    )
    negative_pair_rows = export_negative_pair_summary(conn)
    copy_modules()
    stats = {
        "positive_pair_index_rows": positive_pairs,
        "negative": negative_stats,
        "baseline_compounds": baseline_compounds,
        "new_compounds": new_compounds,
        "baseline_forms": baseline_forms,
        "new_forms": new_forms,
        "protein_rows": protein_rows,
        "positive_evidence_rows": positive_rows,
        "disease_rows": disease_rows,
        "positive_pair_rows": pair_rows,
        "binding_site_rows": site_rows,
        "negative_pair_rows": negative_pair_rows,
    }
    validation = validate_release(conn, stats)
    write_json(REPORT, validation)
    if validation["status"] != "PASS":
        write_json(
            PROGRESS,
            {
                "status": "blocked",
                "stage": "validation_failed",
                "started_utc": started,
                "updated_utc": now(),
                "validation": validation,
                "candidate": str(CANDIDATE),
            },
        )
        conn.close()
        raise RuntimeError(f"V6.2 validation failed: {validation}")
    write_release_docs(validation, stats)
    manifest_and_hashes()
    freeze = {
        "release_version": "v6.2",
        "status": "frozen",
        "frozen_utc": now(),
        "baseline_release": str(V61),
        "validation_report_sha256": sha256(
            CANDIDATE / "V62_VALIDATION_REPORT.json"
        ),
        "sha256sums_sha256": sha256(CANDIDATE / "SHA256SUMS_v6_2.txt"),
    }
    write_json(CANDIDATE / "FREEZE_v6_2.json", freeze)
    CANDIDATE.rename(FINAL)
    write_json(
        PROGRESS,
        {
            "status": "complete",
            "stage": "frozen",
            "started_utc": started,
            "completed_utc": now(),
            "release": str(FINAL),
            "stats": stats,
            "validation": validation,
        },
    )
    conn.close()
    print(json.dumps(json.loads(PROGRESS.read_text(encoding="utf-8")), indent=2))


if __name__ == "__main__":
    main()
