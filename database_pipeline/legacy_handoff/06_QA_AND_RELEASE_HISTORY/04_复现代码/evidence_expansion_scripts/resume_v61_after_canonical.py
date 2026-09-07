#!/usr/bin/env python3
"""Resume the formal V6.1 build from the completed canonical SQLite checkpoint."""

from __future__ import annotations

import csv
import json
import shutil
import sqlite3
from collections import Counter

import finalize_v61_release as v


def load_slim_master() -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    path = v.CANDIDATE / "small_molecule_master_v1_2.tsv"
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["compound_internal_id"]] = {
                "preferred_name": row["preferred_name"],
                "compound_scope_status": row["compound_scope_status"],
                "standard_inchikey": row["standard_inchikey"],
            }
    return result


def rebuild_form_hierarchy(
    con: sqlite3.Connection,
    all_master: dict[str, dict[str, str]],
    baseline_form_ids: set[str],
) -> int:
    with v.FORM0.open("rt", encoding="utf-8-sig", newline="") as handle:
        form_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    for field in (
        "parent_inchikey_v12",
        "parent_relation_v12",
        "standardization_status_v12",
        "release_version_v12",
    ):
        if field not in form_fields:
            form_fields.append(field)

    con.executescript(
        """
        DROP TABLE IF EXISTS candidate_source_agg;
        CREATE TABLE candidate_source_agg AS
        SELECT exact_key,
               GROUP_CONCAT(DISTINCT source_database) AS sources,
               GROUP_CONCAT(source_database || ':' || source_id) AS source_ids,
               COUNT(*) AS source_record_count
        FROM candidate_source
        GROUP BY exact_key;
        CREATE UNIQUE INDEX idx_candidate_source_agg_key
        ON candidate_source_agg(exact_key);
        """
    )
    con.commit()
    path = v.CANDIDATE / "compound_form_hierarchy_v1_2.tsv"
    appended = 0
    with v.atomic_text(path, form_fields) as writer:
        with v.FORM0.open("rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row.update(
                    {
                        "parent_inchikey_v12": all_master[
                            row["compound_internal_id"]
                        ]["standard_inchikey"],
                        "parent_relation_v12": "baseline_form",
                        "standardization_status_v12": "baseline_preserved",
                        "release_version_v12": "small-molecule-v1.2",
                    }
                )
                writer.writerow({field: row.get(field, "") for field in form_fields})
        for row in con.execute(
            """
            SELECT exact_key,parent_key,compound_internal_id,compound_form_id,
                   parent_relation,exact_smiles,exact_inchi,exact_formula,
                   exact_mw,exact_charge,standardization_status,
                   standardization_notes,COALESCE(a.sources,''),
                   COALESCE(a.source_ids,''),COALESCE(a.source_record_count,0)
            FROM canonical_form c
            LEFT JOIN candidate_source_agg a USING(exact_key)
            ORDER BY compound_form_id
            """
        ):
            (
                exact_key,
                parent_key,
                compound_id,
                form_id,
                relation,
                smiles,
                inchi,
                formula,
                mw,
                charge,
                status,
                notes,
                source_text,
                source_id_text,
                source_record_count,
            ) = row
            if form_id in baseline_form_ids:
                continue
            writer.writerow(
                {
                    "compound_form_id": form_id,
                    "compound_internal_id": compound_id,
                    "form_type": relation,
                    "exact_smiles": smiles,
                    "exact_inchi": inchi,
                    "exact_inchikey": exact_key,
                    "exact_identity_key": f"IK:{exact_key}",
                    "molecular_formula": formula,
                    "molecular_weight": mw,
                    "formal_charge": charge,
                    "source_record_count": source_record_count,
                    "source_databases": ";".join(
                        sorted(v.split_values(source_text))
                    ),
                    "source_compound_ids": ";".join(
                        sorted(v.split_values(source_id_text))
                    ),
                    "record_qc_status": "ok",
                    "qc_notes": notes,
                    "release_version_v11": "small-molecule-v1.2",
                    "parent_inchikey_v12": parent_key,
                    "parent_relation_v12": relation,
                    "standardization_status_v12": status,
                    "release_version_v12": "small-molecule-v1.2",
                }
            )
            appended += 1
    return appended


def reconstruct_audit_counts(
    con: sqlite3.Connection,
    max_compound: int,
    max_form: int,
) -> tuple[dict, dict, dict]:
    prefixes = {
        "PubChem BioAssay": "PubChem",
        "PDBe": "PDBe",
        "PDBbind": "PDBbind",
        "BRENDA": "BRENDA",
    }
    identity_counts = {
        f"{prefixes.get(source, source)}:{status}": count
        for source, status, count in con.execute(
            """
            SELECT source_database,safe_status,COUNT(*)
            FROM identity GROUP BY source_database,safe_status
            """
        )
    }
    disposition_counts = {
        action: count
        for action, count in con.execute(
            "SELECT action,COUNT(*) FROM disposition GROUP BY action"
        )
    }
    type_to_source = {
        "PUBCHEM_CID": "PubChem BioAssay",
        "PDB_CCD": "PDBe",
        "PDBBIND_COMPLEX": "PDBbind",
        "BRENDA_NAME": "BRENDA",
    }
    source_disposition_counts = {
        f"{type_to_source.get(source_type, source_type)}:{action}": count
        for source_type, action, count in con.execute(
            """
            SELECT source_type,action,COUNT(*)
            FROM disposition GROUP BY source_type,action
            """
        )
    }
    disposition_unique_rows = sum(disposition_counts.values())
    review_rows_total = v.count_rows(v.REVIEW0)
    reconciliation = {
        "review_rows_total": review_rows_total,
        "disposition_counts": disposition_counts,
        "source_disposition_counts": source_disposition_counts,
        "disposition_unique_source_evidence_ids": disposition_unique_rows,
        "duplicate_source_evidence_id_rows": (
            review_rows_total - disposition_unique_rows
        ),
    }
    canonicalized = con.execute(
        "SELECT COUNT(*) FROM canonical_form"
    ).fetchone()[0]
    failed = con.execute(
        """
        SELECT COUNT(DISTINCT exact_key) FROM disposition
        WHERE reason LIKE 'canonicalization_failed:%'
        """
    ).fetchone()[0]
    new_parents = con.execute(
        """
        SELECT COUNT(DISTINCT compound_internal_id) FROM canonical_form
        WHERE CAST(SUBSTR(compound_internal_id,11) AS INTEGER)>?
        """,
        (max_compound,),
    ).fetchone()[0]
    new_forms = con.execute(
        """
        SELECT COUNT(*) FROM canonical_form
        WHERE CAST(SUBSTR(compound_form_id,11) AS INTEGER)>?
        """,
        (max_form,),
    ).fetchone()[0]
    canonical = {
        "referenced_new_exact_structures": canonicalized + failed,
        "canonicalized_exact_structures": canonicalized,
        "canonicalization_failures": {"total": failed},
        "new_parent_compounds": new_parents,
        "new_exact_forms": new_forms,
        "resume_checkpoint": "canonical_form_sqlite",
    }
    return identity_counts, reconciliation, canonical


def main() -> None:
    v.require_inputs()
    if not v.DB.exists():
        raise FileNotFoundError(f"Canonical checkpoint is missing: {v.DB}")
    if not (v.CANDIDATE / "small_molecule_master_v1_2.tsv").exists():
        raise FileNotFoundError("Completed V6.1 small-molecule master is missing")
    con = sqlite3.connect(v.DB)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA temp_store=FILE")
    master_by_id, _, form_by_key, max_compound, max_form = (
        v.load_baseline_identity()
    )

    v.update_progress("resume_load_master")
    all_master = load_slim_master()
    baseline_form_ids = set(form_by_key.values())
    v.update_progress(
        "resume_rebuild_form_hierarchy",
        master_rows=len(all_master),
        baseline_form_rows=len(baseline_form_ids),
    )
    appended_forms = rebuild_form_hierarchy(
        con, all_master, baseline_form_ids
    )

    identity_counts, reconciliation, canonical = reconstruct_audit_counts(
        con, max_compound, max_form
    )
    master_stats = {
        "master_rows": len(all_master),
        "new_master_rows": len(all_master) - len(master_by_id),
        "new_form_rows": appended_forms,
    }

    for path in v.CANDIDATE.iterdir():
        if path.is_file() and path.name not in {
            "small_molecule_master_v1_2.tsv",
            "compound_form_hierarchy_v1_2.tsv",
        }:
            path.unlink()
    con.executescript(
        """
        DELETE FROM pair_raw;
        DELETE FROM quant_raw;
        DELETE FROM site_new;
        """
    )
    con.commit()

    v.update_progress("rebuild_evidence_layers", master_stats=master_stats)
    evidence_stats = v.build_evidence_layers(con, all_master, form_by_key)
    v.update_progress("build_pair_site_protein_outputs")
    output_stats = v.build_pair_site_protein_outputs(con, all_master)
    v.update_progress("validate_referential_integrity")
    validation = v.duplicate_and_fk_checks(con, all_master, reconciliation)
    validation.update(
        {
            "status": (
                "PASS" if validation["blocking_error_total"] == 0 else "FAIL"
            ),
            "release_version": v.RELEASE_VERSION,
            "release_date": v.RELEASE_DATE,
            "generated_utc": v.now(),
            "identity_counts": identity_counts,
            "reconciliation": reconciliation,
            "canonicalization": canonical,
            "master_stats": master_stats,
            "evidence_stats": evidence_stats,
            "output_stats": output_stats,
            "resumed_from_checkpoint": True,
        }
    )
    v.write_json(v.CANDIDATE / "V61_VALIDATION_REPORT.json", validation)
    v.write_release_docs(
        identity_counts,
        reconciliation,
        canonical,
        master_stats,
        evidence_stats,
        output_stats,
        validation,
    )
    v.manifest_and_hashes()
    if validation["blocking_error_total"]:
        v.update_progress(
            "validation_failed",
            status="failed",
            validation_report=str(
                v.CANDIDATE / "V61_VALIDATION_REPORT.json"
            ),
            blocking_errors=validation["blocking_errors"],
        )
        raise RuntimeError(
            f"V6.1 not frozen; blocking errors: {validation['blocking_errors']}"
        )

    shutil.copytree(v.CANDIDATE, v.FINAL)
    freeze = {
        "status": "frozen",
        "release_version": v.RELEASE_VERSION,
        "release_path": str(v.FINAL),
        "frozen_utc": v.now(),
        "validation_report_sha256": v.sha256_file(
            v.FINAL / "V61_VALIDATION_REPORT.json"
        ),
        "release_manifest_sha256": v.sha256_file(
            v.FINAL / "RELEASE_MANIFEST_v6_1.tsv"
        ),
        "resumed_from_checkpoint": True,
    }
    v.write_json(v.FINAL / "FREEZE_v6_1.json", freeze)
    v.update_progress(
        "complete",
        status="complete",
        release_path=str(v.FINAL),
        validation_status="PASS",
        freeze_manifest=str(v.FINAL / "FREEZE_v6_1.json"),
    )
    con.close()
    print(json.dumps(freeze, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
