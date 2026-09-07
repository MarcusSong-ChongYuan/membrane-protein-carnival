#!/usr/bin/env python3
"""Create experiment-lineage keys and database-count semantics for V6.3.

This is an additive audit layer over frozen V6.2.  It intentionally keeps
database contribution count separate from structure, PubChem, literature and
evidence-modality counts; none of these are treated as interchangeable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd


def arguments():
    p = argparse.ArgumentParser()
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--sites", type=Path, required=True)
    p.add_argument("--pair-summary", type=Path, required=True)
    p.add_argument("--docking-full", type=Path, required=True)
    p.add_argument("--docking-conflict", type=Path, required=True)
    p.add_argument("--protein-classification", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--qa-dir", type=Path, required=True)
    p.add_argument("--work-db", type=Path, required=True)
    return p.parse_args()


def h(path: Path) -> str:
    out = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            out.update(block)
    return out.hexdigest()


def q(path: Path) -> str:
    return path.as_posix().replace("'", "''")


def main():
    a = arguments(); a.output_dir.mkdir(parents=True, exist_ok=True); a.qa_dir.mkdir(parents=True, exist_ok=True); a.work_db.parent.mkdir(parents=True, exist_ok=True)
    if a.work_db.exists():
        a.work_db.unlink()
    con = duckdb.connect(str(a.work_db))
    con.execute("PRAGMA threads=4")
    con.execute("PRAGMA memory_limit='8GB'")
    con.execute(f"CREATE VIEW evidence_raw AS SELECT * FROM read_csv('{q(a.evidence)}', delim='\\t', header=true, all_varchar=true, compression='gzip')")
    con.execute(f"CREATE VIEW sites_raw AS SELECT * FROM read_csv('{q(a.sites)}', delim='\\t', header=true, all_varchar=true, compression='gzip')")
    con.execute(f"CREATE VIEW pair_raw AS SELECT * FROM read_csv('{q(a.pair_summary)}', delim='\\t', header=true, all_varchar=true, compression='gzip')")
    con.execute(f"CREATE VIEW docking_full_raw AS SELECT * FROM read_csv('{q(a.docking_full)}', delim='\\t', header=true, all_varchar=true, compression='gzip')")
    con.execute(f"CREATE VIEW docking_conflict_raw AS SELECT * FROM read_csv('{q(a.docking_conflict)}', delim='\\t', header=true, all_varchar=true, compression='gzip')")
    con.execute(f"CREATE VIEW protein_class_raw AS SELECT * FROM read_csv('{q(a.protein_classification)}', delim='\\t', header=true, all_varchar=true)")

    con.execute(r"""
    CREATE TABLE evidence_lineage AS
    SELECT
      evidence_id, target_uniprot_id, compound_internal_id, compound_form_id,
      source_database, source_record_id, relationship_context_id,
      evidence_tier, evidence_type, evidence_directness, activity_type,
      activity_relation, activity_value, activity_unit, standard_value_nM,
      assay_or_mechanism, pdb_ids, ligand_het_id, pubmed_ids, doi,
      conflict_status_v60, record_qc_status, default_release_inclusion,
      'DBR-' || upper(substr(md5(concat_ws('|',
        upper(trim(coalesce(source_database,''))), trim(coalesce(source_record_id,'')),
        trim(coalesce(relationship_context_id,'')), trim(coalesce(target_uniprot_id,'')),
        trim(coalesce(compound_internal_id,'')))),1,20)) AS database_record_key,
      CASE WHEN trim(coalesce(pdb_ids,'')) <> '' THEN
        'STR-' || upper(substr(md5(concat_ws('|', upper(trim(target_uniprot_id)),
        upper(regexp_replace(trim(pdb_ids),'[;, ]+',';', 'g')),
        upper(trim(coalesce(ligand_het_id,''))), trim(coalesce(compound_internal_id,'')))),1,20))
      ELSE NULL END AS structure_experiment_key,
      CASE WHEN lower(coalesce(source_database,'')) LIKE '%pubchem%' THEN
        'PCA-' || upper(substr(md5(concat_ws('|', trim(coalesce(source_record_id,'')),
        trim(coalesce(relationship_context_id,'')), trim(coalesce(compound_id,'')),
        trim(coalesce(target_uniprot_id,'')), trim(coalesce(assay_or_mechanism,'')))),1,20))
      ELSE NULL END AS pubchem_experiment_key,
      CASE WHEN trim(coalesce(pubmed_ids,'')) <> '' OR trim(coalesce(doi,'')) <> '' THEN
        'LIT-' || upper(substr(md5(concat_ws('|', lower(trim(coalesce(pubmed_ids,''))),
        lower(trim(coalesce(doi,''))), trim(coalesce(target_uniprot_id,'')),
        trim(coalesce(compound_internal_id,'')), lower(trim(coalesce(activity_type,''))),
        trim(coalesce(standard_value_nM,'')), lower(trim(coalesce(assay_or_mechanism,''))))),1,20))
      ELSE NULL END AS literature_experiment_key,
      CASE
        WHEN trim(coalesce(pdb_ids,'')) <> '' OR lower(coalesce(source_database,'')) IN ('pdbe','pdbbind','biolip','sc-pdb') THEN 'experimental_structure'
        WHEN lower(coalesce(source_database,'')) LIKE '%pubchem%' THEN 'bioassay'
        WHEN lower(coalesce(source_database,''))='brenda' THEN 'enzyme_kinetics'
        WHEN lower(coalesce(source_database,'')) IN ('chembl','bindingdb') THEN 'quantitative_bioactivity'
        WHEN lower(coalesce(evidence_type,'')) LIKE '%mechanism%' THEN 'mechanism_annotation'
        WHEN lower(coalesce(evidence_type,'')) LIKE '%site%' THEN 'binding_site_annotation'
        ELSE 'curated_or_other'
      END AS evidence_modality,
      CASE
        WHEN trim(coalesce(pdb_ids,'')) <> '' AND trim(coalesce(target_uniprot_id,'')) <> '' AND trim(coalesce(compound_internal_id,'')) <> '' THEN 'structure_key_resolved'
        WHEN lower(coalesce(source_database,'')) LIKE '%pubchem%' AND trim(coalesce(source_record_id,'')) <> '' THEN 'pubchem_key_resolved'
        WHEN trim(coalesce(pubmed_ids,'')) <> '' OR trim(coalesce(doi,'')) <> '' THEN 'literature_key_resolved'
        ELSE 'database_record_only'
      END AS lineage_resolution_status,
      'evidence_lineage_v0.1.0' AS lineage_rule_version
    FROM evidence_raw
    """)

    con.execute(r"""
    CREATE TABLE structure_site_lineage AS
    SELECT binding_site_instance_id, evidence_id, target_uniprot_id, compound_internal_id,
      source_database, pdb_ids, pdb_chain_ids_v60, residue_or_site_description,
      'SIT-' || upper(substr(md5(concat_ws('|', upper(trim(coalesce(target_uniprot_id,''))),
      upper(regexp_replace(trim(coalesce(pdb_ids,'')),'[;, ]+',';','g')),
      upper(regexp_replace(trim(coalesce(pdb_chain_ids_v60,'')),'[;, ]+',';','g')),
      trim(coalesce(compound_internal_id,'')), trim(coalesce(residue_or_site_description,'')))),1,20)) AS structure_site_lineage_key,
      CASE WHEN trim(coalesce(pdb_ids,'')) <> '' THEN 'resolved' ELSE 'no_pdb_identifier' END AS lineage_status,
      'structure_site_lineage_v0.1.0' AS lineage_rule_version
    FROM sites_raw
    """)

    con.execute(r"""
    CREATE TABLE pair_lineage AS
    SELECT target_uniprot_id, compound_internal_id,
      count(*) AS evidence_row_count_v63,
      count(DISTINCT upper(trim(source_database))) AS distinct_contributing_database_count_v63,
      string_agg(DISTINCT trim(source_database), ';' ORDER BY trim(source_database)) AS contributing_databases_v63,
      count(DISTINCT structure_experiment_key) FILTER (WHERE structure_experiment_key IS NOT NULL) AS independent_structure_count_v63,
      count(DISTINCT pubchem_experiment_key) FILTER (WHERE pubchem_experiment_key IS NOT NULL) AS independent_pubchem_assay_count_v63,
      count(DISTINCT literature_experiment_key) FILTER (WHERE literature_experiment_key IS NOT NULL) AS independent_literature_experiment_count_v63,
      count(DISTINCT evidence_modality) AS independent_evidence_modality_count_v63,
      string_agg(DISTINCT evidence_modality, ';' ORDER BY evidence_modality) AS evidence_modalities_v63,
      sum(CASE WHEN lineage_resolution_status='database_record_only' THEN 1 ELSE 0 END) AS database_only_lineage_rows_v63
    FROM evidence_lineage
    WHERE trim(coalesce(compound_internal_id,'')) <> ''
    GROUP BY target_uniprot_id, compound_internal_id
    """)

    con.execute(r"""
    CREATE TABLE pair_candidate AS
    SELECT p.*,
      coalesce(l.distinct_contributing_database_count_v63,0) AS distinct_contributing_database_count_v63,
      coalesce(l.contributing_databases_v63,'') AS contributing_databases_v63,
      coalesce(l.independent_structure_count_v63,0) AS independent_structure_count_v63,
      coalesce(l.independent_pubchem_assay_count_v63,0) AS independent_pubchem_assay_count_v63,
      coalesce(l.independent_literature_experiment_count_v63,0) AS independent_literature_experiment_count_v63,
      coalesce(l.independent_evidence_modality_count_v63,0) AS independent_evidence_modality_count_v63,
      coalesce(l.evidence_modalities_v63,'') AS evidence_modalities_v63,
      coalesce(l.database_only_lineage_rows_v63,0) AS database_only_lineage_rows_v63,
      CASE WHEN try_cast(p.independent_source_count AS BIGINT)=coalesce(l.distinct_contributing_database_count_v63,0)
        THEN 'legacy_count_matches_database_count' ELSE 'legacy_count_requires_review' END AS legacy_source_count_semantics_audit_v63,
      'distinct database count is not an independent experiment count' AS source_count_semantics_v63,
      'v6.3_candidate' AS release_version_v63
    FROM pair_raw p
    LEFT JOIN pair_lineage l USING(target_uniprot_id, compound_internal_id)
    """)

    con.execute(r"""
    CREATE TABLE review_queue AS
    WITH docking AS (
      SELECT target_uniprot_id, compound_internal_id, 'docking_nonconflict_candidate' AS docking_scope FROM docking_full_raw
      UNION ALL
      SELECT target_uniprot_id, compound_internal_id, 'docking_conflict_candidate' AS docking_scope FROM docking_conflict_raw
    ), issue AS (
      SELECT e.*,
        CASE
          WHEN lower(coalesce(conflict_status_v60,'')) NOT IN ('','none','no_conflict','baseline_policy_preserved') THEN 'existing_conflict_flag'
          WHEN evidence_tier='BE1' AND evidence_modality='experimental_structure' AND structure_experiment_key IS NULL THEN 'BE1_structure_lineage_unresolved'
          WHEN evidence_tier IN ('BE1','BE2') AND lineage_resolution_status='database_record_only' THEN 'BE1_BE2_database_only_lineage'
          ELSE NULL END AS issue_reason
      FROM evidence_lineage e
      WHERE evidence_tier IN ('BE1','BE2') OR lower(coalesce(conflict_status_v60,'')) NOT IN ('','none','no_conflict','baseline_policy_preserved')
    )
    SELECT i.evidence_id, i.target_uniprot_id, i.compound_internal_id, i.source_database,
      i.evidence_tier, i.evidence_modality, i.lineage_resolution_status,
      i.database_record_key, i.structure_experiment_key, i.pubchem_experiment_key,
      i.literature_experiment_key, i.conflict_status_v60,
      coalesce(i.issue_reason, CASE WHEN d.docking_scope IS NOT NULL THEN 'docking_candidate_lineage_audit' END) AS review_reason,
      coalesce(d.docking_scope,'') AS docking_scope,
      CASE WHEN i.issue_reason='existing_conflict_flag' OR d.docking_scope='docking_conflict_candidate' THEN 'P1'
           WHEN i.evidence_tier='BE1' OR d.docking_scope IS NOT NULL THEN 'P2' ELSE 'P3' END AS review_priority,
      'not_manually_adjudicated' AS manual_review_status
    FROM issue i LEFT JOIN docking d USING(target_uniprot_id, compound_internal_id)
    WHERE i.issue_reason IS NOT NULL OR d.docking_scope IS NOT NULL
    """)

    # A bounded, reproducible validation sample is preferable to pretending that
    # hundreds of thousands of BE1/BE2 rows were manually checked.
    con.execute(r"""
    CREATE TABLE review_sample AS
    SELECT * EXCLUDE(rn) FROM (
      SELECT *, row_number() OVER(PARTITION BY source_database, evidence_tier, review_reason ORDER BY md5(evidence_id)) AS rn
      FROM review_queue
    ) WHERE rn <= 25
    """)

    con.execute(r"""
    CREATE TABLE docking_v63 AS
    WITH combined AS (
      SELECT *, 'nonconflict' AS v62_source_set FROM docking_full_raw
      UNION ALL BY NAME
      SELECT *, 'conflict_review' AS v62_source_set FROM docking_conflict_raw
    )
    SELECT d.*,
      coalesce(l.distinct_contributing_database_count_v63,0) AS distinct_contributing_database_count_v63,
      coalesce(l.independent_structure_count_v63,0) AS independent_structure_count_v63,
      coalesce(l.independent_pubchem_assay_count_v63,0) AS independent_pubchem_assay_count_v63,
      coalesce(l.independent_literature_experiment_count_v63,0) AS independent_literature_experiment_count_v63,
      coalesce(l.independent_evidence_modality_count_v63,0) AS independent_evidence_modality_count_v63,
      p.membrane_role_primary_v63, p.structural_family_primary_v63,
      p.molecular_function_primary_v63, p.biological_process_primary_v63,
      (least(coalesce(l.independent_structure_count_v63,0),3)*2
       + least(coalesce(l.independent_literature_experiment_count_v63,0),3)
       + least(coalesce(l.independent_pubchem_assay_count_v63,0),3)
       + least(coalesce(l.independent_evidence_modality_count_v63,0),3)) AS lineage_bonus_v63,
      CASE WHEN try_cast(d.positive_negative_conflict_flag_v62 AS INTEGER)=1 THEN 15 ELSE 0 END AS conflict_penalty_v63,
      coalesce(try_cast(d.ranking_score_v2 AS DOUBLE),0)
       + (least(coalesce(l.independent_structure_count_v63,0),3)*2
       + least(coalesce(l.independent_literature_experiment_count_v63,0),3)
       + least(coalesce(l.independent_pubchem_assay_count_v63,0),3)
       + least(coalesce(l.independent_evidence_modality_count_v63,0),3))
       - CASE WHEN try_cast(d.positive_negative_conflict_flag_v62 AS INTEGER)=1 THEN 15 ELSE 0 END AS ranking_score_v63,
      'v6.2 shortlist membership preserved; V6.3 adds lineage bonus and conflict penalty only' AS ranking_policy_v63
    FROM combined d
    LEFT JOIN pair_lineage l USING(target_uniprot_id, compound_internal_id)
    LEFT JOIN protein_class_raw p USING(target_uniprot_id)
    """)

    outputs = {
        "binding_evidence_lineage_v0_1.parquet": "evidence_lineage",
        "structure_site_lineage_v0_1.parquet": "structure_site_lineage",
        "protein_compound_lineage_summary_v0_1.tsv.gz": "pair_lineage",
        "protein_compound_summary_v6_3_candidate.tsv.gz": "pair_candidate",
        "targeted_evidence_review_queue_v0_1.tsv.gz": "review_queue",
        "targeted_evidence_validation_sample_v0_1.tsv": "review_sample",
        "docking_priority_v6_3_candidate.tsv.gz": "docking_v63",
    }
    for name, table in outputs.items():
        path = a.output_dir / name
        if name.endswith(".parquet"):
            con.execute(f"COPY {table} TO '{q(path)}' (FORMAT PARQUET, COMPRESSION ZSTD)")
        elif name.endswith(".gz"):
            con.execute(f"COPY {table} TO '{q(path)}' (FORMAT CSV, HEADER, DELIMITER '\\t', COMPRESSION GZIP)")
        else:
            con.execute(f"COPY {table} TO '{q(path)}' (FORMAT CSV, HEADER, DELIMITER '\\t')")

    dictionary = pd.DataFrame([
        ["distinct_contributing_database_count_v63", "COUNT(DISTINCT source_database)", "Database labels contributing records; not independent experiments"],
        ["structure_experiment_key", "target + normalized PDB set + ligand HET + canonical compound", "Structure lineage; site table provides chain/residue-granular companion key"],
        ["pubchem_experiment_key", "source record + relationship context/AID + CID + target + assay", "PubChem assay lineage"],
        ["literature_experiment_key", "PubMed/DOI + target + compound + measurement + assay", "Literature experiment proxy; same paper can still contain multiple experiments"],
        ["independent_evidence_modality_count_v63", "distinct controlled evidence modality", "Orthogonal evidence modes, not source count"],
    ], columns=["field", "definition", "interpretation"])
    dictionary.to_csv(a.output_dir / "evidence_lineage_key_dictionary_v0_1.tsv", sep="\t", index=False)

    def scalar(sql): return con.execute(sql).fetchone()[0]
    checks = {
        "evidence_rows_preserved": scalar("select count(*) from evidence_lineage") == scalar("select count(*) from evidence_raw"),
        "evidence_ids_unique": scalar("select count(*)=count(distinct evidence_id) from evidence_lineage"),
        "pair_candidate_rows_preserved": scalar("select count(*) from pair_candidate") == scalar("select count(*) from pair_raw"),
        "pair_candidate_key_unique": scalar("select count(*)=count(distinct target_uniprot_id||'|'||compound_internal_id) from pair_candidate"),
        "database_count_named_correctly": scalar("select count(*) from pair_candidate where distinct_contributing_database_count_v63 < 1") == 0,
        "review_scope_limited": scalar("select count(*) from review_queue where evidence_tier not in ('BE1','BE2') and review_reason <> 'existing_conflict_flag' and docking_scope='' ") == 0,
        "docking_membership_preserved": scalar("select count(*) from docking_v63") == scalar("select count(*) from docking_full_raw") + scalar("select count(*) from docking_conflict_raw"),
    }
    counts = {
        "evidence_rows": scalar("select count(*) from evidence_lineage"),
        "pair_rows": scalar("select count(*) from pair_candidate"),
        "structure_lineage_keys": scalar("select count(distinct structure_experiment_key) from evidence_lineage"),
        "structure_site_lineage_keys": scalar("select count(distinct structure_site_lineage_key) from structure_site_lineage where lineage_status='resolved'"),
        "pubchem_lineage_keys": scalar("select count(distinct pubchem_experiment_key) from evidence_lineage"),
        "literature_lineage_keys": scalar("select count(distinct literature_experiment_key) from evidence_lineage"),
        "review_queue_rows": scalar("select count(*) from review_queue"),
        "review_sample_rows": scalar("select count(*) from review_sample"),
        "docking_candidate_rows": scalar("select count(*) from docking_v63"),
        "legacy_database_count_mismatches": scalar("select count(*) from pair_candidate where legacy_source_count_semantics_audit_v63='legacy_count_requires_review'"),
    }
    validation = {
        "module": "evidence_lineage_and_docking_rerank", "version": "0.1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL", "checks": checks, "counts": counts,
        "inputs_sha256": {"evidence": h(a.evidence), "sites": h(a.sites), "pair_summary": h(a.pair_summary), "docking_full": h(a.docking_full), "docking_conflict": h(a.docking_conflict), "protein_classification": h(a.protein_classification)},
        "semantic_warning": "Database contribution counts are not independent experiment counts. Structure, PubChem, literature and modality counts are separate fields.",
    }
    (a.qa_dir / "EVIDENCE_LINEAGE_V0_1_VALIDATION.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    con.close()
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
