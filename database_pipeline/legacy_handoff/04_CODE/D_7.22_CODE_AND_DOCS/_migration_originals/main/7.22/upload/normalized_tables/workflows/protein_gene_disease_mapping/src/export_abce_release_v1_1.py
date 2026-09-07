from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml
from openpyxl import Workbook

from .common import clean_cell, load_config, utc_now
from .compare_evidence_filter_rules import KEY, RELATION_KEY, sha256
from .export_balanced_two_source_v1_1 import (
    RELEASE_TAG,
    SCHEMA_VERSION,
    asset,
    entity_exclusion_reason,
    replace_section,
    source_features,
    update_checksums,
    write_json,
    write_parquet,
)


SCORING_VERSION = "pgd_evidence_level_abce_v1"
EXPECTED = {"A": 63269, "B": 6978, "C": 3860}


def relation_metrics(frame: pd.DataFrame) -> dict[str, Any]:
    per_protein = frame.groupby("target_uniprot_id")["disease_id"].nunique()
    return {
        "relationship_count": int(len(frame)),
        "protein_count": int(frame["target_uniprot_id"].nunique()),
        "disease_id_count": int(frame["disease_id"].nunique()),
        "median_diseases_per_protein": float(per_protein.median()) if len(per_protein) else 0.0,
        "mean_diseases_per_protein": float(per_protein.mean()) if len(per_protein) else 0.0,
        "duplicate_key_count": int(frame.duplicated(RELATION_KEY).sum()),
    }


def build_relations(release: Path, source_policy: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, str]]:
    inputs = {
        "core": release / "protein_gene_disease_direct_core_v1.parquet",
        "evidence": release / "protein_disease_evidence_summary_v1.parquet",
        "mapping": release / "protein_gene_mapping_v1.parquet",
        "disease": release / "disease_dictionary_v1.parquet",
        "primary": release / "protein_primary_disease_v1.tsv",
    }
    hashes = {name: sha256(path) for name, path in inputs.items()}
    core = pd.read_parquet(inputs["core"])
    evidence = pd.read_parquet(inputs["evidence"], columns=["ensembl_gene_id", "disease_id", "datasource_id"])
    features, unknown_ids = source_features(evidence, source_policy)
    qualifying = evidence.loc[evidence["rule_a_qualifies"]].drop_duplicates(KEY + ["datasource_id"])
    features = features.join(qualifying.groupby(KEY).size().rename("rule_a_qualifying_datasource_count"))
    features["rule_a_qualifying_datasource_count"] = features["rule_a_qualifying_datasource_count"].fillna(0).astype("int32")
    disease = pd.read_parquet(inputs["disease"], columns=[
        "disease_id", "disease_name", "entity_type", "pathology_category", "specificity_status",
    ])
    mapping = pd.read_parquet(inputs["mapping"])
    if unknown_ids:
        raise RuntimeError(f"Unknown datasource IDs prevent formal scoring: {unknown_ids}")
    if core.duplicated(RELATION_KEY).any() or disease.duplicated(["disease_id"]).any() or mapping.duplicated(["target_uniprot_id", "ensembl_gene_id"]).any():
        raise RuntimeError("Input unique-key validation failed")
    relations = core.merge(disease, on="disease_id", how="left", validate="many_to_one").join(
        features, on=KEY, how="left", validate="many_to_one"
    )
    numeric = [
        "raw_datasource_count", "rule_a_qualifying_datasource_count", "independent_source_family_count",
        "core_source_count", "supporting_source_count", "supporting_only_source_count",
        "non_qualifying_source_count", "clinical_or_curated_genetic_source_count",
        "rule_c_eligible_source_count", "second_core_or_human_genetic_source_count", "high_authority_source_count",
    ]
    text = ["all_source_families", "independent_source_families", "core_source_families",
            "supporting_source_families", "raw_datasource_ids"]
    relations[numeric] = relations[numeric].fillna(0).astype("int32")
    relations[text] = relations[text].fillna("")
    relations[["weak_sources_only_flag", "unknown_source_flag"]] = relations[
        ["weak_sources_only_flag", "unknown_source_flag"]
    ].fillna(False).astype(bool)
    relations["entity_exclusion_reason"] = entity_exclusion_reason(relations)
    return relations, mapping, pd.read_csv(inputs["primary"], sep="\t", dtype=str, keep_default_na=False), hashes


def apply_rules(relations: pd.DataFrame) -> pd.DataFrame:
    direct = relations["association_scope"].eq("direct")
    disease = relations["entity_type"].eq("disease")
    known = ~relations["unknown_source_flag"]
    not_weak = ~relations["weak_sources_only_flag"]
    rule_a = (relations["rule_a_qualifying_datasource_count"] >= 2) & known
    rule_b = (
        direct & disease & known & not_weak
        & (relations["independent_source_family_count"] >= 2)
        & (relations["core_source_count"] >= 1)
    )
    rule_c = (
        rule_b
        & (relations["clinical_or_curated_genetic_source_count"] >= 1)
        & (relations["second_core_or_human_genetic_source_count"] >= 1)
    )
    original_strong = (
        ((relations["independent_source_family_count"] >= 3) & (relations["core_source_count"] >= 2))
        | ((relations["independent_source_family_count"] >= 2) & (relations["high_authority_source_count"] >= 2))
    )
    # Rule E explicitly intersects Rule C so the formal hierarchy is guaranteed;
    # it is a corrected recomputation, not a rename of historical Rule D.
    rule_e = rule_c & original_strong & direct & disease & known & not_weak & relations["entity_exclusion_reason"].eq("")
    relations["rule_a_pass"] = rule_a
    relations["rule_b_pass"] = rule_b
    relations["rule_c_pass"] = rule_c
    relations["rule_e_pass"] = rule_e
    relations["default_release_inclusion"] = rule_b
    conditions = [rule_e, rule_c, rule_b, rule_a]
    relations["evidence_level"] = np.select(conditions, ["very_high", "high", "medium", "low"], default="not_included")
    relations["evidence_level_rank"] = np.select(conditions, [4, 3, 2, 1], default=0).astype("int8")
    relations["highest_rule_passed"] = np.select(conditions, ["E", "C", "B", "A"], default="")
    relations["scoring_policy_version"] = SCORING_VERSION
    return relations


def hierarchy_validation(relations: pd.DataFrame, reports: Path) -> tuple[dict[str, int], pd.DataFrame]:
    a, b, c, e = (relations[f"rule_{rule.lower()}_pass"] for rule in "ABCE")
    exceptions = pd.concat([
        relations.loc[b & ~a].assign(hierarchy_exception="B_not_in_A"),
        relations.loc[c & ~b].assign(hierarchy_exception="C_not_in_B"),
        relations.loc[e & ~c].assign(hierarchy_exception="E_not_in_C"),
    ], ignore_index=True)
    fields = RELATION_KEY + ["disease_name", "entity_type", "specificity_status", "independent_source_families",
                             "rule_a_pass", "rule_b_pass", "rule_c_pass", "rule_e_pass", "hierarchy_exception"]
    exceptions[fields].sort_values(RELATION_KEY).to_csv(
        reports / "rule_hierarchy_exceptions.tsv", sep="\t", index=False, lineterminator="\n"
    )
    counts = {
        "A_count": int(a.sum()), "B_count": int(b.sum()), "C_count": int(c.sum()), "E_count": int(e.sum()),
        "B_not_in_A_count": int((b & ~a).sum()), "C_not_in_B_count": int((c & ~b).sum()),
        "E_not_in_C_count": int((e & ~c).sum()),
    }
    status = "passed" if not len(exceptions) and all(counts[rule] == EXPECTED[rule[0]] for rule in ["A_count", "B_count", "C_count"]) else "failed"
    report = f'''# A/B/C/E rule hierarchy validation

- Status: **{status}**
- Rule A: {counts['A_count']:,}
- Rule B: {counts['B_count']:,}
- Rule C: {counts['C_count']:,}
- Rule E (corrected recomputation): {counts['E_count']:,}
- B not in A: {counts['B_not_in_A_count']:,}
- C not in B: {counts['C_not_in_B_count']:,}
- E not in C: {counts['E_not_in_C_count']:,}

Rule E is computed as the intersection of Rule C, the original strong multi-source predicate, direct disease scope, explicit specificity exclusions, known sources and non-weak evidence. Historical Rule D is retained only in the trial comparison report and is not used as a formal evidence level.
'''
    (reports / "RULE_HIERARCHY_VALIDATION.md").write_text(report, encoding="utf-8")
    if status != "passed":
        raise RuntimeError(f"Rule hierarchy validation failed: {counts}")
    return counts, exceptions


def scored_columns(core_columns: list[str]) -> list[str]:
    extras = [
        "disease_name", "entity_type", "pathology_category", "specificity_status",
        "raw_datasource_count", "independent_source_family_count", "core_source_count",
        "independent_source_families", "core_source_families", "supporting_source_families",
        "weak_sources_only_flag", "unknown_source_flag", "rule_a_pass", "rule_b_pass", "rule_c_pass", "rule_e_pass",
        "highest_rule_passed", "evidence_level", "evidence_level_rank", "default_release_inclusion",
        "scoring_policy_version",
    ]
    return core_columns + [column for column in extras if column not in core_columns]


def rank_rule_b(rule_b: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    mapping_fields = ["target_uniprot_id", "ensembl_gene_id", "protein_name", "approved_symbol", "hgnc_id", "ncbi_gene_id"]
    ranked = rule_b.merge(mapping[mapping_fields], on=["target_uniprot_id", "ensembl_gene_id"], how="left", validate="many_to_one")
    ranked["_clinical_rank"] = (~ranked["clinical_genetic_flag"]).astype(int)
    ranked["_human_rank"] = (~ranked["human_genetic_flag"]).astype(int)
    ranked["_evidence_count_sort"] = ranked["evidence_count"].fillna(-1)
    ranked = ranked.sort_values(
        ["target_uniprot_id", "evidence_level_rank", "_clinical_rank", "_human_rank", "core_source_count",
         "independent_source_family_count", "ot_overall_score", "_evidence_count_sort", "disease_id"],
        ascending=[True, False, True, True, False, False, False, False, True], kind="mergesort",
    )
    ranked["rule_b_rank"] = ranked.groupby("target_uniprot_id").cumcount() + 1
    return ranked


def summary_outputs(rule_b: pd.DataFrame, mapping: pd.DataFrame, primary_template: pd.DataFrame,
                    release: Path, qc: dict[str, Any]) -> tuple[Path, Path, Path, int, int]:
    ranked = rank_rule_b(rule_b, mapping)
    detail_fields = [
        "target_uniprot_id", "protein_name", "approved_symbol", "ensembl_gene_id", "hgnc_id", "ncbi_gene_id",
        "disease_id", "disease_name", "entity_type", "pathology_category", "project_evidence_tier", "ot_overall_score",
        "evidence_level", "evidence_level_rank", "highest_rule_passed", "clinical_genetic_flag", "human_genetic_flag",
        "core_source_count", "independent_source_family_count", "independent_source_families", "rule_b_rank",
    ]
    top10 = ranked.loc[ranked["rule_b_rank"] <= 10, detail_fields].copy()
    top10_path = release / "protein_gene_disease_rule_b_top10_v1.tsv"
    top10.to_csv(top10_path, sep="\t", index=False, lineterminator="\n")
    selected = ranked.loc[ranked["rule_b_rank"] == 1, detail_fields].copy()
    selected = selected.rename(columns={
        "disease_id": "primary_disease_id", "disease_name": "primary_disease_name",
        "entity_type": "primary_entity_type", "pathology_category": "primary_pathology_category",
        "project_evidence_tier": "primary_project_evidence_tier", "ot_overall_score": "primary_ot_overall_score",
        "evidence_level": "primary_evidence_level", "evidence_level_rank": "primary_evidence_level_rank",
        "highest_rule_passed": "primary_highest_rule_passed",
    })
    protein_fields = ["target_uniprot_id", "protein_name", "approved_symbol", "protein_class", "membrane_protein_type", "record_qc_status"]
    primary = primary_template[protein_fields].merge(
        selected.drop(columns=["protein_name", "approved_symbol"]), on="target_uniprot_id", how="left", validate="one_to_one"
    )
    primary["primary_status"] = np.where(primary["primary_disease_id"].notna(), "selected_from_rule_b", "no_rule_b_supported_disease")
    for column in primary.columns:
        if primary[column].dtype == object:
            primary[column] = primary[column].fillna("")
    primary_path = release / "protein_primary_disease_rule_b_v1.tsv"
    primary.to_csv(primary_path, sep="\t", index=False, lineterminator="\n")

    workbook_path = release / "protein_gene_disease_rule_b_summary_v1.xlsx"
    workbook = Workbook(write_only=True)

    def add_frame(name: str, frame: pd.DataFrame) -> None:
        ws = workbook.create_sheet(name)
        ws.append(list(frame.columns))
        for row in frame.itertuples(index=False, name=None):
            ws.append([clean_cell(value) for value in row])

    add_frame("1_Rule_B_Direct", rule_b)
    add_frame("2_Primary_Disease", primary)
    add_frame("3_Top10_Diseases", top10)
    level_summary = rule_b.groupby(["evidence_level", "evidence_level_rank", "highest_rule_passed"], as_index=False).agg(
        relationship_count=("disease_id", "size"), protein_count=("target_uniprot_id", "nunique"), disease_count=("disease_id", "nunique")
    ).sort_values("evidence_level_rank")
    add_frame("4_Evidence_Level_Summary", level_summary)
    exploded = rule_b[["target_uniprot_id", "disease_id", "independent_source_families"]].copy()
    exploded["source_family"] = exploded["independent_source_families"].str.split(";")
    exploded = exploded.explode("source_family")
    source_summary = exploded.groupby("source_family", as_index=False).agg(
        relationship_count=("disease_id", "size"), protein_count=("target_uniprot_id", "nunique"), disease_count=("disease_id", "nunique")
    ).sort_values(["relationship_count", "source_family"], ascending=[False, True])
    add_frame("5_Source_Family_Summary", source_summary)
    add_frame("6_Protein_Gene_Mapping", mapping)
    disease_summary = rule_b.groupby(["disease_id", "disease_name", "pathology_category"], as_index=False).agg(
        relationship_count=("target_uniprot_id", "size"), protein_count=("target_uniprot_id", "nunique")
    ).sort_values(["relationship_count", "disease_id"], ascending=[False, True])
    add_frame("7_Disease_Summary", disease_summary)
    add_frame("8_No_Rule_B_Disease", primary.loc[primary["primary_status"] == "no_rule_b_supported_disease"])
    add_frame("9_QC_Report", pd.DataFrame([{"metric": key, "value": clean_cell(value)} for key, value in qc.items()]))
    metadata = pd.DataFrame([
        {"field": "scoring_policy_version", "value": SCORING_VERSION},
        {"field": "default_release_rule", "value": "Rule B"},
        {"field": "archive_note", "value": "A-only low relationships are stored only on the archive branch."},
        {"field": "historical_rule_d", "value": "Historical trial only; not a formal evidence level."},
        {"field": "generated_at_utc", "value": utc_now()},
    ])
    add_frame("10_Metadata", metadata)
    workbook.save(workbook_path)
    return primary_path, top10_path, workbook_path, len(primary), len(top10)


def make_asset(path: Path, root: Path, frame: pd.DataFrame, role: str, unique_key: str = "+".join(RELATION_KEY)) -> dict[str, Any]:
    first = str(frame.iloc[0]["target_uniprot_id"]) if len(frame) and "target_uniprot_id" in frame else ""
    last = str(frame.iloc[-1]["target_uniprot_id"]) if len(frame) and "target_uniprot_id" in frame else ""
    return asset(path, root, len(frame), len(frame.columns), unique_key, role, first, last)


def make_excel_asset(path: Path, root: Path) -> dict[str, Any]:
    record = asset(path, root, 0, 10, "workbook", "rule_b_summary_excel")
    record.update({
        "format": "xlsx",
        "compression": "zip/container",
        "git_tracking_status": "ordinary_git_ok",
    })
    return record


def report_artifact(path: Path, root: Path) -> dict[str, Any]:
    return {
        "filename": path.name,
        "relative_local_path": path.relative_to(root).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def update_documentation(release: Path, counts: dict[str, int], level_counts: dict[str, int]) -> None:
    readme = f'''## v1.1 A/B/C/E evidence scoring

The default `main` release is Rule B and above. Rule A-only relationships are `low` and archive-only; Rule B is the `medium` inclusion threshold; Rule C is `high`; corrected Rule E is `very_high`. The unified Rule B table contains {counts['B_count']:,} relationships and retains `rule_c_pass`, `rule_e_pass`, `highest_rule_passed`, `evidence_level` and `evidence_level_rank` so higher-confidence subsets remain visible without duplicate main-branch datasets.

Historical Rule D had an entity-type defect and is not a formal evidence level. Rule E is recomputed with direct disease scope, specificity exclusions, known sources, non-weak evidence and the Rule C hierarchy. EVA remains `EVA` and is not inferred to be ClinVar. The four levels are evidence-policy categories, not disease-causality probabilities; `medium` is the database inclusion threshold and does not mean poor evidence quality.

Full A/B/C/E scored outputs and independent A/C/E files are stored on `archive/protein-gene-disease-abce-v1.1`.'''
    methods = f'''## v1.1 scoring method

Source channels are collapsed to independent families before counting: UniProt channels to `UNIPROT`, CRISPR channels to `CRISPR`, and EVA channels to `EVA`. Europe PMC and Clinical Precedence do not qualify for Rule B; Expression Atlas and IMPC cannot pass without a core family.

- Rule A (`low` when A-only): at least two qualifying raw datasource IDs.
- Rule B (`medium`, default inclusion): direct disease, at least two independent qualifying families, at least one core family, known source and not weak-only.
- Rule C (`high`): Rule B plus a clinical/curated genetic family and a different eligible core/human-genetic family.
- Rule E (`very_high`): Rule C plus the original strong multi-source predicate and corrected direct/disease/specificity/known/non-weak constraints.

Counts: A={counts['A_count']:,}, B={counts['B_count']:,}, C={counts['C_count']:,}, E={counts['E_count']:,}. Formal levels: low={level_counts['low']:,}, medium={level_counts['medium']:,}, high={level_counts['high']:,}, very_high={level_counts['very_high']:,}. Single high-authority sources remain outside Rule B and require raw evidence validation.'''
    dictionary = '''## v1.1 scoring fields and outputs

`rule_a_pass`, `rule_b_pass`, `rule_c_pass` and `rule_e_pass` are boolean rule membership fields. `highest_rule_passed` is A/B/C/E. `evidence_level` is `low`, `medium`, `high` or `very_high`; `evidence_level_rank` is 1 through 4. `default_release_inclusion` equals `rule_b_pass`. `scoring_policy_version` is `pgd_evidence_level_abce_v1`.

`protein_gene_disease_abce_scored_v1.*` is archive-only and contains all Rule A relationships. `protein_gene_disease_rule_b_scored_v1.*` is the default main release and contains only medium/high/very-high rows. `protein_primary_disease_rule_b_v1.tsv`, `protein_gene_disease_rule_b_top10_v1.tsv` and `protein_gene_disease_rule_b_summary_v1.xlsx` are derived exclusively from Rule B relationships. Independent Rule A/C/E and level-split Parquets are archive-only.'''
    replace_section(release / "README.md", "ABCE_SCORING_V1_1", readme)
    replace_section(release / "METHODS.md", "ABCE_SCORING_V1_1", methods)
    replace_section(release / "DATA_DICTIONARY.md", "ABCE_SCORING_V1_1", dictionary)


def run_export(config_path: str, source_policy_path: str, scoring_policy_path: str) -> dict[str, Any]:
    root, config = load_config(config_path)
    release = root / config["outputs"]["release_dir"]
    reports = root / config["outputs"]["workflow_dir"] / "reports"
    source_policy_file = root / source_policy_path
    scoring_policy_file = root / scoring_policy_path
    source_policy = yaml.safe_load(source_policy_file.read_text(encoding="utf-8"))
    scoring_policy = yaml.safe_load(scoring_policy_file.read_text(encoding="utf-8"))
    if scoring_policy["scoring_policy_version"] != SCORING_VERSION:
        raise RuntimeError("Unexpected scoring policy version")
    relations, mapping, primary_template, protected_before = build_relations(release, source_policy)
    core_columns = list(pq.ParquetFile(release / "protein_gene_disease_direct_core_v1.parquet").schema_arrow.names)
    relations = apply_rules(relations)
    hierarchy_counts, _ = hierarchy_validation(relations, reports)
    scored = relations.loc[relations["rule_a_pass"], scored_columns(core_columns)].sort_values(RELATION_KEY, kind="mergesort").reset_index(drop=True)
    rule_b = scored.loc[scored["rule_b_pass"]].reset_index(drop=True)
    rule_c = scored.loc[scored["rule_c_pass"]].reset_index(drop=True)
    rule_e = scored.loc[scored["rule_e_pass"]].reset_index(drop=True)
    rule_a = scored.copy()
    level_frames = {level: scored.loc[scored["evidence_level"] == level].reset_index(drop=True) for level in ["low", "medium", "high", "very_high"]}
    level_counts = {level: len(frame) for level, frame in level_frames.items()}
    b_metrics = relation_metrics(rule_b)
    expected_b_metrics = {
        "relationship_count": 6978,
        "protein_count": 1119,
        "disease_id_count": 2310,
        "median_diseases_per_protein": 3.0,
        "duplicate_key_count": 0,
    }
    if any(b_metrics[key] != value for key, value in expected_b_metrics.items()):
        raise RuntimeError(f"Rule B validation failed: {b_metrics}")
    if rule_b["evidence_level"].isin({"low", "not_included"}).any() or rule_b["unknown_source_flag"].any() or rule_b["weak_sources_only_flag"].any() or (~rule_b["entity_type"].eq("disease")).any():
        raise RuntimeError("Rule B release contains a forbidden row")

    paths = {
        "abce_parquet": release / "protein_gene_disease_abce_scored_v1.parquet",
        "abce_gzip": release / "protein_gene_disease_abce_scored_v1.tsv.gz",
        "rule_a": release / "protein_gene_disease_rule_a_v1.parquet",
        "rule_b_parquet": release / "protein_gene_disease_rule_b_scored_v1.parquet",
        "rule_b_gzip": release / "protein_gene_disease_rule_b_scored_v1.tsv.gz",
        "rule_c": release / "protein_gene_disease_rule_c_v1.parquet",
        "rule_e": release / "protein_gene_disease_rule_e_v1.parquet",
        "medium": release / "protein_gene_disease_medium_v1.parquet",
        "high": release / "protein_gene_disease_high_v1.parquet",
        "very_high": release / "protein_gene_disease_very_high_v1.parquet",
    }
    for name, frame in [("abce_parquet", scored), ("rule_a", rule_a), ("rule_b_parquet", rule_b),
                        ("rule_c", rule_c), ("rule_e", rule_e), ("medium", level_frames["medium"]),
                        ("high", level_frames["high"]), ("very_high", level_frames["very_high"])]:
        write_parquet(frame, paths[name])
    scored.to_csv(paths["abce_gzip"], sep="\t", index=False, compression={"method": "gzip", "compresslevel": 9, "mtime": 0}, lineterminator="\n")
    rule_b.to_csv(paths["rule_b_gzip"], sep="\t", index=False, compression={"method": "gzip", "compresslevel": 9, "mtime": 0}, lineterminator="\n")

    qc = {
        "status": "passed", "scoring_policy_version": SCORING_VERSION, **hierarchy_counts,
        "low_count": level_counts["low"], "medium_count": level_counts["medium"],
        "high_count": level_counts["high"], "very_high_count": level_counts["very_high"],
        "rule_b_relationship_count": b_metrics["relationship_count"], "rule_b_protein_count": b_metrics["protein_count"],
        "rule_b_disease_id_count": b_metrics["disease_id_count"], "rule_b_median_diseases_per_protein": b_metrics["median_diseases_per_protein"],
        "rule_b_duplicate_key_count": b_metrics["duplicate_key_count"], "low_count_in_rule_b_release": 0,
        "not_included_count_in_rule_b_release": 0, "unknown_source_rows_passed": 0,
        "weak_source_only_rows_passed": 0, "non_disease_rows_passed": 0,
        "historical_rule_d_used_for_scoring": False, "formal_input_hashes_unchanged": True,
    }
    primary_path, top10_path, excel_path, primary_rows, top10_rows = summary_outputs(rule_b, mapping, primary_template, release, qc)
    qc["primary_summary_rows"] = primary_rows
    qc["top10_rows"] = top10_rows
    qc["summary_excel_bytes"] = excel_path.stat().st_size
    if excel_path.stat().st_size >= 100 * 1024 * 1024:
        raise RuntimeError("Rule B summary Excel exceeds 100 MiB")
    write_json(release / "RULE_B_QC_REPORT.json", qc)
    write_json(release / "ABCE_QC_REPORT.json", {
        **qc,
        "scope": "archive A/B/C/E scoring outputs",
        "rule_a_archive_only": True,
        "formal_levels": ["low", "medium", "high", "very_high"],
        "hierarchy_validation_report": "reports/RULE_HIERARCHY_VALIDATION.md",
    })

    frames = {"abce_parquet": scored, "abce_gzip": scored, "rule_a": rule_a, "rule_b_parquet": rule_b,
              "rule_b_gzip": rule_b, "rule_c": rule_c, "rule_e": rule_e,
              "medium": level_frames["medium"], "high": level_frames["high"], "very_high": level_frames["very_high"]}
    archive_assets = [make_asset(paths[name], root, frames[name], name) for name in paths]
    primary_frame = pd.read_csv(primary_path, sep="\t", dtype=str, keep_default_na=False)
    top10_frame = pd.read_csv(top10_path, sep="\t", dtype=str, keep_default_na=False)
    excel_asset = make_excel_asset(excel_path, root)
    summary_assets = [
        make_asset(primary_path, root, primary_frame, "rule_b_primary", "target_uniprot_id"),
        make_asset(top10_path, root, top10_frame, "rule_b_top10", "target_uniprot_id+rule_b_rank"),
        excel_asset,
    ]
    all_assets = archive_assets + summary_assets
    comparison_paths = [
        reports / "EVIDENCE_FILTER_RULE_COMPARISON.md",
        reports / "evidence_filter_rule_counts.tsv",
        reports / "source_combination_frequency.tsv",
        reports / "source_family_mapping.tsv",
        reports / "unmapped_datasources.tsv",
        reports / "single_source_high_confidence_candidates.tsv",
        reports / "RULE_B_TIER_DEPENDENCY_AUDIT.md",
        reports / "rule_b_with_vs_without_tier_diff.tsv",
        reports / "rule_b_formal_vs_trial_diff.tsv",
        reports / "RULE_HIERARCHY_VALIDATION.md",
        reports / "rule_hierarchy_exceptions.tsv",
    ]
    missing_comparison = [str(path) for path in comparison_paths if not path.exists()]
    if missing_comparison:
        raise RuntimeError(f"Missing comparison material: {missing_comparison}")
    write_json(release / "ABCE_SCORING_MANIFEST.json", {
        "generated_at_utc": utc_now(), "scoring_policy_version": SCORING_VERSION,
        "recommended_archive_branch": "archive/protein-gene-disease-abce-v1.1",
        "historical_rule_d_used_for_scoring": False,
        "assets": all_assets,
        "comparison_materials": [report_artifact(path, root) for path in comparison_paths],
    })
    rule_b_names = ["rule_b_parquet", "rule_b_gzip"]
    rule_b_assets = [make_asset(paths[name], root, frames[name], name) for name in rule_b_names] + summary_assets
    write_json(release / "RULE_B_RELEASE_MANIFEST.json", {
        "generated_at_utc": utc_now(), "scoring_policy_version": SCORING_VERSION,
        "default_release_rule": "B", "recommended_release_tag": RELEASE_TAG,
        "archive_branch": "archive/protein-gene-disease-abce-v1.1", "assets": rule_b_assets,
    })
    update_documentation(release, hierarchy_counts, level_counts)

    rule_b_checksum_paths = [root / item["relative_local_path"] for item in rule_b_assets] + [
        release / "README.md", release / "METHODS.md", release / "DATA_DICTIONARY.md",
        release / "RULE_B_QC_REPORT.json", release / "RULE_B_RELEASE_MANIFEST.json",
        source_policy_file, scoring_policy_file,
    ]
    (release / "rule_b_checksums.sha256").write_text(
        "\n".join(f"{sha256(path)}  {path.relative_to(root).as_posix()}" for path in sorted(set(rule_b_checksum_paths))) + "\n",
        encoding="utf-8",
    )
    additions = [root / item["relative_local_path"] for item in all_assets] + comparison_paths + [
        release / "ABCE_SCORING_MANIFEST.json", release / "RULE_B_RELEASE_MANIFEST.json",
        release / "ABCE_QC_REPORT.json", release / "RULE_B_QC_REPORT.json", release / "rule_b_checksums.sha256",
        source_policy_file, scoring_policy_file, release / "README.md", release / "METHODS.md", release / "DATA_DICTIONARY.md",
    ]
    checksum_count = update_checksums(root, release, additions)
    protected_after = {
        "core": sha256(release / "protein_gene_disease_direct_core_v1.parquet"),
        "evidence": sha256(release / "protein_disease_evidence_summary_v1.parquet"),
        "mapping": sha256(release / "protein_gene_mapping_v1.parquet"),
        "disease": sha256(release / "disease_dictionary_v1.parquet"),
        "primary": sha256(release / "protein_primary_disease_v1.tsv"),
    }
    if protected_after != protected_before:
        raise RuntimeError("A protected v1 input changed during A/B/C/E export")
    return {
        "status": "passed", "hierarchy": hierarchy_counts, "level_counts": level_counts,
        "rule_b_metrics": b_metrics, "asset_count": len(all_assets), "checksum_entry_count": checksum_count,
        "rule_b_excel_bytes": excel_path.stat().st_size, "formal_input_hashes_unchanged": True,
    }
