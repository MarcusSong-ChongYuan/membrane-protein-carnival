from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .common import load_config


KEY = ["ensembl_gene_id", "disease_id"]
RELATION_KEY = ["target_uniprot_id", "ensembl_gene_id", "disease_id"]
CORE_LEVELS = {"core", "core_curated", "core_cancer"}
LITERATURE_DATASOURCES = {"europepmc", "uniprot_literature"}
EXPRESSION_ANIMAL_DATASOURCES = {"expression_atlas", "impc"}
POLICY_FIELDS = [
    "source_family", "source_level", "counts_toward_two_source_rule",
    "allowed_as_only_core_source", "rule_a_qualifies", "clinical_or_curated_genetic",
    "second_core_or_human_genetic", "high_authority_source", "source_note",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def joined_values(frame: pd.DataFrame, column: str, output_name: str) -> pd.Series:
    if frame.empty:
        return pd.Series(dtype="string", name=output_name)
    ordered = frame[KEY + [column]].drop_duplicates().sort_values(KEY + [column])
    return ordered.groupby(KEY, sort=False)[column].agg(";".join).rename(output_name)


def grouped_count(frame: pd.DataFrame, output_name: str) -> pd.Series:
    if frame.empty:
        return pd.Series(dtype="int64", name=output_name)
    return frame.drop_duplicates(KEY + ["source_family"]).groupby(KEY, sort=False).size().rename(output_name)


def run_comparison(config_path: str, policy_path: str) -> dict[str, Any]:
    root, config = load_config(config_path)
    policy_file = root / policy_path
    policy = yaml.safe_load(policy_file.read_text(encoding="utf-8"))
    release = root / config["outputs"]["release_dir"]
    reports = root / config["outputs"]["workflow_dir"] / "reports"
    reports.mkdir(parents=True, exist_ok=True)

    inputs = {
        "core": release / "protein_gene_disease_direct_core_v1.parquet",
        "evidence": release / "protein_disease_evidence_summary_v1.parquet",
        "mapping": release / "protein_gene_mapping_v1.parquet",
        "disease": release / "disease_dictionary_v1.parquet",
    }
    missing = [str(path) for path in inputs.values() if not path.exists()]
    if missing:
        raise RuntimeError(f"Missing validated input assets: {missing}")
    hashes_before = {name: sha256(path) for name, path in inputs.items()}

    evidence = pd.read_parquet(inputs["evidence"], columns=["ensembl_gene_id", "disease_id", "datasource_id"])
    evidence["datasource_id"] = evidence["datasource_id"].astype("string")
    observed_ids = sorted(evidence["datasource_id"].dropna().unique().tolist())
    configured = policy.get("source_families", {})
    defaults = policy["unknown_source_defaults"]
    unknown_ids = sorted(set(observed_ids) - set(configured))

    def policy_value(datasource_id: str, field: str):
        return configured.get(datasource_id, defaults).get(field, defaults.get(field))

    for field in POLICY_FIELDS:
        evidence[field] = evidence["datasource_id"].map(lambda value, f=field: policy_value(str(value), f))
    evidence.loc[evidence["datasource_id"].isin(unknown_ids), "source_family"] = (
        "UNKNOWN__" + evidence.loc[evidence["datasource_id"].isin(unknown_ids), "datasource_id"].astype(str).str.upper()
    )

    for field in [
        "counts_toward_two_source_rule", "allowed_as_only_core_source", "rule_a_qualifies",
        "clinical_or_curated_genetic", "second_core_or_human_genetic", "high_authority_source",
    ]:
        evidence[field] = evidence[field].fillna(False).astype(bool)

    # Confirm the policy's Rule C categories can be interpreted as two distinct
    # eligible families, one of which is clinical/curated.
    invalid_c = evidence.loc[
        evidence["clinical_or_curated_genetic"] & ~evidence["second_core_or_human_genetic"], "datasource_id"
    ].unique().tolist()
    if invalid_c:
        raise RuntimeError(f"Rule C policy error; clinical sources must also be second-source eligible: {invalid_c}")

    observed_counts = evidence.groupby("datasource_id", sort=True).agg(
        evidence_row_count=("datasource_id", "size"),
        gene_disease_count=("disease_id", "size"),
    ).reset_index()
    # The input evidence key is already unique, so row count equals the number
    # of gene-disease pairs for each datasource.
    mapping_rows = []
    for row in observed_counts.to_dict("records"):
        datasource_id = str(row["datasource_id"])
        selected = configured.get(datasource_id, defaults)
        mapping_rows.append({
            "datasource_id": datasource_id,
            **{field: selected.get(field, defaults.get(field)) for field in POLICY_FIELDS},
            "evidence_row_count": int(row["evidence_row_count"]),
            "gene_disease_count": int(row["gene_disease_count"]),
            "mapping_status": "unknown_requires_manual_review" if datasource_id in unknown_ids else "mapped",
        })
    source_mapping_fields = ["datasource_id"] + POLICY_FIELDS + ["evidence_row_count", "gene_disease_count", "mapping_status"]
    write_tsv(reports / "source_family_mapping.tsv", mapping_rows, source_mapping_fields)
    write_tsv(
        reports / "unmapped_datasources.tsv",
        [row for row in mapping_rows if row["mapping_status"] != "mapped"],
        source_mapping_fields,
    )

    family_rows = evidence[KEY + ["source_family", "source_level", "counts_toward_two_source_rule",
                                  "clinical_or_curated_genetic", "second_core_or_human_genetic",
                                  "high_authority_source"]].drop_duplicates(KEY + ["source_family"])
    datasource_rows = evidence[KEY + ["datasource_id", "rule_a_qualifies"]].drop_duplicates()

    aggregate = pd.DataFrame(index=pd.MultiIndex.from_frame(evidence[KEY].drop_duplicates()))
    aggregate.index.names = KEY
    aggregate = aggregate.join(evidence.groupby(KEY, sort=False)["datasource_id"].nunique().rename("raw_datasource_count"))
    aggregate = aggregate.join(datasource_rows.loc[datasource_rows["rule_a_qualifies"]].groupby(KEY).size().rename("rule_a_qualifying_datasource_count"))
    aggregate = aggregate.join(grouped_count(family_rows.loc[family_rows["counts_toward_two_source_rule"]], "independent_source_family_count"))
    aggregate = aggregate.join(grouped_count(family_rows.loc[family_rows["source_level"].isin(CORE_LEVELS)], "core_source_count"))
    aggregate = aggregate.join(grouped_count(family_rows.loc[family_rows["source_level"] == "supporting"], "supporting_source_count"))
    aggregate = aggregate.join(grouped_count(family_rows.loc[family_rows["source_level"] == "supporting_only"], "supporting_only_source_count"))
    aggregate = aggregate.join(grouped_count(family_rows.loc[family_rows["source_level"].isin({"non_qualifying", "drug_derived", "unknown"})], "non_qualifying_source_count"))
    aggregate = aggregate.join(grouped_count(family_rows.loc[family_rows["clinical_or_curated_genetic"]], "clinical_or_curated_genetic_source_count"))
    aggregate = aggregate.join(grouped_count(family_rows.loc[family_rows["second_core_or_human_genetic"]], "rule_c_eligible_source_count"))
    aggregate = aggregate.join(grouped_count(family_rows.loc[family_rows["high_authority_source"]], "high_authority_source_count"))
    aggregate = aggregate.join(joined_values(family_rows, "source_family", "all_source_families"))
    aggregate = aggregate.join(joined_values(family_rows.loc[family_rows["counts_toward_two_source_rule"]], "source_family", "independent_source_families"))
    aggregate = aggregate.join(joined_values(family_rows.loc[family_rows["source_level"].isin(CORE_LEVELS)], "source_family", "core_source_families"))
    aggregate = aggregate.join(joined_values(family_rows.loc[family_rows["source_level"] == "supporting"], "source_family", "supporting_source_families"))
    aggregate = aggregate.join(joined_values(datasource_rows, "datasource_id", "raw_datasource_ids"))

    numeric = [
        "raw_datasource_count", "rule_a_qualifying_datasource_count", "independent_source_family_count",
        "core_source_count", "supporting_source_count", "supporting_only_source_count",
        "non_qualifying_source_count", "clinical_or_curated_genetic_source_count",
        "rule_c_eligible_source_count", "high_authority_source_count",
    ]
    aggregate[numeric] = aggregate[numeric].fillna(0).astype("int32")
    strings = ["all_source_families", "independent_source_families", "core_source_families",
               "supporting_source_families", "raw_datasource_ids"]
    aggregate[strings] = aggregate[strings].fillna("")
    aggregate["second_core_or_human_genetic_source_count"] = np.where(
        aggregate["clinical_or_curated_genetic_source_count"] > 0,
        np.maximum(aggregate["rule_c_eligible_source_count"] - 1, 0),
        0,
    ).astype("int32")
    aggregate["weak_sources_only_flag"] = aggregate["core_source_count"] == 0
    aggregate["unknown_source_flag"] = aggregate["all_source_families"].str.contains("UNKNOWN__", regex=False)

    raw_sets = aggregate["raw_datasource_ids"].str.split(";").map(set)
    aggregate["literature_only_flag"] = raw_sets.map(
        lambda values: bool(values) and values.issubset(LITERATURE_DATASOURCES)
    )
    aggregate["expression_plus_animal_only_flag"] = raw_sets.map(
        lambda values: {"expression_atlas", "impc"}.issubset(values) and values.issubset(EXPRESSION_ANIMAL_DATASOURCES)
    )

    core_columns = RELATION_KEY + [
        "association_scope", "project_evidence_tier", "ot_overall_score", "clinical_genetic_flag",
        "human_genetic_flag", "literature_flag", "expression_flag", "animal_model_flag",
    ]
    core = pd.read_parquet(inputs["core"], columns=core_columns)
    diseases = pd.read_parquet(inputs["disease"], columns=["disease_id", "entity_type"])
    mapping = pd.read_parquet(inputs["mapping"], columns=["target_uniprot_id", "ensembl_gene_id"])
    if mapping.duplicated(["target_uniprot_id", "ensembl_gene_id"]).any():
        raise RuntimeError("Protein-gene mapping key is not unique")
    if diseases.duplicated(["disease_id"]).any():
        raise RuntimeError("Disease dictionary key is not unique")
    if core.duplicated(RELATION_KEY).any():
        raise RuntimeError("Core relation key is not unique")

    relations = core.merge(diseases, on="disease_id", how="left", validate="many_to_one")
    relations = relations.join(aggregate, on=KEY, how="left", validate="many_to_one")
    relations[numeric + ["second_core_or_human_genetic_source_count"]] = relations[
        numeric + ["second_core_or_human_genetic_source_count"]
    ].fillna(0).astype("int32")
    relations[strings] = relations[strings].fillna("")
    relations[["weak_sources_only_flag", "unknown_source_flag", "literature_only_flag",
               "expression_plus_animal_only_flag"]] = relations[
        ["weak_sources_only_flag", "unknown_source_flag", "literature_only_flag",
         "expression_plus_animal_only_flag"]
    ].fillna(False).astype(bool)

    disease_direct = (relations["entity_type"] == "disease") & (relations["association_scope"] == "direct")
    relations["rule_A_pass"] = relations["rule_a_qualifying_datasource_count"] >= 2
    relations["rule_B_pass"] = (
        (relations["independent_source_family_count"] >= 2)
        & (relations["core_source_count"] >= 1)
        & ~relations["weak_sources_only_flag"]
        & disease_direct
    )
    relations["rule_C_pass"] = (
        (relations["independent_source_family_count"] >= 2)
        & (relations["clinical_or_curated_genetic_source_count"] >= 1)
        & (relations["second_core_or_human_genetic_source_count"] >= 1)
        & (relations["entity_type"] == "disease")
    )
    relations["rule_D_pass"] = (
        ((relations["independent_source_family_count"] >= 3) & (relations["core_source_count"] >= 2))
        | ((relations["independent_source_family_count"] >= 2) & (relations["high_authority_source_count"] >= 2))
    )

    all_count = len(relations)
    count_rows = []
    for rule in ["A", "B", "C", "D"]:
        retained = relations.loc[relations[f"rule_{rule}_pass"]]
        disease_per_protein = retained.groupby("target_uniprot_id")["disease_id"].nunique()
        tiers = retained["project_evidence_tier"].value_counts()
        count_rows.append({
            "rule": f"Rule {rule}",
            "recommendation": "recommended_default_trial_only" if rule == "B" else "comparison_only",
            "retained_relationship_count": len(retained),
            "retained_protein_count": retained["target_uniprot_id"].nunique(),
            "retained_disease_count": retained["disease_id"].nunique(),
            "median_diseases_per_protein": round(float(disease_per_protein.median()), 4) if len(disease_per_protein) else 0,
            "mean_diseases_per_protein": round(float(disease_per_protein.mean()), 4) if len(disease_per_protein) else 0,
            "tier_A_count": int(tiers.get("Tier A", 0)),
            "tier_B_count": int(tiers.get("Tier B", 0)),
            "tier_C_count": int(tiers.get("Tier C", 0)),
            "tier_D_count": int(tiers.get("Tier D", 0)),
            "clinical_genetic_supported_count": int(retained["clinical_genetic_flag"].sum()),
            "human_genetic_supported_count": int(retained["human_genetic_flag"].sum()),
            "literature_only_count": int(retained["literature_only_flag"].sum()),
            "expression_plus_animal_only_count": int(retained["expression_plus_animal_only_flag"].sum()),
            "unknown_source_count": int(retained["unknown_source_flag"].sum()),
            "percentage_of_all_direct_retained": round(100 * len(retained) / all_count, 6),
        })
    count_fields = list(count_rows[0])
    write_tsv(reports / "evidence_filter_rule_counts.tsv", count_rows, count_fields)

    candidate = relations.loc[
        (relations["independent_source_family_count"] == 1)
        & (relations["high_authority_source_count"] == 1)
        & disease_direct
    ].copy()
    candidate["candidate_status"] = "single_source_high_authority_candidate"
    candidate["validation_status"] = "requires_raw_evidence_validation"
    candidate_fields = RELATION_KEY + [
        "entity_type", "association_scope", "independent_source_families", "raw_datasource_ids",
        "ot_overall_score", "project_evidence_tier", "candidate_status", "validation_status",
    ]
    candidate.sort_values(RELATION_KEY)[candidate_fields].to_csv(
        reports / "single_source_high_confidence_candidates.tsv", sep="\t", index=False, lineterminator="\n"
    )

    rule_columns = [f"rule_{rule}_pass" for rule in "ABCD"]
    combination_rows = []
    for combination, frame in relations.groupby("all_source_families", sort=False):
        row = {
            "source_combination": combination or "NO_EVIDENCE_SUMMARY_SOURCE",
            "relationship_count": len(frame),
            "protein_count": frame["target_uniprot_id"].nunique(),
            "disease_count": frame["disease_id"].nunique(),
        }
        for rule, column in zip("ABCD", rule_columns):
            passed = int(frame[column].sum())
            row[f"rule_{rule}_pass_relationship_count"] = passed
            row[f"rule_{rule}_outcome"] = "all_pass" if passed == len(frame) else "none_pass" if passed == 0 else "mixed"
        combination_rows.append(row)
    combination_rows.sort(key=lambda row: (-row["relationship_count"], row["source_combination"]))
    combination_rows = combination_rows[:50]
    combination_fields = list(combination_rows[0]) if combination_rows else ["source_combination", "relationship_count"]
    write_tsv(reports / "source_combination_frequency.tsv", combination_rows, combination_fields)

    hashes_after = {name: sha256(path) for name, path in inputs.items()}
    if hashes_before != hashes_after:
        raise RuntimeError("One or more formal validated input assets changed during read-only comparison")
    report = build_report(
        policy_file, policy, observed_ids, unknown_ids, mapping_rows, count_rows, combination_rows,
        len(relations), len(candidate), hashes_before,
    )
    (reports / "EVIDENCE_FILTER_RULE_COMPARISON.md").write_text(report, encoding="utf-8")
    return {
        "status": "passed",
        "policy_status": policy["status"],
        "all_direct_relationship_count": len(relations),
        "observed_datasource_count": len(observed_ids),
        "unknown_datasource_count": len(unknown_ids),
        "single_source_high_authority_candidate_count": len(candidate),
        "rule_counts": count_rows,
        "formal_input_hashes_unchanged": True,
        "report": (reports / "EVIDENCE_FILTER_RULE_COMPARISON.md").relative_to(root).as_posix(),
    }


def build_report(policy_file: Path, policy: dict[str, Any], observed_ids: list[str], unknown_ids: list[str],
                 mapping_rows: list[dict[str, Any]], count_rows: list[dict[str, Any]],
                 combinations: list[dict[str, Any]], all_count: int, candidate_count: int,
                 input_hashes: dict[str, str]) -> str:
    source_table = [
        "| datasource_id | family | level | counts | evidence rows |",
        "|---|---|---|:---:|---:|",
    ]
    for row in mapping_rows:
        source_table.append(
            f"| {row['datasource_id']} | {row['source_family']} | {row['source_level']} | "
            f"{str(row['counts_toward_two_source_rule']).lower()} | {row['evidence_row_count']:,} |"
        )
    rule_table = [
        "| rule | relationships | proteins | diseases | retained | Tier A/B/C/D | clinical | human genetic | literature only | expression+animal only | unknown |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in count_rows:
        rule_table.append(
            f"| {row['rule']} | {row['retained_relationship_count']:,} | {row['retained_protein_count']:,} | "
            f"{row['retained_disease_count']:,} | {row['percentage_of_all_direct_retained']:.3f}% | "
            f"{row['tier_A_count']:,}/{row['tier_B_count']:,}/{row['tier_C_count']:,}/{row['tier_D_count']:,} | "
            f"{row['clinical_genetic_supported_count']:,} | {row['human_genetic_supported_count']:,} | "
            f"{row['literature_only_count']:,} | {row['expression_plus_animal_only_count']:,} | {row['unknown_source_count']:,} |"
        )
    combo_table = [
        "| source families | relationships | Rule A | Rule B | Rule C | Rule D |",
        "|---|---:|---|---|---|---|",
    ]
    for row in combinations:
        combo_table.append(
            f"| {row['source_combination'].replace(';', ' + ')} | {row['relationship_count']:,} | "
            f"{row['rule_A_outcome']} | {row['rule_B_outcome']} | {row['rule_C_outcome']} | {row['rule_D_outcome']} |"
        )
    hashes = "\n".join(f"- `{name}`: `{digest}`" for name, digest in input_hashes.items())
    unknown_text = ", ".join(unknown_ids) if unknown_ids else "None"
    return f'''# Evidence filter rule comparison

Status: trial comparison only. No formal filtered dataset was created or replaced.

## Inputs and safety

- All-direct relationships: {all_count:,}
- Observed datasource IDs: {len(observed_ids):,}
- Unmapped datasource IDs: {len(unknown_ids):,} ({unknown_text})
- Policy: `{policy_file.as_posix()}` (`{policy['policy_version']}`)
- The four formal Parquet input hashes were identical before and after this read-only analysis.

{hashes}

## Actual datasource mapping

{chr(10).join(source_table)}

`eva` is retained as family `EVA`, rather than being guessed to mean ClinVar. `crispr` and `crispr_screen` count once as `CRISPR`; `eva` and `eva_somatic` count once as `EVA`; both UniProt channels count once as `UNIPROT`. Unknown sources never count automatically.

## Operational rule definitions

- **Rule A — simple two-source:** at least two raw datasource IDs after excluding only `clinical_precedence`. Europe PMC therefore counts in Rule A; unknown sources do not.
- **Rule B — balanced two-source (recommended_default, trial only):** at least two qualifying independent families, at least one core/core-curated/core-cancer family, not weak-only, `entity_type=disease`, and direct scope. Europe PMC and Clinical Precedence do not count; Expression Atlas and IMPC may only supplement a core family.
- **Rule C — strict two-source:** at least two qualifying independent families, at least one clinical/curated genetic family, and at least one *different* Rule-C-eligible core/human-genetic family, with `entity_type=disease`. Functional screens, Expression Atlas, IMPC, Europe PMC and Clinical Precedence do not satisfy either source slot.
- **Rule D — historical multi-source trial:** either at least three qualifying families with at least two core families, or at least two qualifying families that are both marked high-authority. It lacks the final `entity_type=disease` and specificity correction, so it is retained only for historical comparison and is not a formal evidence level. Corrected Rule E is computed separately by the A/B/C/E exporter.

The `high_authority_source` flag is a source-level trial classification. It is not a Strong/Definitive/Green/expert-panel assertion about any individual record.

## Comparison

{chr(10).join(rule_table)}

The full TSV also reports median and mean diseases per retained protein. Rule B is pre-marked `recommended_default_trial_only` because it requires two independent families and a human/curated/core anchor, preventing text-mining, expression, animal-model, functional-screen or drug-derived evidence from qualifying on their own. It has not been promoted to a formal rule.

## Single high-authority source candidates

There are {candidate_count:,} disease/direct relationships with one qualifying independent family that is marked high-authority at source level. They are kept separately in `single_source_high_confidence_candidates.tsv`, labelled `single_source_high_authority_candidate` and `requires_raw_evidence_validation`; they are not included in Rule B by exception.

## Top 50 source-family combinations

`mixed` means the same source-family combination spans both disease and non-disease entities (or another relation-level condition), so not every relationship has the same rule outcome.

{chr(10).join(combo_table)}

## Outputs

- `reports/evidence_filter_rule_counts.tsv`
- `reports/source_family_mapping.tsv`
- `reports/unmapped_datasources.tsv`
- `reports/source_combination_frequency.tsv`
- `reports/single_source_high_confidence_candidates.tsv`

No `trial_rule_*` relationship exports were generated at this stage; the comparison is statistics-only apart from the explicitly requested single-source candidate review table.
'''
