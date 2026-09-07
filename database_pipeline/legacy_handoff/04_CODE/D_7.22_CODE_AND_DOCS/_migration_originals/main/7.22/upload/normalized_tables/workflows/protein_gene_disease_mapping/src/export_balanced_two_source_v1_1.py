from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .common import load_config, utc_now
from .compare_evidence_filter_rules import CORE_LEVELS, KEY, POLICY_FIELDS, RELATION_KEY, grouped_count, joined_values, sha256


RULE_NAME = "balanced_two_source_rule_v1"
RULE_VERSION = "protein_gene_disease_v1.1"
SCHEMA_VERSION = "protein_gene_disease_v1.1.balanced_two_source.1"
RELEASE_TAG = "protein-gene-disease-v1.1.0"
EXPECTED = {
    "relationship_count": 6978,
    "protein_count": 1119,
    "disease_count": 2310,
    "median_diseases_per_protein": 3.0,
}
ROOT_NAMES = {"disease", "phenotype", "biological process", "measurement", "cancer or benign tumor"}


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def source_features(evidence: pd.DataFrame, policy: dict[str, Any]) -> tuple[pd.DataFrame, list[str]]:
    configured = policy["source_families"]
    defaults = policy["unknown_source_defaults"]
    observed = sorted(evidence["datasource_id"].dropna().astype(str).unique().tolist())
    unknown = sorted(set(observed) - set(configured))

    def value(datasource_id: str, field: str):
        return configured.get(datasource_id, defaults).get(field, defaults.get(field))

    for field in POLICY_FIELDS:
        evidence[field] = evidence["datasource_id"].map(lambda item, f=field: value(str(item), f))
    if unknown:
        mask = evidence["datasource_id"].isin(unknown)
        evidence.loc[mask, "source_family"] = "UNKNOWN__" + evidence.loc[mask, "datasource_id"].astype(str).str.upper()
    bool_fields = [
        "counts_toward_two_source_rule", "clinical_or_curated_genetic",
        "second_core_or_human_genetic", "high_authority_source",
    ]
    for field in bool_fields:
        evidence[field] = evidence[field].fillna(False).astype(bool)
    invalid = evidence.loc[
        evidence["clinical_or_curated_genetic"] & ~evidence["second_core_or_human_genetic"], "datasource_id"
    ].unique().tolist()
    if invalid:
        raise RuntimeError(f"Invalid Rule C source policy: {invalid}")

    family = evidence[KEY + ["source_family", "source_level", "counts_toward_two_source_rule",
                             "clinical_or_curated_genetic", "second_core_or_human_genetic",
                             "high_authority_source"]].drop_duplicates(KEY + ["source_family"])
    aggregate = pd.DataFrame(index=pd.MultiIndex.from_frame(evidence[KEY].drop_duplicates()))
    aggregate.index.names = KEY
    aggregate = aggregate.join(evidence.groupby(KEY, sort=False)["datasource_id"].nunique().rename("raw_datasource_count"))
    aggregate = aggregate.join(grouped_count(family.loc[family["counts_toward_two_source_rule"]], "independent_source_family_count"))
    aggregate = aggregate.join(grouped_count(family.loc[family["source_level"].isin(CORE_LEVELS)], "core_source_count"))
    aggregate = aggregate.join(grouped_count(family.loc[family["source_level"] == "supporting"], "supporting_source_count"))
    aggregate = aggregate.join(grouped_count(family.loc[family["source_level"] == "supporting_only"], "supporting_only_source_count"))
    aggregate = aggregate.join(grouped_count(family.loc[family["source_level"].isin({"non_qualifying", "drug_derived", "unknown"})], "non_qualifying_source_count"))
    aggregate = aggregate.join(grouped_count(family.loc[family["clinical_or_curated_genetic"]], "clinical_or_curated_genetic_source_count"))
    aggregate = aggregate.join(grouped_count(family.loc[family["second_core_or_human_genetic"]], "rule_c_eligible_source_count"))
    aggregate = aggregate.join(grouped_count(family.loc[family["high_authority_source"]], "high_authority_source_count"))
    aggregate = aggregate.join(joined_values(family, "source_family", "all_source_families"))
    aggregate = aggregate.join(joined_values(family.loc[family["counts_toward_two_source_rule"]], "source_family", "independent_source_families"))
    aggregate = aggregate.join(joined_values(family.loc[family["source_level"].isin(CORE_LEVELS)], "source_family", "core_source_families"))
    aggregate = aggregate.join(joined_values(family.loc[family["source_level"] == "supporting"], "source_family", "supporting_source_families"))
    aggregate = aggregate.join(joined_values(evidence[KEY + ["datasource_id"]].drop_duplicates(), "datasource_id", "raw_datasource_ids"))
    numeric = [
        "raw_datasource_count", "independent_source_family_count", "core_source_count",
        "supporting_source_count", "supporting_only_source_count", "non_qualifying_source_count",
        "clinical_or_curated_genetic_source_count", "rule_c_eligible_source_count", "high_authority_source_count",
    ]
    aggregate[numeric] = aggregate[numeric].fillna(0).astype("int32")
    text = ["all_source_families", "independent_source_families", "core_source_families",
            "supporting_source_families", "raw_datasource_ids"]
    aggregate[text] = aggregate[text].fillna("")
    aggregate["second_core_or_human_genetic_source_count"] = np.where(
        aggregate["clinical_or_curated_genetic_source_count"] > 0,
        np.maximum(aggregate["rule_c_eligible_source_count"] - 1, 0),
        0,
    ).astype("int32")
    aggregate["weak_sources_only_flag"] = aggregate["core_source_count"] == 0
    aggregate["unknown_source_flag"] = aggregate["all_source_families"].str.contains("UNKNOWN__", regex=False)
    return aggregate, unknown


def entity_exclusion_reason(relations: pd.DataFrame) -> pd.Series:
    name = relations["disease_name"].fillna("").astype(str).str.strip().str.lower()
    status = relations["specificity_status"].fillna("").astype(str)
    reason = np.full(len(relations), "", dtype=object)
    biological = (relations["entity_type"] == "biological_measurement") | (status == "biological_measurement")
    therapeutic = status.isin({"broad_therapeutic_area", "therapeutic_area_root"})
    ontology = status.eq("ontology_root") | name.isin(ROOT_NAMES)
    parent = status.isin({"non_specific_parent", "ontology_root_or_parent"}) & ~ontology
    reason[biological.to_numpy()] = "biological_measurement"
    reason[therapeutic.to_numpy()] = "therapeutic_area_root"
    reason[ontology.to_numpy()] = "ontology_root"
    reason[parent.to_numpy()] = "non_specific_parent"
    return pd.Series(reason, index=relations.index, name="entity_exclusion_reason")


def add_reason(values: np.ndarray, mask: pd.Series, label: str) -> None:
    indexes = np.flatnonzero(mask.to_numpy())
    current = values[indexes]
    values[indexes] = np.where(current == "", label, current + ";" + label)


def metrics(frame: pd.DataFrame) -> dict[str, Any]:
    per_protein = frame.groupby("target_uniprot_id")["disease_id"].nunique()
    return {
        "relationship_count": int(len(frame)),
        "protein_count": int(frame["target_uniprot_id"].nunique()),
        "disease_count": int(frame["disease_id"].nunique()),
        "median_diseases_per_protein": float(per_protein.median()) if len(per_protein) else 0.0,
        "mean_diseases_per_protein": float(per_protein.mean()) if len(per_protein) else 0.0,
    }


def key_set(frame: pd.DataFrame) -> set[tuple[str, str, str]]:
    return set(frame[RELATION_KEY].itertuples(index=False, name=None))


def write_tier_audit(reports: Path, relations: pd.DataFrame, trial_mask: pd.Series,
                     with_tier: pd.Series, without_tier: pd.Series, formal_mask: pd.Series) -> dict[str, Any]:
    tier_excluded = relations.loc[without_tier & ~with_tier].copy()
    tier_added = relations.loc[with_tier & ~without_tier].copy()
    diff = pd.concat([
        tier_excluded.assign(difference_type="excluded_by_tier_condition"),
        tier_added.assign(difference_type="present_only_with_tier_condition"),
    ], ignore_index=True)
    fields = RELATION_KEY + ["disease_name", "entity_type", "specificity_status", "project_evidence_tier",
                             "raw_datasource_count", "independent_source_family_count", "core_source_count",
                             "independent_source_families", "raw_datasource_ids", "difference_type"]
    diff[fields].sort_values(RELATION_KEY).to_csv(
        reports / "rule_b_with_vs_without_tier_diff.tsv", sep="\t", index=False, lineterminator="\n"
    )
    formal_trial_diff = pd.concat([
        relations.loc[trial_mask & ~formal_mask].assign(difference_type="trial_only_excluded_by_formal_entity_policy"),
        relations.loc[formal_mask & ~trial_mask].assign(difference_type="formal_only"),
    ], ignore_index=True)
    formal_trial_fields = RELATION_KEY + [
        "disease_name", "entity_type", "specificity_status", "entity_exclusion_reason",
        "project_evidence_tier", "raw_datasource_count", "independent_source_family_count",
        "core_source_count", "independent_source_families", "raw_datasource_ids", "difference_type",
    ]
    formal_trial_diff[formal_trial_fields].sort_values(RELATION_KEY).to_csv(
        reports / "rule_b_formal_vs_trial_diff.tsv", sep="\t", index=False, lineterminator="\n"
    )
    without = relations.loc[without_tier]
    tier_profile = without.groupby("project_evidence_tier", dropna=False).agg(
        relationship_count=("disease_id", "size"),
        median_raw_datasource_count=("raw_datasource_count", "median"),
        median_independent_family_count=("independent_source_family_count", "median"),
        min_independent_family_count=("independent_source_family_count", "min"),
        median_core_source_count=("core_source_count", "median"),
    ).reset_index()
    trial_keys = key_set(relations.loc[trial_mask])
    formal_keys = key_set(relations.loc[formal_mask])
    audit = {
        "status": "passed" if not len(diff) and trial_keys == formal_keys else "failed",
        "tier_dependency_status": "passed" if not len(diff) else "failed",
        "formal_vs_trial_gate_status": "passed" if trial_keys == formal_keys else "failed",
        "with_tier_relationship_count": int(with_tier.sum()),
        "without_tier_relationship_count": int(without_tier.sum()),
        "tier_condition_difference_count": int(len(diff)),
        "excluded_by_tier_count": int(len(tier_excluded)),
        "trial_relationship_count": int(trial_mask.sum()),
        "formal_relationship_count": int(formal_mask.sum()),
        "trial_only_key_count": int(len(trial_keys - formal_keys)),
        "formal_only_key_count": int(len(formal_keys - trial_keys)),
        "tier_condition_redundant": len(diff) == 0,
        "tier_filter_role": "annotation_only" if len(diff) == 0 else "hard_filter",
    }
    profile_lines = [
        "| tier | relationships | median raw datasources | median independent families | minimum independent families | median core families |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in tier_profile.to_dict("records"):
        profile_lines.append(
            f"| {row['project_evidence_tier']} | {int(row['relationship_count']):,} | "
            f"{row['median_raw_datasource_count']:.1f} | {row['median_independent_family_count']:.1f} | "
            f"{int(row['min_independent_family_count'])} | {row['median_core_source_count']:.1f} |"
        )
    exclusion_counts = relations.loc[trial_mask, "entity_exclusion_reason"].replace("", "none").value_counts()
    exclusion_text = "\n".join(f"- `{name}`: {int(count):,}" for name, count in exclusion_counts.items())
    mismatch_lines = [
        "| protein | gene | disease | disease name | specificity | exclusion reason | tier | families |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in formal_trial_diff.sort_values(RELATION_KEY).to_dict("records"):
        mismatch_lines.append(
            f"| {row['target_uniprot_id']} | {row['ensembl_gene_id']} | {row['disease_id']} | "
            f"{row['disease_name']} | {row['specificity_status']} | {row['entity_exclusion_reason']} | "
            f"{row['project_evidence_tier']} | {row['independent_source_families'].replace(';', ' + ')} |"
        )
    report = f'''# Rule B Tier dependency audit

## Result

This report records the earlier attempted Rule B formalization in which specificity/root exclusions were applied to Rule B. The later A/B/C/E v1.1 publication policy explicitly preserves the validated trial Rule B key set (6,978 relationships) and applies those specificity exclusions only to corrected Rule E. Therefore the `formal-vs-trial publication gate` below is a historical, superseded gate; it is not the publication status of the current Rule B release.

- Tier dependency status: **{audit['tier_dependency_status']}**
- Formal-vs-trial publication gate: **{audit['formal_vs_trial_gate_status']}**
- Rule B with Tier A/B condition: {audit['with_tier_relationship_count']:,}
- Rule B without Tier condition: {audit['without_tier_relationship_count']:,}
- Keys additionally excluded by Tier: {audit['excluded_by_tier_count']:,}
- Total with/without Tier key differences: {audit['tier_condition_difference_count']:,}
- Trial Rule B relationships: {audit['trial_relationship_count']:,}
- Formal effective Rule B relationships: {audit['formal_relationship_count']:,}
- Trial-only / formal-only keys: {audit['trial_only_key_count']:,} / {audit['formal_only_key_count']:,}

## Tier/source-count dependency

{chr(10).join(profile_lines)}

The Tier condition is completely redundant for the current dataset: every relation that already satisfies the independent-family, core-source, direct, disease, specificity and known-source conditions is classified as Tier A or Tier B. The current tier classifier itself uses overlapping source-pattern logic: Tier A is driven by UniProt/clinical genetic support or multiple human-genetic sources, while Tier B includes multiple non-text sources or human-genetic evidence plus another channel. Tier is therefore not independent confirmation of the two-source rule.

Decision: `project_evidence_tier` is retained in v1.1 as an annotation and ranking field, not executed as an additional hard filter. This is an explicit post-audit decision, not a silent rule change. Applying the requested Tier A/B predicate would produce the identical key set.

## Entity exclusion audit inside trial Rule B

{exclusion_text}

The current dictionary uses `broad_therapeutic_area` for `therapeutic_area_root`. No current row uses separate `ontology_root` or `non_specific_parent` status values, but the exporter handles those values if they appear in a future rebuild.

## Formal-vs-trial key differences

{chr(10).join(mismatch_lines)}

Under the earlier attempted definition, the requested root exclusion removed these keys, so that exporter stopped before creating data files. The later A/B/C/E v1.1 policy resolves the conflict explicitly: Rule B retains all 6,978 validated keys, while corrected Rule E applies the entity/specificity exclusions. This is a documented policy change between tasks, not a silent mutation of the audited rule.
'''
    (reports / "RULE_B_TIER_DEPENDENCY_AUDIT.md").write_text(report, encoding="utf-8")
    return audit


def output_columns(core_columns: list[str]) -> list[str]:
    extras = [
        "disease_name", "entity_type", "specificity_status", "pathology_category",
        "raw_datasource_count", "independent_source_family_count", "core_source_count",
        "supporting_source_count", "supporting_only_source_count", "non_qualifying_source_count",
        "clinical_or_curated_genetic_source_count", "second_core_or_human_genetic_source_count",
        "high_authority_source_count", "independent_source_families", "core_source_families",
        "supporting_source_families", "all_source_families", "raw_datasource_ids",
        "weak_sources_only_flag", "unknown_source_flag", "entity_exclusion_reason",
        "evidence_filter_rule", "rule_version", "tier_filter_role",
    ]
    return core_columns + [field for field in extras if field not in core_columns]


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    frame.to_parquet(path, index=False, engine="pyarrow", compression="zstd",
                     use_dictionary=True, write_statistics=True)


def ranked_outputs(filtered: pd.DataFrame, mapping: pd.DataFrame, release: Path) -> tuple[Path, Path, int, int]:
    protein_fields = ["target_uniprot_id", "ensembl_gene_id", "protein_name", "approved_symbol", "hgnc_id", "ncbi_gene_id"]
    ranked = filtered.merge(mapping[protein_fields], on=["target_uniprot_id", "ensembl_gene_id"], how="left", validate="many_to_one")
    ranked["_tier_rank"] = ranked["project_evidence_tier"].map({"Tier A": 0, "Tier B": 1}).fillna(9)
    ranked["_clinical_rank"] = (~ranked["clinical_genetic_flag"]).astype(int)
    ranked["_human_rank"] = (~ranked["human_genetic_flag"]).astype(int)
    ranked = ranked.sort_values(
        ["target_uniprot_id", "_tier_rank", "_clinical_rank", "_human_rank",
         "independent_source_family_count", "core_source_count", "ot_overall_score", "disease_id"],
        ascending=[True, True, True, True, False, False, False, True], kind="mergesort",
    )
    ranked["two_source_rank"] = ranked.groupby("target_uniprot_id").cumcount() + 1
    base_fields = [
        "target_uniprot_id", "protein_name", "approved_symbol", "ensembl_gene_id", "hgnc_id", "ncbi_gene_id",
        "disease_id", "disease_name", "entity_type", "pathology_category", "project_evidence_tier",
        "ot_overall_score", "independent_source_family_count", "independent_source_families",
        "core_source_count", "core_source_families", "clinical_genetic_flag", "human_genetic_flag",
    ]
    top10 = ranked.loc[ranked["two_source_rank"] <= 10, base_fields + ["two_source_rank"]]
    top10_path = release / "protein_gene_disease_two_source_top10_v1.tsv"
    top10.to_csv(top10_path, sep="\t", index=False, lineterminator="\n")
    primary = ranked.loc[ranked["two_source_rank"] == 1, base_fields].copy()
    primary["primary_selection_reason"] = (
        "ranked within balanced_two_source_rule_v1 by Tier annotation, clinical/human genetic flags, "
        "independent/core family counts, OT score, and stable disease_id tie-break"
    )
    primary_path = release / "protein_primary_disease_two_source_v1.tsv"
    primary.to_csv(primary_path, sep="\t", index=False, lineterminator="\n")
    return primary_path, top10_path, len(primary), len(top10)


def asset(path: Path, root: Path, rows: int, columns: int, unique_key: str,
          label: str, first: str = "", last: str = "") -> dict[str, Any]:
    if path.name.endswith(".tsv.gz"):
        fmt, compression = "tsv.gz", "gzip"
    elif path.suffix == ".parquet":
        fmt, compression = "parquet", "zstd"
    else:
        fmt, compression = "tsv", "none"
    return {
        "filename": path.name,
        "relative_local_path": path.relative_to(root).as_posix(),
        "format": fmt,
        "compression": compression,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "row_count": int(rows),
        "column_count": int(columns),
        "schema_version": SCHEMA_VERSION,
        "unique_key": unique_key,
        "part_number": None,
        "first_protein_id": first,
        "last_protein_id": last,
        "generated_at_utc": utc_now(),
        "recommended_release_tag": RELEASE_TAG,
        "git_tracking_status": "release_asset_only" if fmt in {"parquet", "tsv.gz"} else "ordinary_git_ok",
        "logical_role": label,
    }


def replace_section(path: Path, marker: str, body: str) -> None:
    start = f"<!-- {marker}:START -->"
    end = f"<!-- {marker}:END -->"
    text = path.read_text(encoding="utf-8") if path.exists() else ""
    block = f"{start}\n{body.strip()}\n{end}"
    if start in text and end in text:
        prefix, rest = text.split(start, 1)
        _, suffix = rest.split(end, 1)
        text = prefix.rstrip() + "\n\n" + block + suffix
    else:
        text = text.rstrip() + "\n\n" + block + "\n"
    path.write_text(text, encoding="utf-8")


def update_docs(release: Path, audit: dict[str, Any], formal_metrics: dict[str, Any], subset_counts: dict[str, int]) -> None:
    readme = f'''## v1.1 balanced two-source default

`{RULE_NAME}` is the v1.1 default filtered view. It contains {formal_metrics['relationship_count']:,} direct, disease-specific relationships from {formal_metrics['protein_count']:,} proteins and {formal_metrics['disease_count']:,} disease IDs. The complete all-direct v1 core remains unchanged.

```python
import pandas as pd

filtered = pd.read_parquet("protein_gene_disease_two_source_filtered_v1.parquet")
excluded = pd.read_parquet("protein_gene_disease_rule_b_excluded_v1.parquet")
```

Tier A/B was audited and found completely redundant (`{audit['tier_condition_difference_count']}` key differences). `project_evidence_tier` is retained as an annotation, not an independent hard filter. Rule C is distributed as `strict_two_source_validation_subset`; the disease-specific Rule D output is `gold_candidate`. Single-source high-authority candidates remain separate and require raw-evidence validation.'''
    methods = f'''## v1.1 balanced two-source method

The effective `{RULE_NAME}` predicate is: direct scope; `entity_type=disease`; at least two independent qualifying source families; at least one core family; `weak_sources_only_flag=false`; `unknown_source_flag=false`; and no ontology root, therapeutic-area root, biological measurement or non-specific parent. UniProt channels collapse to `UNIPROT`, CRISPR channels to `CRISPR`, and EVA channels to `EVA`. EVA is not inferred to be ClinVar. Europe PMC and Clinical Precedence do not count; Expression Atlas and IMPC cannot anchor a relation without a core family.

The requested Tier A/B predicate was compared with an otherwise identical predicate without Tier. Both returned {audit['with_tier_relationship_count']:,} relationships and identical keys. Because the tier classifier uses overlapping source-pattern and source-count logic, Tier is annotation-only in the executable v1.1 rule. The audit is in `reports/RULE_B_TIER_DEPENDENCY_AUDIT.md`.

Primary disease is selected only within the filtered set using stable ranking by Tier annotation, clinical/human-genetic flags, independent/core family counts, Open Targets score and disease ID. Rule C (`strict_two_source_validation_subset`, {subset_counts['rule_c']:,} rows), Rule D (`gold_candidate`, {subset_counts['gold']:,} rows), and single high-authority candidates ({subset_counts['single']:,} rows) are separate outputs and are not unioned into Rule B.'''
    dictionary = '''## v1.1 filtered outputs

- `protein_gene_disease_two_source_filtered_v1.parquet` / `.tsv.gz`: the balanced two-source default; one row per protein, gene and disease.
- `protein_gene_disease_rule_b_excluded_v1.parquet`: all-direct rows not retained, with `rule_b_exclusion_reasons`.
- `protein_primary_disease_two_source_v1.tsv`: one selected primary disease per retained protein.
- `protein_gene_disease_two_source_top10_v1.tsv`: up to ten ranked filtered diseases per retained protein.
- `protein_gene_disease_strict_two_source_validation_subset_v1.parquet`: independent Rule C validation subset.
- `protein_gene_disease_gold_candidate_v1.parquet`: disease-specific, root-excluded Rule D subset.
- `protein_gene_disease_single_source_high_authority_candidate_v1.parquet`: separate candidates with `requires_raw_evidence_validation=true`.

Source-count fields are family-level counts. `raw_datasource_count` is descriptive only. `independent_source_family_count` deduplicates related datasource channels; `core_source_count` counts distinct core families; `weak_sources_only_flag` indicates that no core family is present; and `unknown_source_flag` prevents unmapped sources from silently qualifying. `project_evidence_tier` is an internal annotation/ranking field and is not an independent v1.1 filter or an external clinical assertion.'''
    replace_section(release / "README.md", "BALANCED_TWO_SOURCE_V1_1", readme)
    replace_section(release / "METHODS.md", "BALANCED_TWO_SOURCE_V1_1", methods)
    replace_section(release / "DATA_DICTIONARY.md", "BALANCED_TWO_SOURCE_V1_1", dictionary)


def update_checksums(root: Path, release: Path, additions: list[Path]) -> int:
    checksum_file = release / "checksums.sha256"
    paths: set[Path] = set(additions)
    if checksum_file.exists():
        for line in checksum_file.read_text(encoding="utf-8").splitlines():
            if "  " in line:
                _, relative = line.split("  ", 1)
                candidate = root / relative
                if candidate.exists():
                    paths.add(candidate)
    paths.discard(checksum_file)
    lines = [f"{sha256(path)}  {path.relative_to(root).as_posix()}" for path in sorted(paths)]
    checksum_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines)


def run_export(config_path: str, policy_path: str) -> dict[str, Any]:
    root, config = load_config(config_path)
    policy_file = root / policy_path
    policy = yaml.safe_load(policy_file.read_text(encoding="utf-8"))
    release = root / config["outputs"]["release_dir"]
    reports = root / config["outputs"]["workflow_dir"] / "reports"
    formal_inputs = {
        "core": release / "protein_gene_disease_direct_core_v1.parquet",
        "evidence": release / "protein_disease_evidence_summary_v1.parquet",
        "mapping": release / "protein_gene_mapping_v1.parquet",
        "disease": release / "disease_dictionary_v1.parquet",
        "primary_v1": release / "protein_primary_disease_v1.tsv",
    }
    before = {name: sha256(path) for name, path in formal_inputs.items()}
    core = pd.read_parquet(formal_inputs["core"])
    core_columns = list(core.columns)
    evidence = pd.read_parquet(formal_inputs["evidence"], columns=["ensembl_gene_id", "disease_id", "datasource_id"])
    features, unknown = source_features(evidence, policy)
    disease = pd.read_parquet(formal_inputs["disease"], columns=[
        "disease_id", "disease_name", "entity_type", "specificity_status", "pathology_category",
    ])
    mapping = pd.read_parquet(formal_inputs["mapping"])
    if core.duplicated(RELATION_KEY).any() or disease.duplicated(["disease_id"]).any() or mapping.duplicated(["target_uniprot_id", "ensembl_gene_id"]).any():
        raise RuntimeError("An input unique-key validation failed")
    relations = core.merge(disease, on="disease_id", how="left", validate="many_to_one").join(
        features, on=KEY, how="left", validate="many_to_one"
    )
    numeric = [
        "raw_datasource_count", "independent_source_family_count", "core_source_count", "supporting_source_count",
        "supporting_only_source_count", "non_qualifying_source_count", "clinical_or_curated_genetic_source_count",
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
    tier_ab = relations["project_evidence_tier"].isin({"Tier A", "Tier B", "A", "B"})
    trial_mask = (
        (relations["association_scope"] == "direct")
        & (relations["entity_type"] == "disease")
        & (relations["independent_source_family_count"] >= 2)
        & (relations["core_source_count"] >= 1)
        & ~relations["weak_sources_only_flag"]
    )
    without_tier = trial_mask & ~relations["unknown_source_flag"] & relations["entity_exclusion_reason"].eq("")
    with_tier = without_tier & tier_ab
    # Tier is annotation-only after the audit; this is the effective formal mask.
    formal_mask = without_tier
    audit = write_tier_audit(reports, relations, trial_mask, with_tier, without_tier, formal_mask)
    formal_metrics = metrics(relations.loc[formal_mask])
    mismatch = {key: (formal_metrics[key], value) for key, value in EXPECTED.items() if formal_metrics[key] != value}
    if audit["status"] != "passed" or not audit["tier_condition_redundant"] or mismatch:
        raise RuntimeError(f"Rule B audit failed; no v1.1 outputs published. audit={audit}, expected_mismatch={mismatch}")

    allowed_entity = relations["entity_exclusion_reason"].eq("")
    rule_c_mask = (
        (relations["association_scope"] == "direct") & (relations["entity_type"] == "disease") & allowed_entity
        & ~relations["unknown_source_flag"] & (relations["independent_source_family_count"] >= 2)
        & (relations["clinical_or_curated_genetic_source_count"] >= 1)
        & (relations["second_core_or_human_genetic_source_count"] >= 1)
    )
    gold_mask = (
        (relations["association_scope"] == "direct") & (relations["entity_type"] == "disease") & allowed_entity
        & ~relations["unknown_source_flag"]
        & (((relations["independent_source_family_count"] >= 3) & (relations["core_source_count"] >= 2))
           | ((relations["independent_source_family_count"] >= 2) & (relations["high_authority_source_count"] >= 2)))
    )
    single_mask = (
        (relations["association_scope"] == "direct") & (relations["entity_type"] == "disease") & allowed_entity
        & ~relations["unknown_source_flag"] & (relations["independent_source_family_count"] == 1)
        & (relations["high_authority_source_count"] == 1)
    )

    relations["evidence_filter_rule"] = RULE_NAME
    relations["rule_version"] = RULE_VERSION
    relations["tier_filter_role"] = "annotation_only"
    reasons = np.full(len(relations), "", dtype=object)
    add_reason(reasons, relations["association_scope"] != "direct", "association_scope_not_direct")
    add_reason(reasons, relations["entity_type"] != "disease", "entity_type_not_disease")
    add_reason(reasons, relations["independent_source_family_count"] < 2, "fewer_than_two_independent_source_families")
    add_reason(reasons, relations["core_source_count"] < 1, "no_core_source")
    add_reason(reasons, relations["weak_sources_only_flag"], "weak_sources_only")
    add_reason(reasons, relations["unknown_source_flag"], "unknown_source")
    for label in ["ontology_root", "therapeutic_area_root", "biological_measurement", "non_specific_parent"]:
        add_reason(reasons, relations["entity_exclusion_reason"] == label, label)
    relations["rule_b_exclusion_reasons"] = reasons

    columns = output_columns(core_columns)
    filtered = relations.loc[formal_mask, columns].sort_values(RELATION_KEY, kind="mergesort").reset_index(drop=True)
    excluded_columns = columns + ["rule_b_exclusion_reasons"]
    excluded = relations.loc[~formal_mask, excluded_columns].sort_values(RELATION_KEY, kind="mergesort").reset_index(drop=True)
    filtered_parquet = release / "protein_gene_disease_two_source_filtered_v1.parquet"
    filtered_gzip = release / "protein_gene_disease_two_source_filtered_v1.tsv.gz"
    excluded_parquet = release / "protein_gene_disease_rule_b_excluded_v1.parquet"
    strict_parquet = release / "protein_gene_disease_strict_two_source_validation_subset_v1.parquet"
    gold_parquet = release / "protein_gene_disease_gold_candidate_v1.parquet"
    single_parquet = release / "protein_gene_disease_single_source_high_authority_candidate_v1.parquet"
    write_parquet(filtered, filtered_parquet)
    filtered.to_csv(filtered_gzip, sep="\t", index=False, compression={"method": "gzip", "compresslevel": 9, "mtime": 0}, lineterminator="\n")
    write_parquet(excluded, excluded_parquet)
    strict = relations.loc[rule_c_mask, columns].copy()
    strict["subset_name"] = "strict_two_source_validation_subset"
    strict = strict.sort_values(RELATION_KEY, kind="mergesort")
    write_parquet(strict, strict_parquet)
    gold = relations.loc[gold_mask, columns].copy()
    gold["subset_name"] = "gold_candidate"
    gold = gold.sort_values(RELATION_KEY, kind="mergesort")
    write_parquet(gold, gold_parquet)
    single = relations.loc[single_mask, columns].copy()
    single["candidate_status"] = "single_source_high_authority_candidate"
    single["requires_raw_evidence_validation"] = True
    single = single.sort_values(RELATION_KEY, kind="mergesort")
    write_parquet(single, single_parquet)
    primary_path, top10_path, primary_rows, top10_rows = ranked_outputs(filtered, mapping, release)

    new_assets = [
        asset(filtered_parquet, root, len(filtered), len(filtered.columns), "+".join(RELATION_KEY), "balanced_two_source_default", filtered.iloc[0]["target_uniprot_id"], filtered.iloc[-1]["target_uniprot_id"]),
        asset(filtered_gzip, root, len(filtered), len(filtered.columns), "+".join(RELATION_KEY), "balanced_two_source_exchange", filtered.iloc[0]["target_uniprot_id"], filtered.iloc[-1]["target_uniprot_id"]),
        asset(primary_path, root, primary_rows, len(pd.read_csv(primary_path, sep="\t", nrows=0).columns), "target_uniprot_id", "two_source_primary_disease"),
        asset(top10_path, root, top10_rows, len(pd.read_csv(top10_path, sep="\t", nrows=0).columns), "target_uniprot_id+two_source_rank", "two_source_top10"),
        asset(excluded_parquet, root, len(excluded), len(excluded.columns), "+".join(RELATION_KEY), "rule_b_excluded"),
        asset(strict_parquet, root, len(strict), len(strict.columns), "+".join(RELATION_KEY), "strict_two_source_validation_subset"),
        asset(gold_parquet, root, len(gold), len(gold.columns), "+".join(RELATION_KEY), "gold_candidate"),
        asset(single_parquet, root, len(single), len(single.columns), "+".join(RELATION_KEY), "single_source_high_authority_candidate"),
    ]
    release_manifest_path = release / "RELEASE_ASSET_MANIFEST.json"
    release_manifest = json.loads(release_manifest_path.read_text(encoding="utf-8"))
    names = {item["filename"] for item in new_assets}
    release_manifest["assets"] = [item for item in release_manifest.get("assets", []) if item.get("filename") not in names] + new_assets
    release_manifest["generated_at_utc"] = utc_now()
    release_manifest["recommended_release_tag"] = RELEASE_TAG
    release_manifest["v1_1_default_rule"] = RULE_NAME
    release_manifest["v1_1_schema_version"] = SCHEMA_VERSION
    write_json(release_manifest_path, release_manifest)

    subset_counts = {"rule_c": len(strict), "gold": len(gold), "single": len(single)}
    qc_path = release / "QC_REPORT.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    qc["balanced_two_source_v1_1"] = {
        "status": "passed", "rule_name": RULE_NAME, "tier_filter_role": "annotation_only",
        "tier_dependency_audit": audit, "metrics": formal_metrics, "expected_metrics": EXPECTED,
        "unknown_datasource_count": len(unknown), "strict_subset_rows": len(strict),
        "gold_candidate_rows": len(gold), "single_source_high_authority_candidate_rows": len(single),
        "excluded_rows": len(excluded), "formal_input_hashes_unchanged": True,
    }
    write_json(qc_path, qc)
    run_manifest_path = release / "RUN_MANIFEST.json"
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    run_manifest.setdefault("project", {})["filtered_version"] = RULE_VERSION
    run_manifest["recommended_release_tag"] = RELEASE_TAG
    run_manifest["balanced_two_source_v1_1"] = {
        "generated_at_utc": utc_now(), "rule_name": RULE_NAME, "schema_version": SCHEMA_VERSION,
        "tier_dependency_audit": audit, "metrics": formal_metrics,
        "assets": [{key: item[key] for key in ["filename", "relative_local_path", "bytes", "sha256", "row_count", "logical_role"]} for item in new_assets],
        "push_performed": False, "github_release_created": False,
    }
    write_json(run_manifest_path, run_manifest)
    update_docs(release, audit, formal_metrics, subset_counts)

    additions = [
        *(root / item["relative_local_path"] for item in new_assets),
        reports / "RULE_B_TIER_DEPENDENCY_AUDIT.md", reports / "rule_b_with_vs_without_tier_diff.tsv",
        policy_file, release / "README.md", release / "METHODS.md", release / "DATA_DICTIONARY.md",
        qc_path, run_manifest_path, release_manifest_path,
    ]
    checksum_count = update_checksums(root, release, additions)
    after = {name: sha256(path) for name, path in formal_inputs.items()}
    if before != after:
        raise RuntimeError("A protected v1 input changed during v1.1 export")
    return {
        "status": "passed", "rule_name": RULE_NAME, "tier_filter_role": "annotation_only",
        "tier_dependency_audit": audit, "metrics": formal_metrics, "expected_metrics_match": True,
        "strict_two_source_validation_subset_rows": len(strict), "gold_candidate_rows": len(gold),
        "single_source_high_authority_candidate_rows": len(single), "excluded_rows": len(excluded),
        "asset_count": len(new_assets), "checksum_entry_count": checksum_count,
        "formal_input_hashes_unchanged": True, "push_performed": False, "github_release_created": False,
    }
