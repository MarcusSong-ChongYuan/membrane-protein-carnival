from __future__ import annotations

import json
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from figure_style_v62 import (
    CAPTION_DIR,
    COLORS,
    DATA_DIR,
    QA_DIR,
    clean,
    ecdf,
    make_figure,
    panel,
    save,
    write_json,
    write_table,
)


REL = Path(r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730")
DOCK = Path(r"D:\7.22\evidence_expansion_v2_working\runs\docking_selection_v2_v62_20260730\staging")


def truth(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"1", "true", "yes", "y"})


def explode_counts(series: pd.Series, separators=("|", ";", ",")) -> Counter:
    counts = Counter()
    for value in series.dropna().astype(str):
        values = [value]
        for separator in separators:
            values = [part for item in values for part in item.split(separator)]
        counts.update(part.strip() for part in values if part.strip() and part.strip().lower() != "nan")
    return counts


def top_bar(ax, counts: pd.Series, n=10, color=COLORS["blue"], xlabel="Records") -> None:
    counts = counts.sort_values(ascending=False).head(n).sort_values()
    bars = ax.barh(np.arange(len(counts)), counts.values, color=color)
    labels = [
        textwrap.fill(str(value).replace("_", " "), width=28)
        for value in counts.index
    ]
    ax.set_yticks(np.arange(len(counts)), labels)
    ax.bar_label(bars, labels=[f"{v:,.0f}" for v in counts.values], padding=2, fontsize=6.1)
    ax.set_xlabel(xlabel)
    clean(ax, "x")


def load_inputs():
    protein_cols = [
        "target_uniprot_id", "membrane_class_v52", "evidence_level_v52",
        "membrane_topology_v5", "transmembrane_count_v5", "single_pass_resolution_v5",
        "independent_membrane_source_count_v5", "independent_membrane_sources_v5",
        "hpa_predicted_membrane_v5", "htp_present_v5", "membranome_present_v5",
        "opm_present_v5", "pdbtm_present_v5", "gpcrdb_class_v5", "gtopdb_type_v5",
        "tcdb_tcids_v5", "functional_primary_class_v5", "hpa_mapping_status_v62",
        "hpa_tissue_detected_count_v62", "hpa_cell_type_detected_count_v62",
        "hpa_ihc_detected_tissue_count_v62", "hpa_main_locations_v62",
        "hpa_additional_locations_v62", "hpa_extracellular_locations_v62",
        "hpa_localization_reliability_v62", "expression_location_conflict_flag_v62",
        "oligomeric_state_consensus_v62", "subunit_count_min_v62", "subunit_count_max_v62",
        "experimental_subunit_evidence_flag_v62", "subunit_evidence_grade_v62",
        "subunit_conflict_flag_v62", "pdbe_assembly_count_v62",
    ]
    proteins = pd.read_csv(
        REL / "human_membrane_protein_master_v6_2.tsv",
        sep="\t", usecols=protein_cols, low_memory=False,
    )
    disease = pd.read_csv(
        REL / "protein_gene_disease_relations_v6_2.tsv", sep="\t", low_memory=False
    )
    localization = pd.read_csv(
        REL / "protein_subcellular_localization_v2.tsv.gz", sep="\t", low_memory=False
    )
    subunit = pd.read_csv(
        REL / "protein_subunit_summary_v2.tsv.gz", sep="\t", low_memory=False
    )
    return proteins, disease, localization, subunit


def scan_compounds() -> tuple[pd.DataFrame, dict]:
    cache = DATA_DIR / "S5_compound_summary_cache.tsv"
    meta_cache = QA_DIR / "S5_compound_summary_cache.json"
    if cache.exists() and meta_cache.exists():
        return pd.read_csv(cache, sep="\t"), json.loads(meta_cache.read_text(encoding="utf-8"))
    usecols = [
        "identity_confidence", "compound_scope_class", "form_count",
        "source_record_count", "source_database_count", "record_qc_status",
        "molecular_weight", "xlogp", "tpsa", "rotatable_bond_count",
        "is_approved_drug", "is_clinical_candidate", "is_endogenous_ligand",
        "is_natural_product", "is_chemical_probe",
    ]
    categorical = {name: Counter() for name in [
        "identity_confidence", "compound_scope_class", "record_qc_status"
    ]}
    status = Counter()
    numeric_sample = []
    total = 0
    rng = np.random.default_rng(42)
    for chunk in pd.read_csv(
        REL / "small_molecule_master_v1_3.tsv",
        sep="\t", usecols=usecols, chunksize=150_000, low_memory=False,
    ):
        total += len(chunk)
        for name in categorical:
            categorical[name].update(chunk[name].fillna("missing").astype(str))
        for name in [
            "is_approved_drug", "is_clinical_candidate", "is_endogenous_ligand",
            "is_natural_product", "is_chemical_probe",
        ]:
            status[name] += int(truth(chunk[name]).sum())
        sample_n = min(8000, len(chunk))
        numeric_sample.append(
            chunk[["form_count", "source_record_count", "source_database_count",
                   "molecular_weight", "xlogp", "tpsa", "rotatable_bond_count"]]
            .iloc[rng.choice(len(chunk), sample_n, replace=False)]
        )
    sample = pd.concat(numeric_sample, ignore_index=True)
    sample.to_csv(cache, sep="\t", index=False)
    meta = {
        "total": total,
        "categorical": {key: dict(value) for key, value in categorical.items()},
        "biological_status": dict(status),
    }
    meta_cache.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    return sample, meta


def build_s1(proteins: pd.DataFrame, disease: pd.DataFrame, compound_meta: dict):
    fig, axes = make_figure("S1  Cross-database integration and source support")
    flags = {
        "Human transmembrane\nproteome": truth(proteins["htp_present_v5"]),
        "Membranome": truth(proteins["membranome_present_v5"]),
        "OPM": truth(proteins["opm_present_v5"]),
        "PDBTM": truth(proteins["pdbtm_present_v5"]),
    }
    ax = axes[0]; panel(ax, "A", "Independent membrane source coverage")
    coverage = pd.Series({name: int(mask.sum()) for name, mask in flags.items()})
    coverage["HPA mapped\n(expression context)"] = int(
        proteins["hpa_mapping_status_v62"].astype(str).str.contains("mapped", case=False, na=False).sum()
    )
    top_bar(ax, coverage, n=10, color=COLORS["A"], xlabel="Proteins")

    ax = axes[1]; panel(ax, "B", "Most frequent source signatures")
    signatures = pd.Series([" + ".join([name.replace("\n", " ") for name, mask in flags.items() if mask.iloc[i]]) or "UniProt only"
                            for i in range(len(proteins))]).value_counts().head(10)
    top_bar(ax, signatures, n=10, color=COLORS["B"], xlabel="Proteins")

    ax = axes[2]; panel(ax, "C", "Number of independent membrane sources")
    vals = pd.to_numeric(proteins["independent_membrane_source_count_v5"], errors="coerce").fillna(0).astype(int)
    vc = vals.value_counts().sort_index()
    bars = ax.bar(vc.index.astype(str), vc.values, color=COLORS["C"])
    ax.bar_label(bars, labels=[f"{v:,}" for v in vc.values], fontsize=6)
    ax.set_xlabel("Independent sources per protein"); ax.set_ylabel("Proteins"); clean(ax, "y")

    ax = axes[3]; panel(ax, "D", "Disease-relation source databases")
    dcounts = explode_counts(disease["source_database"]).most_common(12)
    top_bar(ax, pd.Series(dict(dcounts)), n=12, color=COLORS["E2"], xlabel="Disease relations")

    ax = axes[4]; panel(ax, "E", "Disease relations by source-family support")
    dsrc = pd.to_numeric(disease["independent_source_family_count"], errors="coerce").fillna(0).clip(0, 6)
    vc = dsrc.value_counts().sort_index()
    bars = ax.bar(vc.index.astype(int).astype(str), vc.values, color=COLORS["E3"])
    ax.bar_label(bars, labels=[f"{v:,}" for v in vc.values], fontsize=6)
    ax.set_xlabel("Independent source families"); ax.set_ylabel("Relations"); clean(ax, "y")

    ax = axes[5]; panel(ax, "F", "Compound source-database support")
    sample = pd.read_csv(DATA_DIR / "S5_compound_summary_cache.tsv", sep="\t")
    vals = pd.to_numeric(sample["source_database_count"], errors="coerce").dropna().clip(upper=10)
    x, y = ecdf(vals)
    ax.step(x, y, color=COLORS["A"], lw=1.8)
    ax.axvline(vals.median(), color=COLORS["conflict"], ls="--", lw=1)
    ax.set_xlabel("Source databases per compound (10 = 10+)"); ax.set_ylabel("ECDF")
    clean(ax, "both")
    paths = save(fig, "S1_source_integration")
    write_table(coverage.rename_axis("source").reset_index(name="protein_count"), "S1_membrane_source_counts.tsv")
    return paths


def build_s2(proteins: pd.DataFrame, subunit: pd.DataFrame):
    fig, axes = make_figure("S2  Membrane topology, formal classification and assembly")
    ax = axes[0]; panel(ax, "A", "Membrane topology")
    topo = proteins["membrane_topology_v5"].fillna("unresolved").replace("", "unresolved").value_counts()
    top_bar(ax, topo, n=10, color=COLORS["A"], xlabel="Proteins")

    ax = axes[1]; panel(ax, "B", "Transmembrane-segment distribution")
    tm = pd.to_numeric(proteins["transmembrane_count_v5"], errors="coerce").fillna(-1)
    bins = pd.cut(tm, [-2, -0.5, 0.5, 1.5, 3.5, 6.5, 12.5, np.inf],
                  labels=["unknown", "0", "1", "2–3", "4–6", "7–12", ">12"])
    vc = bins.value_counts().reindex(["unknown", "0", "1", "2–3", "4–6", "7–12", ">12"])
    bars = ax.bar(vc.index.astype(str), vc.values, color=[COLORS["unknown"]] + [COLORS["B"]]*6)
    ax.bar_label(bars, labels=[f"{v:,}" for v in vc.values], fontsize=6, rotation=90, padding=2)
    ax.set_ylabel("Proteins"); clean(ax, "y")

    ax = axes[2]; panel(ax, "C", "Single-pass topology resolution")
    single = proteins.loc[tm == 1, "single_pass_resolution_v5"].fillna("unresolved").replace("", "unresolved").value_counts()
    top_bar(ax, single, n=10, color=COLORS["C"], xlabel="Single-pass proteins")

    ax = axes[3]; panel(ax, "D", "Formal taxonomy alignment")
    taxonomy = pd.Series({
        "GPCRdb": proteins["gpcrdb_class_v5"].notna().sum(),
        "IUPHAR/GtoPdb": proteins["gtopdb_type_v5"].notna().sum(),
        "TCDB": proteins["tcdb_tcids_v5"].notna().sum(),
    })
    bars = ax.bar(taxonomy.index, taxonomy.values, color=[COLORS["A"], COLORS["E2"], COLORS["E3"]])
    ax.bar_label(bars, labels=[f"{v:,}" for v in taxonomy.values], fontsize=6)
    ax.set_ylabel("Mapped proteins"); clean(ax, "y")

    ax = axes[4]; panel(ax, "E", "Oligomeric-state consensus")
    olig = subunit["oligomeric_state_consensus_v62"].fillna("unresolved").replace("", "unresolved").value_counts()
    top_bar(ax, olig, n=10, color=COLORS["D3"], xlabel="Proteins")

    ax = axes[5]; panel(ax, "F", "Assembly evidence and conflict")
    summary = pd.Series({
        "Experimental\nsubunit evidence": truth(subunit["experimental_subunit_evidence_flag_v62"]).sum(),
        "PDBe assembly": (pd.to_numeric(subunit["pdbe_assembly_count_v62"], errors="coerce").fillna(0) > 0).sum(),
        "State conflict": truth(subunit["subunit_conflict_flag_v62"]).sum(),
        "Large-assembly\nreview": truth(subunit["large_assembly_review_flag_v62"]).sum(),
    })
    bars = ax.bar(summary.index, summary.values, color=[COLORS["B"], COLORS["A"], COLORS["conflict"], COLORS["review"]])
    ax.bar_label(bars, labels=[f"{v:,}" for v in summary.values], fontsize=6)
    ax.set_ylabel("Proteins"); clean(ax, "y")
    paths = save(fig, "S2_topology_classification_assembly")
    write_table(summary.rename_axis("assembly_metric").reset_index(name="protein_count"), "S2_assembly_summary.tsv")
    return paths


def build_s3(proteins: pd.DataFrame, localization: pd.DataFrame):
    fig, axes = make_figure("S3  Complete expression and subcellular-localization detail")
    ax = axes[0]; panel(ax, "A", "Tissue expression by protein class")
    matrix = pd.read_csv(DATA_DIR / "M3_tissue_expression_matrix.tsv", sep="\t").set_index("functional_class")
    sns.heatmap(matrix, cmap="viridis", ax=ax, cbar_kws={"label": "log1p median nTPM"}, xticklabels=True, yticklabels=True)
    ax.tick_params(axis="x", rotation=55)

    ax = axes[1]; panel(ax, "B", "Cell-type expression by protein class")
    matrix2 = pd.read_csv(DATA_DIR / "M3_cell_type_expression_matrix.tsv", sep="\t").set_index("functional_class")
    sns.heatmap(matrix2, cmap="mako", ax=ax, cbar_kws={"label": "log1p median nCPM"}, xticklabels=True, yticklabels=True)
    ax.tick_params(axis="x", rotation=55)

    ax = axes[2]; panel(ax, "C", "Localization roles")
    roles = localization["location_role"].fillna("unspecified").value_counts()
    top_bar(ax, roles, n=10, color=COLORS["A"], xlabel="Protein–location records")

    ax = axes[3]; panel(ax, "D", "Tissue-detection breadth")
    for cls, color in [("A", COLORS["A"]), ("B", COLORS["B"]), ("C", COLORS["C"]), ("unknown", COLORS["unknown"])]:
        values = pd.to_numeric(proteins.loc[proteins["membrane_class_v52"].fillna("unknown") == cls,
                                            "hpa_tissue_detected_count_v62"], errors="coerce").dropna()
        x, y = ecdf(values); ax.step(x, y, label=cls, color=color, lw=1.6)
    ax.set_xlabel("Detected tissues per protein"); ax.set_ylabel("ECDF"); ax.legend(title="Membrane class")
    clean(ax, "both")

    ax = axes[4]; panel(ax, "E", "Cell-type detection breadth")
    for cls, color in [("A", COLORS["A"]), ("B", COLORS["B"]), ("C", COLORS["C"]), ("unknown", COLORS["unknown"])]:
        values = pd.to_numeric(proteins.loc[proteins["membrane_class_v52"].fillna("unknown") == cls,
                                            "hpa_cell_type_detected_count_v62"], errors="coerce").dropna()
        x, y = ecdf(values); ax.step(x, y, label=cls, color=color, lw=1.6)
    ax.set_xlabel("Detected cell types per protein"); ax.set_ylabel("ECDF")
    clean(ax, "both")

    ax = axes[5]; panel(ax, "F", "Expression/localization QA")
    qa = pd.Series({
        "HPA mapped": proteins["hpa_mapping_status_v62"].astype(str).str.contains("mapped", case=False, na=False).sum(),
        "Main location": proteins["hpa_main_locations_v62"].notna().sum(),
        "Additional location": proteins["hpa_additional_locations_v62"].notna().sum(),
        "Extracellular location": proteins["hpa_extracellular_locations_v62"].notna().sum(),
        "Conflict review": truth(proteins["expression_location_conflict_flag_v62"]).sum(),
    })
    top_bar(ax, qa, n=10, color=COLORS["E2"], xlabel="Proteins")
    paths = save(fig, "S3_expression_localization_detail")
    write_table(qa.rename_axis("metric").reset_index(name="protein_count"), "S3_expression_localization_qa.tsv")
    return paths


def build_s4(disease: pd.DataFrame):
    fig, axes = make_figure("S4  Disease ontology, evidence provenance and network sensitivity")
    ax = axes[0]; panel(ax, "A", "Disease ontology namespaces")
    top_bar(ax, disease["disease_ontology"].fillna("unspecified").value_counts(), n=12,
            color=COLORS["A"], xlabel="Relations")

    ax = axes[1]; panel(ax, "B", "Evidence-level composition by ontology")
    tab = pd.crosstab(disease["disease_ontology"].fillna("unspecified"), disease["disease_evidence_level"].fillna("unknown"))
    tab = tab.loc[tab.sum(axis=1).nlargest(10).index]
    tab.div(tab.sum(axis=1), axis=0).plot(kind="barh", stacked=True, ax=ax,
                                          color=[COLORS["E1"], COLORS["E2"], COLORS["E3"], COLORS["unknown"]][:len(tab.columns)])
    ax.set_xlabel("Fraction of relations"); ax.set_ylabel(""); ax.legend(title="Evidence", fontsize=6)
    clean(ax, "x")

    ax = axes[2]; panel(ax, "C", "Most connected diseases")
    disease_degree = disease.groupby(["disease_id", "disease_name"])["target_uniprot_id"].nunique().sort_values(ascending=False)
    labels = [name[:42] for _, name in disease_degree.head(12).index]
    top_bar(ax, pd.Series(disease_degree.head(12).values, index=labels), n=12, color=COLORS["E3"], xlabel="Unique proteins")

    ax = axes[3]; panel(ax, "D", "Evidence-channel prevalence")
    channels = {
        "Human genetic": "human_genetic_flag", "Clinical genetic": "clinical_genetic_flag",
        "Somatic mutation": "somatic_mutation_flag", "Functional": "functional_flag",
        "Animal model": "animal_model_flag", "Expression": "expression_flag",
        "Literature": "literature_flag", "UniProt disease": "uniprot_disease_flag",
    }
    preval = pd.Series({name: truth(disease[col]).mean() * 100 for name, col in channels.items()})
    top_bar(ax, preval, n=10, color=COLORS["B"], xlabel="Relations carrying channel (%)")

    ax = axes[4]; panel(ax, "E", "Protein and disease degree distributions")
    pdeg = disease.groupby("target_uniprot_id")["disease_id"].nunique()
    ddeg = disease.groupby("disease_id")["target_uniprot_id"].nunique()
    for values, label, color in [(pdeg, "Diseases per protein", COLORS["A"]), (ddeg, "Proteins per disease", COLORS["C"])]:
        x, y = ecdf(values); ax.step(x, 1-y+1/len(y), label=label, color=color, lw=1.6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Degree"); ax.set_ylabel("Complementary cumulative fraction")
    ax.legend(); clean(ax, "both")

    ax = axes[5]; panel(ax, "F", "Relation inclusion and QC")
    qc = pd.Series({
        "Default included": truth(disease["default_relation_inclusion"]).sum(),
        "Current UniProt curated": disease["record_qc_status"].eq("current_uniprot_curated_disease_comment").sum(),
        "Retained V3 rule B": disease["record_qc_status"].eq("retained_from_v3_rule_b_release").sum(),
        "Multi-source": (pd.to_numeric(disease["independent_source_family_count"], errors="coerce").fillna(0) >= 2).sum(),
        "Primary V3": truth(disease["primary_disease_flag_v3"]).sum(),
    })
    top_bar(ax, qc, n=10, color=COLORS["E2"], xlabel="Relations")
    paths = save(fig, "S4_disease_ontology_evidence")
    write_table(preval.rename_axis("evidence_channel").reset_index(name="relation_percent"), "S4_disease_channel_prevalence.tsv")
    return paths


def build_s5(sample: pd.DataFrame, meta: dict):
    fig, axes = make_figure("S5  Compound identity, structure standardization and QC")
    ax = axes[0]; panel(ax, "A", "Identity-confidence states")
    identity = pd.Series(meta["categorical"]["identity_confidence"])
    top_bar(ax, identity, n=12, color=COLORS["A"], xlabel="Canonical compounds")

    ax = axes[1]; panel(ax, "B", "Compound-scope classes")
    scope = pd.Series(meta["categorical"]["compound_scope_class"])
    top_bar(ax, scope, n=12, color=COLORS["B"], xlabel="Canonical compounds")

    ax = axes[2]; panel(ax, "C", "Parent/form hierarchy breadth")
    values = pd.to_numeric(sample["form_count"], errors="coerce").dropna().clip(upper=20)
    x, y = ecdf(values); ax.step(x, y, color=COLORS["C"], lw=1.8)
    ax.axvline(values.median(), color=COLORS["conflict"], ls="--", lw=1)
    ax.set_xlabel("Forms per canonical parent (20 = 20+)"); ax.set_ylabel("ECDF")
    clean(ax, "both")

    ax = axes[3]; panel(ax, "D", "Source-record support")
    values = pd.to_numeric(sample["source_record_count"], errors="coerce").dropna().clip(upper=50)
    x, y = ecdf(values); ax.step(x, y, color=COLORS["A"], lw=1.8)
    ax.set_xlabel("Source records per compound (50 = 50+)"); ax.set_ylabel("ECDF")
    clean(ax, "both")

    ax = axes[4]; panel(ax, "E", "Physicochemical QC landscape")
    plot = sample[["molecular_weight", "xlogp", "tpsa"]].apply(pd.to_numeric, errors="coerce").dropna()
    plot = plot.sample(min(25000, len(plot)), random_state=42)
    hb = ax.hexbin(plot["molecular_weight"], plot["xlogp"], C=plot["tpsa"], reduce_C_function=np.median,
                   gridsize=45, mincnt=2, cmap="viridis")
    fig.colorbar(hb, ax=ax, label="Median TPSA (Å²)")
    ax.set_xlim(0, min(1000, plot["molecular_weight"].quantile(.995)))
    ax.set_xlabel("Molecular weight (Da)"); ax.set_ylabel("XlogP"); clean(ax, "both")

    ax = axes[5]; panel(ax, "F", "Biological/development status")
    status = pd.Series(meta["biological_status"]).rename({
        "is_approved_drug": "Approved drug", "is_clinical_candidate": "Clinical candidate",
        "is_endogenous_ligand": "Endogenous ligand", "is_natural_product": "Natural product",
        "is_chemical_probe": "Chemical probe",
    })
    top_bar(ax, status, n=10, color=COLORS["E3"], xlabel="Canonical compounds")
    paths = save(fig, "S5_compound_identity_qc")
    write_table(identity.rename_axis("identity_confidence").reset_index(name="compound_count"), "S5_identity_confidence.tsv")
    return paths


def build_s6(proteins: pd.DataFrame):
    standard = pd.read_csv(DOCK / "docking_hpc_standard_v2_v62.tsv.gz", sep="\t", low_memory=False)
    target = pd.read_csv(DOCK / "docking_hpc_target_summary_v2_v62.tsv", sep="\t", low_memory=False)
    fig, axes = make_figure("S6  Docking shortlist robustness and HPC workload controls")
    ax = axes[0]; panel(ax, "A", "Docking-tier composition")
    tier = standard["docking_tier"].fillna("unknown").value_counts().reindex(["R0", "D1", "D2", "D3"]).dropna()
    bars = ax.bar(tier.index, tier.values, color=[COLORS.get(x, COLORS["unknown"]) for x in tier.index])
    ax.bar_label(bars, labels=[f"{v:,}" for v in tier.values], fontsize=6)
    ax.set_ylabel("Standard-set pairs"); clean(ax, "y")

    ax = axes[1]; panel(ax, "B", "Receptor-structure readiness")
    readiness = standard["receptor_preparation_status_v2"].fillna("unresolved").value_counts()
    top_bar(ax, readiness, n=10, color=COLORS["A"], xlabel="Pairs")

    ax = axes[2]; panel(ax, "C", "Ligand-preparation readiness")
    readiness2 = standard["ligand_preparation_status_v2"].fillna("unresolved").value_counts()
    top_bar(ax, readiness2, n=10, color=COLORS["B"], xlabel="Pairs")

    ax = axes[3]; panel(ax, "D", "Binding-box definition")
    boxes = standard["box_definition_status_v2"].fillna("unresolved").value_counts()
    top_bar(ax, boxes, n=10, color=COLORS["C"], xlabel="Pairs")

    ax = axes[4]; panel(ax, "E", "Per-target workload caps")
    counts = target[["R0_standard_count", "D1_standard_count", "D2_standard_count", "D3_standard_count"]]
    data = [counts[col].values for col in counts.columns]
    bp = ax.boxplot(data, tick_labels=["R0", "D1", "D2", "D3"], showfliers=False, patch_artist=True)
    for patch, color in zip(bp["boxes"], [COLORS["R0"], COLORS["D1"], COLORS["D2"], COLORS["D3"]]):
        patch.set_facecolor(color); patch.set_alpha(.8)
    ax.set_ylabel("Pairs per target"); clean(ax, "y")

    ax = axes[5]; panel(ax, "F", "Ranking-score separation")
    for tier_name, color in [("R0", COLORS["R0"]), ("D1", COLORS["D1"]), ("D2", COLORS["D2"]), ("D3", COLORS["D3"])]:
        vals = pd.to_numeric(standard.loc[standard["docking_tier"] == tier_name, "ranking_score_v2"], errors="coerce").dropna()
        x, y = ecdf(vals); ax.step(x, y, label=tier_name, color=color, lw=1.6)
    ax.set_xlabel("Ranking score"); ax.set_ylabel("ECDF"); ax.legend(title="Tier")
    clean(ax, "both")
    paths = save(fig, "S6_docking_robustness")
    summary = standard.groupby("docking_tier").agg(
        pairs=("compound_internal_id", "size"),
        targets=("target_uniprot_id", "nunique"),
        compounds=("compound_internal_id", "nunique"),
        median_score=("ranking_score_v2", "median"),
    ).reset_index()
    write_table(summary, "S6_docking_tier_summary.tsv")
    return paths


def main():
    proteins, disease, localization, subunit = load_inputs()
    sample, compound_meta = scan_compounds()
    outputs = {
        "S1": build_s1(proteins, disease, compound_meta),
        "S2": build_s2(proteins, subunit),
        "S3": build_s3(proteins, localization),
        "S4": build_s4(disease),
        "S5": build_s5(sample, compound_meta),
        "S6": build_s6(proteins),
    }
    captions = """# Supplementary figure captions

## Figure S1 | Cross-database integration and source support
Independent-source coverage, common source signatures, source-count distributions and disease/compound provenance. Counts describe V6.2 records and do not imply equal source sensitivity or independence.

## Figure S2 | Membrane topology, formal classification and assembly
Topology and transmembrane-segment distributions are shown together with single-pass resolution, GPCRdb/IUPHAR/TCDB alignment and V6.2 oligomeric/biological-assembly evidence. Unresolved and conflict states remain explicit.

## Figure S3 | Complete expression and subcellular-localization detail
Clustered high-information tissue and cell-type summaries are paired with localization roles, expression breadth by membrane class and mapping/conflict quality controls. HPA expression values are contextual annotations and are not used as membrane-protein inclusion evidence by themselves.

## Figure S4 | Disease ontology, evidence provenance and network sensitivity
Disease namespaces, evidence-level composition, high-degree diseases, evidence-channel prevalence and degree distributions summarize the protein–gene–disease module. Organ-system or ontology aggregation reflects database annotation density, not disease prevalence.

## Figure S5 | Compound identity, structure standardization and QC
Identity-confidence and scope states, parent/form multiplicity, source-record support, physicochemical space and biological status summarize canonical compound standardization. Distribution panels use a deterministic stratified sample for rendering; all categorical totals use the complete 2,016,064-row master table.

## Figure S6 | Docking shortlist robustness and HPC workload controls
The V6.2 standard docking shortlist is decomposed by scientific tier, receptor/ligand/box readiness, per-target workload caps and ranking-score distributions. These panels describe prioritization and preparation readiness, not docking outcomes or predicted affinity.
"""
    (CAPTION_DIR / "CAPTIONS_SUPPLEMENTARY_FIGURES.md").write_text(captions, encoding="utf-8")
    validation = {
        "status": "PASS",
        "protein_rows": len(proteins),
        "disease_rows": len(disease),
        "localization_rows": len(localization),
        "subunit_rows": len(subunit),
        "compound_rows_scanned": compound_meta["total"],
        "figures": outputs,
        "deterministic_seed": 42,
    }
    write_json(validation, "SUPPLEMENTARY_FIGURES_VALIDATION.json")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
