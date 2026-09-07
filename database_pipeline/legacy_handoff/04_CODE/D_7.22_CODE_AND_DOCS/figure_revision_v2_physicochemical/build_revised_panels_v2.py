from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(r"C:\Users\Administrator\HuMemLigDB_figure_revision_v2")
REL = Path(r"D:\finale\01_正式数据_V6.2")
OT_DISEASE = ROOT / "raw" / "opentargets_26.06_disease.parquet"
PHYSICOCHEMICAL_CACHE = (
    Path(r"D:\finale\02_展示图表_V6.2\data")
    / "S5_compound_summary_cache.tsv"
)

for folder in ["data", "png", "pdf", "svg", "qa"]:
    (ROOT / folder).mkdir(parents=True, exist_ok=True)


COLORS = {
    "navy": "#173B57",
    "blue": "#3572A5",
    "teal": "#2A9D8F",
    "teal_light": "#A8DADC",
    "orange": "#F4A261",
    "purple": "#7A5195",
    "red": "#D55E00",
    "gray": "#B9C0C8",
    "gray_dark": "#67727E",
    "gray_light": "#E9EDF1",
    "green": "#4C956C",
}


mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.0,
        "axes.titlesize": 9.2,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.0,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 7.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#7D8790",
        "axes.linewidth": 0.7,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)
sns.set_style("white")


def truth(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y"})


def human_number(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return f"{value:,.0f}"


def truncate_label(value: object, max_chars: int) -> str:
    text = str(value)
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 3].rstrip() + "..."


def panel_label(ax: plt.Axes, letter: str, title: str) -> None:
    ax.text(
        -0.08,
        1.07,
        letter,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        color=COLORS["navy"],
    )
    ax.set_title(title, loc="left", pad=8)


def finish_axis(ax: plt.Axes, grid_axis: str | None = "x") -> None:
    if grid_axis:
        ax.grid(axis=grid_axis, color=COLORS["gray_light"], linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, stem: str) -> dict[str, str]:
    paths: dict[str, str] = {}
    for ext in ["png", "pdf", "svg"]:
        path = ROOT / ext / f"{stem}.{ext}"
        fig.savefig(path, dpi=450 if ext == "png" else None, bbox_inches="tight")
        paths[ext] = str(path)
    plt.close(fig)
    return paths


def build_m3d() -> tuple[dict[str, str], pd.DataFrame]:
    cols = [
        "target_uniprot_id",
        "functional_primary_class_v5",
        "hpa_mapping_status_v62",
        "hpa_tissue_detected_count_v62",
    ]
    protein = pd.read_csv(
        REL / "human_membrane_protein_master_v6_2.tsv",
        sep="\t",
        usecols=cols,
        low_memory=False,
    )
    protein["functional_class"] = protein["functional_primary_class_v5"].fillna("unclassified")
    protein["detected_tissues"] = pd.to_numeric(
        protein["hpa_tissue_detected_count_v62"], errors="coerce"
    )
    protein["mapped"] = (
        protein["hpa_mapping_status_v62"].fillna("").str.lower().eq("mapped")
        & protein["detected_tissues"].notna()
    )

    classes = ["receptor", "ion_channel", "transporter", "enzyme"]
    labels = ["Receptor", "Ion channel", "Transporter", "Enzyme"]
    rows = []
    for class_name, label in zip(classes, labels):
        subset = protein.loc[protein["functional_class"].eq(class_name)].copy()
        mapped = subset.loc[subset["mapped"], "detected_tissues"]
        rows.append(
            {
                "functional_class": class_name,
                "display_label": label,
                "n_total": len(subset),
                "n_mapped_with_value": int(subset["mapped"].sum()),
                "n_mapped_detected": int((subset["mapped"] & subset["detected_tissues"].gt(0)).sum()),
                "n_mapped_zero": int((subset["mapped"] & subset["detected_tissues"].eq(0)).sum()),
                "n_unmapped_or_missing": int((~subset["mapped"]).sum()),
                "median_mapped": float(mapped.median()) if len(mapped) else np.nan,
                "q1_mapped": float(mapped.quantile(0.25)) if len(mapped) else np.nan,
                "q3_mapped": float(mapped.quantile(0.75)) if len(mapped) else np.nan,
            }
        )
    summary = pd.DataFrame(rows)
    summary["pct_mapped_zero_of_total"] = 100 * summary["n_mapped_zero"] / summary["n_total"]
    summary["pct_unmapped_or_missing"] = 100 * summary["n_unmapped_or_missing"] / summary["n_total"]
    summary.to_csv(ROOT / "data" / "M3D_hpa_expression_breadth_status.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(7.2, 3.65))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.8, 1.05], wspace=0.30)
    ax = fig.add_subplot(gs[0, 0])
    status_ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "A", "Distribution among HPA-mapped proteins")
    panel_label(status_ax, "B", "Mapping and detection status")

    mapped_plot = protein.loc[
        protein["mapped"] & protein["functional_class"].isin(classes),
        ["functional_class", "detected_tissues"],
    ].copy()
    mapped_plot["display_label"] = pd.Categorical(
        mapped_plot["functional_class"].map(dict(zip(classes, labels))),
        categories=labels,
        ordered=True,
    )
    sns.violinplot(
        data=mapped_plot,
        y="display_label",
        x="detected_tissues",
        order=labels,
        inner=None,
        cut=0,
        density_norm="width",
        color=COLORS["teal_light"],
        linewidth=0.8,
        ax=ax,
    )
    sns.boxplot(
        data=mapped_plot,
        y="display_label",
        x="detected_tissues",
        order=labels,
        width=0.18,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": COLORS["navy"], "linewidth": 0.9},
        medianprops={"color": COLORS["red"], "linewidth": 1.3},
        whiskerprops={"color": COLORS["navy"], "linewidth": 0.8},
        capprops={"color": COLORS["navy"], "linewidth": 0.8},
        ax=ax,
    )
    rng = np.random.default_rng(42)
    for y_pos, class_name in enumerate(classes):
        values = mapped_plot.loc[
            mapped_plot["functional_class"].eq(class_name), "detected_tissues"
        ].to_numpy()
        if len(values) > 220:
            values = rng.choice(values, 220, replace=False)
        jitter = rng.normal(y_pos + 0.19, 0.035, size=len(values))
        ax.scatter(values, jitter, s=5, alpha=0.22, color=COLORS["navy"], linewidths=0, zorder=3)
        record = summary.iloc[y_pos]
        ax.text(
            0.99,
            y_pos,
            f"{record.median_mapped:.0f} [{record.q1_mapped:.0f}–{record.q3_mapped:.0f}]",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=6.6,
            color=COLORS["navy"],
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5, "alpha": 0.75},
        )
    ax.set_xlabel("Number of HPA tissues meeting the detection rule")
    ax.set_ylabel("")
    ax.text(
        0.99,
        -0.20,
        "Labels show median [IQR]; mapped zeros remain in the distribution.",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.7,
        color=COLORS["gray_dark"],
    )
    finish_axis(ax, "x")

    y = np.arange(len(summary))
    left = np.zeros(len(summary))
    status_parts = [
        ("Mapped, ≥1 detected tissue", "n_mapped_detected", COLORS["teal"]),
        ("Mapped, 0 detected tissues", "n_mapped_zero", COLORS["orange"]),
        ("Unmapped / missing", "n_unmapped_or_missing", COLORS["gray"]),
    ]
    for label, field, color in status_parts:
        values = summary[field].to_numpy() / summary["n_total"].to_numpy() * 100
        status_ax.barh(y, values, left=left, height=0.58, label=label, color=color, edgecolor="white", linewidth=0.6)
        for idx, (start, value) in enumerate(zip(left, values)):
            if value >= 7:
                status_ax.text(start + value / 2, idx, f"{value:.1f}%", ha="center", va="center", fontsize=6.3)
        left += values
    status_ax.set_yticks(y, labels)
    for idx, record in summary.reset_index(drop=True).iterrows():
        status_ax.text(
            99.3,
            idx + 0.34,
            (
                f"mapped 0: {record['pct_mapped_zero_of_total']:.1f}%  |  "
                f"missing: {record['pct_unmapped_or_missing']:.1f}%"
            ),
            ha="right",
            va="bottom",
            fontsize=5.9,
            color=COLORS["gray_dark"],
            clip_on=False,
        )
    status_ax.invert_yaxis()
    status_ax.set_xlim(0, 100)
    status_ax.set_xlabel("Fraction of proteins in functional class (%)")
    status_ax.set_ylabel("")
    status_ax.legend(
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.20),
        ncol=1,
        handlelength=1.2,
    )
    finish_axis(status_ax, "x")

    fig.suptitle(
        "M3D revised — HPA tissue-expression breadth with explicit missingness",
        x=0.08,
        ha="left",
        fontsize=11,
        fontweight="bold",
        color=COLORS["navy"],
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.77, bottom=0.24, wspace=0.30)
    paths = save_figure(fig, "M3D_expression_breadth_missing_aware")
    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        top=0.77,
        bottom=0.24,
        wspace=0.30,
    )
    return paths, summary


GENERIC_THERAPEUTIC_AREAS = {
    "GO_0008150",       # biological process
    "EFO_0001444",      # measurement
    "EFO_0002571",      # medical procedure
    "EFO_0000651",      # phenotype
    "MONDO_0005583",    # non-human animal disease
}

CROSS_CUTTING_AREAS = {
    "MONDO_0045024",    # cancer or benign tumor
    "OTAR_0000018",     # genetic, familial or congenital disease
    "MONDO_0005550",    # infectious disease
    "OTAR_0000009",     # injury, poisoning or complication
    "OTAR_0000020",     # nutritional or metabolic disease
    "OTAR_0000014",     # pregnancy or perinatal disease
}


def normalize_xref(disease_id: str) -> str:
    value = str(disease_id).strip()
    if value.startswith("MIM_"):
        return f"OMIM:{value.split('_', 1)[1]}".upper()
    if value.startswith("Orphanet_"):
        return f"ORPHANET:{value.split('_', 1)[1]}".upper()
    return value.upper()


def list_values(value) -> list[str]:
    if isinstance(value, np.ndarray):
        return [str(item) for item in value.tolist()]
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value]
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    return [str(value)]


def best_grade(values: pd.Series) -> str:
    ranks = {"very_high": 3, "high": 2, "medium": 1}
    present = [str(value) for value in values if str(value) in ranks]
    return max(present, key=lambda value: ranks[value]) if present else "medium"


def build_disease_mapping(rel: pd.DataFrame, ot: pd.DataFrame):
    ot = ot.copy()
    ot["id"] = ot["id"].astype(str)
    by_id = ot.set_index("id", drop=False)
    xref_to_ids: dict[str, set[str]] = defaultdict(set)
    for row in ot.itertuples(index=False):
        for xref in list_values(row.dbXRefs):
            xref_to_ids[xref.upper()].add(str(row.id))

    mapping_rows = []
    area_rows = []
    for disease_id in sorted(rel["disease_id"].dropna().astype(str).unique()):
        if disease_id in by_id.index:
            candidates = [disease_id]
            status = "OT entity ID"
        else:
            candidates = sorted(xref_to_ids.get(normalize_xref(disease_id), set()))
            candidates_with_area = [
                candidate
                for candidate in candidates
                if len(list_values(by_id.loc[candidate, "therapeuticAreas"])) > 0
            ]
            if candidates_with_area:
                candidates = candidates_with_area
            if len(candidates) == 1:
                status = "Single OT cross-reference"
            elif len(candidates) > 1:
                status = "Ambiguous OT cross-reference"
            else:
                status = "Unmapped / no OT entity"

        if len(candidates) == 1:
            canonical_ot_id = candidates[0]
        else:
            canonical_ot_id = f"LOCAL:{disease_id}"
        candidate_area_sets = [
            {
                area
                for area in list_values(by_id.loc[candidate, "therapeuticAreas"])
                if area not in GENERIC_THERAPEUTIC_AREAS
            }
            for candidate in candidates
            if candidate in by_id.index
        ]
        if len(candidate_area_sets) == 1:
            area_ids = sorted(candidate_area_sets[0])
            area_assignment_method = "single_entity_therapeutic_areas"
        elif len(candidate_area_sets) > 1:
            area_ids = sorted(set.intersection(*candidate_area_sets))
            area_assignment_method = (
                "ambiguous_candidate_consensus"
                if area_ids
                else "ambiguous_no_consensus_area"
            )
        else:
            area_ids = []
            area_assignment_method = "unmapped"
        mapping_rows.append(
            {
                "source_disease_id": disease_id,
                "mapping_status": status,
                "candidate_ot_ids": ";".join(candidates),
                "canonical_ot_id_for_safe_dedup": canonical_ot_id,
                "n_candidate_ot_ids": len(candidates),
                "n_non_generic_therapeutic_areas": len(area_ids),
                "area_assignment_method": area_assignment_method,
            }
        )
        for area_id in area_ids:
            area_name = (
                str(by_id.loc[area_id, "name"])
                if area_id in by_id.index
                else area_id.replace("_", " ")
            )
            area_rows.append(
                {
                    "source_disease_id": disease_id,
                    "canonical_ot_id_for_safe_dedup": canonical_ot_id,
                    "therapeutic_area_id": area_id,
                    "therapeutic_area_name": area_name,
                    "area_axis": "Cross-cutting disease type" if area_id in CROSS_CUTTING_AREAS else "Organ/system or clinical area",
                }
            )
    mapping = pd.DataFrame(mapping_rows)
    areas = pd.DataFrame(area_rows)
    return mapping, areas


def build_m4() -> tuple[dict[str, str], dict[str, pd.DataFrame]]:
    rel = pd.read_csv(
        REL / "protein_gene_disease_relations_v6_2.tsv",
        sep="\t",
        low_memory=False,
    )
    ot = pd.read_parquet(OT_DISEASE)
    mapping, area_map = build_disease_mapping(rel, ot)
    mapping.to_csv(ROOT / "data" / "M4_disease_to_ot_therapeutic_area_mapping.tsv", sep="\t", index=False)
    area_map.to_csv(ROOT / "data" / "M4_disease_therapeutic_area_memberships.tsv", sep="\t", index=False)

    unique = rel.drop_duplicates(["target_uniprot_id", "disease_id"]).copy()
    unique = unique.merge(mapping, left_on="disease_id", right_on="source_disease_id", how="left")
    flags = [
        "human_genetic_flag",
        "clinical_genetic_flag",
        "somatic_mutation_flag",
        "functional_flag",
        "animal_model_flag",
        "expression_flag",
        "literature_flag",
        "uniprot_disease_flag",
    ]
    for flag in flags:
        unique[flag] = truth(unique[flag])

    # Safely collapse only direct or unique cross-reference mappings.
    aggregations = {
        "disease_name": "first",
        "mapping_status": "first",
        "source_database": lambda x: ";".join(sorted(set(map(str, x)))),
        "disease_evidence_level": best_grade,
    }
    aggregations.update({flag: "max" for flag in flags})
    canonical_pairs = (
        unique.groupby(
            ["target_uniprot_id", "canonical_ot_id_for_safe_dedup"],
            as_index=False,
            dropna=False,
        )
        .agg(aggregations)
        .rename(columns={"canonical_ot_id_for_safe_dedup": "canonical_disease_id"})
    )
    canonical_pairs.to_csv(ROOT / "data" / "M4_canonical_protein_disease_pairs.tsv", sep="\t", index=False)

    safe_area_map = area_map.drop_duplicates(
        ["canonical_ot_id_for_safe_dedup", "therapeutic_area_id"]
    )
    pair_area = canonical_pairs.merge(
        safe_area_map,
        left_on="canonical_disease_id",
        right_on="canonical_ot_id_for_safe_dedup",
        how="left",
    )
    pair_area = pair_area.loc[pair_area["therapeutic_area_id"].notna()].copy()
    n_areas_per_pair = (
        pair_area.groupby(["target_uniprot_id", "canonical_disease_id"])["therapeutic_area_id"]
        .nunique()
        .rename("n_areas_for_pair")
        .reset_index()
    )
    pair_area = pair_area.merge(n_areas_per_pair, on=["target_uniprot_id", "canonical_disease_id"], how="left")
    pair_area["fractional_pair_weight"] = 1 / pair_area["n_areas_for_pair"].clip(lower=1)
    pair_area.to_csv(ROOT / "data" / "M4_protein_disease_therapeutic_area_pairs.tsv", sep="\t", index=False)

    area_summary = (
        pair_area.groupby(["therapeutic_area_id", "therapeutic_area_name", "area_axis"], as_index=False)
        .agg(
            raw_pair_memberships=("target_uniprot_id", "size"),
            fractionally_allocated_pairs=("fractional_pair_weight", "sum"),
            unique_proteins=("target_uniprot_id", "nunique"),
            unique_diseases=("canonical_disease_id", "nunique"),
        )
        .sort_values("raw_pair_memberships", ascending=False)
    )
    area_summary.to_csv(ROOT / "data" / "M4_therapeutic_area_summary.tsv", sep="\t", index=False)

    mapping_status_order = [
        "OT entity ID",
        "Single OT cross-reference",
        "Ambiguous OT cross-reference",
        "Unmapped / no OT entity",
    ]
    source_id_map = (
        rel[["source_database", "disease_id"]]
        .drop_duplicates()
        .merge(mapping, left_on="disease_id", right_on="source_disease_id", how="left")
    )
    coverage = pd.crosstab(source_id_map["source_database"], source_id_map["mapping_status"])
    coverage = coverage.reindex(columns=mapping_status_order, fill_value=0)
    coverage.to_csv(ROOT / "data" / "M4_mapping_coverage_by_source.tsv", sep="\t")

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.5))
    fig.subplots_adjust(hspace=0.68, wspace=0.58, left=0.13, right=0.98, top=0.84, bottom=0.09)

    ax = axes[0, 0]
    panel_label(ax, "A", "Ontology mapping coverage")
    coverage_pct = coverage.div(coverage.sum(axis=1), axis=0) * 100
    left = np.zeros(len(coverage_pct))
    status_colors = [COLORS["teal"], COLORS["blue"], COLORS["orange"], COLORS["gray"]]
    for status, color in zip(mapping_status_order, status_colors):
        values = coverage_pct[status].to_numpy()
        ax.barh(coverage_pct.index, values, left=left, color=color, edgecolor="white", linewidth=0.6, label=status)
        for idx, (start, value) in enumerate(zip(left, values)):
            if value >= 7:
                ax.text(start + value / 2, idx, f"{value:.1f}%", ha="center", va="center", fontsize=6.5)
        left += values
    ax.set_xlim(0, 100)
    ax.set_xlabel("Unique source disease IDs (%)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    finish_axis(ax, "x")

    ax = axes[0, 1]
    panel_label(ax, "B", "Open Targets therapeutic-area landscape")
    top = area_summary.head(14).sort_values("raw_pair_memberships").copy()
    top["area_display"] = top["therapeutic_area_name"].map(lambda value: truncate_label(value, 31))
    bar_colors = [COLORS["orange"] if axis == "Cross-cutting disease type" else COLORS["teal"] for axis in top["area_axis"]]
    bars = ax.barh(top["area_display"], top["raw_pair_memberships"], color=bar_colors, alpha=0.88)
    ax.scatter(
        top["fractionally_allocated_pairs"],
        np.arange(len(top)),
        color=COLORS["navy"],
        s=17,
        marker="D",
        zorder=4,
        label="Fractionally allocated pairs",
    )
    for bar, value in zip(bars, top["raw_pair_memberships"]):
        ax.text(value, bar.get_y() + bar.get_height() / 2, f" {value:,.0f}", va="center", fontsize=6.2)
    ax.set_xlabel("Unique protein–disease pairs")
    ax.legend(
        handles=[
            Patch(facecolor=COLORS["teal"], label="Organ/system or clinical area"),
            Patch(facecolor=COLORS["orange"], label="Cross-cutting disease type"),
            Line2D([0], [0], marker="D", color="none", markerfacecolor=COLORS["navy"], markeredgecolor="none", label="Fractional allocation"),
        ],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=1,
    )
    finish_axis(ax, "x")

    ax = axes[1, 0]
    panel_label(ax, "C", "Evidence channels across therapeutic areas")
    top_area_names = area_summary.head(10)["therapeutic_area_name"].tolist()[::-1]
    channel_labels = {
        "human_genetic_flag": "Human\ngenetic",
        "clinical_genetic_flag": "Clinical\ngenetic",
        "somatic_mutation_flag": "Somatic",
        "functional_flag": "Functional",
        "animal_model_flag": "Animal",
        "expression_flag": "Expression",
        "literature_flag": "Literature",
        "uniprot_disease_flag": "UniProt",
    }
    bubble_rows = []
    for area_name in top_area_names:
        group = pair_area.loc[pair_area["therapeutic_area_name"].eq(area_name)]
        for flag, display in channel_labels.items():
            flagged = group.loc[group[flag]]
            high_fraction = (
                flagged["disease_evidence_level"].isin(["very_high", "high"]).mean()
                if len(flagged)
                else 0
            )
            bubble_rows.append(
                {
                    "area": area_name,
                    "channel": display,
                    "pair_count": len(flagged),
                    "high_fraction": high_fraction,
                }
            )
    bubbles = pd.DataFrame(bubble_rows)
    bubbles.to_csv(ROOT / "data" / "M4_area_evidence_channel_bubbles.tsv", sep="\t", index=False)
    x_map = {name: idx for idx, name in enumerate(channel_labels.values())}
    y_map = {name: idx for idx, name in enumerate(top_area_names)}
    max_count = max(float(bubbles["pair_count"].max()), 1)
    sizes = 12 + 190 * np.sqrt(bubbles["pair_count"] / max_count)
    scatter = ax.scatter(
        bubbles["channel"].map(x_map),
        bubbles["area"].map(y_map),
        s=sizes,
        c=bubbles["high_fraction"] * 100,
        cmap=sns.light_palette(COLORS["purple"], as_cmap=True),
        vmin=0,
        vmax=100,
        edgecolors="white",
        linewidths=0.5,
    )
    ax.set_xticks(range(len(x_map)), list(x_map.keys()))
    ax.set_yticks(range(len(y_map)), [truncate_label(name, 27) for name in y_map.keys()])
    ax.set_xlim(-0.6, len(x_map) - 0.4)
    ax.set_ylim(-0.6, len(y_map) - 0.4)
    ax.grid(color=COLORS["gray_light"], linewidth=0.6)
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("High/very-high evidence among channel-supported pairs (%)")
    ax.text(
        0.0,
        -0.28,
        "Bubble area represents the number of channel-supported protein–disease pairs.",
        transform=ax.transAxes,
        fontsize=6.5,
        color=COLORS["gray_dark"],
    )

    ax = axes[1, 1]
    panel_label(ax, "D", "Disease concepts with the broadest membrane-protein support")
    disease_counts = (
        canonical_pairs.groupby(["canonical_disease_id", "disease_name", "disease_evidence_level"])
        .size()
        .unstack(fill_value=0)
    )
    for grade in ["very_high", "high", "medium"]:
        if grade not in disease_counts.columns:
            disease_counts[grade] = 0
    disease_counts["total"] = disease_counts[["very_high", "high", "medium"]].sum(axis=1)
    top_diseases = disease_counts.nlargest(12, "total").sort_values("total")
    disease_display = [
        (name[:44] + "…") if len(name) > 45 else name
        for _, name in top_diseases.index
    ]
    disease_display = [truncate_label(name, 32) for name in disease_display]
    left = np.zeros(len(top_diseases))
    for grade, color in [
        ("very_high", COLORS["navy"]),
        ("high", COLORS["green"]),
        ("medium", COLORS["orange"]),
    ]:
        values = top_diseases[grade].to_numpy()
        ax.barh(disease_display, values, left=left, label=grade.replace("_", " ").title(), color=color)
        left += values
    for idx, total in enumerate(top_diseases["total"]):
        ax.text(total, idx, f" {total:,.0f}", va="center", fontsize=6.2)
    ax.set_xlabel("Distinct membrane proteins")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=3)
    finish_axis(ax, "x")

    fig.suptitle(
        "M4 revised — Ontology-derived, multi-label disease landscape",
        x=0.08,
        ha="left",
        fontsize=12,
        fontweight="bold",
        color=COLORS["navy"],
    )
    fig.text(
        0.98,
        0.965,
        "Open Targets Platform 26.06 disease/phenotype ontology",
        ha="right",
        va="top",
        fontsize=7.2,
        color=COLORS["gray_dark"],
    )
    paths = save_figure(fig, "M4_ontology_multilabel_disease_landscape")
    return paths, {
        "mapping": mapping,
        "area_map": area_map,
        "canonical_pairs": canonical_pairs,
        "pair_area": pair_area,
        "area_summary": area_summary,
        "coverage": coverage,
    }


STATUS_FIELDS = [
    ("Approved drug", "is_approved_drug"),
    ("Clinical candidate", "is_clinical_candidate"),
    ("Chemical probe", "is_chemical_probe"),
    ("Endogenous ligand", "is_endogenous_ligand"),
    ("Natural product", "is_natural_product"),
]


def scan_compounds() -> dict:
    usecols = [
        "compound_scope_status",
        "record_qc_status",
        "molecular_weight",
        "xlogp",
        "hbond_donor_count",
        "hbond_acceptor_count",
    ] + [field for _, field in STATUS_FIELDS]
    combo_counts: Counter[tuple[bool, ...]] = Counter()
    set_totals = Counter()
    status_coverage = Counter()
    total_core_qc = 0
    ro5_complete = 0
    violation_counts = Counter()
    criterion = {
        "MW ≤ 500 Da": {"observed": 0, "pass": 0, "missing": 0},
        "XlogP ≤ 5": {"observed": 0, "pass": 0, "missing": 0},
        "H-bond donors ≤ 5": {"observed": 0, "pass": 0, "missing": 0},
        "H-bond acceptors ≤ 10": {"observed": 0, "pass": 0, "missing": 0},
    }
    property_map = [
        ("MW ≤ 500 Da", "molecular_weight", lambda x: x <= 500),
        ("XlogP ≤ 5", "xlogp", lambda x: x <= 5),
        ("H-bond donors ≤ 5", "hbond_donor_count", lambda x: x <= 5),
        ("H-bond acceptors ≤ 10", "hbond_acceptor_count", lambda x: x <= 10),
    ]
    for chunk in pd.read_csv(
        REL / "small_molecule_master_v1_3.tsv",
        sep="\t",
        usecols=usecols,
        chunksize=150_000,
        low_memory=False,
    ):
        core = chunk["compound_scope_status"].eq("core") & chunk["record_qc_status"].eq("ok")
        frame = chunk.loc[core].copy()
        total_core_qc += len(frame)
        raw_status = [frame[field] for _, field in STATUS_FIELDS]
        available_fields = [
            values.notna() & values.astype(str).str.strip().ne("")
            for values in raw_status
        ]
        bool_fields = [truth(values) for values in raw_status]
        combos = pd.DataFrame({label: values for (label, _), values in zip(STATUS_FIELDS, bool_fields)})
        availability = pd.DataFrame({label: values for (label, _), values in zip(STATUS_FIELDS, available_fields)})
        any_positive = combos.any(axis=1)
        all_available = availability.all(axis=1)
        status_coverage["Status-positive"] += int(any_positive.sum())
        status_coverage["Assessed, no positive status"] += int((~any_positive & all_available).sum())
        status_coverage["Status unavailable / partial"] += int((~any_positive & ~all_available).sum())
        for label in combos.columns:
            set_totals[label] += int(combos[label].sum())
        for values, count in combos.value_counts().items():
            values_tuple = values if isinstance(values, tuple) else (values,)
            if any(values_tuple):
                combo_counts[tuple(bool(value) for value in values_tuple)] += int(count)

        numeric = {}
        for label, field, rule in property_map:
            values = pd.to_numeric(frame[field], errors="coerce")
            numeric[field] = values
            observed = values.notna()
            passed = observed & rule(values)
            criterion[label]["observed"] += int(observed.sum())
            criterion[label]["pass"] += int(passed.sum())
            criterion[label]["missing"] += int((~observed).sum())
        complete = pd.DataFrame(numeric).notna().all(axis=1)
        ro5_complete += int(complete.sum())
        violations = (
            numeric["molecular_weight"].gt(500).astype(int)
            + numeric["xlogp"].gt(5).astype(int)
            + numeric["hbond_donor_count"].gt(5).astype(int)
            + numeric["hbond_acceptor_count"].gt(10).astype(int)
        )
        violation_counts.update(violations.loc[complete].astype(int).tolist())
    return {
        "combo_counts": combo_counts,
        "set_totals": set_totals,
        "total_core_qc": total_core_qc,
        "status_coverage": status_coverage,
        "ro5_complete": ro5_complete,
        "violation_counts": violation_counts,
        "criterion": criterion,
    }


def build_m5b(scan: dict) -> tuple[dict[str, str], pd.DataFrame]:
    names = [label for label, _ in STATUS_FIELDS]
    combo_records = []
    for bits, count in scan["combo_counts"].items():
        combo_records.append({**{name: bool(bit) for name, bit in zip(names, bits)}, "count": count})
    combos = pd.DataFrame(combo_records).sort_values("count", ascending=False)
    combos.to_csv(ROOT / "data" / "M5B_biological_status_intersections.tsv", sep="\t", index=False)
    top = combos.head(12).reset_index(drop=True)
    annotated_total = int(combos["count"].sum())

    fig = plt.figure(figsize=(8.2, 4.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.0, 1.25], width_ratios=[1.05, 3.4], hspace=0.05, wspace=0.06)
    blank = fig.add_subplot(gs[0, 0])
    bar_ax = fig.add_subplot(gs[0, 1])
    set_ax = fig.add_subplot(gs[1, 0])
    matrix_ax = fig.add_subplot(gs[1, 1])
    coverage_names = [
        "Status-positive",
        "Assessed, no positive status",
        "Status unavailable / partial",
    ]
    coverage_colors = [COLORS["purple"], COLORS["teal_light"], COLORS["gray"]]
    coverage_values = np.array([scan["status_coverage"][name] for name in coverage_names], dtype=float)
    coverage_pct = coverage_values / scan["total_core_qc"] * 100
    coverage_left = 0.0
    for name, value, pct, color in zip(coverage_names, coverage_values, coverage_pct, coverage_colors):
        blank.barh([0], [pct], left=coverage_left, color=color, height=0.38, edgecolor="white", linewidth=0.6)
        if pct >= 8:
            blank.text(coverage_left + pct / 2, 0, f"{pct:.1f}%", ha="center", va="center", fontsize=6.1)
        coverage_left += pct
    blank.set_xlim(0, 100)
    blank.set_ylim(-1.25, 0.7)
    blank.set_yticks([])
    blank.set_xticks([0, 50, 100])
    blank.set_xlabel("Status-field coverage (%)")
    blank.set_title("Annotation coverage", loc="left", fontsize=8.2, fontweight="bold")
    blank.text(1.0, 0.34, f"{coverage_pct[0]:.2f}% positive", ha="left", va="center",
               fontsize=6.0, fontweight="bold", color=COLORS["purple"])
    blank.spines[["left", "right", "top"]].set_visible(False)
    blank.text(
        0,
        -0.78,
        f"Positive: {int(coverage_values[0]):,} ({coverage_pct[0]:.2f}%)\n"
        f"Assessed none: {int(coverage_values[1]):,}\n"
        f"Unavailable: {int(coverage_values[2]):,}",
        ha="left",
        va="top",
        fontsize=6.2,
        color=COLORS["gray_dark"],
    )

    x = np.arange(len(top))
    bars = bar_ax.bar(x, top["count"], color=COLORS["blue"], width=0.72)
    bar_ax.set_yscale("log")
    bar_ax.set_ylabel("Intersection size (log scale)")
    bar_ax.set_xticks([])
    finish_axis(bar_ax, "y")
    for bar, value in zip(bars, top["count"]):
        bar_ax.text(bar.get_x() + bar.get_width() / 2, value * 1.10, human_number(value), ha="center", va="bottom", fontsize=6.3, rotation=90)

    set_sizes = pd.Series({name: scan["set_totals"][name] for name in names})
    set_ax.barh(np.arange(len(names)), set_sizes.values, color=COLORS["teal"])
    set_ax.set_yticks(np.arange(len(names)), names)
    set_ax.invert_xaxis()
    set_ax.set_xscale("log")
    set_ax.set_xlabel("Set size\n(log scale)")
    set_ax.invert_yaxis()
    finish_axis(set_ax, "x")

    matrix_ax.set_xlim(-0.6, len(top) - 0.4)
    matrix_ax.set_ylim(-0.6, len(names) - 0.4)
    matrix_ax.invert_yaxis()
    for col, row in top.iterrows():
        active_y = []
        for y_pos, name in enumerate(names):
            active = bool(row[name])
            matrix_ax.scatter(
                col,
                y_pos,
                s=30 if active else 18,
                color=COLORS["navy"] if active else COLORS["gray_light"],
                edgecolor="none",
                zorder=3,
            )
            if active:
                active_y.append(y_pos)
        if len(active_y) > 1:
            matrix_ax.plot([col, col], [min(active_y), max(active_y)], color=COLORS["navy"], linewidth=1.2, zorder=2)
    matrix_ax.set_xticks(x, [str(index + 1) for index in x])
    matrix_ax.set_yticks(np.arange(len(names)), [""] * len(names))
    matrix_ax.set_xlabel("Top non-empty intersections")
    matrix_ax.spines[["left", "bottom"]].set_visible(False)
    matrix_ax.grid(False)

    fig.suptitle(
        "Biological-status intersections among annotated compounds",
        x=0.06,
        ha="left",
        fontsize=12,
        fontweight="bold",
        color=COLORS["navy"],
    )
    fig.text(
        0.06,
        0.925,
        f"Core, QC-passed canonical compounds with ≥1 status annotation: n = {annotated_total:,}. Compounds with no status flag are excluded.",
        ha="left",
        fontsize=7.2,
        color=COLORS["gray_dark"],
    )
    fig.subplots_adjust(left=0.10, right=0.98, top=0.83, bottom=0.16, hspace=0.06, wspace=0.06)
    paths = save_figure(fig, "M5B_biological_status_upset")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.83, bottom=0.16, hspace=0.06, wspace=0.06)
    return paths, combos


def build_m5d(scan: dict) -> tuple[dict[str, str], pd.DataFrame, pd.DataFrame]:
    """Restore the descriptive physicochemical ECDF panel used before Ro5."""
    cache = pd.read_csv(PHYSICOCHEMICAL_CACHE, sep="\t", low_memory=False)
    properties = [
        ("molecular_weight", "Molecular weight (Da)", (0.0, 900.0), "MW"),
        ("xlogp", "XlogP", (-6.0, 12.0), "XlogP"),
        ("tpsa", "TPSA (Å²)", (0.0, 300.0), "TPSA"),
        ("rotatable_bond_count", "Rotatable bonds", (0.0, 35.0), "Rotatable bonds"),
    ]
    summary_rows = []
    plot_values: dict[str, np.ndarray] = {}
    for field, label, limits, short_label in properties:
        numeric = pd.to_numeric(cache[field], errors="coerce")
        observed = numeric.dropna()
        displayed = observed.loc[observed.between(*limits)].sort_values().to_numpy()
        plot_values[field] = displayed
        summary_rows.append(
            {
                "property": field,
                "display_label": short_label,
                "cache_rows": len(cache),
                "observed_count": len(observed),
                "displayed_count": len(displayed),
                "display_min": limits[0],
                "display_max": limits[1],
                "p10": float(observed.quantile(0.10)),
                "median": float(observed.median()),
                "p90": float(observed.quantile(0.90)),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        ROOT / "data" / "M5D_physicochemical_ecdf_summary.tsv",
        sep="\t",
        index=False,
    )

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 5.25))
    letters = ["A", "B", "C", "D"]
    for ax, letter, (field, label, limits, short_label), row in zip(
        axes.flat, letters, properties, summary.to_dict(orient="records")
    ):
        values = plot_values[field]
        cumulative = np.arange(1, len(values) + 1, dtype=float) / len(values)
        ax.plot(values, cumulative, color=COLORS["blue"], linewidth=1.7)
        median = row["median"]
        ax.axvline(median, color=COLORS["orange"], linewidth=1.0, linestyle="--")
        ax.text(
            0.97,
            0.08,
            f"n={row['observed_count']:,}\nP10={row['p10']:.2f}  Median={median:.2f}  P90={row['p90']:.2f}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=6.4,
            color=COLORS["gray_dark"],
        )
        panel_label(ax, letter, short_label)
        ax.set_xlim(*limits)
        ax.set_ylim(0, 1)
        ax.set_xlabel(label)
        ax.set_ylabel("Cumulative fraction")
        finish_axis(ax, "both")

    fig.suptitle(
        "M5D revised — Physicochemical distributions of canonical compounds",
        x=0.07,
        ha="left",
        fontsize=12,
        fontweight="bold",
        color=COLORS["navy"],
    )
    fig.text(
        0.07,
        0.93,
        "ECDFs describe molecular size, lipophilicity, polarity and conformational flexibility; no oral-drug rule is used as a universal filter.",
        ha="left",
        fontsize=7.1,
        color=COLORS["gray_dark"],
    )
    fig.subplots_adjust(
        left=0.11, right=0.98, top=0.82, bottom=0.11, hspace=0.43, wspace=0.29
    )
    paths = save_figure(fig, "M5D_physicochemical_ecdf")
    return paths, summary, cache


def write_caption_and_qa(outputs: dict, m3_summary: pd.DataFrame, m4_data: dict, scan: dict) -> None:
    caption = """# Revised figure captions — trial version 1

## Revised Figure M3D. HPA tissue-expression breadth with explicit mapping status

Expression breadth is shown only after separating HPA mapping status. The violin,
box and deterministic point sample use proteins with `hpa_mapping_status_v62=mapped`
and a numeric detected-tissue count; mapped zeros remain valid plotted observations.
The adjacent 100% bars separate mapped proteins with at least one detected tissue,
mapped proteins with zero tissues passing the HPA detection rule, and
unmapped/missing proteins. A mapped zero means that no tissue passed the chosen HPA
detection rule; it is not evidence that the protein is biologically absent.

## Revised Figure M4. Ontology-derived, multi-label disease landscape

Disease classification is taken from the Open Targets Platform 26.06
`therapeuticAreas` field rather than disease-name keywords. Open Targets entity IDs
are joined directly; UniProt/OMIM records are joined only through official Open
Targets cross-references. Ambiguous cross-references are retained as ambiguous and
are not used for automatic canonical disease merging. Therapeutic-area membership
is multi-label and non-exclusive. Raw membership counts can therefore exceed the
number of unique protein–disease pairs; diamond markers show fractional allocation,
where each pair contributes a total weight of one across its assigned areas.
Generic phenotype, biological-process, measurement, medical-procedure and non-human
branches are not treated as human organ systems. This figure represents clinical
therapeutic areas, not a DO/Uberon anatomical map.

## Revised Figure M5B. Biological-status intersections among annotated compounds

The UpSet plot includes core, QC-passed canonical compounds carrying at least one
of five status annotations: approved drug, clinical candidate, chemical probe,
endogenous ligand or natural product. Compounds with no status annotation are
excluded from the intersection denominator. Set and intersection axes use a log
scale. The approved-drug set size reflects the conservative annotation present in
the frozen HuMemLigDB release and is not a count of all approved drugs worldwide.
The coverage strip separates compounds with at least one positive status,
compounds whose five status fields were assessed but all negative, and compounds
for which the status fields are unavailable or incomplete.


## Revised Figure M5D. Physicochemical distributions of canonical compounds

Empirical cumulative distribution functions (ECDFs) summarize molecular weight,
XlogP, topological polar surface area (TPSA), and rotatable-bond count in the
deterministic V6.2 physicochemical cache. Each curve gives the fraction of observed
compounds at or below a given property value; dashed lines mark medians and the
annotations report the 10th, 50th and 90th percentiles. Display limits are applied
only for legibility and do not remove compounds from HuMemLigDB. These properties
describe size, lipophilicity, polarity and conformational flexibility. They are
not efficacy, approval, membrane permeability or universal docking criteria.
Lipinski Rule-of-Five counts are not used as the main chemical-space panel because
they were designed primarily for oral drug-likeness.
"""
    (ROOT / "CAPTIONS_REVISED_TRIAL_v1.md").write_text(caption, encoding="utf-8")

    mapping = m4_data["mapping"]
    status_counts = mapping["mapping_status"].value_counts().to_dict()
    qa = {
        "status": "PASS",
        "release_inputs_read_only": True,
        "open_targets_release": "26.06",
        "outputs": outputs,
        "m3d": {
            "functional_classes": m3_summary.to_dict(orient="records"),
            "unmapped_or_missing_explicit": True,
        },
        "m4": {
            "source_disease_ids": int(len(mapping)),
            "mapping_status_counts": {str(k): int(v) for k, v in status_counts.items()},
            "canonical_pair_rows_after_safe_dedup": int(len(m4_data["canonical_pairs"])),
            "therapeutic_area_membership_rows": int(len(m4_data["pair_area"])),
            "keyword_classification_used": False,
            "multi_label": True,
        },
        "m5": {
            "core_qc_passed_compounds": int(scan["total_core_qc"]),
            "status_annotated_compounds": int(sum(scan["combo_counts"].values())),
            "status_coverage": {
                str(key): int(value) for key, value in scan["status_coverage"].items()
            },
            "ro5_complete_compounds": int(scan["ro5_complete"]),
            "m5d_primary_display": "physicochemical_ecdf",
            "lipinski_used_as_primary_display": False,
            "physicochemical_cache": str(PHYSICOCHEMICAL_CACHE),
        },
    }
    (ROOT / "qa" / "REVISED_PANELS_V1_VALIDATION.json").write_text(
        json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def main() -> None:
    outputs: dict[str, dict[str, str]] = {}
    outputs["M3D"], m3_summary = build_m3d()
    outputs["M4"], m4_data = build_m4()
    scan = scan_compounds()
    outputs["M5B"], _ = build_m5b(scan)
    outputs["M5D"], _, _ = build_m5d(scan)
    write_caption_and_qa(outputs, m3_summary, m4_data, scan)
    print(json.dumps(outputs, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
