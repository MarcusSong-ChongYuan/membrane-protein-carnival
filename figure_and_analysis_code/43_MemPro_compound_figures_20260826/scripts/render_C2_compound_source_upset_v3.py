"""Source-aware visual-only C2 UpSet redraw using frozen, precomputed tables."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compound_figure_pool"
MEMBERSHIP = OUT / "C2_compound_source_membership.tsv"
INTERSECTIONS = OUT / "C2_compound_source_intersections.tsv"
SUMMARY = OUT / "C2_source_overlap_summary.tsv"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.linewidth": .8,
})

SOURCES = [
    "PubChem BioAssay", "ChEMBL", "BindingDB", "PDBe",
    "IUPHAR/BPS Guide to PHARMACOLOGY", "PDBbind", "BRENDA", "sc-PDB",
]
DISPLAY = {"IUPHAR/BPS Guide to PHARMACOLOGY": "GtoPdb (IUPHAR/BPS)"}
NAVY, TEAL, LIGHT, INK, SLATE = "#173B6C", "#168C8C", "#DCE3E8", "#24364B", "#5E7182"


def tokens(combination):
    return str(combination).split("|")


def contains(combination, source):
    return source in tokens(combination)


def fmt(value):
    value = int(round(value))
    return f"{value / 1_000_000:.1f}M" if value >= 1_000_000 else (f"{value / 1_000:.1f}k" if value >= 1_000 else str(value))


def main():
    # Read the three already-calculated C2 tables only. No source membership or
    # intersection count is recomputed from the formal database.
    membership = pd.read_csv(MEMBERSHIP, sep="\t")
    intersections = pd.read_csv(INTERSECTIONS, sep="\t")
    _ = pd.read_csv(SUMMARY, sep="\t")
    analysis_n = len(membership)

    # Omit the aggregated "Other sources" display construct and combinations
    # containing hidden sources; the displayed dot pattern then exactly matches
    # the complete source combination printed in the source data.
    eligible = intersections.loc[
        ~intersections["source_combination"].str.contains("Other sources", regex=False)
    ].copy()
    eligible = eligible.loc[eligible["source_combination"].map(lambda combo: set(tokens(combo)).issubset(SOURCES))]
    eligible = eligible.sort_values("canonical_compound_count", ascending=False).reset_index(drop=True)

    selected = eligible.head(8).copy()
    # Guarantee that every displayed external source has at least one solid dot.
    for source in SOURCES:
        if not selected["source_combination"].map(lambda combo: contains(combo, source)).any():
            representative = eligible.loc[eligible["source_combination"].map(lambda combo: contains(combo, source))].head(1)
            selected = pd.concat([selected, representative], ignore_index=True)
    selected = selected.drop_duplicates("source_combination")

    # Retain meaningful singleton combinations when absent, then fill to a
    # compact 12 columns with the next largest exact intersections.
    for source in SOURCES:
        singleton = eligible.loc[eligible["source_combination"].eq(source)]
        if len(singleton) and singleton.iloc[0]["source_combination"] not in set(selected["source_combination"]):
            selected = pd.concat([selected, singleton.head(1)], ignore_index=True)
    selected = selected.drop_duplicates("source_combination")
    remaining = eligible.loc[~eligible["source_combination"].isin(selected["source_combination"])]
    if len(selected) < 12:
        selected = pd.concat([selected, remaining.head(12 - len(selected))], ignore_index=True)
    selected = selected.drop_duplicates("source_combination").sort_values("canonical_compound_count", ascending=False).head(14).reset_index(drop=True)
    selected.insert(0, "display_order", np.arange(1, len(selected) + 1))
    selected.to_csv(OUT / "C2_displayed_intersections.tsv", sep="\t", index=False)

    source_coverage = membership.drop(columns="compound_internal_id").sum().reindex(SOURCES).fillna(0).astype(int)
    qc_rows, contribution_rows = [], []
    for source in SOURCES:
        source_rows = intersections.loc[intersections["source_combination"].map(lambda combo: contains(combo, source))]
        largest = source_rows.sort_values("canonical_compound_count", ascending=False).iloc[0]
        source_only = intersections.loc[intersections["source_combination"].eq(source), "canonical_compound_count"]
        source_only_n = int(source_only.iloc[0]) if len(source_only) else 0
        total = int(source_coverage.loc[source])
        displayed = bool(selected["source_combination"].map(lambda combo: contains(combo, source)).any())
        qc_rows.append({
            "source": source, "total_compound_coverage": total, "source_only_count": source_only_n,
            "largest_intersection_containing_source": largest["source_combination"],
            "largest_intersection_count": int(largest["canonical_compound_count"]),
            "displayed_in_main_figure": displayed,
        })
        contribution_rows.append({
            "source": source, "total_compound_coverage": total, "source_only_compounds": source_only_n,
            "source_only_fraction": source_only_n / total if total else np.nan,
            "compounds_shared_with_at_least_1_other_source": total - source_only_n,
            "shared_fraction": (total - source_only_n) / total if total else np.nan,
            "largest_partner_or_combination": largest["source_combination"],
            "largest_intersection_count": int(largest["canonical_compound_count"]),
        })
    qc = pd.DataFrame(qc_rows)
    contributions = pd.DataFrame(contribution_rows)
    qc.to_csv(OUT / "C2_source_intersection_visibility_qc.tsv", sep="\t", index=False)
    contributions.to_csv(OUT / "C2_source_contribution_summary.tsv", sep="\t", index=False)
    assert qc["displayed_in_main_figure"].all(), "A displayed source is missing from the matrix"

    x = np.arange(len(selected))
    singleton = selected["source_count_display"].eq(1).to_numpy()
    fig = plt.figure(figsize=(1700 / 150, 1050 / 150), dpi=600)
    ax_bar = fig.add_axes([.43, .63, .49, .22])
    ax_matrix = fig.add_axes([.43, .21, .49, .34])
    ax_sets = fig.add_axes([.07, .21, .22, .34])

    bar_colors = np.where(singleton, TEAL, NAVY)
    ax_bar.bar(x, selected["canonical_compound_count"], color=bar_colors, width=.66)
    ax_bar.set_yscale("log")
    ax_bar.set_ylim(10, 700_000)
    ticks = [10, 100, 1_000, 10_000, 100_000, 500_000]
    ax_bar.set_yticks(ticks)
    ax_bar.yaxis.set_major_formatter(FuncFormatter(lambda value, _: fmt(value)))
    ax_bar.set(xticks=[], ylabel="Intersection size\n(log scale)")
    ax_bar.tick_params(axis="y", labelsize=8, colors=INK)
    ax_bar.spines[["top", "right", "bottom"]].set_visible(False)
    ax_bar.grid(axis="y", which="major", color="#E4EAED", lw=.65)
    ax_bar.set_axisbelow(True)
    top6 = set(selected.nlargest(6, "canonical_compound_count").index)
    for i, row in selected.iterrows():
        if i in top6 or row["source_count_display"] == 1:
            ax_bar.text(i, row["canonical_compound_count"] * 1.25, fmt(row["canonical_compound_count"]), ha="center", va="bottom", fontsize=7.6, color=INK)

    y = np.arange(len(SOURCES))
    ax_matrix.set(xlim=(-.6, len(selected) - .4), ylim=(len(SOURCES) - .5, -.5), xticks=x, yticks=y, yticklabels=[DISPLAY.get(source, source) for source in SOURCES])
    ax_matrix.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_matrix.tick_params(axis="y", labelsize=8.5, length=0, pad=8, colors=INK)
    ax_matrix.spines[:].set_visible(False)
    for i, row in selected.iterrows():
        included = [SOURCES.index(source) for source in tokens(row["source_combination"]) if source in SOURCES]
        ax_matrix.scatter([i] * len(SOURCES), y, s=37, color=LIGHT, zorder=1)
        mark_color = TEAL if row["source_count_display"] == 1 else NAVY
        ax_matrix.scatter([i] * len(included), included, s=78, color=mark_color, zorder=3)
        if len(included) > 1:
            ax_matrix.plot([i, i], [min(included), max(included)], color=NAVY, lw=1.65, zorder=2, solid_capstyle="round")

    ax_sets.barh(y, source_coverage.values, color=NAVY, height=.56)
    ax_sets.set_xscale("log")
    ax_sets.set(xlim=(800, 900_000), yticks=y, yticklabels=[""] * len(SOURCES), xlabel="Source coverage\n(log scale)")
    ax_sets.set_xticks([1_000, 10_000, 100_000, 500_000])
    ax_sets.xaxis.set_major_formatter(FuncFormatter(lambda value, _: fmt(value)))
    ax_sets.invert_yaxis()
    ax_sets.tick_params(axis="x", labelsize=7.5, colors=INK)
    ax_sets.tick_params(axis="y", length=0)
    ax_sets.spines[["top", "right", "left"]].set_visible(False)
    for i, value in enumerate(source_coverage.values):
        ax_sets.text(value * 1.12, i, fmt(value), va="center", fontsize=7.7, color=INK)

    fig.text(.43, .955, "Compound-source overlap and complementarity", fontsize=14.5, weight="bold", color=INK, ha="left")
    fig.text(.43, .910, "Canonical compounds are counted once per direct external source.", fontsize=8.8, color=SLATE, ha="left")
    fig.text(.43, .875, f"Analysis N = {analysis_n:,} canonical compounds with direct source provenance", fontsize=8.2, color=SLATE, ha="left")
    fig.text(.43, .842, "Source-aware intersections are displayed; complete intersections are available in Source Data.", fontsize=7.7, color=SLATE, ha="left")
    fig.text(.43, .085, "Teal denotes source-only intersections; navy denotes multi-source intersections.", fontsize=7.5, color=SLATE, ha="left")
    for suffix, kwargs in {".svg": {}, ".pdf": {}, ".png": {"dpi": 600}}.items():
        fig.savefig(OUT / f"C2_compound_source_upset_v3{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
