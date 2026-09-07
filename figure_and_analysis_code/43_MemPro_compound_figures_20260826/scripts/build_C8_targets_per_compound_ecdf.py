"""C8: Unique membrane-protein targets per canonical compound.

This is a one-table, one-figure analysis. It reads the frozen canonical
protein–compound pair table, does not use evidence-record multiplicity, and
does not alter pair inclusion or compound identities.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compound_figure_pool"
PAIR = Path(r"D:\finale\FORMAL\01_core_v72\01_release_tables\protein_compound_pair_v72.tsv.gz")

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.linewidth": 0.8,
})

INK = "#24364B"
SLATE = "#5E7182"
LINE = "#173B6C"
ACCENT = "#D6A23A"


def nearest_quantile(values: np.ndarray, q: float) -> int:
    return int(np.quantile(values, q, method="nearest"))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(
        PAIR,
        sep="\t",
        compression="gzip",
        usecols=["pair_id", "target_uniprot_id", "compound_internal_id", "default_web_inclusion_v72", "searchable_inclusion_v72"],
        low_memory=False,
    )
    # The frozen pair table is already canonical. Duplicate protection is
    # retained so the analysis cannot be inflated by future source joins.
    deduplicated = pairs.drop_duplicates(["compound_internal_id", "target_uniprot_id"]).copy()
    target_counts = (
        deduplicated.groupby("compound_internal_id", as_index=False)
        .agg(
            unique_target_count=("target_uniprot_id", "nunique"),
            unique_pair_count=("pair_id", "nunique"),
            default_web_pair_count=("default_web_inclusion_v72", "sum"),
            searchable_pair_count=("searchable_inclusion_v72", "sum"),
        )
        .sort_values(["unique_target_count", "compound_internal_id"], kind="stable")
        .reset_index(drop=True)
    )
    target_counts["is_single_target"] = target_counts["unique_target_count"].eq(1)
    target_counts.to_csv(OUT / "C8_targets_per_compound.tsv", sep="\t", index=False)

    values = target_counts["unique_target_count"].to_numpy()
    total = int(values.size)
    summary = pd.DataFrame([{
        "analysis_unit": "one canonical compound",
        "frozen_pair_rows": int(len(pairs)),
        "duplicate_compound_target_rows_removed": int(len(pairs) - len(deduplicated)),
        "analysis_compounds": total,
        "unique_membrane_proteins_represented": int(deduplicated["target_uniprot_id"].nunique()),
        "median_targets_per_compound": nearest_quantile(values, 0.50),
        "p75_targets_per_compound": nearest_quantile(values, 0.75),
        "p90_targets_per_compound": nearest_quantile(values, 0.90),
        "p95_targets_per_compound": nearest_quantile(values, 0.95),
        "maximum_targets_per_compound": int(values.max()),
        "single_target_compounds": int((values == 1).sum()),
        "single_target_fraction": float((values == 1).mean()),
        "pair_inclusion_rule": "All rows in frozen canonical protein_compound_pair_v72 table; default-web status retained in source data but not used as a filter.",
    }])
    summary.to_csv(OUT / "C8_targets_per_compound_summary.tsv", sep="\t", index=False)

    x = np.sort(values)
    ecdf = np.arange(1, total + 1) / total
    fig = plt.figure(figsize=(178 / 25.4, 94 / 25.4), dpi=600)
    ax = fig.add_axes([0.14, 0.20, 0.74, 0.68])
    ax.step(x, ecdf, where="post", color=LINE, linewidth=1.55)
    ax.fill_between(x, ecdf, step="post", color=LINE, alpha=0.09)
    ax.set_xscale("log")
    ax.set_xlim(0.85, values.max() * 1.25)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel("Unique membrane proteins per canonical compound", fontsize=8, color=INK)
    ax.set_ylabel("Cumulative fraction of canonical compounds", fontsize=8, color=INK)
    ax.tick_params(axis="both", labelsize=8, colors=INK)
    ax.set_yticks(np.linspace(0, 1, 6), [f"{v:.0%}" for v in np.linspace(0, 1, 6)])
    ax.grid(which="major", color="#E4EAED", linewidth=0.65)
    ax.grid(which="minor", axis="x", color="#F1F4F5", linewidth=0.35)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    median = int(summary.at[0, "median_targets_per_compound"])
    p75 = int(summary.at[0, "p75_targets_per_compound"])
    p90 = int(summary.at[0, "p90_targets_per_compound"])
    p95 = int(summary.at[0, "p95_targets_per_compound"])
    annotation_positions = {
        median: (1.16, 0.54),
        p75: (2.18, 0.73),
        p90: (3.35, 0.84),
        p95: (4.85, 0.975),
    }
    # Collapsed quantiles (if present) share one biologically meaningful line;
    # the fixed positions avoid label collisions at the steep ECDF head.
    shown_values = []
    for value, label in [(median, "Median"), (p75, "75th percentile"), (p90, "90th percentile"), (p95, "95th percentile")]:
        if value not in shown_values:
            shown_values.append(value)
        else:
            continue
        frac = float((values <= value).mean())
        ax.axvline(value, color="#AAB7C0", linewidth=0.75, linestyle="--", zorder=1)
        ax.scatter(value, frac, s=26, color=ACCENT, edgecolor="white", linewidth=0.7, zorder=4)
        ax.annotate(
            f"{label}: {value}",
            xy=(value, frac),
            xytext=annotation_positions[value],
            fontsize=7,
            color=INK,
            arrowprops={"arrowstyle": "-", "color": SLATE, "linewidth": 0.65},
        )
    single_fraction = float(summary.at[0, "single_target_fraction"])
    ax.text(
        0.02,
        0.43,
        f"Single-target compounds: {single_fraction:.1%}",
        transform=ax.transAxes,
        fontsize=7.2,
        color=INK,
        va="top",
    )
    ax.text(
        0,
        1.03,
        f"Analysis N = {total:,} canonical compounds; {deduplicated['target_uniprot_id'].nunique():,} membrane proteins",
        transform=ax.transAxes,
        fontsize=7.2,
        color=SLATE,
    )
    ax.text(
        0,
        -0.23,
        "Each compound–protein pair is counted once; x-axis is logarithmic. A compound connected to several roles remains one compound in this distribution.",
        transform=ax.transAxes,
        fontsize=7.0,
        color=SLATE,
    )
    for suffix, kwargs in {".svg": {}, ".pdf": {}, ".png": {"dpi": 600}}.items():
        fig.savefig(OUT / f"C8_targets_per_compound_ecdf{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)

    (OUT / "C8_targets_per_compound_statistics.txt").write_text(
        "\n".join([
            "Figure\tC8 — Target breadth per canonical compound",
            f"Analysis N\t{total:,} canonical compounds",
            f"Input pair rows\t{len(pairs):,}",
            f"Unique canonical pairs\t{len(deduplicated):,}",
            f"Single-target compounds\t{int(summary.at[0, 'single_target_compounds']):,} ({single_fraction:.4%})",
            f"Median / P75 / P90 / P95\t{median} / {p75} / {p90} / {p95}",
            f"Maximum\t{int(summary.at[0, 'maximum_targets_per_compound'])}",
            "Missingness\tNo pair-table compounds were excluded; compounds without formal positive pairs are outside the defined analysis denominator.",
        ]) + "\n",
        encoding="utf-8",
    )
    (OUT / "C8_reproducibility_manifest.json").write_text(json.dumps({
        "figure": "C8_targets_per_compound_ecdf",
        "analysis_unit": "one canonical compound",
        "input": str(PAIR),
        "deduplication_key": ["compound_internal_id", "target_uniprot_id"],
        "output": ["C8_targets_per_compound_ecdf.svg", "C8_targets_per_compound_ecdf.pdf", "C8_targets_per_compound_ecdf.png", "C8_targets_per_compound.tsv", "C8_targets_per_compound_summary.tsv"],
        "x_scale": "logarithmic native target count",
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
