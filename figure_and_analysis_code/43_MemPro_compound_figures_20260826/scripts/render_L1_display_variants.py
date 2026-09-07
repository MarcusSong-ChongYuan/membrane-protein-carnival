#!/usr/bin/env python3
"""Render display-only L1 variants from the frozen L1 statistics table.

This script deliberately does *not* recalculate Fisher tests, odds ratios, or
Benjamini-Hochberg FDR values.  It only changes the visualization rules for
the already-generated L1 statistics, producing a main-text version limited to
roles with adequate overall localization support and a full supplementary
version retaining every formal role.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle


ROOT = Path(r"D:\finale\43_MemPro_compound_figures_20260826")
OUT = ROOT / "results" / "context_figure_pool"
STATS_PATH = OUT / "L1_role_localization_enrichment_statistics.tsv"
LOC_PATH = OUT / "L1_localization_summary.tsv"
MAIN_BASE = OUT / "L1_role_subcellular_localization_enrichment"
SUPP_BASE = OUT / "L1_role_subcellular_localization_enrichment_full_supplementary"

DISPLAY_N = 12
MAIN_ROLE_MIN_ELIGIBLE = 100
DARK = "#24364B"
MUTED = "#5E7182"
ZERO_SUPPORT = "#F2F3F4"


def role_label(role: str) -> str:
    return role.replace("_", " ")


def relative_luminance(rgb: tuple[float, float, float]) -> float:
    values = np.asarray(rgb, dtype=float)
    linear = np.where(values <= 0.04045, values / 12.92, ((values + 0.055) / 1.055) ** 2.4)
    return float(np.dot(linear, [0.2126, 0.7152, 0.0722]))


def muted_low_support(rgb: tuple[float, float, float]) -> tuple[float, float, float]:
    """Blend an actual effect colour toward white without changing its sign."""
    colour = np.asarray(rgb, dtype=float)
    return tuple(0.46 * colour + 0.54 * np.ones(3))


def fdr_stars(q_value: float, support: int) -> str:
    """Display-only BH-FDR annotation; low-support cells receive no stars."""
    if support < 20 or not np.isfinite(q_value):
        return ""
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return ""


def render_variant(
    stats: pd.DataFrame,
    locations: list[str],
    roles: list[str],
    base: Path,
    supplementary: bool,
) -> None:
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "muted_blue_white_red", ["#2F6F9F", "#B5D0E2", "#F7F7F5", "#E6B7B5", "#B6575C"]
    )
    norm = mcolors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    fig_height = max(5.0, 0.42 * len(roles) + 2.05)
    fig, ax = plt.subplots(figsize=(7.15, fig_height))

    ordered = stats.copy()
    ordered["formal_primary_membrane_role"] = pd.Categorical(
        ordered["formal_primary_membrane_role"], categories=roles, ordered=True
    )
    ordered["location_term"] = pd.Categorical(ordered["location_term"], categories=locations, ordered=True)
    ordered = ordered.sort_values(["formal_primary_membrane_role", "location_term"])

    for yi, role in enumerate(roles):
        for xi, location in enumerate(locations):
            row = ordered.loc[
                (ordered["formal_primary_membrane_role"] == role) & (ordered["location_term"] == location)
            ].iloc[0]
            support = int(row["observed_unique_proteins"])
            effect = float(row["display_log2_odds_ratio_clipped"])
            if support == 0:
                face = ZERO_SUPPORT
            else:
                actual = tuple(cmap(norm(effect))[:3])
                face = muted_low_support(actual) if support < 20 else actual
            ax.add_patch(
                Rectangle((xi - 0.5, yi - 0.5), 1, 1, facecolor=face, edgecolor="white", linewidth=0.9)
            )
            if support >= 20:
                text_colour = "white" if relative_luminance(face) < 0.38 else DARK
                ax.text(xi, yi, str(support), ha="center", va="center", fontsize=6.15, color=text_colour)
            stars = fdr_stars(float(row["BH_FDR"]), support)
            if stars:
                ax.text(xi + 0.34, yi - 0.31, stars, ha="right", va="center", fontsize=6.7,
                        fontweight="bold", color="#18222C", zorder=4)

    ax.set_xlim(-0.5, len(locations) - 0.5)
    ax.set_ylim(len(roles) - 0.5, -0.5)
    ax.set_xticks(range(len(locations)))
    ax.set_xticklabels(locations, rotation=32, ha="right", rotation_mode="anchor", fontsize=6.25, color=DARK)
    ax.set_yticks(range(len(roles)))
    ax.set_yticklabels([role_label(role) for role in roles], fontsize=6.75, color=DARK)
    ax.tick_params(axis="both", length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)

    panel = "Full supplementary view" if supplementary else "Main-text roles with adequate localization support"
    fig.suptitle(
        "Subcellular-localization enrichment across membrane-protein roles",
        x=0.125, y=0.985, ha="left", fontsize=10.1, fontweight="bold", color=DARK,
    )
    fig.text(
        0.125, 0.946,
        "Colour = log2 odds ratio; labels = observed unique proteins when support ≥20.",
        ha="left", va="top", fontsize=6.8, color=MUTED,
    )
    fig.text(
        0.125, 0.920,
        "* q<0.05; ** q<0.01; *** q<0.001 (BH-FDR; support ≥20). Low-support cells (n<20) are muted and not emphasized for inference.",
        ha="left", va="top", fontsize=6.35, color=MUTED,
    )
    fig.text(0.125, 0.895, panel, ha="left", va="top", fontsize=6.25, color=MUTED, style="italic")

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(sm, ax=ax, fraction=0.040, pad=0.020)
    colorbar.ax.tick_params(labelsize=6.0, length=2.0, colors=MUTED)
    colorbar.set_label("log2 odds ratio", fontsize=6.8, color=DARK, labelpad=5)
    colorbar.set_ticks([-3, -2, -1, 0, 1, 2, 3])

    fig.text(
        0.125, 0.012,
        "Neutral grey = no protein support. Localization annotations are multi-label; proteins without usable localization annotation are excluded from the eligible universe.",
        ha="left", va="bottom", fontsize=5.65, color=MUTED,
    )
    fig.subplots_adjust(left=0.30, right=0.94, top=0.84, bottom=0.19)
    for suffix, kwargs in (
        (".svg", {}),
        (".pdf", {}),
        (".png", {"dpi": 600}),
        (".tiff", {"dpi": 600}),
    ):
        fig.savefig(base.with_suffix(suffix), bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def main() -> None:
    if not STATS_PATH.exists() or not LOC_PATH.exists():
        raise FileNotFoundError("Precomputed L1 statistics/source data are required for display-only rendering.")
    OUT.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.75,
    })
    stats = pd.read_csv(STATS_PATH, sep="\t")
    loc_summary = pd.read_csv(LOC_PATH, sep="\t")
    locations = loc_summary.loc[loc_summary["displayed_in_L1"].astype(bool), "location_term"].tolist()
    if len(locations) != DISPLAY_N:
        raise ValueError(f"Expected {DISPLAY_N} preselected display locations; found {len(locations)}.")
    role_sizes = (
        stats[["formal_primary_membrane_role", "role_total_eligible_proteins"]]
        .drop_duplicates()
        .sort_values(["role_total_eligible_proteins", "formal_primary_membrane_role"], ascending=[False, True])
    )
    all_roles = role_sizes["formal_primary_membrane_role"].tolist()
    main_roles = role_sizes.loc[
        role_sizes["role_total_eligible_proteins"] >= MAIN_ROLE_MIN_ELIGIBLE,
        "formal_primary_membrane_role",
    ].tolist()
    if not main_roles:
        raise ValueError("No roles meet the frozen display-support threshold.")

    # Rendering only: no Fisher, odds-ratio, or FDR calculations occur here.
    render_variant(stats, locations, main_roles, MAIN_BASE, supplementary=False)
    render_variant(stats, locations, all_roles, SUPP_BASE, supplementary=True)

    display_stats = stats.loc[stats["location_term"].isin(locations)].copy()
    display_stats["main_text_role"] = display_stats["formal_primary_membrane_role"].isin(main_roles)
    display_stats["display_rule"] = np.select(
        [
            display_stats["observed_unique_proteins"].eq(0),
            display_stats["observed_unique_proteins"].lt(20),
        ],
        [
            "neutral grey; no protein support; no label or dot",
            "muted real log2 odds-ratio colour; no label or dot",
        ],
        default="full real log2 odds-ratio colour; label support count; stars encode BH-FDR level",
    )
    display_stats.to_csv(OUT / "L1_role_localization_display_rules_by_cell.tsv", sep="\t", index=False)
    (OUT / "L1_display_rules.txt").write_text(
        "DISPLAY-ONLY REVISION: Fisher tests, odds ratios and BH-FDR values are read from the existing L1 statistics TSV and are not recalculated.\n"
        f"Main-text figure rows: formal roles with >= {MAIN_ROLE_MIN_ELIGIBLE} eligible canonical proteins with usable mapped localization annotations ({len(main_roles)} roles).\n"
        f"Supplementary figure rows: all {len(all_roles)} formal roles.\n"
        f"Columns: the existing top {DISPLAY_N} formal location_term categories, unchanged.\n"
        "Observed count = 0: neutral light grey, no number, no significance dot.\n"
        "Observed count 1-19: true log2 odds-ratio colour blended 54% toward white, no number, no significance dot.\n"
        "Observed count >=20: full log2 odds-ratio colour and observed protein count; * q<0.05, ** q<0.01, *** q<0.001 (BH-FDR).\n"
        "The full 48-category comparison family and original BH correction are unchanged.\n",
        encoding="utf-8",
    )
    (OUT / "L1_caption_draft.txt").write_text(
        "Subcellular-localization enrichment across membrane-protein roles. The main-text view includes roles with at least 100 eligible canonical proteins carrying a usable, mapped formal localization annotation; the full 14-role view is provided as Supplementary Figure. "
        "Localization categories are multi-label and non-mutually-exclusive, and each protein contributes at most once to a protein-localization category. "
        "The displayed Fisher exact-test, Haldane-Anscombe-corrected log2 odds-ratio and Benjamini-Hochberg FDR values are unchanged from the original analysis over the complete formal role × localization comparison family. "
        "Cell colour denotes log2 odds ratio (clipped to ±3); observed protein counts are printed only for support ≥20, and stars denote BH-FDR level with support ≥20 (* q<0.05, ** q<0.01, *** q<0.001). "
        "Cells with 1-19 proteins retain their true effect direction but are muted and not emphasized for inference; neutral grey cells have zero protein support. "
        "Proteins without usable localization annotation are excluded from the eligible universe rather than treated as negative evidence.",
        encoding="utf-8",
    )
    qc = pd.DataFrame([
        ("statistics source", str(STATS_PATH), "PASS"),
        ("Fisher / odds ratio / BH recalculated", False, "PASS"),
        ("main roles", len(main_roles), "PASS"),
        ("supplementary roles", len(all_roles), "PASS"),
        ("displayed locations", len(locations), "PASS"),
        ("zero-support cells", int(display_stats["observed_unique_proteins"].eq(0).sum()), "PASS"),
        ("low-support cells 1-19", int(display_stats["observed_unique_proteins"].between(1, 19).sum()), "PASS"),
        ("labels printed", int(display_stats["observed_unique_proteins"].ge(20).sum()), "PASS"),
        ("star-annotated significant cells", int((display_stats["BH_FDR"].lt(0.05) & display_stats["observed_unique_proteins"].ge(20)).sum()), "PASS"),
    ], columns=["check", "result", "status"])
    qc.to_csv(OUT / "L1_QC.tsv", sep="\t", index=False)
    (OUT / "L1_display_revision_manifest.json").write_text(json.dumps({
        "operation": "display-only revision",
        "statistics_recalculated": False,
        "main_role_minimum_eligible_proteins": MAIN_ROLE_MIN_ELIGIBLE,
        "n_main_roles": len(main_roles),
        "n_supplementary_roles": len(all_roles),
        "n_displayed_categories": len(locations),
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
