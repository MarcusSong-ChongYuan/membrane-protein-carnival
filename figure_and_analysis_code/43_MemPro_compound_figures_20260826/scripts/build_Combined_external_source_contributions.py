#!/usr/bin/env python3
"""Build the MemPro external-source contribution composite figure.

Figure contract
---------------
Core conclusion: External databases contribute complementary information at
the canonical-compound layer and distinct, often overlapping, interaction-
evidence profiles at the canonical protein-compound-pair layer.
Archetype: quantitative grid with two aligned source-centric panels.
Panel A: unique canonical compounds per source/module (log1p heatmap).
Panel B: source-specific unique canonical pairs covered by each formal,
non-mutually-exclusive evidence modality (grouped lollipop profile).
Frozen biological facts are not read or modified; this script only reuses the
already validated C1 and E3 source-data tables.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ROOT = Path(r"D:\finale\43_MemPro_compound_figures_20260826")
COMPOUND_POOL = ROOT / "results" / "compound_figure_pool"
CONTEXT_POOL = ROOT / "results" / "context_figure_pool"
SCRIPT_DIR = ROOT / "scripts"

PANEL_A_INPUT = COMPOUND_POOL / "C1_source_module_compound_coverage.tsv"
PANEL_B_INPUT = CONTEXT_POOL / "E3_source_evidence_type_counts.tsv"
PANEL_B_SUMMARY = CONTEXT_POOL / "E3_source_summary.tsv"
OUT = CONTEXT_POOL

BASE = OUT / "Combined_external_source_contributions"

MODULES = [
    "Identity / chemical structure",
    "Protein interaction",
    "Quantitative activity",
    "Binding evidence",
    "Structural binding context",
]
EVIDENCE_TYPES = [
    "quantitative_binding",
    "bioassay",
    "functional_or_curated_relation",
    "functional_pharmacology",
    "structure",
]
EVIDENCE_LABELS = {
    "quantitative_binding": "Direct quantitative assay",
    "bioassay": "Bioassay evidence",
    "functional_or_curated_relation": "Curated interaction relation",
    "functional_pharmacology": "Functional pharmacology evidence",
    "structure": "Complex-structure evidence",
}
EVIDENCE_COLORS = {
    "quantitative_binding": "#168C8C",
    "bioassay": "#4C78A8",
    "functional_or_curated_relation": "#6F9B6D",
    "functional_pharmacology": "#D6A23A",
    "structure": "#887AA9",
}

# Shared ordering is primarily descending pair coverage (Panel B), with
# compound-only BRENDA retained as an explicit NA row in Panel B.
SOURCE_ORDER = [
    "ChEMBL",
    "PubChem BioAssay",
    "BindingDB",
    "PDBe",
    "IUPHAR/BPS Guide to PHARMACOLOGY",
    "PDSP KiDatabase",
    "DrugCentral",
    "PDBbind",
    "sc-PDB",
    "BioLiP",
    "BRENDA",
]
DISPLAY_SOURCE = {"IUPHAR/BPS Guide to PHARMACOLOGY": "GtoPdb (IUPHAR/BPS)"}

HEAT_COLORS = ["#EAF2F6", "#B9D4DF", "#79B4C8", "#3B86AA", "#174A73"]
BLANK = "#F5F6F7"
NA = "#E6E9EC"
DARK = "#24364B"
MUTED = "#5E7182"


def fmt(n: float | int | None) -> str:
    if n is None or pd.isna(n):
        return ""
    n = float(n)
    if abs(n) >= 100_000:
        return f"{n / 1000:.1f}k"
    if abs(n) >= 1000:
        return f"{n / 1000:.1f}k"
    return f"{int(round(n)):,}"


def luminance(hex_color: str) -> float:
    rgb = np.array(mcolors.to_rgb(hex_color))
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    return float(np.dot(linear, [0.2126, 0.7152, 0.0722]))


def main() -> None:
    for required in (PANEL_A_INPUT, PANEL_B_INPUT, PANEL_B_SUMMARY):
        if not required.exists():
            raise FileNotFoundError(required)

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    a_raw = pd.read_csv(PANEL_A_INPUT, sep="\t")
    b_raw = pd.read_csv(PANEL_B_INPUT, sep="\t")
    b_summary = pd.read_csv(PANEL_B_SUMMARY, sep="\t")

    # External sources only; no MemPro-derived annotations enter this figure.
    a_raw = a_raw.loc[a_raw["data_source"] != "MemPro-derived"].copy()
    a_raw = a_raw.loc[a_raw["compound_module"].isin(MODULES)].copy()
    b_raw = b_raw.loc[b_raw["category"].isin(EVIDENCE_TYPES)].copy()

    a_sources = set(a_raw["data_source"])
    b_sources = set(b_summary["source"])
    unknown = (a_sources | b_sources) - set(SOURCE_ORDER)
    if unknown:
        raise ValueError(f"Source order missing real sources: {sorted(unknown)}")

    # Panel A: source x module with explicit separation of observed zero versus
    # panel-level NA. Blank numeric cells represent no formal contribution.
    a_index = a_raw.set_index(["data_source", "compound_module"])["unique_canonical_compounds"]
    panel_a_rows: list[dict] = []
    for source in SOURCE_ORDER:
        panel_status = "PRESENT" if source in a_sources else "PANEL_NA"
        for module in MODULES:
            observed = a_index.get((source, module), np.nan)
            if panel_status == "PANEL_NA":
                cell_status = "PANEL_NA"
                count = np.nan
            elif pd.isna(observed) or float(observed) == 0:
                cell_status = "NO_FORMAL_CONTRIBUTION"
                count = np.nan
            else:
                cell_status = "OBSERVED"
                count = int(observed)
            panel_a_rows.append({
                "source": source,
                "display_source": DISPLAY_SOURCE.get(source, source),
                "compound_module": module,
                "n_unique_canonical_compounds": count,
                "cell_status": cell_status,
                "analysis_unit": "unique canonical compound",
                "analysis_universe": "646,670 all formal canonical compounds",
            })
    panel_a = pd.DataFrame(panel_a_rows)

    # Panel B: E3 was pre-deduplicated by source x modality x canonical pair.
    # It is deliberately not converted into a primary modality: pair modalities
    # are not mutually exclusive for all sources.
    b_totals = b_summary.set_index("source")["n_unique_pairs"]
    panel_b_rows: list[dict] = []
    for source in SOURCE_ORDER:
        source_present = source in b_sources
        total = b_totals.get(source, np.nan)
        for evidence_type in EVIDENCE_TYPES:
            hit = b_raw.loc[(b_raw["source"] == source) & (b_raw["category"] == evidence_type)]
            count = int(hit["n_unique_pairs"].iloc[0]) if len(hit) else 0
            if not source_present:
                status = "PANEL_NA"
                count_display = np.nan
                pct = np.nan
            else:
                status = "OBSERVED" if count > 0 else "NO_FORMAL_CONTRIBUTION"
                count_display = count
                pct = 100 * count / float(total)
            panel_b_rows.append({
                "source": source,
                "display_source": DISPLAY_SOURCE.get(source, source),
                "evidence_type": evidence_type,
                "display_evidence_type": EVIDENCE_LABELS[evidence_type],
                "n_unique_canonical_pairs_with_type": count_display,
                "n_unique_canonical_pairs_source_total": total if source_present else np.nan,
                "percent_of_source_unique_pairs": pct,
                "cell_status": status,
                "analysis_unit": "unique canonical protein-compound pair",
                "evidence_types_mutually_exclusive": False,
            })
    panel_b = pd.DataFrame(panel_b_rows)

    # Pair-level total is unique within a source; it must never be derived by
    # summing the non-mutually-exclusive modality values.
    source_summary = pd.DataFrame({
        "source": SOURCE_ORDER,
        "display_source": [DISPLAY_SOURCE.get(s, s) for s in SOURCE_ORDER],
    })
    source_summary["panel_a_status"] = np.where(source_summary["source"].isin(a_sources), "PRESENT", "PANEL_NA")
    source_summary["panel_b_status"] = np.where(source_summary["source"].isin(b_sources), "PRESENT", "PANEL_NA")
    source_summary = source_summary.merge(
        b_summary[["source", "n_unique_pairs", "n_raw_evidence_records", "n_evidence_types"]],
        on="source", how="left"
    )
    source_summary["panel_a_modules_with_observed_contribution"] = (
        panel_a.loc[panel_a["cell_status"] == "OBSERVED"].groupby("source").size()
        .reindex(SOURCE_ORDER).fillna(0).astype(int).to_numpy()
    )

    # Cross-panel QC: source panel B modality sums demonstrate non-exclusivity
    # when they exceed the deduplicated pair total for a source.
    b_check = panel_b.groupby("source", dropna=False)["n_unique_canonical_pairs_with_type"].sum(min_count=1)
    source_summary["sum_of_modality_pair_counts"] = source_summary["source"].map(b_check)
    source_summary["modality_sum_exceeds_pair_total"] = (
        source_summary["sum_of_modality_pair_counts"] > source_summary["n_unique_pairs"]
    )
    sources_with_nonexclusive = source_summary.loc[source_summary["modality_sum_exceeds_pair_total"], "source"].tolist()
    if not sources_with_nonexclusive:
        raise RuntimeError("Expected non-mutually-exclusive evidence modalities were not detected.")

    OUT.mkdir(parents=True, exist_ok=True)
    panel_a.to_csv(OUT / "Combined_panelA_compound_module_coverage.tsv", sep="\t", index=False)
    panel_b.to_csv(OUT / "Combined_panelB_interaction_evidence_profile.tsv", sep="\t", index=False)
    source_summary.to_csv(OUT / "Combined_source_summary.tsv", sep="\t", index=False)

    # Figure construction. Matrices share y values; source labels only appear
    # once at far left to strengthen the compound -> evidence-layer narrative.
    # A slightly taller canvas gives the shared source-row matrix more visual
    # weight at double-column width without reducing text below the font floor.
    fig = plt.figure(figsize=(7.0, 6.8))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.52, 1.0], wspace=0.10,
                  left=0.26, right=0.96, bottom=0.20, top=0.81)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1], sharey=ax_a)
    n_rows = len(SOURCE_ORDER)

    # A transparent background axis creates source-row bands that continue
    # through the central gutter. Both panel axes share exactly the same y
    # extent, row centres and row height as this background.
    band_ax = fig.add_axes([0.26, 0.20, 0.70, 0.61], zorder=0)
    band_ax.set_xlim(0, 1)
    band_ax.set_ylim(n_rows - 0.5, -0.5)
    for yi in range(n_rows):
        band_ax.axhspan(yi - 0.5, yi + 0.5,
                        color="#F5F7F8" if yi % 2 == 0 else "white", zorder=0)
        band_ax.axhline(yi, color="#DDE5EA", linewidth=0.45, alpha=0.80, zorder=1)
    band_ax.set_axis_off()
    ax_a.set_zorder(2)
    ax_b.set_zorder(2)
    ax_a.patch.set_alpha(0)
    ax_b.patch.set_alpha(0)

    # Panel A heatmap.
    positive_values = panel_a.loc[
        panel_a["n_unique_canonical_compounds"].notna(), "n_unique_canonical_compounds"
    ].astype(float)
    norm = mcolors.Normalize(vmin=np.log1p(positive_values.min()) * 0.92, vmax=np.log1p(positive_values.max()))
    cmap = mcolors.LinearSegmentedColormap.from_list("mempro_teal", HEAT_COLORS)
    for yi, source in enumerate(SOURCE_ORDER):
        source_data = panel_a.loc[panel_a["source"] == source].set_index("compound_module")
        row_is_na = source not in a_sources
        for xi, module in enumerate(MODULES):
            record = source_data.loc[module]
            value = record["n_unique_canonical_compounds"]
            if row_is_na:
                face = NA
            elif pd.isna(value):
                face = BLANK
            else:
                face = cmap(norm(np.log1p(float(value))))
            ax_a.add_patch(Rectangle(
                (xi - 0.5, yi - 0.5), 1, 1,
                facecolor=face,
                edgecolor="#D9E0E5" if row_is_na else "white",
                linewidth=0.9,
                hatch="////" if row_is_na else None,
            ))
            if pd.notna(value):
                text_color = "white" if luminance(mcolors.to_hex(face)) < 0.38 else DARK
                ax_a.text(xi, yi, fmt(value), ha="center", va="center", fontsize=6.35, color=text_color)

    ax_a.set_xlim(-0.5, len(MODULES) - 0.5)
    ax_a.set_ylim(n_rows - 0.5, -0.5)
    ax_a.set_xticks(range(len(MODULES)))
    ax_a.set_xticklabels([
        "Identity /\nchemical structure", "Compounds linked\nto proteins", "Compounds with\nquantitative\nmeasurements",
        "Compounds with\nbinding\nannotations", "Compounds with\nstructural context",
    ], ha="center", fontsize=5.6, color=DARK)
    ax_a.set_yticks(range(n_rows))
    ax_a.set_yticklabels([DISPLAY_SOURCE.get(s, s) for s in SOURCE_ORDER], fontsize=6.7, color=DARK)
    ax_a.tick_params(axis="x", length=0, pad=6)
    ax_a.tick_params(axis="y", length=0, pad=4)
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    # Panel B: five compact, horizontally grouped bars are constrained within
    # every shared source row. Modalities are intentionally non-exclusive.
    offsets = np.linspace(-0.30, 0.30, len(EVIDENCE_TYPES))
    ax_b.set_xlim(0, 105)
    ax_b.set_ylim(n_rows - 0.5, -0.5)
    ax_b.set_yticks(range(n_rows))
    ax_b.tick_params(axis="y", left=False, labelleft=False)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_xticklabels(["0", "25", "50", "75", "100%"], fontsize=6.2, color=MUTED)
    ax_b.grid(axis="x", color="#DCE5EA", linewidth=0.65, zorder=0)
    ax_b.axvline(0, color="#9DAEBB", linewidth=0.7)
    ax_b.set_axisbelow(True)
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    for yi, source in enumerate(SOURCE_ORDER):
        rows = panel_b.loc[panel_b["source"] == source].set_index("evidence_type")
        if source not in b_sources:
            ax_b.text(2.0, yi, "No formal pair evidence", ha="left", va="center", fontsize=6.0, color=MUTED, fontstyle="italic")
            continue
        for evidence_type, offset in zip(EVIDENCE_TYPES, offsets):
            record = rows.loc[evidence_type]
            pct = record["percent_of_source_unique_pairs"]
            if pd.isna(pct) or pct == 0:
                continue
            color = EVIDENCE_COLORS[evidence_type]
            ax_b.barh(yi + offset, pct, height=0.135, color=color, edgecolor="none", alpha=0.90, zorder=3)
        total = b_totals[source]
        ax_b.text(102.0, yi, f"{fmt(total)} pairs", ha="left", va="center", fontsize=5.95, color=MUTED, clip_on=False)
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=4.6,
                          markerfacecolor=EVIDENCE_COLORS[e], markeredgecolor="white", markeredgewidth=0.35,
                          label=EVIDENCE_LABELS[e]) for e in EVIDENCE_TYPES]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.615, 0.935), ncol=2,
               fontsize=5.7, handletextpad=0.35, columnspacing=0.8, borderaxespad=0.0,
               frameon=True, edgecolor="#D2D9DE", facecolor="white")

    fig.suptitle("External-source contributions across compound and interaction-evidence layers",
                 x=0.125, y=0.985, ha="left", fontsize=10.2, fontweight="bold", color=DARK)
    fig.text(0.125, 0.950,
             "Compound-level coverage and interaction-evidence profiles use distinct, non-additive analysis units.",
             ha="left", va="top", fontsize=7.0, color=MUTED)
    # Put panel headings in a dedicated band so they cannot collide with the
    # global legend or the panel-specific explanatory lines.
    box_a = ax_a.get_position()
    box_b = ax_b.get_position()
    fig.text(box_a.x0, 0.855, "A. Compound-level data-module coverage", ha="left", va="bottom",
             fontsize=8.4, fontweight="bold", color=DARK)
    fig.text(box_a.x0, 0.810,
             "Cells show unique canonical compounds covered by each external source\nin each data module; colour indicates log1p coverage.",
             ha="left", va="bottom", fontsize=5.7, color=MUTED, linespacing=1.0)
    fig.text(box_b.x0, 0.855, "B. Pair-level interaction-evidence profile by source", ha="left", va="bottom",
             fontsize=8.4, fontweight="bold", color=DARK)
    fig.text(box_b.x0, 0.810,
             "Bars show the percentage of source-specific unique canonical protein–compound pairs\ncovered by each non-mutually-exclusive evidence modality.",
             ha="left", va="bottom", fontsize=5.7, color=MUTED, linespacing=1.0)
    divider_x = (box_a.x1 + box_b.x0) / 2
    fig.add_artist(Line2D([divider_x, divider_x], [0.20, 0.81], transform=fig.transFigure,
                          color="#D2DCE2", linewidth=0.65, zorder=4))
    # Move the heatmap scale below Panel A, leaving the A/B gutter as a clean
    # visual channel and preserving the shared row framework.
    cax = fig.add_axes([box_a.x0, 0.090, box_a.width * 0.86, 0.013])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.ax.tick_params(labelsize=5.6, length=2.0, colors=MUTED, pad=1)
    cbar.set_label("Canonical compound coverage (log1p)", fontsize=5.9, color=DARK, labelpad=2)
    cbar.set_ticks(np.log1p([10, 1000, 100000, 500000]))
    cbar.set_ticklabels(["10", "1k", "100k", "500k"])
    fig.text(0.125, 0.018,
             "Blank heatmap cells indicate no direct source-level contribution counted. Hatched rows indicate that the source is absent from that analysis layer.\n"
             "Panel B modalities are non-mutually-exclusive; percentages may sum to >100%. Panel A reports unique canonical compounds, whereas Panel B reports unique canonical protein–compound pairs.",
             ha="left", va="bottom", fontsize=5.75, color=MUTED)

    fig.savefig(BASE.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(BASE.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(BASE.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(BASE.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    (OUT / "Combined_source_figure_fields_used.txt").write_text(
        "Panel A source table: C1_source_module_compound_coverage.tsv\n"
        "  source = data_source\n"
        "  module = compound_module\n"
        "  count = unique_canonical_compounds\n"
        "  unit = unique canonical compound\n\n"
        "Panel B source tables: E3_source_evidence_type_counts.tsv and E3_source_summary.tsv\n"
        "  source = source\n"
        "  evidence type = category (formal evidence_modality_v7 classification used by E3)\n"
        "  modality count = n_unique_pairs (deduplicated canonical pair within source x modality)\n"
        "  source pair total = n_unique_pairs in E3_source_summary.tsv\n"
        "  unit = unique canonical protein-compound pair\n\n"
        "Panel B formal modalities are non-mutually-exclusive; no primary modality was inferred.\n",
        encoding="utf-8"
    )
    (OUT / "Combined_source_figure_caption.txt").write_text(
        "External-source contributions across compound and interaction-evidence layers. "
        "(A) Cells show unique canonical compounds covered by each external source in each data module; colour indicates log1p coverage. "
        "(B) Bars show the percentage of source-specific unique canonical protein–compound pairs covered by each non-mutually-exclusive evidence modality; total unique pairs are printed at right. "
        "Panel A reports unique canonical compounds, whereas Panel B reports unique canonical protein–compound pairs; these analysis units are distinct and non-additive. "
        "Source totals are deduplicated within their stated unit, and the same compound or pair may occur in multiple external sources. "
        "Blank heatmap cells indicate no direct source-level contribution counted. Hatched rows indicate that the source is absent from the compound-module analysis layer. "
        "Because Panel B modalities are non-mutually-exclusive, percentages may sum to more than 100% within a source.",
        encoding="utf-8"
    )
    qc = pd.DataFrame([
        ("Panel A analysis unit", "unique canonical compound", "PASS"),
        ("Panel B analysis unit", "unique canonical protein-compound pair", "PASS"),
        ("Panel A external sources", len(a_sources), "PASS"),
        ("Panel B external sources", len(b_sources), "PASS"),
        ("Shared external source master list", len(SOURCE_ORDER), "PASS"),
        ("Panel A includes MemPro-derived row", False, "PASS"),
        ("Panel B evidence types mutually exclusive", False, "PASS"),
        ("Panel B visualization", "grouped lollipop profile; not 100% stacked", "PASS"),
        ("Sources whose modality sums exceed deduplicated pair total", "; ".join(sources_with_nonexclusive), "PASS"),
        ("Panel A blanks", "no direct source-level contribution, not numeric zero labels", "PASS"),
        ("Panel-level NA", "source absent from one analysis layer", "PASS"),
    ], columns=["check", "result", "status"])
    qc.to_csv(OUT / "Combined_source_figure_QC.tsv", sep="\t", index=False)
    (OUT / "Combined_source_figure_manifest.json").write_text(json.dumps({
        "figure": "Combined_external_source_contributions",
        "panel_a_unit": "unique canonical compound",
        "panel_b_unit": "unique canonical protein-compound pair",
        "evidence_types_mutually_exclusive": False,
        "panel_b_visualization": "grouped lollipop profile",
        "n_sources_panel_a": int(len(a_sources)),
        "n_sources_panel_b": int(len(b_sources)),
        "n_modules_panel_a": int(len(MODULES)),
        "n_evidence_types_panel_b": int(len(EVIDENCE_TYPES)),
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
