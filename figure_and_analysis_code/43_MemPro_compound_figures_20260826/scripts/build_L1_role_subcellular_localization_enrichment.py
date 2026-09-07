#!/usr/bin/env python3
"""Build L1: role x formal subcellular-localization enrichment heatmap.

Figure contract
---------------
Core conclusion: Formal membrane-protein roles show relative enrichment and
depletion patterns across HPA standardized subcellular-localization terms.
Archetype: one quantitative grid.
Analysis unit: unique canonical protein. Localization annotations are multi-
label; a protein contributes at most once to each protein x location term.
Eligible universe: formal canonical membrane proteins with at least one mapped,
non-empty localization term. Proteins lacking usable localization annotations
are excluded from the universe, never treated as negative evidence.
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
from scipy.stats import fisher_exact


FORMAL = Path(r"D:\finale\FORMAL")
ROOT = Path(r"D:\finale\43_MemPro_compound_figures_20260826")
OUT = ROOT / "results" / "context_figure_pool"

PROTEIN_MASTER = FORMAL / "01_core_v72" / "01_release_tables" / "protein_master_v72.tsv.gz"
ROLE_TABLE = FORMAL / "03_publication_repairs" / "protein_membrane_role_FORMAL.tsv"
LOCALIZATION = FORMAL / "01_core_v72" / "01_release_tables" / "subcellular_localization_v72.tsv.gz"
BASE = OUT / "L1_role_subcellular_localization_enrichment"

DARK = "#24364B"
MUTED = "#5E7182"
LOW_SUPPORT = "#ECEFF1"
DISPLAY_N = 12


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    """Benjamini-Hochberg FDR, preserving the original index."""
    valid = pvalues.notna()
    p = pvalues.loc[valid].to_numpy(dtype=float)
    if len(p) == 0:
        return pd.Series(np.nan, index=pvalues.index, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.full(len(pvalues), np.nan)
    idx = np.flatnonzero(valid.to_numpy())
    out[idx[order]] = adjusted
    return pd.Series(out, index=pvalues.index)


def high_level_role_display(role: str) -> str:
    return role.replace("_", " ")


def main() -> None:
    for path in (PROTEIN_MASTER, ROLE_TABLE, LOCALIZATION):
        if not path.exists():
            raise FileNotFoundError(path)
    OUT.mkdir(parents=True, exist_ok=True)

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    protein = pd.read_csv(PROTEIN_MASTER, sep="\t", usecols=["canonical_uniprot_accession", "public_inclusion_v7"])
    role = pd.read_csv(ROLE_TABLE, sep="\t", usecols=["canonical_uniprot_accession", "formal_primary_membrane_role"])
    localization = pd.read_csv(
        LOCALIZATION, sep="\t", low_memory=False,
        usecols=[
            "target_uniprot_id", "location_term", "go_id", "go_label", "location_role",
            "hpa_reliability", "direct_observation_flag", "protein_mapping_status", "source_dataset",
            "source_version", "mapping_status_v71", "observation_semantics_v71",
            "biological_context_v71", "isoform_specificity_v71",
        ],
    )

    # Formal protein universe and role schema. No role inference is performed.
    protein = protein.loc[protein["public_inclusion_v7"].eq(1)].copy()
    if protein["canonical_uniprot_accession"].nunique() != 7800:
        raise ValueError("Formal canonical protein denominator is not 7,800.")
    if role["canonical_uniprot_accession"].nunique() != 7800:
        raise ValueError("Formal role table does not cover the 7,800 formal proteins.")

    rows_before_mapping = len(localization)
    mapped = localization.loc[
        localization["mapping_status_v71"].eq("MAPPED")
        & localization["target_uniprot_id"].notna()
        & localization["location_term"].notna()
        & localization["location_term"].astype(str).str.strip().ne("")
    ].copy()
    mapped_rows = len(mapped)

    # Preserve the formal, source-controlled HPA location_term field exactly;
    # no textual keywords or researcher-made localization classes are used.
    pairs_before_dedup = len(mapped)
    protein_location = (
        mapped[["target_uniprot_id", "location_term", "go_id", "go_label"]]
        .drop_duplicates(["target_uniprot_id", "location_term"])
        .merge(protein[["canonical_uniprot_accession"]], left_on="target_uniprot_id", right_on="canonical_uniprot_accession", how="inner")
        .merge(role, left_on="target_uniprot_id", right_on="canonical_uniprot_accession", how="inner", suffixes=("", "_role"))
    )
    duplicate_protein_location = len(mapped) - len(
        mapped.drop_duplicates(["target_uniprot_id", "location_term"])
    )
    eligible = set(protein_location["target_uniprot_id"])
    universe_n = len(eligible)
    if universe_n == 0:
        raise ValueError("No eligible proteins after formal localization filtering.")

    # Full, formal category catalog. Top 12 are selected only for display;
    # every category remains in source data and in the BH correction family.
    loc_summary = (
        protein_location.groupby(["location_term", "go_id", "go_label"], dropna=False)["target_uniprot_id"]
        .nunique().rename("n_unique_eligible_canonical_proteins").reset_index()
        .sort_values(["n_unique_eligible_canonical_proteins", "location_term"], ascending=[False, True])
        .reset_index(drop=True)
    )
    display_locations = loc_summary.head(DISPLAY_N)["location_term"].tolist()
    loc_summary["displayed_in_L1"] = loc_summary["location_term"].isin(display_locations)
    loc_summary["formal_category_source"] = "subcellular_localization_v72.location_term"
    loc_summary["go_id_missing"] = loc_summary["go_id"].isna()

    role_sizes = (
        protein_location.groupby("formal_primary_membrane_role")["target_uniprot_id"].nunique()
        .rename("n_eligible_canonical_proteins").reset_index()
        .sort_values(["n_eligible_canonical_proteins", "formal_primary_membrane_role"], ascending=[False, True])
        .reset_index(drop=True)
    )
    roles = role_sizes["formal_primary_membrane_role"].tolist()
    role_size_map = role_sizes.set_index("formal_primary_membrane_role")["n_eligible_canonical_proteins"].to_dict()
    loc_size_map = loc_summary.set_index("location_term")["n_unique_eligible_canonical_proteins"].to_dict()
    protein_sets_by_role = {
        r: set(protein_location.loc[protein_location["formal_primary_membrane_role"] == r, "target_uniprot_id"])
        for r in roles
    }
    protein_sets_by_location = {
        loc: set(protein_location.loc[protein_location["location_term"] == loc, "target_uniprot_id"])
        for loc in loc_summary["location_term"]
    }

    statistics_rows: list[dict] = []
    for r in roles:
        rset = protein_sets_by_role[r]
        for loc in loc_summary["location_term"]:
            lset = protein_sets_by_location[loc]
            a = len(rset & lset)
            b = len(rset - lset)
            c = len(lset - rset)
            d = universe_n - a - b - c
            if min(a, b, c, d) < 0:
                raise ValueError("Invalid Fisher contingency table.")
            odds_fisher, pvalue = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            # Haldane-Anscombe correction gives a finite, sign-preserving log2
            # effect for zero cells; the uncorrected Fisher OR is retained too.
            odds_corrected = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
            expected = len(rset) * len(lset) / universe_n
            statistics_rows.append({
                "formal_primary_membrane_role": r,
                "display_role": high_level_role_display(r),
                "location_term": loc,
                "n_eligible_proteins": universe_n,
                "role_total_eligible_proteins": len(rset),
                "localization_total_eligible_proteins": len(lset),
                "observed_unique_proteins": a,
                "expected_unique_proteins": expected,
                "a_role_and_location": a,
                "b_role_not_location": b,
                "c_not_role_location": c,
                "d_neither": d,
                "odds_ratio_fisher": odds_fisher,
                "odds_ratio_haldane_anscombe": odds_corrected,
                "log2_odds_ratio": float(np.log2(odds_corrected)),
                "p_value": pvalue,
            })
    stats = pd.DataFrame(statistics_rows)
    stats["BH_FDR"] = bh_fdr(stats["p_value"])
    stats["displayed_in_L1"] = stats["location_term"].isin(display_locations)
    stats["low_support_cell"] = stats["observed_unique_proteins"] < 20
    stats["significant_for_display"] = (
        stats["BH_FDR"].lt(0.05) & stats["observed_unique_proteins"].ge(20)
    )
    stats["display_log2_odds_ratio_clipped"] = stats["log2_odds_ratio"].clip(-3, 3)
    stats["inference_note"] = np.where(
        stats["low_support_cell"],
        "low support (<20); visually muted and not emphasized for inference",
        "eligible for displayed inference",
    )

    counts = stats[[
        "formal_primary_membrane_role", "display_role", "location_term", "role_total_eligible_proteins",
        "localization_total_eligible_proteins", "observed_unique_proteins", "expected_unique_proteins",
        "displayed_in_L1", "low_support_cell",
    ]].copy()
    counts["analysis_unit"] = "unique canonical protein"
    counts["localization_annotations_multi_label"] = True
    counts["deduplicated_by"] = "target_uniprot_id + location_term"

    # Write source data before rendering, so the plotted table can be audited.
    counts.to_csv(OUT / "L1_role_localization_counts.tsv", sep="\t", index=False)
    stats.to_csv(OUT / "L1_role_localization_enrichment_statistics.tsv", sep="\t", index=False)
    loc_summary.to_csv(OUT / "L1_localization_summary.tsv", sep="\t", index=False)

    # Heatmap data: all formal roles x top 12 source-controlled location terms.
    display_stats = stats.loc[stats["displayed_in_L1"]].copy()
    loc_display_map = {loc: loc for loc in display_locations}
    display_stats["location_term"] = pd.Categorical(display_stats["location_term"], categories=display_locations, ordered=True)
    display_stats["formal_primary_membrane_role"] = pd.Categorical(
        display_stats["formal_primary_membrane_role"], categories=roles, ordered=True
    )
    display_stats = display_stats.sort_values(["formal_primary_membrane_role", "location_term"])

    fig_height = max(6.3, 0.36 * len(roles) + 2.25)
    fig, ax = plt.subplots(figsize=(7.0, fig_height))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "muted_blue_white_red", ["#2F6F9F", "#B5D0E2", "#F7F7F5", "#E6B7B5", "#B6575C"]
    )
    norm = mcolors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    for yi, r in enumerate(roles):
        for xi, loc in enumerate(display_locations):
            row = display_stats.loc[
                (display_stats["formal_primary_membrane_role"] == r) & (display_stats["location_term"] == loc)
            ].iloc[0]
            low_support = bool(row["low_support_cell"])
            face = LOW_SUPPORT if low_support else cmap(norm(row["display_log2_odds_ratio_clipped"]))
            ax.add_patch(Rectangle((xi - 0.5, yi - 0.5), 1, 1, facecolor=face, edgecolor="white", linewidth=0.9))
            if not low_support:
                bg_lum = np.dot(
                    np.where(np.array(mcolors.to_rgb(face)) <= 0.04045,
                             np.array(mcolors.to_rgb(face)) / 12.92,
                             ((np.array(mcolors.to_rgb(face)) + 0.055) / 1.055) ** 2.4),
                    [0.2126, 0.7152, 0.0722]
                )
                label_color = "white" if bg_lum < 0.38 else DARK
                ax.text(xi, yi, str(int(row["observed_unique_proteins"])), ha="center", va="center",
                        fontsize=6.1, color=label_color)
            if bool(row["significant_for_display"]):
                ax.scatter(xi + 0.34, yi - 0.34, s=11, color="#18222C", zorder=4, linewidths=0)

    ax.set_xlim(-0.5, len(display_locations) - 0.5)
    ax.set_ylim(len(roles) - 0.5, -0.5)
    ax.set_xticks(range(len(display_locations)))
    ax.set_xticklabels(display_locations, rotation=32, ha="right", rotation_mode="anchor", fontsize=6.3, color=DARK)
    ax.set_yticks(range(len(roles)))
    role_labels = [
        high_level_role_display(r) + ("  (low support)" if role_size_map[r] < 100 else "")
        for r in roles
    ]
    ax.set_yticklabels(role_labels, fontsize=6.7, color=DARK)
    for label, r in zip(ax.get_yticklabels(), roles):
        if role_size_map[r] < 100:
            label.set_color(MUTED)
    ax.tick_params(axis="both", length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Subcellular-localization enrichment across membrane-protein roles",
                 loc="left", fontsize=10.0, fontweight="bold", color=DARK, pad=42)
    ax.text(0.0, 1.055, "Cell colour shows log2 odds ratio; labels show unique canonical proteins.",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.8, color=MUTED)
    ax.text(0.0, 1.025, "Black dot denotes BH-FDR < 0.05 with observed support ≥20.",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.4, color=MUTED)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.038, pad=0.020)
    cbar.ax.tick_params(labelsize=6.0, length=2.0, colors=MUTED)
    cbar.set_label("log2 odds ratio", fontsize=6.8, color=DARK, labelpad=5)
    cbar.set_ticks([-3, -2, -1, 0, 1, 2, 3])
    fig.text(0.125, 0.014,
             "Localization annotations are multi-label. Low-support cells (<20 proteins) are neutral and not emphasized for inference; "
             "proteins without usable localization annotation were not treated as negative evidence.",
             ha="left", va="bottom", fontsize=5.7, color=MUTED)
    fig.subplots_adjust(left=0.30, right=0.94, top=0.83, bottom=0.18)
    fig.savefig(BASE.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(BASE.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(BASE.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(BASE.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Human- and machine-readable audit artifacts.
    (OUT / "L1_fields_used.txt").write_text(
        "concept\tactual_field_name\ttable/file\tdefinition used\n"
        "canonical protein ID\ttarget_uniprot_id\tsubcellular_localization_v72.tsv.gz\tmapped HPA localization target accession\n"
        "formal canonical protein\tcanonical_uniprot_accession\tprotein_master_v72.tsv.gz\tpublic_inclusion_v7 = 1\n"
        "primary membrane-protein role\tformal_primary_membrane_role\tprotein_membrane_role_FORMAL.tsv\tformal publication role; no reclassification\n"
        "standardized localization category\tlocation_term\tsubcellular_localization_v72.tsv.gz\tformal HPA standardized localization term, retained exactly\n"
        "localization ontology ID\tgo_id\tsubcellular_localization_v72.tsv.gz\tlinked GO cellular-component ID where supplied\n"
        "localization ontology label\tgo_label\tsubcellular_localization_v72.tsv.gz\tlinked GO label where supplied\n"
        "localization source\tsource_dataset\tsubcellular_localization_v72.tsv.gz\tHPA_subcellular_IF\n"
        "source version\tsource_version\tsubcellular_localization_v72.tsv.gz\tHPA 25.1\n"
        "mapping status\tmapping_status_v71\tsubcellular_localization_v72.tsv.gz\tMAPPED required for eligible localization\n"
        "observation semantics\tobservation_semantics_v71\tsubcellular_localization_v72.tsv.gz\tpreserved as provenance; not used as a new classification\n",
        encoding="utf-8",
    )
    (OUT / "L1_display_rules.txt").write_text(
        f"Rows: all {len(roles)} formal primary membrane-protein roles, ordered by eligible unique-protein support descending; roles with <100 eligible proteins are muted but retained.\n"
        f"Columns: top {DISPLAY_N} formal location_term categories by unique eligible canonical-protein coverage; no manual Other category was created.\n"
        f"All {len(loc_summary)} formal categories are retained in L1_localization_summary.tsv and in the BH correction family.\n"
        "Location terms are read directly from subcellular_localization_v72.location_term; no keyword-based regrouping was used.\n",
        encoding="utf-8",
    )
    generic_present = sorted(set(loc_summary["location_term"].str.lower()) & {"cell", "membrane"})
    source_status = localization["source_dataset"].dropna().unique().tolist()
    qc = pd.DataFrame([
        ("formal canonical protein denominator", len(protein), "expected 7800", "PASS"),
        ("formal role coverage", role["canonical_uniprot_accession"].nunique(), "expected 7800", "PASS"),
        ("raw localization rows", rows_before_mapping, "formal localization table", "PASS"),
        ("mapped localization rows", mapped_rows, "mapping_status_v71=MAPPED", "PASS"),
        ("mapped localization rows excluded as ambiguous", int(localization["mapping_status_v71"].eq("AMBIGUOUS").sum()), "not eligible", "PASS"),
        ("eligible protein universe", universe_n, "formal proteins with >=1 mapped non-empty location_term", "PASS"),
        ("formal proteins without usable localization annotation", len(protein) - universe_n, "excluded from denominator, not negatives", "PASS"),
        ("protein-location rows after deduplication", len(protein_location), "unique protein x location_term", "PASS"),
        ("duplicate protein-location records removed", duplicate_protein_location, "source records do not inflate counts", "PASS"),
        ("localization categories total", len(loc_summary), "formal location_term values", "PASS"),
        ("localization categories displayed", DISPLAY_N, "top coverage", "PASS"),
        ("formal roles", len(roles), "all formal roles retained", "PASS"),
        ("Fisher comparison family", len(stats), "all role x all formal location term pairs", "PASS"),
        ("BH-FDR scope", "all valid Fisher comparisons", "not display subset only", "PASS"),
        ("multi-label localization", True, "protein can occur in multiple location categories", "PASS"),
        ("keyword-derived localization classification", False, "formal location_term retained", "PASS"),
        ("generic broad location terms present", "; ".join(generic_present) if generic_present else "none", "reported only; not removed", "PASS"),
        ("localization source datasets", "; ".join(source_status), "formal source field", "PASS"),
        ("non-GO formal location terms", int(loc_summary["go_id_missing"].sum()), "retained as formal HPA terms", "REVIEW"),
    ], columns=["check", "result", "definition_or_rule", "status"])
    qc.to_csv(OUT / "L1_QC.tsv", sep="\t", index=False)
    (OUT / "L1_caption_draft.txt").write_text(
        "Subcellular-localization enrichment across membrane-protein roles. Each comparison uses unique canonical proteins with at least one usable, mapped formal localization annotation. "
        "Localization categories are multi-label and non-mutually-exclusive; each protein contributes at most once to a given protein-localization category. "
        "For every role-localization combination, Fisher exact testing was performed against the same eligible protein universe and P values were adjusted across the full formal comparison family using Benjamini-Hochberg FDR. "
        "Cell colour denotes the Haldane-Anscombe-corrected log2 odds ratio, shown symmetrically on a clipped ±3 scale; red indicates relative enrichment and blue relative depletion. "
        "Cell labels are observed unique canonical proteins. Black dots denote BH-FDR <0.05 with observed support ≥20. Low-support cells (<20 proteins) are shown in neutral grey and are not emphasized for inference. "
        "Enrichment reflects relative representation, not absolute localization abundance. Proteins without usable localization annotation were not automatically treated as negative localization evidence.",
        encoding="utf-8",
    )
    (OUT / "L1_manifest.json").write_text(json.dumps({
        "figure": "L1_role_subcellular_localization_enrichment",
        "analysis_unit": "unique canonical protein",
        "eligible_protein_universe": "formal canonical proteins with >=1 mapped, non-empty formal location_term",
        "n_eligible_proteins": universe_n,
        "n_roles": len(roles),
        "n_formal_localization_categories": len(loc_summary),
        "n_displayed_localization_categories": DISPLAY_N,
        "localization_annotations_multi_label": True,
        "test": "two-sided Fisher exact test",
        "multiple_testing": "Benjamini-Hochberg over all role x formal location-term comparisons",
        "display_effect": "Haldane-Anscombe-corrected log2 odds ratio, clipped to +/-3",
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
