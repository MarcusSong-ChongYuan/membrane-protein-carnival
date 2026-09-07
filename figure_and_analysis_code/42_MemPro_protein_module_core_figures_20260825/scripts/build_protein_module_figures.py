"""Build standalone MemPro protein-module figures from frozen FORMAL tables.

Figure contract
---------------
Core conclusion: MemPro's 7,800 formal membrane proteins have an explicit
membrane-connection/evidence composition, traceable source coverage, and a
non-random functional organization.
Archetype: three independent quantitative figures (not a dashboard).
Backend: Python/matplotlib only.
Input: FORMAL frozen tables plus the final role overlay; no biological facts
are changed by this script.
"""
from __future__ import annotations

import hashlib
import json
import math
import textwrap
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch
from scipy.stats import fisher_exact


ROOT = Path(r"D:\finale\42_MemPro_protein_module_core_figures_20260825")
FORMAL = Path(r"D:\finale\FORMAL")
OUT = ROOT / "figures"
SRC = ROOT / "source_data"
QA = ROOT / "qa"
for p in [OUT / "png", OUT / "pdf", OUT / "svg", OUT / "tiff", SRC, QA]:
    p.mkdir(parents=True, exist_ok=True)

MASTER = FORMAL / "01_core_v72" / "01_release_tables" / "protein_master_v72.tsv.gz"
CROSS = FORMAL / "01_core_v72" / "01_release_tables" / "protein_cross_classification_summary_v72.tsv.gz"
ROLE = FORMAL / "03_publication_repairs" / "protein_membrane_role_FORMAL.tsv"
EXPRESSION = FORMAL / "01_core_v72" / "01_release_tables" / "expression_measurement_v72.tsv.gz"
LOCALIZATION = FORMAL / "01_core_v72" / "01_release_tables" / "subcellular_localization_v72.tsv.gz"
SITES = FORMAL / "01_core_v72" / "01_release_tables" / "interaction_site_v72.tsv.gz"
MODELS = FORMAL / "02_companion_v721" / "02_companion_tables" / "protein_structure_model_v721.tsv.gz"
DISEASE = FORMAL / "01_core_v72" / "01_release_tables" / "protein_disease_relation_v72.tsv.gz"
COMPLEX = FORMAL / "01_core_v72" / "01_release_tables" / "protein_complex_v72.tsv.gz"

# User-selected restrained Seaborn palette.  Dark text remains neutral for
# legibility; colour encodes categories or signed effect direction only.
SEABORN = ["#7b95c6", "#49c2d9", "#a1d8e8", "#67a583", "#a2c986", "#d0e2c0", "#fded95", "#ffc1a6", "#f59c7c", "#f47254", "#c85e62"]
NAVY = "#2D405C"
TEAL = "#49c2d9"
SAGE = "#67a583"
GOLD = "#fded95"
PURPLE = "#7b95c6"
GREY = "#A9B1B8"
TEXT = "#24364B"
LIGHT = "#E9EEF2"
WHITE = "#FFFFFF"
TIER_COLORS = {"E1": "#7b95c6", "E2": "#49c2d9", "E3": "#a1d8e8"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.labelcolor": TEXT,
    "xtick.color": TEXT,
    "ytick.color": TEXT,
    "text.color": TEXT,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": WHITE,
    "axes.facecolor": WHITE,
})


def clean_string(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def humanize(s: str) -> str:
    return s.replace("_", " ").replace("/", " / ").strip().title()


def bh_fdr(pvalues: list[float]) -> np.ndarray:
    """Benjamini-Hochberg correction without a non-core dependency."""
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranks = np.arange(1, len(p) + 1)
    adjusted_sorted = p[order] * len(p) / ranks
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = np.clip(adjusted_sorted, 0, 1)
    return adjusted


def save_all(fig: mpl.figure.Figure, stem: str) -> None:
    # Explicit extensions are kept here so export intent is auditable as well
    # as reproducible. SVG/PDF retain editable text; PNG/TIFF are 600 dpi.
    fig.savefig(OUT / "svg" / f"{stem}.svg", bbox_inches="tight", facecolor=WHITE)
    fig.savefig(OUT / "pdf" / f"{stem}.pdf", bbox_inches="tight", facecolor=WHITE)
    fig.savefig(OUT / "png" / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor=WHITE)
    fig.savefig(OUT / "tiff" / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_proteins() -> pd.DataFrame:
    p = pd.read_csv(MASTER, sep="\t", compression="gzip", usecols=[
        "canonical_uniprot_accession", "membrane_class_v7", "membrane_evidence_level_v7",
        "default_web_inclusion_v72", "public_inclusion_v7",
    ])
    # The FORMAL protein table itself is the 7,800-protein public set. The
    # `default_web_inclusion_v72` flag separates the E1/E2 default layer from
    # retained E3 records and must not silently remove E3 from this census.
    p = p.loc[p["public_inclusion_v7"].eq(1)].copy()
    p["membrane_class_v7"] = clean_string(p["membrane_class_v7"])
    p["membrane_evidence_level_v7"] = clean_string(p["membrane_evidence_level_v7"])
    p = p.rename(columns={"canonical_uniprot_accession": "target_uniprot_id"})
    if p.target_uniprot_id.nunique() != 7800:
        raise ValueError(f"Expected 7,800 formal proteins, observed {p.target_uniprot_id.nunique():,}.")
    if not set(p.membrane_class_v7).issubset({"A", "B", "C"}):
        raise ValueError("Unexpected membrane class in formal proteins.")
    if not set(p.membrane_evidence_level_v7).issubset({"E1", "E2", "E3"}):
        raise ValueError("Unexpected evidence tier in formal proteins.")
    return p


def load_current_roles(proteins: pd.DataFrame) -> pd.DataFrame:
    r = pd.read_csv(ROLE, sep="\t")
    idcol = next(c for c in r.columns if c in {"target_uniprot_id", "canonical_uniprot_accession"})
    rolecol = next(c for c in r.columns if "membrane_role" in c and "primary" in c)
    r = r[[idcol, rolecol]].drop_duplicates(idcol).rename(columns={idcol: "target_uniprot_id", rolecol: "membrane_role"})
    out = proteins[["target_uniprot_id"]].merge(r, on="target_uniprot_id", how="left")
    out["membrane_role"] = clean_string(out["membrane_role"]).replace("", "family_defined_membrane_role_unresolved")
    return out


def figure_01_radial(proteins: pd.DataFrame) -> pd.DataFrame:
    classes = ["A", "B", "C"]
    tiers = ["E1", "E2", "E3"]
    total = len(proteins)
    table = (
        proteins.groupby(["membrane_class_v7", "membrane_evidence_level_v7"]).size()
        .reindex(pd.MultiIndex.from_product([classes, tiers], names=["membrane_class", "evidence_tier"]), fill_value=0)
        .rename("protein_count").reset_index()
    )
    table["class_total"] = table.groupby("membrane_class")["protein_count"].transform("sum")
    table["class_percentage"] = table["class_total"] / total * 100
    table["within_class_percentage"] = np.where(table["class_total"] > 0, table["protein_count"] / table["class_total"] * 100, np.nan)
    table.to_csv(SRC / "Figure1A_membrane_class_by_evidence_tier.tsv", sep="\t", index=False)

    # Shared angular scale: 0–6,000 canonical proteins. Each class is an annular bar;
    # tiers are contiguous arc segments, so ring length encodes the class total.
    fig = plt.figure(figsize=(7.01, 7.01), constrained_layout=False)
    ax = fig.add_axes([0.05, 0.06, 0.72, 0.86], projection="polar")
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_ylim(0, 3.65)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)
    max_count = 6000
    span = np.deg2rad(292)
    theta0 = np.deg2rad(-146)
    ring = {"A": (2.42, 0.47), "B": (1.77, 0.47), "C": (1.12, 0.47)}
    for cls in classes:
        radius, height = ring[cls]
        class_total = int(table.loc[table.membrane_class.eq(cls), "class_total"].iloc[0])
        arc = span * class_total / max_count
        # Light reference track makes the common count scale visible.
        ax.bar(theta0 + span / 2, height, width=span, bottom=radius, color="#F1F4F6", edgecolor="none", align="center", zorder=0)
        start = theta0
        for tier in tiers:
            n = int(table.loc[(table.membrane_class.eq(cls)) & (table.evidence_tier.eq(tier)), "protein_count"].iloc[0])
            w = span * n / max_count
            if n:
                ax.bar(start + w / 2, height, width=w, bottom=radius, color=TIER_COLORS[tier], edgecolor=WHITE, linewidth=0.9, align="center", zorder=3)
                if w > np.deg2rad(18):
                    angle = start + w / 2
                    ax.text(angle, radius + height / 2, f"{tier}\n{n:,}", fontsize=7.2, weight="bold", color=WHITE,
                            ha="center", va="center", rotation=0, rotation_mode="anchor", zorder=4)
            start += w

    # Minimal shared scale, positioned above the arcs rather than a complete polar axis.
    for n in [0, 2000, 4000, 6000]:
        ang = theta0 + span * n / max_count
        ax.plot([ang, ang], [2.92, 3.05], color="#8090A0", lw=0.7, zorder=2)
        ax.text(ang, 3.16, f"{n:,}", fontsize=7, ha="center", va="center", color="#5E7182")
    ax.text(theta0 + span / 2, 3.38, "Canonical proteins (shared arc scale)", fontsize=8, ha="center", va="center", weight="bold")
    ax.text(np.deg2rad(0), 0.30, "MemPro\nformal protein set", fontsize=10.5, weight="bold", ha="center", va="center", color=NAVY)
    ax.text(np.deg2rad(0), 0.02, f"n = {total:,}", fontsize=9, ha="center", va="center", color=TEXT)
    handles = [Patch(facecolor=TIER_COLORS[t], edgecolor="none", label=f"{t} membrane-evidence tier") for t in tiers]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.79, 0.67), fontsize=8.2, labelspacing=1.1, frameon=False)
    # External, vertically ordered class census: direct labels without forcing
    # text into narrow B/C rings. The swatches follow A/B/C ring order.
    class_description = {"A": "integral / transmembrane", "B": "membrane-embedded / lipid-anchored", "C": "peripheral membrane-associated"}
    class_y = {"A": .47, "B": .37, "C": .27}
    for cls in classes:
        n = int(table.loc[table.membrane_class.eq(cls), "class_total"].iloc[0])
        fig.text(.80, class_y[cls], f"{cls}  {class_description[cls]}\n    {n:,}  ({n / total:.1%})", fontsize=8.3, color=TEXT, weight="bold", va="top")
    fig.text(0.05, 0.965, "Membrane connection mechanism and evidence strength", fontsize=11, weight="bold", color=NAVY)
    fig.text(0.05, 0.935, "Each annulus is one non-overlapping A/B/C class; arc length is the canonical-protein count.", fontsize=7.5, color="#5E7182")
    save_all(fig, "Figure1A_membrane_class_evidence_radial")
    return table


def _contributor_sets(cross: pd.DataFrame) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {k: set() for k in ["HTP", "Membranome", "PDBTM", "OPM", "HPA membrane evidence"]}
    for row in cross.itertuples(index=False):
        uid = getattr(row, "target_uniprot_id")
        values = str(getattr(row, "contributing_membrane_databases_v63") or "").split(";")
        values = {v.strip() for v in values if v.strip()}
        if any(v == "UniTmp_HTP" for v in values): sets["HTP"].add(uid)
        if "Membranome" in values: sets["Membranome"].add(uid)
        if "PDBTM" in values: sets["PDBTM"].add(uid)
        if "OPM" in values: sets["OPM"].add(uid)
        if any(v.startswith("HPA_") for v in values): sets["HPA membrane evidence"].add(uid)
    return sets


def figure_01b_source_heatmap(proteins: pd.DataFrame) -> pd.DataFrame:
    formal_ids = set(proteins.target_uniprot_id)
    cross = pd.read_csv(CROSS, sep="\t", compression="gzip", usecols=["target_uniprot_id", "contributing_membrane_databases_v63"])
    cross = cross.loc[cross.target_uniprot_id.isin(formal_ids)].copy()
    contrib = _contributor_sets(cross)
    expr = pd.read_csv(EXPRESSION, sep="\t", compression="gzip", usecols=["target_uniprot_id"]) 
    expr_ids = set(expr.target_uniprot_id.dropna().unique()) & formal_ids
    loc = pd.read_csv(LOCALIZATION, sep="\t", compression="gzip", usecols=["target_uniprot_id"])
    loc_ids = set(loc.target_uniprot_id.dropna().unique()) & formal_ids
    site = pd.read_csv(SITES, sep="\t", compression="gzip", usecols=["target_uniprot_id", "source_database"])
    site = site.loc[site.target_uniprot_id.isin(formal_ids)].copy()
    pdbe_ids = set(site.loc[site.source_database.eq("PDBe"), "target_uniprot_id"])
    pdbbind_ids = set(site.loc[site.source_database.eq("PDBbind"), "target_uniprot_id"])
    biolip_ids = set(site.loc[site.source_database.eq("BioLiP"), "target_uniprot_id"])
    models = pd.read_csv(MODELS, sep="\t", compression="gzip", usecols=["target_uniprot_id", "model_availability_status"])
    af_ids = set(models.loc[models.model_availability_status.eq("MODEL_AVAILABLE"), "target_uniprot_id"]) & formal_ids
    disease = pd.read_csv(DISEASE, sep="\t", compression="gzip", usecols=["target_uniprot_id", "source_databases"])
    disease = disease.loc[disease.target_uniprot_id.isin(formal_ids)].copy()
    ot_ids = set(disease.loc[disease.source_databases.fillna("").str.contains("Open Targets", regex=False), "target_uniprot_id"])
    uniprot_disease_ids = set(disease.loc[disease.source_databases.fillna("").str.contains("UniProt", regex=False), "target_uniprot_id"])
    complex = pd.read_csv(COMPLEX, sep="\t", compression="gzip", usecols=["default_membrane_component_ids"])
    complex_ids = set()
    for values in complex.default_membrane_component_ids.dropna().astype(str):
        complex_ids.update(v.strip() for v in values.split(";") if v.strip())
    complex_ids &= formal_ids

    classification = pd.read_csv(
        FORMAL / "01_core_v72" / "01_release_tables" / "protein_cross_classification_v72.tsv.gz",
        sep="\t", compression="gzip", usecols=["target_uniprot_id", "classification_source"],
    )
    classification = classification.loc[classification.target_uniprot_id.isin(formal_ids)].copy()
    source_text = classification.classification_source.fillna("").astype(str)
    go_ids = set(classification.loc[source_text.str.contains("Gene Ontology|/GO/", regex=True), "target_uniprot_id"])
    interpro_ids = set(classification.loc[source_text.str.contains("InterPro|Pfam", regex=True), "target_uniprot_id"])
    reactome_ids = set(classification.loc[source_text.str.contains("Reactome", regex=False), "target_uniprot_id"])
    uniprot_function_ids = set(classification.loc[source_text.str.contains("UniProt", regex=False), "target_uniprot_id"])

    columns = ["Identity /\nsequence", "Membrane\nevidence", "Functional\nannotation", "Expression", "Localization", "Membrane\nstructure", "PDB\nstructure", "Binding\nsite", "Predicted\nmodel", "Disease"]
    coverage = {
        "UniProtKB": [formal_ids, formal_ids, uniprot_function_ids, set(), set(), set(), set(), set(), set(), uniprot_disease_ids],
        "HPA": [set(), contrib["HPA membrane evidence"], set(), expr_ids, loc_ids, set(), set(), set(), set(), set()],
        "Gene Ontology": [set(), set(), go_ids, set(), set(), set(), set(), set(), set(), set()],
        "InterPro / Pfam": [set(), set(), interpro_ids, set(), set(), set(), set(), set(), set(), set()],
        "Reactome": [set(), set(), reactome_ids, set(), set(), set(), set(), set(), set(), set()],
        "HTP": [set(), contrib["HTP"], set(), set(), set(), set(), set(), set(), set(), set()],
        "Membranome": [set(), contrib["Membranome"], set(), set(), set(), set(), set(), set(), set(), set()],
        "PDBTM": [set(), contrib["PDBTM"], set(), set(), set(), contrib["PDBTM"], set(), set(), set(), set()],
        "OPM": [set(), contrib["OPM"], set(), set(), set(), contrib["OPM"], set(), set(), set(), set()],
        "PDBe": [set(), set(), set(), set(), set(), set(), pdbe_ids, pdbe_ids, set(), set()],
        "PDBbind": [set(), set(), set(), set(), set(), set(), pdbbind_ids, pdbbind_ids, set(), set()],
        "BioLiP": [set(), set(), set(), set(), set(), set(), biolip_ids, biolip_ids, set(), set()],
        "AlphaFold DB": [set(), set(), set(), set(), set(), set(), set(), set(), af_ids, set()],
        "Open Targets": [set(), set(), set(), set(), set(), set(), set(), set(), set(), ot_ids],
    }
    rows = []
    for source, sets in coverage.items():
        for module, ids in zip(columns, sets):
            rows.append({"source": source, "protein_module": module, "canonical_protein_count": len(ids), "formal_protein_denominator": len(formal_ids), "coverage_percentage": len(ids) / len(formal_ids) * 100})
    data = pd.DataFrame(rows)
    data.to_csv(SRC / "Figure1B_source_by_protein_module_coverage.tsv", sep="\t", index=False)

    display_sources = ["UniProtKB", "HPA", "Gene Ontology", "InterPro / Pfam", "Reactome", "HTP", "Membranome", "PDBTM", "OPM", "PDBe", "PDBbind", "BioLiP", "AlphaFold DB", "Open Targets"]
    matrix = data.pivot(index="source", columns="protein_module", values="canonical_protein_count").reindex(index=display_sources, columns=columns)
    logmat = np.log1p(matrix.values.astype(float))
    cmap = LinearSegmentedColormap.from_list("mempro_source", ["#F8FAFB", "#d0e2c0", "#a2c986", "#67a583"])
    fig = plt.figure(figsize=(7.01, 6.35), constrained_layout=False)
    ax = fig.add_axes([0.22, 0.28, 0.64, 0.57])
    im = ax.imshow(logmat, cmap=cmap, vmin=0, vmax=np.log1p(7800), aspect="auto")
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(columns, fontsize=6.8, rotation=28, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(len(display_sources)))
    ax.set_yticklabels(display_sources, fontsize=8)
    ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False, pad=6)
    ax.tick_params(axis="y", left=False, pad=4)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            n = int(matrix.iloc[i, j])
            if n:
                label = f"{n / 1000:.1f}k" if n >= 1000 else f"{n:,}"
                ax.text(j, i, label, ha="center", va="center", fontsize=6.9, color=WHITE if logmat[i, j] > 5.4 else TEXT, weight="bold" if n >= 1000 else "normal")
    for x in np.arange(-.5, len(columns), 1): ax.axvline(x, color=WHITE, lw=0.9)
    for y in np.arange(-.5, len(display_sources), 1): ax.axhline(y, color=WHITE, lw=0.9)
    for spine in ax.spines.values(): spine.set_visible(False)
    cax = fig.add_axes([0.89, 0.33, 0.018, 0.45])
    cb = fig.colorbar(im, cax=cax, ticks=np.log1p([0, 100, 1000, 7800]))
    cb.ax.set_yticklabels(["0", "100", "1k", "7.8k"])
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("log1p(unique canonical proteins)", fontsize=7.2, labelpad=6)
    fig.text(0.05, 0.95, "Source-by-module protein coverage", fontsize=11, weight="bold", color=NAVY)
    fig.text(0.05, 0.915, "Cell value: unique formal canonical proteins covered by that source for the indicated data module (n = 7,800).", fontsize=7.2, color="#5E7182")
    fig.text(0.22, 0.06, "Blank = no direct protein-level coverage counted in this module; colour is log1p(unique canonical proteins).", fontsize=7.1, color="#5E7182")
    save_all(fig, "Figure1B_source_by_module_coverage_heatmap")
    return data


def figure_01c_role_function(proteins: pd.DataFrame) -> pd.DataFrame:
    roles = load_current_roles(proteins)
    cross = pd.read_csv(CROSS, sep="\t", compression="gzip", usecols=["target_uniprot_id", "molecular_function_primary_v63"])
    data = proteins[["target_uniprot_id"]].merge(roles, on="target_uniprot_id", how="left").merge(cross, on="target_uniprot_id", how="left")
    data["molecular_function"] = clean_string(data["molecular_function_primary_v63"]).replace("", "unclassified")
    # Keep the functional census complete: only the three extremely sparse role labels are collapsed.
    role_n = data.membrane_role.value_counts()
    data["role_display"] = data.membrane_role.where(data.membrane_role.map(role_n).ge(20), "other low-frequency roles")
    role_order = data.role_display.value_counts().index.tolist()
    fun_order = data.molecular_function.value_counts().index.tolist()
    rows = []
    pvals = []
    for role in role_order:
        for fun in fun_order:
            a = int(((data.role_display == role) & (data.molecular_function == fun)).sum())
            b = int((data.role_display == role).sum() - a)
            c = int((data.molecular_function == fun).sum() - a)
            d = int(len(data) - a - b - c)
            _, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            # Haldane-Anscombe correction guarantees finite log2 odds for display.
            log2or = math.log2(((a + .5) * (d + .5)) / ((b + .5) * (c + .5)))
            expected = (a + b) * (a + c) / len(data)
            rows.append({"membrane_role": role, "molecular_function": fun, "observed_protein_count": a, "expected_count": expected, "log2_odds_ratio": log2or, "fisher_p": p, "support_n": a})
            pvals.append(p)
    stats = pd.DataFrame(rows)
    stats["bh_fdr_q"] = bh_fdr(pvals)
    stats["low_support"] = stats.support_n < 20
    stats.to_csv(SRC / "Figure1C_membrane_role_by_molecular_function_enrichment.tsv", sep="\t", index=False)

    mat = stats.pivot(index="membrane_role", columns="molecular_function", values="log2_odds_ratio").reindex(index=role_order, columns=fun_order)
    qmat = stats.pivot(index="membrane_role", columns="molecular_function", values="bh_fdr_q").reindex(index=role_order, columns=fun_order)
    nmat = stats.pivot(index="membrane_role", columns="molecular_function", values="support_n").reindex(index=role_order, columns=fun_order)
    fig = plt.figure(figsize=(7.01, 7.01), constrained_layout=False)
    ax = fig.add_axes([0.31, 0.50, 0.57, 0.28])
    vmax = max(2.0, np.nanquantile(np.abs(mat.values), .97))
    cmap = LinearSegmentedColormap.from_list("mempro_div", ["#7b95c6", "#F7F8F8", "#f47254"])
    im = ax.imshow(mat.values, cmap=cmap, norm=Normalize(vmin=-vmax, vmax=vmax), aspect="auto")
    ax.set_xticks(np.arange(len(fun_order)))
    function_codes = [f"F{i+1}" for i in range(len(fun_order))]
    ax.set_xticklabels(function_codes, fontsize=7.2, rotation=0, ha="center", weight="bold")
    ax.set_yticks(np.arange(len(role_order)))
    ax.set_yticklabels([textwrap.fill(humanize(x), 23) for x in role_order], fontsize=7.35)
    ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False, pad=8)
    ax.tick_params(axis="y", left=False, pad=3)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if qmat.iloc[i, j] < .05 and nmat.iloc[i, j] >= 20:
                ax.plot(j, i, marker="o", markersize=3.2, color="#111111", markeredgewidth=0, zorder=3)
    for x in np.arange(-.5, len(fun_order), 1): ax.axvline(x, color=WHITE, lw=.7)
    for y in np.arange(-.5, len(role_order), 1): ax.axhline(y, color=WHITE, lw=.7)
    for spine in ax.spines.values(): spine.set_visible(False)
    cax = fig.add_axes([0.91, 0.54, .018, .20])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("log2 odds ratio", fontsize=7.2, labelpad=4)
    fig.text(.05, .945, "Membrane role and molecular-function organization", fontsize=11, weight="bold", color=NAVY)
    fig.text(.05, .912, "Values are Fisher exact-test log2 odds ratios across the 7,800-protein formal set; dot: BH-FDR < 0.05 and observed support ≥20.", fontsize=7.1, color="#5E7182")
    # Full, non-truncated function labels use a compact three-column key; F-codes
    # prevent long ontology labels from colliding below 178 mm print width.
    key = [(code, humanize(fun)) for code, fun in zip(function_codes, fun_order)]
    for idx, (code, label) in enumerate(key):
        col, row = idx % 3, idx // 3
        wrapped = textwrap.fill(label, width=21, subsequent_indent="    ")
        fig.text(.31 + col * .195, .39 - row * .105, f"{code}  {wrapped}", fontsize=6.5, color=TEXT, va="top")
    fig.text(.31, .055, "Black dot: BH-FDR < 0.05 with observed support ≥20. Unresolved is retained; it is not imputed as a role.", fontsize=6.8, color="#5E7182")
    save_all(fig, "Figure1C_membrane_role_molecular_function_enrichment")
    return stats


def write_notes(radial: pd.DataFrame, coverage: pd.DataFrame, enrich: pd.DataFrame) -> None:
    notes = f"""# MemPro protein-module figure package

## Figure contract

These three standalone quantitative figures support one claim: MemPro's formal protein layer contains **7,800 canonical membrane proteins** with traceable membrane evidence, multi-source annotation coverage, and functional organization beyond a single family label.

## Output figures

1. `Figure1A_membrane_class_evidence_radial`: A/B/C membrane classes, each split into E1/E2/E3 evidence tiers. Each annular arc uses a common count scale; individual annuli are not mutually stacked.
2. `Figure1B_source_by_module_coverage_heatmap`: source-specific coverage measured as unique formal canonical proteins, rather than evidence rows. Blank cells are no direct protein-level coverage counted in the relevant input tables.
3. `Figure1C_membrane_role_molecular_function_enrichment`: Haldane–Anscombe-corrected Fisher-test log2 odds ratios. A black dot requires observed support >=20 and BH-FDR q<0.05. `family_defined_membrane_role_unresolved` is intentionally visible.

## Data integrity

- Frozen core protein denominator: 7,800.
- No rows were changed, inferred, or reclassified.
- Source identifiers in the heatmap are derived from source-linked protein fields; source rows and evidence records were not used as the statistical unit.
- The heatmap does not claim that contributing database count equals independent experimental count.
"""
    (ROOT / "README.md").write_text(notes, encoding="utf-8")
    summary = {
        "formal_protein_denominator": 7800,
        "radial_rows": int(len(radial)),
        "coverage_rows": int(len(coverage)),
        "enrichment_cells": int(len(enrich)),
        "significant_enrichment_cells": int(((enrich.bh_fdr_q < 0.05) & (enrich.support_n >= 20)).sum()),
        "exports": ["SVG (editable text)", "PDF (editable TrueType text)", "PNG 600 dpi", "TIFF 600 dpi"],
    }
    (QA / "build_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    manifest = []
    for path in sorted(ROOT.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.tsv":
            manifest.append({"path": str(path.relative_to(ROOT)), "sha256": sha256(path), "bytes": path.stat().st_size})
    pd.DataFrame(manifest).to_csv(QA / "SHA256SUMS.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build standalone MemPro protein-module figures.")
    parser.add_argument("--figure", choices=["radial", "source_heatmap", "role_function", "all"], default="all", help="Render only one reviewed figure at a time, or all.")
    args = parser.parse_args()
    proteins = load_proteins()
    radial = pd.DataFrame()
    coverage = pd.DataFrame()
    enrichment = pd.DataFrame()
    if args.figure in {"radial", "all"}:
        radial = figure_01_radial(proteins)
    if args.figure in {"source_heatmap", "all"}:
        coverage = figure_01b_source_heatmap(proteins)
    if args.figure in {"role_function", "all"}:
        enrichment = figure_01c_role_function(proteins)
    if args.figure == "all":
        write_notes(radial, coverage, enrichment)
    print(f"Generated {args.figure} figure(s) from frozen FORMAL tables.")


if __name__ == "__main__":
    main()
