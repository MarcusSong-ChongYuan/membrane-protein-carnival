from __future__ import annotations

import os
import re
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(os.environ.get("MEMPRO_REVISION_ROOT", r"D:\finale\07_M1-M8图表修订_20260806"))
V62 = Path(r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730")
V63 = Path(r"D:\finale\02_V6.3_candidate_20260803")
COL = {"navy":"#334E60", "blue":"#6F8798", "teal":"#4F8F88", "orange":"#C98257", "purple":"#80739E", "gray":"#9AA1A6", "text":"#2F3437"}


def truth(s):
    return s.astype(str).str.lower().isin(["1", "true", "yes"])


def safe_int(v):
    return 0 if pd.isna(v) else int(v)


def sclass(r):
    formula = "" if pd.isna(r.molecular_formula) else str(r.molecular_formula)
    rings = safe_int(r.ring_count)
    maxring = safe_int(r.max_ring_size)
    heavy = safe_int(r.heavy_atom_count)
    amide = safe_int(r.amide_bond_count)
    if not re.search("C", formula): return "inorganic / no carbon"
    if heavy >= 25 and amide >= 4: return "peptide-like"
    if maxring >= 12: return "macrocycle"
    if rings >= 4: return "polycyclic (>=4 rings)"
    if rings >= 2: return "polycyclic (2-3 rings)"
    if rings == 1: return "monocyclic"
    return "acyclic organic"


def panel(ax, label, title, src):
    ax.text(-.08, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top", color=COL["navy"])
    ax.set_title(title, loc="left", pad=8)
    ax.text(0, -.20, "Source: " + src, transform=ax.transAxes, fontsize=6.2, color="#666D72", va="top")


def clean(ax, axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)
    ax.grid(True, axis=axis, alpha=.5, color="#D6D9DC")


def labels(ax, vals):
    m = max(vals) if len(vals) else 1
    for i, v in enumerate(vals):
        ax.text(v + m * .012, i, f"{int(v):,}", va="center", fontsize=6.3)


def main():
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update({
        "font.family":"sans-serif", "font.sans-serif":["Arial", "Microsoft YaHei", "DejaVu Sans"],
        "font.size":8.2, "axes.titlesize":9.6, "axes.titleweight":"bold", "axes.labelsize":8.3,
        "xtick.labelsize":7.1, "ytick.labelsize":7.1, "text.color":COL["text"],
        "pdf.fonttype":42, "svg.fonttype":"none",
    })
    cols = ["compound_scope_status","record_qc_status","identity_confidence","molecular_formula","molecular_weight","xlogp","tpsa","hbond_donor_count","hbond_acceptor_count","rotatable_bond_count","heavy_atom_count","ring_count","max_ring_size","amide_bond_count","is_approved_drug","is_clinical_candidate","is_endogenous_ligand","is_natural_product","is_chemical_probe"]
    nums = ["molecular_weight","xlogp","tpsa","hbond_donor_count","hbond_acceptor_count","rotatable_bond_count","heavy_atom_count","ring_count","max_ring_size","amide_bond_count"]
    ident, classes, status, desc = Counter(), Counter(), Counter(), []
    total = 0
    for ch in pd.read_csv(V62 / "small_molecule_master_v1_3.tsv", sep="\t", usecols=cols, chunksize=140000, low_memory=False):
        d = ch[ch.compound_scope_status.eq("core") & ch.record_qc_status.eq("ok")].copy()
        total += len(d)
        for c in nums: d[c] = pd.to_numeric(d[c], errors="coerce")
        ident.update(d.identity_confidence.fillna("missing"))
        classes.update(d.apply(sclass, axis=1))
        for name, c in [("Approved","is_approved_drug"),("Clinical","is_clinical_candidate"),("Endogenous","is_endogenous_ligand"),("Natural","is_natural_product"),("Probe","is_chemical_probe")]:
            status[name] += int(truth(d[c]).sum())
        v = d[["molecular_weight","xlogp","tpsa","rotatable_bond_count"]].dropna()
        if len(v) and sum(len(x) for x in desc) < 50000:
            desc.append(v.sample(min(4000, len(v)), random_state=42))
    ident = pd.Series(ident)
    classes = pd.Series(classes)
    status = pd.Series(status)
    desc = pd.concat(desc).head(50000)

    fig, ax = plt.subplots(2, 3, figsize=(14.2, 8.25), constrained_layout=True)
    ax = ax.ravel()
    fig.suptitle("M5  Compound identity, deterministic classes and evaluated PAGTN space", x=.015, y=1.015, ha="left", color=COL["navy"], fontweight="bold")

    panel(ax[0], "A", "Canonical identity confidence", "MemPro compound master v1.3")
    x = ident.sort_values(); ax[0].barh(x.index, x.values, color=COL["blue"], alpha=.84); labels(ax[0], x.values); ax[0].set_xlabel("Core QC canonical compounds"); clean(ax[0], "x")

    panel(ax[1], "B", "Biological-status annotations overlap", "ChEMBL/ChEBI/PubChem and source annotations")
    x = status.sort_values(); ax[1].barh(x.index, x.values, color=sns.color_palette("Paired", len(x), desat=.72), alpha=.84); labels(ax[1], x.values); ax[1].set_xlabel("Annotated compounds (non-exclusive)"); clean(ax[1], "x")

    panel(ax[2], "C", "One deterministic structural classification", "Descriptor-based mutually exclusive rules")
    x = classes.sort_values(); ax[2].barh(x.index, x.values, color=COL["teal"], alpha=.84); labels(ax[2], x.values); ax[2].set_xlabel("Core QC canonical compounds"); clean(ax[2], "x")
    classes.rename("compound_count").to_csv(ROOT / "02_data/M5C_recomputed_structural_classes.tsv", sep="\t")

    panel(ax[3], "D", "Descriptor distributions guide docking", "Calculated properties; 10th/50th/90th percentiles")
    q = desc.quantile([.1, .5, .9]).T; q.columns = ["p10", "median", "p90"]
    z = q.div(q["median"].replace(0, 1), axis=0)
    sns.heatmap(z, cmap="flare", annot=q.round(1), fmt=".1f", linewidths=.4, cbar_kws={"label":"relative to median"}, ax=ax[3]); ax[3].set_xlabel(""); ax[3].set_ylabel("")

    panel(ax[4], "E", "Property-trained PAGTN embedding + HDBSCAN", "Fixed QC sample n=2,400; seed=42; Morgan baseline retained in QA")
    coord = pd.read_csv(ROOT / "02_data/M5E_pagtn_morgan_coordinates.tsv", sep="\t")
    labels_p = coord.pagtn_cluster.astype(int)
    clusters = sorted([c for c in labels_p.unique() if c >= 0])
    palette = dict(zip(clusters, sns.color_palette("Paired", max(2, len(clusters)), desat=.68)))
    colors = ["#B8BDC1" if c < 0 else palette[c] for c in labels_p]
    ax[4].scatter(coord.pagtn_umap1, coord.pagtn_umap2, c=colors, s=8, alpha=.32, linewidths=0)
    ax[4].set_xlabel("UMAP 1 (no direct chemical unit)"); ax[4].set_ylabel("UMAP 2 (no direct chemical unit)")
    ax[4].grid(False); ax[4].spines[["top","right"]].set_visible(False)
    metrics = pd.read_csv(ROOT / "02_data/M5E_pagtn_cluster_summary.tsv", sep="\t")
    def metric(method, name):
        return float(metrics[(metrics.method == method) & (metrics.metric == name)].value.iloc[0])
    summary = (
        f"PAGTN: silhouette {metric('PAGTN_property_supervised','silhouette_nonnoise'):.3f}; "
        f"noise {metric('PAGTN_property_supervised','noise_fraction'):.1%}\n"
        f"10-NN scaffold purity: PAGTN {metric('PAGTN_property_supervised','10NN_scaffold_purity'):.1%}; "
        f"Morgan {metric('Morgan_radius2_1024','10NN_scaffold_purity'):.1%}"
    )
    ax[4].text(.02, .98, summary, transform=ax[4].transAxes, va="top", fontsize=6.5, color="#4A5054", bbox={"facecolor":"white","edgecolor":"#D7DADD","alpha":.82,"boxstyle":"round,pad=.28"})

    panel(ax[5], "F", "Scaffold diversity remains directly interpretable", "Bemis-Murcko scaffolds, V6.3 candidate")
    sf = pd.read_csv(V63 / "04_disease/scaffold_frequency_v0_1.tsv", sep="\t").head(1000)
    ax[5].plot(sf.scaffold_rank, sf.compound_count, color=COL["blue"], lw=1.5)
    ax[5].set_xscale("log"); ax[5].set_yscale("log"); ax[5].set_xlabel("Scaffold rank (log)"); ax[5].set_ylabel("Compounds/scaffold (log)"); clean(ax[5], "both")

    for ext, folder, kwargs in [("png","03_figures/png",{"dpi":600}),("svg","03_figures/svg",{}),("pdf","03_figures/pdf",{})]:
        fig.savefig(ROOT / folder / f"M5_chemical_identity_space_revised.{ext}", bbox_inches="tight", **kwargs)
    plt.close(fig)
    print({"status":"PASS", "core_qc":total, "pagtn_sample":len(coord), "classes":classes.to_dict()})


if __name__ == "__main__":
    main()
