from __future__ import annotations

import hashlib
import json
import math
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import squarify
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.ML.Cluster import Butina
from scipy import sparse
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler, StandardScaler


FORMAL = Path(r"D:\finale\FORMAL")
OUT = Path(r"D:\finale\compound_classification_exploration")
ANALYSIS = OUT / "02_classification_tables/INTERACTION_LINKED_COMPOUND_ANALYSIS.tsv.gz"
ROLE_COUNTS = OUT / "04_heatmaps/COMPOUND_ROLE_COUNTS.tsv.gz"
PAIR = FORMAL / "01_core_v72/01_release_tables/protein_compound_pair_v72.tsv.gz"
ROLE = FORMAL / "03_publication_repairs/protein_membrane_role_FORMAL.tsv"
SEED = 42
RNG = np.random.default_rng(SEED)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})
sns.set_theme(style="white", context="paper", font_scale=0.85)

PALETTE = ["#7b95c6", "#49c2d9", "#67a583", "#a2c986", "#fded95", "#ffc1a6", "#f59c7c", "#c85e62", "#9b83b4"]
UNKNOWN = "#b8b8b8"

ROLE_ORDER = [
    "membrane_associated_enzyme", "receptor", "transporter", "ion_channel",
    "signaling_regulator", "membrane_scaffold_or_linker", "adhesion_recognition",
    "immune_or_cell_recognition", "membrane_trafficking", "membrane_organizer",
    "junction_or_adhesion", "membrane_lipid_transfer", "membrane_lipid_translocase",
    "family_defined_membrane_role_unresolved",
]


def stable_hash(s: str) -> int:
    return int(hashlib.sha256(s.encode("utf-8")).hexdigest()[:16], 16)


def export(fig: plt.Figure, base: Path, png_dpi: int = 350, svg: bool = True, pdf: bool = True) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".png"), dpi=png_dpi, bbox_inches="tight", facecolor="white")
    if svg:
        fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    if pdf:
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def deterministic_stratified_sample(df: pd.DataFrame, category: str, n_total: int, min_per_group: int = 250) -> pd.DataFrame:
    d = df.copy()
    d["_hash"] = d["compound_internal_id"].map(stable_hash)
    counts = d[category].value_counts(dropna=False)
    alloc = {}
    for group, count in counts.items():
        alloc[group] = min(count, max(min_per_group, int(round(n_total * count / len(d)))))
    while sum(alloc.values()) > n_total:
        largest = max((g for g in alloc if alloc[g] > min(counts[g], min_per_group)), key=lambda g: alloc[g], default=None)
        if largest is None:
            break
        alloc[largest] -= 1
    while sum(alloc.values()) < n_total:
        candidates = [g for g in alloc if alloc[g] < counts[g]]
        if not candidates:
            break
        g = max(candidates, key=lambda x: counts[x] - alloc[x])
        alloc[g] += 1
    parts = [d[d[category] == g].nsmallest(n, "_hash") for g, n in alloc.items() if n]
    return pd.concat(parts, ignore_index=True).drop(columns="_hash")


def role_display(x: str) -> str:
    return x.replace("family_defined_membrane_role_unresolved", "Unresolved role").replace("_", " ").title()


def gini(values: np.ndarray) -> float:
    x = np.sort(np.asarray(values, dtype=float))
    if x.size == 0 or x.sum() == 0:
        return float("nan")
    return float((2 * np.arange(1, x.size + 1) @ x) / (x.size * x.sum()) - (x.size + 1) / x.size)


def robust_z(df: pd.DataFrame) -> pd.DataFrame:
    med = df.median(axis=0)
    mad = (df - med).abs().median(axis=0).replace(0, np.nan)
    z = (df - med) / (1.4826 * mad)
    return z.fillna(0).clip(-3, 3)


def plot_v1(analysis: pd.DataFrame) -> None:
    out = OUT / "02_classification_tables"
    counts = analysis["chemical_regime_local"].value_counts().rename_axis("chemical_regime_local").reset_index(name="compound_count")
    counts["denominator"] = len(analysis)
    counts["percent"] = counts["compound_count"] / len(analysis) * 100
    tags = analysis.assign(structural_tag=analysis["authoritative_chebi_direct_class"].where(
        analysis["authoritative_chebi_direct_class"] != "NOT_AVAILABLE", "No authoritative direct class"
    )).groupby(["chemical_regime_local", "structural_tag"]).size().rename("compound_count").reset_index()
    counts.to_csv(out / "V1_chemical_taxonomy_source.tsv", sep="\t", index=False)
    tags.to_csv(out / "V1_chemical_taxonomy_hierarchy_source.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    d = counts.sort_values("compound_count")
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(d))]
    ax.barh(d["chemical_regime_local"], d["compound_count"], color=colors, alpha=0.82)
    for y, v in enumerate(d["compound_count"]):
        ax.text(v, y, f" {v:,} ({v/len(analysis)*100:.1f}%)", va="center", fontsize=6.5)
    ax.set_xlabel("Interaction-linked compounds")
    ax.set_ylabel("")
    ax.set_title("Local rule-based chemical regimes", loc="left", fontweight="bold")
    ax.text(1, -0.2, f"n = {len(analysis):,}; not ClassyFire", transform=ax.transAxes, ha="right", fontsize=6, color="#555555")
    export(fig, out / "V1A_chemical_taxonomy_ranked_bar")

    # Nested treemap: regime rectangles subdivided by top ChEBI direct labels; long tail grouped.
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.set_axis_off()
    top_rects = squarify.squarify(squarify.normalize_sizes(counts.compound_count.tolist(), 100, 60), 0, 0, 100, 60)
    for i, (rect, row) in enumerate(zip(top_rects, counts.itertuples())):
        color = PALETTE[i % len(PALETTE)]
        ax.add_patch(plt.Rectangle((rect["x"], rect["y"]), rect["dx"], rect["dy"], facecolor=color, alpha=0.33, edgecolor="white", lw=1.2))
        sub = tags[tags.chemical_regime_local == row.chemical_regime_local].sort_values("compound_count", ascending=False)
        sub_top = sub.head(5).copy()
        other = sub.iloc[5:].compound_count.sum()
        if other:
            sub_top = pd.concat([sub_top, pd.DataFrame({"chemical_regime_local":[row.chemical_regime_local], "structural_tag":["Other / unmapped"], "compound_count":[other]})])
        if rect["dx"] > 8 and rect["dy"] > 6:
            inner = squarify.squarify(squarify.normalize_sizes(sub_top.compound_count.tolist(), rect["dx"], rect["dy"]), rect["x"], rect["y"], rect["dx"], rect["dy"])
            for sr, srow in zip(inner, sub_top.itertuples()):
                ax.add_patch(plt.Rectangle((sr["x"], sr["y"]), sr["dx"], sr["dy"], facecolor=color, alpha=0.33, edgecolor="white", lw=0.35))
                if sr["dx"] > 9 and sr["dy"] > 5:
                    ax.text(sr["x"]+0.5, sr["y"]+0.5, textwrap.fill(str(srow.structural_tag), 18), fontsize=5.5, va="bottom")
        if rect["dx"] > 10 and rect["dy"] > 7:
            ax.text(rect["x"]+0.6, rect["y"]+rect["dy"]-0.8, textwrap.fill(row.chemical_regime_local, 20), fontsize=7, fontweight="bold", va="top")
    ax.set_xlim(0, 100); ax.set_ylim(0, 60)
    ax.set_title("Local regime → available ChEBI direct class", loc="left", fontweight="bold")
    export(fig, out / "V1B_chemical_taxonomy_hierarchy")


def plot_v2(analysis: pd.DataFrame) -> None:
    out = OUT / "03_scaffold"
    for code, col, title in [
        ("V2A_exact", "exact_scaffold_id", "Exact Bemis–Murcko scaffold rank–abundance"),
        ("V2B_generic", "generic_scaffold_id", "Generic Bemis–Murcko scaffold rank–abundance"),
    ]:
        c = analysis[col].value_counts().rename_axis(col).reset_index(name="compound_count")
        c["rank"] = np.arange(1, len(c)+1)
        c.to_csv(out / f"{code}_scaffold_rank_abundance_source.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6.5, 4.0))
        ax.plot(c["rank"], c["compound_count"], color="#7b95c6", lw=1.2)
        ax.fill_between(c["rank"], 1, c["compound_count"], color="#a1d8e8", alpha=0.35)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Scaffold rank (log scale)"); ax.set_ylabel("Compounds per scaffold (log scale)")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.text(0.99, 0.98, f"{len(c):,} scaffolds\nGini = {gini(c.compound_count.values):.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=6.5)
        for row in c.head(7).itertuples():
            ax.annotate(f"{row.compound_count:,}", (row.rank, row.compound_count), xytext=(3, 2), textcoords="offset points", fontsize=5.5)
        export(fig, out / f"{code}_scaffold_rank_abundance")


def plot_v3(analysis: pd.DataFrame, pair_role: pd.DataFrame) -> None:
    out = OUT / "03_scaffold"
    sc = analysis.groupby(["generic_scaffold_id", "generic_murcko_smiles_recalculated"]).agg(
        compound_count=("compound_internal_id", "nunique")
    ).reset_index().sort_values("compound_count", ascending=False)
    top = sc[sc.generic_scaffold_id != "GMS-ACYCLIC"].head(16).copy()
    map_comp = analysis[["compound_internal_id", "generic_scaffold_id"]]
    pr = pair_role.merge(map_comp, on="compound_internal_id", how="inner")
    target_n = pr[pr.generic_scaffold_id.isin(top.generic_scaffold_id)].groupby("generic_scaffold_id").target_uniprot_id.nunique()
    role_n = pr[pr.generic_scaffold_id.isin(top.generic_scaffold_id)].groupby(["generic_scaffold_id", "formal_primary_membrane_role"]).size().rename("n").reset_index()
    dom = role_n.sort_values(["generic_scaffold_id", "n", "formal_primary_membrane_role"], ascending=[True,False,True]).groupby("generic_scaffold_id").first()["formal_primary_membrane_role"]
    top["protein_target_count"] = top.generic_scaffold_id.map(target_n).fillna(0).astype(int)
    top["dominant_role_descriptive"] = top.generic_scaffold_id.map(dom).fillna("not available")
    top.to_csv(out / "V3_top_scaffold_gallery_source.tsv", sep="\t", index=False)

    mols, legends = [], []
    for r in top.itertuples():
        mol = Chem.MolFromSmiles(r.generic_murcko_smiles_recalculated)
        mols.append(mol)
        legends.append(f"{r.generic_scaffold_id}\n{r.compound_count:,} compounds | {r.protein_target_count:,} targets\n{role_display(r.dominant_role_descriptive)}")
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(520, 360), legends=legends, useSVG=False)
    img.save(out / "V3_top_scaffold_gallery.png")
    # Place the raster molecular grid into a PDF review sheet; molecule source SMILES remain editable/reproducible in TSV.
    fig, ax = plt.subplots(figsize=(11.5, 8.0)); ax.imshow(img); ax.axis("off")
    fig.savefig(out / "V3_top_scaffold_gallery.pdf", bbox_inches="tight")
    plt.close(fig)


def heatmap_fig(matrix: pd.DataFrame, title: str, cmap: str, cbar_label: str, base: Path, center=None) -> None:
    h = max(3.5, 0.36 * len(matrix) + 1.5)
    fig, ax = plt.subplots(figsize=(8.1, h))
    sns.heatmap(matrix, cmap=cmap, center=center, ax=ax, cbar_kws={"label": cbar_label}, linewidths=0.25, linecolor="white")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("Formal primary membrane role"); ax.set_ylabel("")
    ax.set_xticklabels([role_display(x) for x in matrix.columns], rotation=42, ha="right")
    export(fig, base)


def plot_v4_v5(analysis: pd.DataFrame, role_counts: pd.DataFrame) -> None:
    out = OUT / "04_heatmaps"
    cr = role_counts.merge(analysis[["compound_internal_id", "chemical_regime_local", "generic_scaffold_id"]], on="compound_internal_id", how="inner")
    cols = [r for r in ROLE_ORDER if r in set(cr.formal_primary_membrane_role)]
    count = cr.pivot_table(index="chemical_regime_local", columns="formal_primary_membrane_role", values="protein_target_count_in_role", aggfunc="sum", fill_value=0).reindex(columns=cols, fill_value=0)
    count = count.loc[count.sum(axis=1).sort_values(ascending=False).index]
    profile = count.div(count.sum(axis=1), axis=0)
    count.to_csv(out / "V4A_chemical_class_role_count_source.tsv", sep="\t")
    profile.to_csv(out / "V4B_chemical_class_role_profile_source.tsv", sep="\t")
    heatmap_fig(np.log1p(count), "Local chemical regime × membrane role: pair counts", "flare", "log1p(pair count)", out / "V4A_chemical_class_role_count_heatmap")
    heatmap_fig(profile, "Local chemical regime × membrane role: within-regime profile", "crest", "Within-row proportion", out / "V4B_chemical_class_role_profile_heatmap")

    sc_count = cr.pivot_table(index="generic_scaffold_id", columns="formal_primary_membrane_role", values="protein_target_count_in_role", aggfunc="sum", fill_value=0).reindex(columns=cols, fill_value=0)
    support = sc_count.sum(axis=1).sort_values(ascending=False)
    top100 = sc_count.loc[support.head(100).index]
    top100.div(top100.sum(axis=1), axis=0).to_csv(out / "V5_scaffold_role_matrix.tsv", sep="\t")
    for n in [30, 50, 100]:
        prof = top100.head(n).div(top100.head(n).sum(axis=1), axis=0)
        prof.to_csv(out / f"V5_scaffold_role_matrix_top{n}.tsv", sep="\t")
    prof = top100.head(50).div(top100.head(50).sum(axis=1), axis=0)
    heatmap_fig(prof, "Top 50 generic scaffolds × membrane role", "crest", "Within-scaffold proportion", out / "V5A_scaffold_role_heatmap_ordered")
    row_order = leaves_list(linkage(prof.values, method="average", metric="cosine"))
    row_clustered = prof.iloc[row_order]
    heatmap_fig(row_clustered, "Top 50 scaffolds: role-profile row clustering", "crest", "Within-scaffold proportion", out / "V5B_scaffold_role_heatmap_rowclustered")
    col_order = leaves_list(linkage(prof.values.T, method="average", metric="cosine"))
    bic = prof.iloc[row_order, col_order]
    heatmap_fig(bic, "Top 50 scaffolds: row + column profile clustering", "crest", "Within-scaffold proportion", out / "V5C_scaffold_role_heatmap_biclustered")
    (out / "V5_clustering_parameters.md").write_text(
        "# V5 clustering parameters\n\nRows are generic Bemis–Murcko scaffolds. Matrix values are within-scaffold proportions across formal primary membrane roles. Top 50 rows are selected by total eligible formal pair support. Average linkage with cosine distance is used. V5A keeps abundance and biological column order; V5B clusters rows only; V5C clusters rows and columns. Clusters are profile-based groups, not chemical classes.\n",
        encoding="utf-8",
    )


def plot_v6(analysis: pd.DataFrame) -> None:
    out = OUT / "04_heatmaps"
    desc = ["molecular_weight", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "ring_count", "aromatic_ring_count", "fraction_csp3", "formal_charge"]
    raw = analysis.groupby("chemical_regime_local")[desc].median()
    raw.to_csv(out / "V6_physchem_class_raw_medians.tsv", sep="\t")
    z = robust_z(raw)
    heatmap_fig(z, "Robust physicochemical profiles by local chemical regime", "vlag", "Robust z-score", out / "V6_physchem_class_heatmap", center=0)


def plot_v7(analysis: pd.DataFrame) -> pd.DataFrame:
    out = OUT / "05_chemical_space"
    desc = ["molecular_weight", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "ring_count", "aromatic_ring_count", "fraction_csp3", "formal_charge", "heavy_atom_count"]
    sample = deterministic_stratified_sample(analysis, "chemical_regime_local", 50_000, 300)
    X = sample[desc].replace([np.inf, -np.inf], np.nan).fillna(sample[desc].median())
    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=5, random_state=SEED).fit(Xs)
    coords = pca.transform(Xs)
    for i in range(5): sample[f"PCA{i+1}"] = coords[:, i]
    sample.to_csv(out / "V7_PCA_coordinates.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame({"component":[f"PCA{i+1}" for i in range(5)], "explained_variance_ratio":pca.explained_variance_ratio_}).to_csv(out / "V7_PCA_explained_variance.tsv", sep="\t", index=False)

    color_specs = [
        ("V7A_PCA_chemical_superclass", "chemical_regime_local", "Local chemical regime"),
        ("V7B_PCA_target_role", "dominant_membrane_role_descriptive", "Dominant membrane role (descriptive)"),
        ("V7C_PCA_scaffold_family", "generic_scaffold_id", "Selected major generic scaffold"),
    ]
    top_scaff = set(analysis.generic_scaffold_id.value_counts().head(8).index)
    sample["selected_scaffold"] = sample.generic_scaffold_id.where(sample.generic_scaffold_id.isin(top_scaff), "Other scaffolds")
    for base, col, title in color_specs:
        if base.startswith("V7C"): col = "selected_scaffold"
        vc = sample[col].value_counts(); top = list(vc.head(10).index)
        plot_df = sample.copy(); plot_df["_group"] = plot_df[col].where(plot_df[col].isin(top), "Other")
        fig, ax = plt.subplots(figsize=(6.6, 5.2))
        for j, (g, d) in enumerate(plot_df.groupby("_group", sort=False)):
            ax.scatter(d.PCA1, d.PCA2, s=5, alpha=0.20, color=(UNKNOWN if g == "Other" else PALETTE[j % len(PALETTE)]), label=str(g), rasterized=True)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title(title, loc="left", fontweight="bold"); ax.legend(markerscale=2.5, fontsize=5.5, bbox_to_anchor=(1.02,1), loc="upper left")
        ax.text(0.99, 0.02, f"Deterministic stratified sample n={len(sample):,}", transform=ax.transAxes, ha="right", fontsize=5.5)
        export(fig, out / base)
    return sample


def plot_v8(analysis: pd.DataFrame) -> pd.DataFrame:
    out = OUT / "05_chemical_space"
    import umap
    sample = deterministic_stratified_sample(analysis, "chemical_regime_local", 20_000, 200)
    gen = GetMorganGenerator(radius=2, fpSize=2048)
    rows, cols = [], []
    for i, smi in enumerate(sample.canonical_smiles_recalculated):
        fp = gen.GetFingerprint(Chem.MolFromSmiles(smi))
        on = list(fp.GetOnBits())
        rows.extend([i]*len(on)); cols.extend(on)
    X = sparse.csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(len(sample), 2048))
    reducer = umap.UMAP(n_neighbors=20, min_dist=0.15, n_components=2, metric="jaccard", random_state=SEED, low_memory=True)
    emb = reducer.fit_transform(X)
    sample["UMAP1"], sample["UMAP2"] = emb[:,0], emb[:,1]
    sample.to_csv(out / "V8_UMAP_coordinates.tsv.gz", sep="\t", index=False, compression="gzip")
    params = {"sample_n":len(sample), "sample":"deterministic stratified by local chemical regime", "morgan_radius":2, "morgan_bits":2048, "n_neighbors":20, "min_dist":0.15, "metric":"jaccard", "random_seed":SEED}
    (out / "V8_UMAP_parameters.json").write_text(json.dumps(params, indent=2), encoding="utf-8")
    top_scaff = set(analysis.generic_scaffold_id.value_counts().head(8).index)
    sample["selected_scaffold"] = sample.generic_scaffold_id.where(sample.generic_scaffold_id.isin(top_scaff), "Other scaffolds")
    specs = [
        ("V8A_UMAP_chemical_superclass", "chemical_regime_local", "Morgan/UMAP by local chemical regime"),
        ("V8B_UMAP_scaffold_family", "selected_scaffold", "Morgan/UMAP by abundant generic scaffold"),
        ("V8C_UMAP_target_role", "dominant_membrane_role_descriptive", "Morgan/UMAP by dominant membrane role"),
    ]
    for base, col, title in specs:
        top = list(sample[col].value_counts().head(10).index)
        d0 = sample.copy(); d0["_group"] = d0[col].where(d0[col].isin(top), "Other")
        fig, ax = plt.subplots(figsize=(6.6, 5.2))
        for j, (g, d) in enumerate(d0.groupby("_group", sort=False)):
            ax.scatter(d.UMAP1, d.UMAP2, s=4, alpha=0.22, color=(UNKNOWN if g=="Other" else PALETTE[j%len(PALETTE)]), label=str(g), rasterized=True)
        ax.set_xlabel("UMAP 1 (no direct chemical meaning)"); ax.set_ylabel("UMAP 2 (no direct chemical meaning)")
        ax.set_title(title, loc="left", fontweight="bold"); ax.legend(markerscale=2.5, fontsize=5.5, bbox_to_anchor=(1.02,1), loc="upper left")
        ax.text(0.99, 0.02, f"Same embedding; n={len(sample):,}", transform=ax.transAxes, ha="right", fontsize=5.5)
        export(fig, out / base)
    return sample


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(values); return x, np.arange(1, len(x)+1)/len(x)


def plot_v9(analysis: pd.DataFrame) -> None:
    out = OUT / "06_target_breadth"
    vals = analysis.unique_protein_target_count.fillna(0).astype(int).values
    x,y = ecdf(vals)
    pd.DataFrame({"target_count":x,"ecdf":y}).to_csv(out / "V9A_target_breadth_distribution_source.tsv",sep="\t",index=False)
    fig, ax=plt.subplots(figsize=(6.2,4)); ax.step(x,y,where="post",color="#7b95c6",lw=1.5); ax.set_xscale("symlog",linthresh=1); ax.set_xlabel("Unique formal membrane-protein targets per compound"); ax.set_ylabel("ECDF"); ax.set_title("Target breadth of interaction-linked compounds",loc="left",fontweight="bold"); export(fig,out/"V9A_target_breadth_distribution")

    regimes=analysis.chemical_regime_local.value_counts().index[:8]
    d=analysis[analysis.chemical_regime_local.isin(regimes)].copy(); d.to_csv(out/"V9B_target_breadth_by_class_source.tsv.gz",sep="\t",index=False,compression="gzip")
    fig,ax=plt.subplots(figsize=(7.4,4.6)); sns.boxenplot(data=d,x="chemical_regime_local",y="unique_protein_target_count",color="#a1d8e8",ax=ax,showfliers=False); ax.set_yscale("log"); ax.set_xlabel(""); ax.set_ylabel("Unique targets (log scale)"); ax.tick_params(axis="x",rotation=35); ax.set_title("Target breadth by local chemical regime",loc="left",fontweight="bold"); export(fig,out/"V9B_target_breadth_by_class")

    top=set(analysis.generic_scaffold_id.value_counts().head(12).index); ds=analysis[analysis.generic_scaffold_id.isin(top)].copy(); ds.to_csv(out/"V9C_target_breadth_by_scaffold_source.tsv.gz",sep="\t",index=False,compression="gzip")
    order=ds.groupby("generic_scaffold_id").unique_protein_target_count.median().sort_values(ascending=False).index
    fig,ax=plt.subplots(figsize=(8.0,4.6)); sns.boxenplot(data=ds,x="generic_scaffold_id",y="unique_protein_target_count",order=order,color="#d0e2c0",ax=ax,showfliers=False); ax.set_yscale("log"); ax.set_xlabel("Abundant generic scaffold ID"); ax.set_ylabel("Unique targets (log scale)"); ax.tick_params(axis="x",rotation=48,labelsize=5); ax.set_title("Target breadth by abundant generic scaffold",loc="left",fontweight="bold"); export(fig,out/"V9C_target_breadth_by_scaffold")


def plot_v10(analysis: pd.DataFrame) -> None:
    out=OUT/"07_optional_3d"; d=deterministic_stratified_sample(analysis,"chemical_regime_local",20_000,200); d.to_csv(out/"V10_3D_physchem_space_source.tsv.gz",sep="\t",index=False,compression="gzip")
    fig=plt.figure(figsize=(7,5.4)); ax=fig.add_subplot(111,projection="3d"); groups=list(d.chemical_regime_local.value_counts().index)
    for i,g in enumerate(groups):
        x=d[d.chemical_regime_local==g]; ax.scatter(x.molecular_weight,x.logp,x.tpsa,s=3+2*np.log10(x.unique_protein_target_count+1),alpha=.18,color=PALETTE[i%len(PALETTE)],label=g,depthshade=False)
    ax.set_xlabel("MW"); ax.set_ylabel("LogP"); ax.set_zlabel("TPSA"); ax.set_title("Exploratory 3D physicochemical space",loc="left",fontweight="bold"); ax.legend(fontsize=5,bbox_to_anchor=(1.08,1)); fig.savefig(out/"V10_3D_physchem_space.png",dpi=350,bbox_inches="tight"); plt.close(fig)
    try:
        import plotly.express as px
        figp=px.scatter_3d(d,x="molecular_weight",y="logp",z="tpsa",color="chemical_regime_local",size=np.log10(d.unique_protein_target_count+1)+0.2,opacity=.35,hover_name="compound_internal_id")
        figp.write_html(out/"V10_3D_physchem_space.html",include_plotlyjs="cdn")
    except Exception as exc:
        (out/"V10_3D_HTML_ERROR.txt").write_text(str(exc),encoding="utf-8")


def butina_benchmark(analysis: pd.DataFrame) -> None:
    out=OUT/"03_scaffold"; d=deterministic_stratified_sample(analysis,"chemical_regime_local",3000,100).reset_index(drop=True); gen=GetMorganGenerator(radius=2,fpSize=2048); fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in d.canonical_smiles_recalculated]
    dists=[]
    for i in range(1,len(fps)):
        sims=DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i]); dists.extend([1-x for x in sims])
    clusters=Butina.ClusterData(dists,len(fps),0.4,isDistData=True,reordering=True)
    assign=np.empty(len(d),dtype=int)
    for ci,members in enumerate(clusters): assign[list(members)]=ci
    d["butina_cluster"]=assign
    summaries=[]
    for ci,x in d.groupby("butina_cluster"):
        summaries.append({"cluster":ci,"size":len(x),"exact_scaffold_purity":x.exact_scaffold_id.value_counts(normalize=True).iloc[0],"chemical_regime_purity":x.chemical_regime_local.value_counts(normalize=True).iloc[0],"dominant_role_purity":x.dominant_membrane_role_descriptive.value_counts(normalize=True).iloc[0]})
    s=pd.DataFrame(summaries).sort_values("size",ascending=False); d.to_csv(out/"OPTIONAL_BUTINA_ASSIGNMENTS.tsv.gz",sep="\t",index=False,compression="gzip"); s.to_csv(out/"OPTIONAL_BUTINA_CLUSTER_SUMMARY.tsv",sep="\t",index=False)
    fig,ax=plt.subplots(figsize=(6.3,4)); ax.plot(np.arange(1,len(s)+1),s["size"],color="#7b95c6"); ax.set_yscale("log"); ax.set_xlabel("Cluster rank"); ax.set_ylabel("Cluster size (log scale)"); ax.set_title("Exploratory Morgan/Butina cluster-size distribution",loc="left",fontweight="bold"); ax.text(.99,.97,f"n={len(d):,}; Tanimoto ≥0.60",transform=ax.transAxes,ha="right",va="top",fontsize=6); export(fig,out/"OPTIONAL_BUTINA_CLUSTER_SIZE")
    (out/"OPTIONAL_BUTINA_PARAMETERS.md").write_text("# Optional Butina benchmark\n\nDeterministic stratified sample n=3,000; Morgan radius 2, 2048 bits; Butina distance threshold 0.40 (Tanimoto similarity >=0.60). These clusters are parameter-dependent exploratory neighbourhoods, not chemical classes.\n",encoding="utf-8")


def main() -> None:
    analysis=pd.read_csv(ANALYSIS,sep="\t",low_memory=False)
    role_counts=pd.read_csv(ROLE_COUNTS,sep="\t")
    pairs=pd.read_csv(PAIR,sep="\t",usecols=["target_uniprot_id","compound_internal_id"]).drop_duplicates()
    roles=pd.read_csv(ROLE,sep="\t",usecols=["canonical_uniprot_accession","formal_primary_membrane_role"])
    pair_role=pairs.merge(roles,left_on="target_uniprot_id",right_on="canonical_uniprot_accession",how="left",validate="many_to_one")
    pair_role=pair_role[pair_role.compound_internal_id.isin(set(analysis.compound_internal_id))]
    plot_v1(analysis); plot_v2(analysis); plot_v3(analysis,pair_role); plot_v4_v5(analysis,role_counts); plot_v6(analysis)
    plot_v7(analysis); plot_v8(analysis); plot_v9(analysis); plot_v10(analysis); butina_benchmark(analysis)
    print(json.dumps({"status":"PASS","analysis_compounds":len(analysis),"eligible_pairs":len(pair_role)},indent=2))


if __name__ == "__main__":
    main()
