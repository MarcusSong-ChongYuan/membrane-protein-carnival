"""Exploratory K-means comparison for the complete MemPro protein landscape.

This does not replace or modify the frozen Leiden graph communities. K-means is
run on the same 50-dimensional representation after L2 normalization, so its
Euclidean distance corresponds to the cosine geometry used by the kNN graph.
All 7,800 proteins are retained in every K-means result and static view.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import linear_sum_assignment
from scipy.stats import chi2
from sklearn.cluster import KMeans
from sklearn.covariance import MinCovDet
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.preprocessing import normalize

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "leiden_structural_results"
DEST = OUT / "kmeans_same_input_comparison"
DEST.mkdir(parents=True, exist_ok=True)
SEED = 42
TARGET_K = 28
LARGE_N = 150

PALETTE = [
    "#7B95C6", "#49C2D9", "#A1D8E8", "#67A583", "#A2C986", "#D0E2C0",
    "#FDED95", "#FFC1A6", "#F59C7C", "#F47254", "#C85E62", "#A67FC5",
    "#B58CCF", "#819CCB", "#89BFA5", "#D89B6E", "#79B5C5", "#B6A8D4",
    "#D7A2A1", "#86A98A", "#D6C480", "#A6B9D7", "#C79D90", "#8EBCD1",
    "#B8CFA2", "#D5A7CC", "#B3B3B3", "#9DA9A0",
]
OVERLAP_CMAP = LinearSegmentedColormap.from_list("mempro_overlap", ["#FFF9F5", "#F7B2A0", "#C85E62", "#55223A"])

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7, "axes.linewidth": .55, "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})


def stable_size_order(labels: np.ndarray) -> np.ndarray:
    """Rename K-means labels by descending cluster size; identities remain arbitrary."""
    sizes = pd.Series(labels).value_counts().sort_values(ascending=False)
    remap = {old: new for new, old in enumerate(sizes.index)}
    return np.asarray([remap[x] for x in labels], dtype=int)


def robust_ellipsoid(ax, points: np.ndarray, color: str):
    if len(points) < 8:
        return False
    try:
        model = MinCovDet(random_state=SEED).fit(points)
        val, vec = np.linalg.eigh(model.covariance_)
        if np.any(val <= 1e-12) or np.max(val) / np.min(val) > 1e7:
            return False
        r = np.sqrt(chi2.ppf(.90, df=3))
        u = np.linspace(0, 2*np.pi, 28); v = np.linspace(0, np.pi, 20)
        sphere = np.vstack([np.outer(np.cos(u),np.sin(v)).ravel(), np.outer(np.sin(u),np.sin(v)).ravel(), np.outer(np.ones_like(u),np.cos(v)).ravel()])
        xyz = vec @ np.diag(np.sqrt(val)*r) @ sphere + model.location_[:,None]
        xx,yy,zz = [a.reshape(len(u),len(v)) for a in xyz]
        ax.plot_surface(xx, yy, zz, color=color, alpha=.043, linewidth=0, antialiased=False, shade=False)
        return True
    except Exception:
        return False


def export(fig, stem: Path):
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")


def main():
    membership = pd.read_csv(OUT / "protein_leiden_structural_umap3d.tsv", sep="\t")
    x = np.load(OUT / "leiden_input_pca50.npy")
    if len(membership) != 7800 or len(x) != 7800:
        raise RuntimeError("Expected exactly 7,800 proteins in the frozen comparison input.")
    if membership["macro_community"].nunique() != 28:
        raise RuntimeError("Expected 28 pre-existing Leiden communities.")

    # Cosine kNN in the Leiden pipeline is based on direction. Unit normalization
    # makes ordinary Euclidean K-means a comparable (not identical) alternative.
    z = normalize(x, norm="l2")
    pca = PCA(n_components=3, random_state=SEED).fit(z)
    coords = pca.transform(z)
    evr = pca.explained_variance_ratio_ * 100
    ks = list(range(20, 41))
    rng = np.random.default_rng(SEED)
    sil_idx = np.sort(rng.choice(len(z), size=2000, replace=False))
    scan_path = DEST / "kmeans_k20_to_k40_scan.tsv"
    cached_membership_path = DEST / "kmeans28_all_proteins_source_data.tsv"
    if scan_path.exists() and cached_membership_path.exists():
        # The expensive deterministic K scan may already have completed before
        # an export-only interruption. Reuse it only when all 7,800 labels exist.
        cached = pd.read_csv(cached_membership_path, sep="\t")
        required = {"uniprot_id", "kmeans28_cluster", "KMeans_PC1", "KMeans_PC2", "KMeans_PC3"}
        if len(cached) == 7800 and required.issubset(cached.columns) and cached.uniprot_id.nunique() == 7800:
            membership = cached
            scan = pd.read_csv(scan_path, sep="\t")
        else:
            scan_path.unlink(missing_ok=True)
            cached_membership_path.unlink(missing_ok=True)
            scan = None
    else:
        scan = None
    if scan is None:
        rows = []
        fits: dict[int, np.ndarray] = {}
        for k in ks:
            km = KMeans(n_clusters=k, n_init=50, max_iter=500, algorithm="lloyd", random_state=SEED)
            label = stable_size_order(km.fit_predict(z))
            fits[k] = label
            rows.append({
                "k": k,
                "inertia": float(km.inertia_),
                "silhouette_cosine_sample_n2000": float(silhouette_score(z[sil_idx], label[sil_idx], metric="cosine")),
                "ARI_vs_Leiden28": float(adjusted_rand_score(membership.macro_community, label)),
                "NMI_vs_Leiden28": float(normalized_mutual_info_score(membership.macro_community, label)),
                "min_cluster_n": int(pd.Series(label).value_counts().min()),
                "max_cluster_n": int(pd.Series(label).value_counts().max()),
            })
        scan = pd.DataFrame(rows)
        scan.to_csv(scan_path, sep="\t", index=False)
        labels = fits[TARGET_K]
        membership = membership.copy()
        membership["kmeans28_cluster"] = labels
        membership[["KMeans_PC1", "KMeans_PC2", "KMeans_PC3"]] = coords
    size = membership.kmeans28_cluster.value_counts().sort_index()
    color = {int(k): PALETTE[int(k) % len(PALETTE)] for k in size.index}
    membership["kmeans_color"] = membership.kmeans28_cluster.map(color)
    membership.to_csv(DEST / "kmeans28_all_proteins_source_data.tsv", sep="\t", index=False)

    # Centroid summaries and an exact best-overlap matching for comparison only.
    centers = membership.groupby("kmeans28_cluster")[["KMeans_PC1","KMeans_PC2","KMeans_PC3"]].mean().reset_index()
    centers["n_proteins"] = centers.kmeans28_cluster.map(size)
    contingency = pd.crosstab(membership.macro_community, membership.kmeans28_cluster)
    rr, cc = linear_sum_assignment(-contingency.to_numpy())
    matched_fraction = int(contingency.to_numpy()[rr,cc].sum()) / len(membership)
    centers["envelope_drawn"] = centers.n_proteins.ge(LARGE_N)
    centers.to_csv(DEST / "kmeans28_centroids.tsv", sep="\t", index=False)
    contingency.to_csv(DEST / "leiden28_by_kmeans28_contingency.tsv", sep="\t")

    lo, hi = coords.min(axis=0), coords.max(axis=0); pad = (hi-lo)*.065
    limits = list(zip(lo-pad, hi+pad))
    fig = plt.figure(figsize=(7.01, 5.15))
    grid = fig.add_gridspec(2, 3, width_ratios=[1.26, 1.05, .95], height_ratios=[1, 1],
                            left=.07, right=.985, top=.92, bottom=.10, wspace=.58, hspace=.58)
    ax3d = fig.add_subplot(grid[:, :2], projection="3d")
    axscan = fig.add_subplot(grid[0, 2])
    axheat = fig.add_subplot(grid[1, 2])
    for cluster in sorted(size.index):
        group = membership[membership.kmeans28_cluster.eq(cluster)]
        p = group[["KMeans_PC1","KMeans_PC2","KMeans_PC3"]].to_numpy()
        ax3d.scatter(p[:,0],p[:,1],p[:,2],s=1.8,c=color[int(cluster)],alpha=.12,linewidths=0,depthshade=False,rasterized=True)
        if len(group) >= LARGE_N:
            robust_ellipsoid(ax3d, p, color[int(cluster)])
    for _, row in centers.iterrows():
        major = row.n_proteins >= LARGE_N
        ax3d.scatter([row.KMeans_PC1],[row.KMeans_PC2],[row.KMeans_PC3],s=32 if major else 12,
                     c=color[int(row.kmeans28_cluster)], edgecolors="#252525", linewidths=.42, depthshade=False, zorder=7)
    for _, row in centers.nlargest(8,"n_proteins").iterrows():
        ax3d.text(row.KMeans_PC1,row.KMeans_PC2,row.KMeans_PC3,f" K{int(row.kmeans28_cluster)}",fontsize=6.2,color="#262626")
    ax3d.set_xlim(*limits[0]); ax3d.set_ylim(*limits[1]); ax3d.set_zlim(*limits[2]); ax3d.set_box_aspect([1,1,.76]); ax3d.view_init(elev=23,azim=42)
    ax3d.set_xlabel(f"PC1 ({evr[0]:.1f}%)",labelpad=3,fontsize=7); ax3d.set_ylabel(f"PC2 ({evr[1]:.1f}%)",labelpad=3,fontsize=7); ax3d.set_zlabel(f"PC3 ({evr[2]:.1f}%)",labelpad=3,fontsize=7); ax3d.tick_params(labelsize=6,pad=-1)
    ax3d.xaxis.pane.fill=False; ax3d.yaxis.pane.fill=False; ax3d.zaxis.pane.fill=False; ax3d.grid(False)
    ax3d.set_title("K-means, K = 28",loc="left",fontsize=8,fontweight="bold",pad=3)
    ax3d.text2D(.01,-.115,"All 7,800 proteins; 90% robust core envelopes are drawn only for clusters with n ≥ 150. K-means labels are arbitrary and ordered by size.",transform=ax3d.transAxes,fontsize=6.2,color="#424242",wrap=True)
    ax3d.text2D(-.06,1.03,"a",transform=ax3d.transAxes,fontsize=8,fontweight="bold")

    axscan.plot(scan.k,scan.silhouette_cosine_sample_n2000,"o-",color="#7B95C6",lw=1,ms=3)
    axscan.axvline(TARGET_K,color="#C85E62",lw=.8,ls="--"); axscan.text(TARGET_K+.4,scan.silhouette_cosine_sample_n2000.max(),"K=28",fontsize=5.8,color="#A2444D",va="top")
    axscan.set_xlabel("K",fontsize=7); axscan.set_ylabel("Cosine silhouette\n(2,000-protein sample)",fontsize=6.5); axscan.tick_params(labelsize=6,length=2); axscan.grid(color="#E7E7E7",lw=.45)
    axscan.set_title("K scan",loc="left",fontsize=7,fontweight="bold",pad=3); axscan.text(-.23,1.04,"b",transform=axscan.transAxes,fontsize=8,fontweight="bold")

    # Row-normalized overlap avoids hiding small Leiden modules; it does not
    # assume either partition is ground truth.
    row_norm = contingency.div(contingency.sum(axis=1),axis=0).to_numpy()
    im=axheat.imshow(row_norm,cmap=OVERLAP_CMAP,vmin=0,vmax=1,aspect="auto",interpolation="nearest")
    axheat.set_xlabel("K-means cluster",fontsize=7); axheat.set_ylabel("Leiden community",fontsize=7)
    ticks=[0,7,14,21,27]; axheat.set_xticks(ticks); axheat.set_yticks(ticks); axheat.tick_params(labelsize=5.7,length=2)
    axheat.set_title("Partition overlap",loc="left",fontsize=7,fontweight="bold",pad=3); axheat.text(-.23,1.04,"c",transform=axheat.transAxes,fontsize=8,fontweight="bold")
    cb=fig.colorbar(im,ax=axheat,fraction=.08,pad=.03); cb.set_label("row fraction",fontsize=5.7); cb.ax.tick_params(labelsize=5.6,length=2)
    fig.text(.07,.97,"K-means sensitivity test on the same MemPro 50-dimensional input",fontsize=9,fontweight="bold",ha="left")
    export(fig, DEST / "kmeans28_same_input_comparison")
    plt.close(fig)

    interactive=go.Figure()
    for cluster in sorted(size.index):
        group=membership[membership.kmeans28_cluster.eq(cluster)]
        interactive.add_trace(go.Scatter3d(x=group.KMeans_PC1,y=group.KMeans_PC2,z=group.KMeans_PC3,mode="markers",name=f"K{cluster} (n={len(group):,})",marker=dict(size=2.5,color=color[int(cluster)],opacity=.58),customdata=group[["uniprot_id","macro_community","kmeans28_cluster","posthoc_membrane_role"]],hovertemplate="%{customdata[0]}<br>Leiden: M%{customdata[1]}<br>K-means: K%{customdata[2]}<br>role: %{customdata[3]}<extra></extra>"))
    interactive.update_layout(template="simple_white",width=1240,height=900,title="K-means K=28 sensitivity map — all 7,800 MemPro proteins",scene=dict(xaxis_title=f"PC1 ({evr[0]:.1f}%)",yaxis_title=f"PC2 ({evr[1]:.1f}%)",zaxis_title=f"PC3 ({evr[2]:.1f})"),legend_title_text="K-means cluster")
    interactive.write_html(DEST / "kmeans28_interactive.html",include_plotlyjs=True)

    target_row=scan.loc[scan.k.eq(TARGET_K)].iloc[0].to_dict()
    summary={"method":"KMeans on L2-normalized frozen 50D graph input","n_proteins":7800,"kmeans_n_init":50,"target_k":TARGET_K,"leiden_communities":28,"ARI_vs_Leiden28":target_row["ARI_vs_Leiden28"],"NMI_vs_Leiden28":target_row["NMI_vs_Leiden28"],"best_one_to_one_overlap_fraction":matched_fraction,"silhouette_cosine_sample_n2000":target_row["silhouette_cosine_sample_n2000"],"note":"K-means is a sensitivity analysis, not a replacement for graph-based Leiden communities. Cluster labels are arbitrary."}
    (DEST / "summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    (DEST / "README.md").write_text("# K-means sensitivity analysis\n\nK-means was applied to the same frozen 50-dimensional representation after L2 normalization. It is compared against the existing 28 Leiden communities but does not modify them. The 3D plot is a PCA display only; clustering was performed in 50 dimensions. Every formal protein remains in the source data and interactive map.\n",encoding="utf-8")

if __name__=="__main__":
    main()
