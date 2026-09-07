"""Reproducible K=35 sensitivity map for the frozen 7,800-protein input.

This script preserves the original Leiden communities and runs K-means only as
an alternative partition on the same L2-normalized 50-dimensional input.
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

HERE = Path(__file__).resolve().parent
SEED = 42
K = 35
N_PROTEINS = 7800
LARGE_N = 150
PALETTE = [
    "#7B95C6", "#49C2D9", "#A1D8E8", "#67A583", "#A2C986", "#D0E2C0",
    "#FDED95", "#FFC1A6", "#F59C7C", "#F47254", "#C85E62", "#A67FC5",
    "#B58CCF", "#819CCB", "#89BFA5", "#D89B6E", "#79B5C5", "#B6A8D4",
    "#D7A2A1", "#86A98A", "#D6C480", "#A6B9D7", "#C79D90", "#8EBCD1",
    "#B8CFA2", "#D5A7CC", "#B3B3B3", "#9DA9A0", "#7AA6B8", "#C69AA0",
    "#9AB883", "#D2AD6E", "#9296C4", "#79B9AE", "#C5A1D4",
]
OVERLAP_CMAP = LinearSegmentedColormap.from_list(
    "mempro_overlap", ["#FFF9F5", "#F7B2A0", "#C85E62", "#55223A"]
)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.linewidth": 0.55,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def stable_size_order(labels: np.ndarray) -> np.ndarray:
    """Rename arbitrary K-means labels by descending cluster size."""
    sizes = pd.Series(labels).value_counts().sort_values(ascending=False)
    remap = {old: new for new, old in enumerate(sizes.index)}
    return np.asarray([remap[x] for x in labels], dtype=int)


def draw_robust_envelope(ax, points: np.ndarray, color: str) -> bool:
    """Draw a conservative robust 90% ellipsoid for only large clusters."""
    if len(points) < LARGE_N:
        return False
    try:
        model = MinCovDet(random_state=SEED).fit(points)
        values, vectors = np.linalg.eigh(model.covariance_)
        if np.any(values <= 1e-12) or np.max(values) / np.min(values) > 1e7:
            return False
        radius = np.sqrt(chi2.ppf(0.90, df=3))
        u = np.linspace(0, 2 * np.pi, 28)
        v = np.linspace(0, np.pi, 20)
        sphere = np.vstack([
            np.outer(np.cos(u), np.sin(v)).ravel(),
            np.outer(np.sin(u), np.sin(v)).ravel(),
            np.outer(np.ones_like(u), np.cos(v)).ravel(),
        ])
        xyz = vectors @ np.diag(np.sqrt(values) * radius) @ sphere + model.location_[:, None]
        xx, yy, zz = [array.reshape(len(u), len(v)) for array in xyz]
        ax.plot_surface(xx, yy, zz, color=color, alpha=0.04, linewidth=0,
                        antialiased=False, shade=False)
        return True
    except Exception:
        return False


def export(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")


def main() -> None:
    membership = pd.read_csv(HERE / "protein_leiden_structural_umap3d.tsv", sep="\t")
    input_50d = np.load(HERE / "leiden_input_pca50.npy")
    if len(membership) != N_PROTEINS or len(input_50d) != N_PROTEINS:
        raise RuntimeError("Expected the frozen input to contain exactly 7,800 proteins.")
    if membership["uniprot_id"].nunique() != N_PROTEINS:
        raise RuntimeError("Protein IDs are not unique in the frozen membership table.")
    if membership["macro_community"].nunique() != 28:
        raise RuntimeError("Expected 28 retained Leiden communities for comparison.")

    normalized_50d = normalize(input_50d, norm="l2")
    labels = stable_size_order(KMeans(
        n_clusters=K, n_init=50, max_iter=500, algorithm="lloyd", random_state=SEED
    ).fit_predict(normalized_50d))
    projection = PCA(n_components=3, random_state=SEED).fit(normalized_50d)
    coordinates = projection.transform(normalized_50d)
    explained = projection.explained_variance_ratio_ * 100

    membership = membership.copy()
    membership["kmeans35_cluster"] = labels
    membership[["KMeans_PC1", "KMeans_PC2", "KMeans_PC3"]] = coordinates
    sizes = membership["kmeans35_cluster"].value_counts().sort_index()
    colors = {cluster: PALETTE[cluster % len(PALETTE)] for cluster in sizes.index}
    membership["kmeans_color"] = membership["kmeans35_cluster"].map(colors)
    membership.to_csv(HERE / "kmeans35_all_proteins_source_data.tsv", sep="\t", index=False)

    centers = membership.groupby("kmeans35_cluster")[["KMeans_PC1", "KMeans_PC2", "KMeans_PC3"]].mean().reset_index()
    centers["n_proteins"] = centers["kmeans35_cluster"].map(sizes)
    centers["robust_envelope_drawn"] = centers["n_proteins"].ge(LARGE_N)
    centers.to_csv(HERE / "kmeans35_centroids.tsv", sep="\t", index=False)

    contingency = pd.crosstab(membership["macro_community"], membership["kmeans35_cluster"])
    contingency.to_csv(HERE / "leiden28_by_kmeans35_contingency.tsv", sep="\t")
    rows, cols = linear_sum_assignment(-contingency.to_numpy())
    best_overlap = float(contingency.to_numpy()[rows, cols].sum() / N_PROTEINS)
    rng = np.random.default_rng(SEED)
    sample_index = np.sort(rng.choice(N_PROTEINS, size=2000, replace=False))
    silhouette = float(silhouette_score(
        normalized_50d[sample_index], labels[sample_index], metric="cosine"
    ))
    ari = float(adjusted_rand_score(membership["macro_community"], labels))
    nmi = float(normalized_mutual_info_score(membership["macro_community"], labels))

    lower, upper = coordinates.min(axis=0), coordinates.max(axis=0)
    padding = (upper - lower) * 0.065
    limits = list(zip(lower - padding, upper + padding))
    fig = plt.figure(figsize=(7.10, 5.15))
    grid = fig.add_gridspec(2, 3, width_ratios=[1.28, 1.05, 0.95], height_ratios=[1, 1],
                            left=0.07, right=0.985, top=0.92, bottom=0.10, wspace=0.58, hspace=0.58)
    ax_3d = fig.add_subplot(grid[:, :2], projection="3d")
    ax_sizes = fig.add_subplot(grid[0, 2])
    ax_overlap = fig.add_subplot(grid[1, 2])

    for cluster in sorted(sizes.index):
        group = membership[membership["kmeans35_cluster"].eq(cluster)]
        points = group[["KMeans_PC1", "KMeans_PC2", "KMeans_PC3"]].to_numpy()
        ax_3d.scatter(points[:, 0], points[:, 1], points[:, 2], s=1.8,
                      c=colors[cluster], alpha=0.12, linewidths=0,
                      depthshade=False, rasterized=True)
        draw_robust_envelope(ax_3d, points, colors[cluster])
    for _, row in centers.iterrows():
        major = row["n_proteins"] >= LARGE_N
        ax_3d.scatter([row["KMeans_PC1"]], [row["KMeans_PC2"]], [row["KMeans_PC3"]],
                      s=32 if major else 12, c=colors[int(row["kmeans35_cluster"])],
                      edgecolors="#252525", linewidths=0.42, depthshade=False, zorder=7)
    for _, row in centers.nlargest(8, "n_proteins").iterrows():
        ax_3d.text(row["KMeans_PC1"], row["KMeans_PC2"], row["KMeans_PC3"],
                   f" K{int(row['kmeans35_cluster'])}", fontsize=6.2, color="#262626")
    ax_3d.set_xlim(*limits[0]); ax_3d.set_ylim(*limits[1]); ax_3d.set_zlim(*limits[2])
    ax_3d.set_box_aspect([1, 1, 0.76]); ax_3d.view_init(elev=23, azim=42)
    ax_3d.set_xlabel(f"PC1 ({explained[0]:.1f}%)", labelpad=3, fontsize=7)
    ax_3d.set_ylabel(f"PC2 ({explained[1]:.1f}%)", labelpad=3, fontsize=7)
    ax_3d.set_zlabel(f"PC3 ({explained[2]:.1f}%)", labelpad=3, fontsize=7)
    ax_3d.tick_params(labelsize=6, pad=-1)
    ax_3d.xaxis.pane.fill = False; ax_3d.yaxis.pane.fill = False; ax_3d.zaxis.pane.fill = False; ax_3d.grid(False)
    ax_3d.set_title("K-means, K = 35", loc="left", fontsize=8, fontweight="bold", pad=3)
    ax_3d.text2D(0.01, -0.115, "All 7,800 proteins; 90% robust core envelopes are shown only for clusters with n ≥ 150. Labels are arbitrary and ordered by cluster size.", transform=ax_3d.transAxes, fontsize=6.2, color="#424242", wrap=True)
    ax_3d.text2D(-0.06, 1.03, "a", transform=ax_3d.transAxes, fontsize=8, fontweight="bold")

    ax_sizes.bar(np.arange(K), sizes.to_numpy(), color=[colors[i] for i in sizes.index], width=0.85, linewidth=0)
    ax_sizes.axhline(LARGE_N, color="#6A6A6A", lw=0.7, ls="--")
    ax_sizes.text(K - 0.25, LARGE_N + 12, "envelope threshold", ha="right", va="bottom", fontsize=5.7, color="#5C5C5C")
    ax_sizes.set_xlabel("K-means cluster (size-ranked)", fontsize=7)
    ax_sizes.set_ylabel("proteins", fontsize=7)
    ax_sizes.set_xlim(-0.8, K - 0.2); ax_sizes.set_xticks([0, 8, 17, 26, 34]); ax_sizes.set_xticklabels(["K0", "K8", "K17", "K26", "K34"])
    ax_sizes.tick_params(labelsize=5.8, length=2); ax_sizes.grid(axis="y", color="#E7E7E7", lw=0.45)
    ax_sizes.set_title("Cluster sizes", loc="left", fontsize=7, fontweight="bold", pad=3)
    ax_sizes.text(-0.23, 1.04, "b", transform=ax_sizes.transAxes, fontsize=8, fontweight="bold")

    normalized_rows = contingency.div(contingency.sum(axis=1), axis=0).to_numpy()
    image = ax_overlap.imshow(normalized_rows, cmap=OVERLAP_CMAP, vmin=0, vmax=1,
                              aspect="auto", interpolation="nearest")
    ax_overlap.set_xlabel("K-means cluster", fontsize=7); ax_overlap.set_ylabel("Leiden community", fontsize=7)
    ax_overlap.set_xticks([0, 8, 17, 26, 34]); ax_overlap.set_yticks([0, 7, 14, 21, 27])
    ax_overlap.tick_params(labelsize=5.7, length=2)
    ax_overlap.set_title("Overlap with Leiden", loc="left", fontsize=7, fontweight="bold", pad=3)
    ax_overlap.text(-0.23, 1.04, "c", transform=ax_overlap.transAxes, fontsize=8, fontweight="bold")
    colorbar = fig.colorbar(image, ax=ax_overlap, fraction=0.08, pad=0.03)
    colorbar.set_label("Leiden-row fraction", fontsize=5.7); colorbar.ax.tick_params(labelsize=5.6, length=2)
    fig.text(0.07, 0.97, "K-means sensitivity test on the same MemPro 50-dimensional input", fontsize=9, fontweight="bold", ha="left")
    export(fig, HERE / "kmeans35_same_input_comparison")
    plt.close(fig)

    interactive = go.Figure()
    for cluster in sorted(sizes.index):
        group = membership[membership["kmeans35_cluster"].eq(cluster)]
        interactive.add_trace(go.Scatter3d(
            x=group["KMeans_PC1"], y=group["KMeans_PC2"], z=group["KMeans_PC3"],
            mode="markers", name=f"K{cluster} (n={len(group):,})",
            marker=dict(size=2.5, color=colors[cluster], opacity=0.58),
            customdata=group[["uniprot_id", "macro_community", "kmeans35_cluster", "posthoc_membrane_role"]],
            hovertemplate="%{customdata[0]}<br>Leiden: M%{customdata[1]}<br>K-means: K%{customdata[2]}<br>role: %{customdata[3]}<extra></extra>",
        ))
    interactive.update_layout(
        template="simple_white", width=1240, height=900,
        title="K-means K=35 sensitivity map — all 7,800 MemPro proteins",
        scene=dict(xaxis_title=f"PC1 ({explained[0]:.1f}%)", yaxis_title=f"PC2 ({explained[1]:.1f}%)", zaxis_title=f"PC3 ({explained[2]:.1f}%)"),
        legend_title_text="K-means cluster",
    )
    interactive.write_html(HERE / "kmeans35_interactive.html", include_plotlyjs=True)

    summary = {
        "method": "K-means on L2-normalized frozen 50D graph input",
        "n_proteins": N_PROTEINS,
        "target_k": K,
        "kmeans_n_init": 50,
        "max_iter": 500,
        "leiden_communities": 28,
        "silhouette_cosine_sample_n2000": silhouette,
        "ARI_vs_Leiden28": ari,
        "NMI_vs_Leiden28": nmi,
        "best_one_to_one_overlap_fraction": best_overlap,
        "min_cluster_n": int(sizes.min()),
        "max_cluster_n": int(sizes.max()),
        "note": "K-means is a sensitivity analysis. The 3D PCA display is not the clustering space.",
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (HERE / "README.md").write_text(
        "# K=35 K-means sensitivity analysis\n\n"
        "K-means is applied to the same frozen 50-dimensional protein representation used for the K=28 comparison, after L2 normalization. "
        "All 7,800 formal proteins are included. It is a sensitivity analysis and does not replace Leiden communities.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
