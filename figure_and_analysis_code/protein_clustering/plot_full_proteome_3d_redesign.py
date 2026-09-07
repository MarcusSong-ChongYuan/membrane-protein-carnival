"""Redesigned complete MemPro protein landscape.

This script changes *only* the presentation of the previously frozen,
all-protein 3D PCA projection. It does not re-cluster, filter, or reassign any
protein. Every formal MemPro protein is plotted in every projection.

Visual design: muted all-member point cloud, robust 90% core envelopes only
for large communities (n >= 150), all-community centroids, and two orthogonal
2D projections. The envelopes are descriptive core summaries, not boundaries.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.lines import Line2D
from scipy.stats import chi2
from sklearn.covariance import MinCovDet
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "leiden_structural_results"
DEST = OUT / "full_proteome_3d_redesign"
DEST.mkdir(parents=True, exist_ok=True)

# User-selected restrained MemPro palette. The colours encode community only;
# community membership itself is not a biological annotation used for clustering.
PALETTE = [
    "#7B95C6", "#49C2D9", "#A1D8E8", "#67A583", "#A2C986", "#D0E2C0",
    "#FDED95", "#FFC1A6", "#F59C7C", "#F47254", "#C85E62", "#A67FC5",
    "#B58CCF", "#819CCB", "#89BFA5", "#D89B6E", "#79B5C5", "#B6A8D4",
    "#D7A2A1", "#86A98A", "#D6C480", "#A6B9D7", "#C79D90", "#8EBCD1",
    "#B8CFA2", "#D5A7CC", "#B3B3B3", "#9DA9A0",
]
LARGE_N = 150
SEED = 42

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.linewidth": 0.55,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "figure.facecolor": "white",
})


def robust_ellipsoid_surface(points: np.ndarray, coverage: float = 0.90):
    """Return a 90% MCD covariance ellipsoid, or None if not stable.

    MCD stops isolated projection points from generating the triangular spikes
    associated with ordinary convex hulls. It is not used to omit those points.
    """
    if len(points) < 8:
        return None
    try:
        fit = MinCovDet(random_state=SEED).fit(points)
        eigenvalues, eigenvectors = np.linalg.eigh(fit.covariance_)
        if np.any(eigenvalues <= 1e-12) or np.max(eigenvalues) / np.min(eigenvalues) > 1e7:
            return None
        radius = np.sqrt(chi2.ppf(coverage, df=3))
        u = np.linspace(0, 2 * np.pi, 28)
        v = np.linspace(0, np.pi, 20)
        sphere = np.vstack([
            np.outer(np.cos(u), np.sin(v)).ravel(),
            np.outer(np.sin(u), np.sin(v)).ravel(),
            np.outer(np.ones_like(u), np.cos(v)).ravel(),
        ])
        xyz = eigenvectors @ np.diag(np.sqrt(eigenvalues) * radius) @ sphere
        xyz += fit.location_[:, None]
        return tuple(part.reshape(len(u), len(v)) for part in xyz), fit.location_
    except Exception:
        return None


def add_core_envelope(ax, points: np.ndarray, color: str):
    result = robust_ellipsoid_surface(points)
    if result is None:
        return False
    (xx, yy, zz), _ = result
    ax.plot_surface(xx, yy, zz, color=color, alpha=0.045, linewidth=0,
                    antialiased=False, shade=False, zorder=0)
    return True


def community_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for community, group in frame.groupby("macro_community", sort=True):
        xyz = group[["PC1", "PC2", "PC3"]].to_numpy()
        center = xyz.mean(axis=0)
        rows.append({
            "community": int(community),
            "n_proteins": int(len(group)),
            "PC1_centroid": center[0], "PC2_centroid": center[1], "PC3_centroid": center[2],
            "envelope_drawn": bool(len(group) >= LARGE_N and robust_ellipsoid_surface(xyz) is not None),
        })
    return pd.DataFrame(rows).sort_values("community")


def set_3d_style(ax, limits, evr):
    ax.set_xlim(*limits[0]); ax.set_ylim(*limits[1]); ax.set_zlim(*limits[2])
    ax.set_box_aspect([1.0, 1.0, 0.76])
    ax.view_init(elev=23, azim=42)
    ax.set_xlabel(f"PC1 ({evr[0]:.1f}%)", labelpad=3, fontsize=7)
    ax.set_ylabel(f"PC2 ({evr[1]:.1f}%)", labelpad=3, fontsize=7)
    ax.set_zlabel(f"PC3 ({evr[2]:.1f}%)", labelpad=3, fontsize=7)
    ax.tick_params(labelsize=6, pad=-1)
    ax.xaxis.pane.fill = False; ax.yaxis.pane.fill = False; ax.zaxis.pane.fill = False
    ax.grid(False)


def style_2d(ax, xlabel: str, ylabel: str):
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.tick_params(labelsize=6, length=2, width=.5)
    ax.grid(color="#E7E7E7", linewidth=.45, zorder=0)


def add_panel_label(ax, label: str):
    ax.text2D(-.06, 1.03, label, transform=ax.transAxes, fontsize=8, fontweight="bold",
              ha="left", va="bottom") if hasattr(ax, "text2D") else ax.text(-.18, 1.04, label,
              transform=ax.transAxes, fontsize=8, fontweight="bold", ha="left", va="bottom")


def export(fig, stem: Path):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")


def main():
    source = OUT / "protein_leiden_structural_umap3d.tsv"
    graph_input = OUT / "leiden_input_pca50.npy"
    frame = pd.read_csv(source, sep="\t")
    x = np.load(graph_input)
    if len(frame) != 7800 or len(x) != 7800:
        raise RuntimeError(f"Expected 7,800 proteins in both inputs; observed rows={len(frame)}, features={len(x)}")
    if frame["uniprot_id"].nunique() != 7800:
        raise RuntimeError("Protein identity is not unique.")
    if frame["macro_community"].nunique() != 28:
        raise RuntimeError("Expected 28 pre-existing Leiden communities.")

    pca = PCA(n_components=3, random_state=SEED).fit(x)
    frame[["PC1", "PC2", "PC3"]] = pca.transform(x)
    evr = pca.explained_variance_ratio_ * 100
    communities = sorted(frame["macro_community"].unique())
    colors = {community: PALETTE[i] for i, community in enumerate(communities)}
    frame["community_color"] = frame["macro_community"].map(colors)
    centers = community_table(frame)
    centers["community_color"] = centers["community"].map(colors)
    large = centers.loc[centers.n_proteins >= LARGE_N].copy()

    # Source data are a full point-level copy plus the rendered summary geometry.
    frame.to_csv(DEST / "full_proteome_3d_redesign_source_data.tsv", sep="\t", index=False)
    centers.to_csv(DEST / "community_centroids_and_envelope_status.tsv", sep="\t", index=False)

    points = frame[["PC1", "PC2", "PC3"]].to_numpy()
    lo, hi = points.min(axis=0), points.max(axis=0)
    pad = (hi - lo) * .065
    limits = list(zip(lo - pad, hi + pad))

    # Asymmetric mixed-modality figure: complete 3D hero plus two transparent
    # orthogonal 2D projections. Points are never sampled or removed.
    fig = plt.figure(figsize=(7.01, 5.15))  # 178 mm manuscript width
    grid = fig.add_gridspec(2, 3, width_ratios=[1.12, 1.12, .93], height_ratios=[1, 1],
                            left=.065, right=.985, top=.925, bottom=.10, wspace=.52, hspace=.52)
    hero = fig.add_subplot(grid[:, :2], projection="3d")
    panel_b = fig.add_subplot(grid[0, 2])
    panel_c = fig.add_subplot(grid[1, 2])

    # Every protein remains visible. Low opacity makes density readable rather
    # than turning the whole landscape into opaque coloured blocks.
    for community in communities:
        subset = frame.loc[frame.macro_community.eq(community), ["PC1", "PC2", "PC3"]].to_numpy()
        hero.scatter(subset[:, 0], subset[:, 1], subset[:, 2], s=1.8, c=colors[community],
                     alpha=.115, linewidths=0, depthshade=False, rasterized=True)
        if len(subset) >= LARGE_N:
            add_core_envelope(hero, subset, colors[community])
    labelled_communities = set(centers.nlargest(8, "n_proteins").community.astype(int))
    for _, row in centers.iterrows():
        is_large = row.n_proteins >= LARGE_N
        hero.scatter([row.PC1_centroid], [row.PC2_centroid], [row.PC3_centroid],
                     s=32 if is_large else 12, c=row.community_color,
                     edgecolors="#262626", linewidths=.42, marker="o", depthshade=False, zorder=8)
        if int(row.community) in labelled_communities:
            hero.text(row.PC1_centroid, row.PC2_centroid, row.PC3_centroid, f" M{int(row.community)}",
                      fontsize=6.2, color="#252525", zorder=9)
    set_3d_style(hero, limits, evr)
    hero.set_title("Complete protein landscape", fontsize=8, loc="left", pad=3, fontweight="bold")
    hero.text2D(.01, -.115, "All 7,800 proteins are shown. Soft envelopes summarize only the 90% robust core of communities with n ≥ 150; they are not membership boundaries.",
                transform=hero.transAxes, fontsize=6.2, color="#424242", wrap=True)
    legend = [
        Line2D([0], [0], marker="o", color="w", label="protein", markerfacecolor="#8D99A8", markeredgecolor="none", markersize=4),
        Line2D([0], [0], marker="o", color="w", label="community centroid", markerfacecolor="#FFFFFF", markeredgecolor="#262626", markersize=5),
        Line2D([0], [0], color="#798A9D", lw=3, alpha=.35, label="robust core envelope (n ≥ 150)"),
    ]
    hero.legend(handles=legend, loc="upper left", bbox_to_anchor=(.01, .98), fontsize=5.8,
                handlelength=1.1, labelspacing=.3, borderpad=.18)
    add_panel_label(hero, "a")

    for community in communities:
        subset = frame.loc[frame.macro_community.eq(community)]
        panel_b.scatter(subset.PC1, subset.PC2, s=1.45, c=colors[community], alpha=.115,
                        linewidths=0, rasterized=True)
        panel_c.scatter(subset.PC1, subset.PC3, s=1.45, c=colors[community], alpha=.115,
                        linewidths=0, rasterized=True)
    panel_b.scatter(centers.PC1_centroid, centers.PC2_centroid, s=10, c=centers.community_color,
                    edgecolor="#2A2A2A", linewidth=.25, zorder=3)
    panel_c.scatter(centers.PC1_centroid, centers.PC3_centroid, s=10, c=centers.community_color,
                    edgecolor="#2A2A2A", linewidth=.25, zorder=3)
    style_2d(panel_b, f"PC1 ({evr[0]:.1f}%)", f"PC2 ({evr[1]:.1f}%)")
    style_2d(panel_c, f"PC1 ({evr[0]:.1f}%)", f"PC3 ({evr[2]:.1f}%)")
    panel_b.set_title("PC1–PC2 projection", fontsize=7, loc="left", pad=3, fontweight="bold")
    panel_c.set_title("PC1–PC3 projection", fontsize=7, loc="left", pad=3, fontweight="bold")
    add_panel_label(panel_b, "b"); add_panel_label(panel_c, "c")

    fig.text(.065, .973, "3D PCA landscape of 7,800 curated membrane proteins across 28 Leiden communities",
             fontsize=9, fontweight="bold", ha="left")
    export(fig, DEST / "full_proteome_leiden_3d_redesign")
    plt.close(fig)

    # Interactive inspection retains full membership and direct hover labels.
    interactive = go.Figure()
    for community in communities:
        subset = frame.loc[frame.macro_community.eq(community)]
        interactive.add_trace(go.Scatter3d(
            x=subset.PC1, y=subset.PC2, z=subset.PC3, mode="markers",
            name=f"M{community} (n={len(subset):,})",
            marker=dict(size=2.5, color=colors[community], opacity=.58),
            customdata=subset[["uniprot_id", "membrane_class_v7", "primary_membrane_mode_v7", "fine_community"]],
            hovertemplate="%{customdata[0]}<br>membrane class: %{customdata[1]}<br>primary mode: %{customdata[2]}<br>fine community: %{customdata[3]}<extra></extra>",
        ))
    interactive.update_layout(
        template="simple_white", width=1240, height=900,
        title="Complete MemPro protein landscape — all 7,800 formal proteins",
        scene=dict(xaxis_title=f"PC1 ({evr[0]:.1f}%)", yaxis_title=f"PC2 ({evr[1]:.1f}%)", zaxis_title=f"PC3 ({evr[2]:.1f}%)"),
        legend_title_text="Pre-existing Leiden community",
    )
    interactive.write_html(DEST / "full_proteome_leiden_3d_redesign_interactive.html", include_plotlyjs=True)

    contract = """# Figure contract — complete MemPro protein landscape

**Core conclusion.** All 7,800 formal MemPro proteins can be displayed in a common low-dimensional structural-feature landscape; the visual summaries show community cores without treating projected envelopes as biological boundaries.

**Archetype.** Asymmetric mixed-modality quantitative figure: one 3D hero panel plus two orthogonal 2D projections.

**Evidence hierarchy.** The coordinate input is the frozen 50-dimensional graph representation (sequence embedding, topology features and structural-family representation). The 28 Leiden communities are pre-existing assignments. Membrane role, disease and ligand annotations were held out from the graph construction.

**Display rule.** Every panel contains all 7,800 proteins. Community-specific robust MCD 90% envelopes are drawn only for communities with n ≥ 150. All 28 centroids are plotted; only large-community centroids are labelled. The envelopes are visual core summaries, not membership filters or strict separations.

**Reviewer boundary.** The first three PCs explain only a limited fraction of the 50-dimensional variance, so this plot is a landscape and navigation view rather than evidence of strict inter-community distances. Quantitative centroid distances and graph topology are supplied separately.
"""
    (DEST / "FIGURE_CONTRACT.md").write_text(contract, encoding="utf-8")
    readme = """# Redesigned complete protein 3D landscape

This revision preserves the prior full dataset and community assignments. It redraws the all-protein landscape with reduced occlusion: low-opacity points for every protein, 90% MCD core envelopes for large communities only, all-community centroids and two orthogonal 2D views. It intentionally does not use convex hulls because one projected outlying member can create a misleading cone-like hull.

Use the interactive HTML to inspect individual protein identities. Use the accompanying topology, centroid-distance and post-hoc enrichment packages for quantitative interpretation; this 3D display alone is not a cluster-separation test.
"""
    (DEST / "README.md").write_text(readme, encoding="utf-8")
    manifest = {
        "n_proteins": int(len(frame)),
        "n_communities": int(len(communities)),
        "all_members_plotted_in_every_static_panel": True,
        "community_assignment": "pre-existing Leiden assignment; not modified",
        "coordinate_method": "PCA projection of frozen 50-dimensional graph input",
        "pca_explained_variance_percent": [float(x) for x in evr],
        "core_envelope": "90% minimum covariance determinant ellipsoid",
        "large_community_threshold": LARGE_N,
        "n_envelopes": int(centers.envelope_drawn.sum()),
        "seed": SEED,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "graph_input_sha256": hashlib.sha256(graph_input.read_bytes()).hexdigest(),
    }
    (DEST / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
