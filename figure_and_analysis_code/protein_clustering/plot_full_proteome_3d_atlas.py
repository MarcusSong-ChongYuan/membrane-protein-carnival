"""Complete 7,800-protein 3D Leiden atlas.

This is deliberately not a selected-community graphic.  Every formal MemPro
protein is included exactly once.  A global panel aggregates complete
communities with translucent convex hulls, while the 28 small multiples retain
the individual protein point clouds without the occlusion of one large panel.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.covariance import MinCovDet
from scipy.stats import chi2
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "leiden_structural_results"
DEST = OUT / "full_proteome_3d_atlas"
DEST.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
})

PALETTE = [
    "#7B95C6", "#49C2D9", "#A1D8E8", "#67A583", "#A2C986", "#D0E2C0",
    "#FDED95", "#FFC1A6", "#F59C7C", "#F47254", "#C85E62", "#A67FC5",
    "#B58CCF", "#819CCB", "#89BFA5", "#D89B6E", "#79B5C5", "#B6A8D4",
    "#D7A2A1", "#86A98A", "#D6C480", "#A6B9D7", "#C79D90", "#8EBCD1",
    "#B8CFA2", "#D5A7CC", "#B3B3B3", "#9DA9A0",
]


def add_robust_ellipsoid(ax, points, color, alpha=.05, coverage=.90):
    """Robust 3D core envelope.  It summarizes, but never excludes, point marks.

    Ordinary convex hulls are dominated by isolated high-leverage points in this
    embedding and form misleading triangular spikes.  The 90% MCD ellipsoid is
    therefore used only as a visual envelope; every member remains plotted.
    """
    if len(points) < 8:
        return
    try:
        mcd = MinCovDet(random_state=42).fit(points)
        eigval, eigvec = np.linalg.eigh(mcd.covariance_)
        if np.any(eigval <= 1e-12):
            return
        radius = np.sqrt(chi2.ppf(coverage, df=3))
        u = np.linspace(0, 2*np.pi, 26)
        v = np.linspace(0, np.pi, 18)
        sphere = np.array([
            np.outer(np.cos(u), np.sin(v)),
            np.outer(np.sin(u), np.sin(v)),
            np.outer(np.ones_like(u), np.cos(v)),
        ]).reshape(3, -1)
        transformed = (eigvec @ np.diag(np.sqrt(eigval) * radius) @ sphere) + mcd.location_[:, None]
        xx, yy, zz = [z.reshape(len(u), len(v)) for z in transformed]
        ax.plot_surface(xx, yy, zz, color=color, alpha=alpha, linewidth=0,
                        antialiased=False, shade=False)
    except Exception:
        return


def add_static_cloud(ax, points, color, marker_size=5, alpha=0.40, envelope_alpha=0.035):
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=marker_size,
               c=color, alpha=alpha, linewidths=0, depthshade=False)
    add_robust_ellipsoid(ax, points, color, alpha=envelope_alpha)


def style_3d(ax, limits, evr, show_ticks=False):
    ax.set_xlim(*limits[0]); ax.set_ylim(*limits[1]); ax.set_zlim(*limits[2])
    ax.view_init(elev=21, azim=42)
    ax.set_box_aspect([1, 1, .72])
    if show_ticks:
        ax.set_xlabel(f"PC1 ({evr[0]:.1f}%)", labelpad=5)
        ax.set_ylabel(f"PC2 ({evr[1]:.1f}%)", labelpad=5)
        ax.set_zlabel(f"PC3 ({evr[2]:.1f}%)", labelpad=5)
        ax.tick_params(labelsize=6, pad=0)
    else:
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.xaxis.pane.fill = False; ax.yaxis.pane.fill = False; ax.zaxis.pane.fill = False
    ax.grid(False)


def main():
    d = pd.read_csv(OUT / "protein_leiden_structural_umap3d.tsv", sep="\t")
    x = np.load(OUT / "leiden_input_pca50.npy")
    if len(d) != len(x):
        raise RuntimeError(f"Embedding/data mismatch: {len(d)} vs {len(x)}")
    pca = PCA(n_components=3, random_state=42).fit(x)
    d[["PC1", "PC2", "PC3"]] = pca.transform(x)
    evr = pca.explained_variance_ratio_ * 100
    comms = sorted(d.macro_community.unique())
    if len(d) != 7800 or len(comms) != 28:
        raise RuntimeError(f"Expected 7,800 proteins in 28 communities; got {len(d)} in {len(comms)}")

    labels = pd.read_csv(OUT / "leiden_macro_module_candidate_labels.tsv", sep="\t")
    label_map = (labels.drop_duplicates("community").set_index("community")
                 .candidate_label.to_dict())
    color_map = {c: PALETTE[i] for i, c in enumerate(comms)}
    d["community_color"] = d.macro_community.map(color_map)
    d["candidate_label"] = d.macro_community.map(label_map)
    d["candidate_label"] = d["candidate_label"].fillna("no independent candidate label")
    d.loc[d.candidate_label.astype(str).str.lower().eq("nan"), "candidate_label"] = "no independent candidate label"
    d.to_csv(DEST / "full_proteome_pca3d_source_data.tsv", sep="\t", index=False)

    pts = d[["PC1", "PC2", "PC3"]].to_numpy()
    mins, maxs = pts.min(0), pts.max(0)
    padding = (maxs - mins) * .06
    limits = list(zip(mins-padding, maxs+padding))

    # Global: all 7,800 proteins, with point marks deliberately faint. Robust
    # ellipsoids represent each community's 90% core, avoiding hull spikes from
    # isolated proteins without removing them from the point layer.
    fig = plt.figure(figsize=(8.0, 6.1))
    ax = fig.add_subplot(projection="3d")
    for c in comms:
        g = d[d.macro_community.eq(c)]
        gpts = g[["PC1", "PC2", "PC3"]].to_numpy()
        add_static_cloud(ax, gpts, color_map[c], marker_size=3.2, alpha=.16, envelope_alpha=.028)
    style_3d(ax, limits, evr, show_ticks=True)
    ax.set_title("Full MemPro protein landscape: 7,800 proteins in 28 structural communities", fontsize=9, pad=12)
    ax.text2D(.02, .02, "Each point = one protein; translucent envelope = robust 90% community core.\nAll 7,800 members, including peripheral points, remain visible in the point layer and interactive view.", transform=ax.transAxes, fontsize=6.5, color="#424242")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(DEST / f"full_proteome_leiden_global_3d.{ext}", dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Full atlas: no downsampling. Each panel shares global limits, so shape,
    # scale and separation are comparable while all proteins remain visible.
    fig = plt.figure(figsize=(13.8, 10.4))
    for i, c in enumerate(comms, 1):
        ax = fig.add_subplot(4, 7, i, projection="3d")
        g = d[d.macro_community.eq(c)]
        gpts = g[["PC1", "PC2", "PC3"]].to_numpy()
        alpha = .16 if len(gpts) > 500 else .48
        msize = 1.2 if len(gpts) > 500 else 5.5
        add_static_cloud(ax, gpts, color_map[c], marker_size=msize, alpha=alpha, envelope_alpha=.075)
        style_3d(ax, limits, evr, show_ticks=False)
        short_label = str(d.loc[d.macro_community.eq(c), "candidate_label"].iloc[0])
        if len(short_label) > 32:
            short_label = short_label[:29] + "..."
        ax.set_title(f"M{c}  n={len(g):,}\n{short_label}", fontsize=6.0, pad=1.5, loc="left")
    fig.text(.5, .98, "Complete MemPro 3D community atlas (all 7,800 formal membrane proteins)", ha="center", va="top", fontsize=10)
    fig.text(.5, .012, f"Shared coordinates: PC1 {evr[0]:.1f}% · PC2 {evr[1]:.1f}% · PC3 {evr[2]:.1f}%. Each panel displays every protein assigned to that community; no sampling.", ha="center", fontsize=7, color="#424242")
    fig.subplots_adjust(left=.015, right=.99, bottom=.035, top=.945, wspace=.00, hspace=.11)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(DEST / f"full_proteome_leiden_community_atlas.{ext}", dpi=400, facecolor="white")
    plt.close(fig)

    # Interactive full dataset: each community is a separate trace, allowing
    # viewers to toggle a community but never replacing full-data defaults.
    interactive = go.Figure()
    for c in comms:
        g = d[d.macro_community.eq(c)]
        interactive.add_trace(go.Scatter3d(
            x=g.PC1, y=g.PC2, z=g.PC3, mode="markers", name=f"M{c} (n={len(g)})",
            marker=dict(size=2.7, color=color_map[c], opacity=.55),
            customdata=g[["uniprot_id", "candidate_label", "posthoc_membrane_role", "fine_community"]],
            hovertemplate="%{customdata[0]}<br>candidate label: %{customdata[1]}<br>role: %{customdata[2]}<br>fine community: %{customdata[3]}<extra></extra>"))
    interactive.update_layout(
        template="simple_white", width=1200, height=850,
        title="Full MemPro protein landscape — all 7,800 proteins",
        scene=dict(xaxis_title=f"PC1 ({evr[0]:.1f}%)", yaxis_title=f"PC2 ({evr[1]:.1f}%)", zaxis_title=f"PC3 ({evr[2]:.1f}%)"),
        legend_title_text="Leiden community (click to filter)")
    interactive.write_html(DEST / "full_proteome_leiden_interactive.html", include_plotlyjs=True)

    manifest = {
        "scope": "all formal MemPro proteins; no sampling or exclusion",
        "n_proteins": int(len(d)), "n_communities": int(len(comms)),
        "coordinate_method": "PCA 3D projection of 50-dimensional graph input representation",
        "graph_input": "ESM sequence embedding + topology features + structural-family annotations; roles, disease and ligand annotations excluded",
        "community_method": "Leiden community assignment on a k=25 cosine-neighbor graph",
        "interactive_default": "all 7,800 points visible; individual communities can be toggled",
        "static_envelope": "robust 90% MCD covariance ellipsoid; it is a visual core summary, not a membership filter",
        "static_atlas": "28 shared-scale small multiples; each protein displayed once in its assigned community",
        "pca_explained_variance_percent": evr.tolist(),
    }
    (DEST / "README.md").write_text(
        "# Complete MemPro 3D protein community atlas\n\n"
        "This package contains the complete formal protein set (7,800 proteins), not a selected subset. "
        "The global view uses low-opacity point marks and robust 90% covariance envelopes to control occlusion without discarding outlying members. The 28-panel atlas displays all proteins once, grouped by Leiden community, while preserving shared PCA coordinate limits. "
        "The interactive view starts with all proteins visible and provides hover/filter interaction. Candidate labels are independent post-hoc annotations, not inputs to clustering.\n",
        encoding="utf-8")
    (DEST / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
