"""Community-topology exploration for the complete 7,800-protein Leiden map.

The graph and centroid calculations use the same role-blind 50-D input that
was used for Leiden.  This script does not re-cluster proteins.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import cosine_distances
from sklearn.neighbors import NearestNeighbors
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import seaborn as sns

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "leiden_structural_results"
DEST = OUT / "community_topology_exploration"
DEST.mkdir(exist_ok=True)
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.top": False, "axes.spines.right": False,
})


def reconstruct_weighted_knn(z, k=25):
    """Exact weighted-edge rule used by the Leiden preprocessing script."""
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine", n_jobs=-1).fit(z)
    distances, indices = nn.kneighbors(z)
    nonzero = distances[:, 1:][distances[:, 1:] > 0]
    scale = float(np.median(nonzero)) if len(nonzero) else 1.0
    weights = {}
    for i in range(len(z)):
        for j, dist in zip(indices[i, 1:], distances[i, 1:]):
            a, b = sorted((int(i), int(j)))
            w = float(np.exp(-((float(dist) / scale) ** 2)))
            weights[(a, b)] = max(weights.get((a, b), 0.0), w)
    return weights, scale


def save_all(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(DEST / f"{stem}.{ext}", dpi=450, bbox_inches="tight", facecolor="white")


def main():
    data = pd.read_csv(OUT / "protein_leiden_structural_umap3d.tsv", sep="\t")
    z = np.load(OUT / "leiden_input_pca50.npy")
    if len(data) != 7800 or len(z) != len(data):
        raise RuntimeError("Expected the complete 7,800-protein embedding")
    communities = np.array(data.macro_community, dtype=int)
    community_ids = np.sort(np.unique(communities))
    if len(community_ids) != 28:
        raise RuntimeError("Expected 28 Leiden communities")
    sizes = pd.Series(communities).value_counts().sort_index()

    candidate = pd.read_csv(OUT / "leiden_macro_module_candidate_labels.tsv", sep="\t")
    label_map = candidate.drop_duplicates("community").set_index("community").candidate_label.to_dict()
    label_map = {int(k): ("unresolved" if pd.isna(v) or str(v).lower() == "nan" else str(v)) for k, v in label_map.items()}

    edges, distance_scale = reconstruct_weighted_knn(z)
    degree = np.zeros(len(z), dtype=float)
    observed = np.zeros((len(community_ids), len(community_ids)), dtype=float)
    edge_counts = np.zeros_like(observed, dtype=int)
    local_index = {c: i for i, c in enumerate(community_ids)}
    for (i, j), w in edges.items():
        degree[i] += w; degree[j] += w
        ci, cj = local_index[communities[i]], local_index[communities[j]]
        if ci != cj:
            observed[ci, cj] += w; observed[cj, ci] += w
            edge_counts[ci, cj] += 1; edge_counts[cj, ci] += 1
    degree_sum = np.array([degree[communities == c].sum() for c in community_ids])
    total_degree = degree.sum()
    expected = np.outer(degree_sum, degree_sum) / total_degree
    np.fill_diagonal(expected, np.nan)
    ratio = observed / expected
    np.fill_diagonal(ratio, np.nan)

    # Community centroids in the original 50-D map; MDS preserves these
    # centroid distances more directly than UMAP does.
    centroids = np.vstack([z[communities == c].mean(axis=0) for c in community_ids])
    cdist = cosine_distances(centroids)
    np.fill_diagonal(cdist, 0)
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42,
              n_init=12, max_iter=1000, normalized_stress="auto")
    xy = mds.fit_transform(cdist)
    linkage_matrix = linkage(squareform(cdist, checks=False), method="average", optimal_ordering=True)

    node = pd.DataFrame({
        "community": community_ids, "protein_count": [int(sizes.loc[c]) for c in community_ids],
        "candidate_label": [label_map.get(int(c), "unresolved") for c in community_ids],
        "centroid_mds1": xy[:, 0], "centroid_mds2": xy[:, 1],
        "weighted_knn_degree": degree_sum,
    })
    edge_rows = []
    for ai, c1 in enumerate(community_ids):
        for bi, c2 in enumerate(community_ids):
            if bi <= ai or observed[ai, bi] <= 0:
                continue
            edge_rows.append({
                "community_a": int(c1), "community_b": int(c2),
                "cross_knn_edge_count": int(edge_counts[ai, bi]),
                "cross_knn_weight": float(observed[ai, bi]),
                "configuration_expected_weight": float(expected[ai, bi]),
                "normalized_cross_community_connectivity": float(ratio[ai, bi]),
                "centroid_cosine_distance": float(cdist[ai, bi]),
            })
    edge_df = pd.DataFrame(edge_rows)
    # PAGA-like visual backbone: union of each node's two strongest
    # normalized cross-community connections, not every possible pair.
    keep = set()
    for c in community_ids:
        subset = edge_df[(edge_df.community_a == c) | (edge_df.community_b == c)]
        for _, row in subset.nlargest(2, "normalized_cross_community_connectivity").iterrows():
            keep.add((int(row.community_a), int(row.community_b)))
    edge_df["shown_in_topology_backbone"] = [
        (int(a), int(b)) in keep for a, b in zip(edge_df.community_a, edge_df.community_b)
    ]
    node.to_csv(DEST / "community_nodes.tsv", sep="\t", index=False)
    edge_df.to_csv(DEST / "community_edges.tsv", sep="\t", index=False)
    pd.DataFrame(cdist, index=[f"M{c}" for c in community_ids], columns=[f"M{c}" for c in community_ids]).to_csv(DEST / "community_centroid_cosine_distance.tsv", sep="\t")

    # A. PAGA-like topology. Force-directed positions keep the topological
    # backbone readable; the separate MDS panel is used to read metric distance.
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    backbone = edge_df[edge_df.shown_in_topology_backbone].copy()
    graph = nx.Graph()
    graph.add_nodes_from(map(int, community_ids))
    for row in backbone.itertuples(index=False):
        graph.add_edge(int(row.community_a), int(row.community_b),
                       weight=float(np.log1p(row.normalized_cross_community_connectivity)))
    initial = {int(c): xy[i] for i, c in enumerate(community_ids)}
    pos = nx.spring_layout(graph, pos=initial, seed=42, weight="weight", iterations=600, k=1.10)
    vals = np.log2(np.maximum(backbone.normalized_cross_community_connectivity.to_numpy(), 1e-8))
    norm = plt.Normalize(np.quantile(vals, .10), np.quantile(vals, .90), clip=True)
    for row, value in zip(backbone.itertuples(index=False), vals):
        a, b = pos[row.community_a], pos[row.community_b]
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#728096", alpha=.18 + .48*norm(value),
                lw=.45 + 1.65*norm(value), zorder=1)
    node_sizes = 24 + 9.5 * np.sqrt(node.protein_count.to_numpy())
    node_xy = np.vstack([pos[int(c)] for c in community_ids])
    ax.scatter(node_xy[:, 0], node_xy[:, 1], s=node_sizes, c="#6B8FBF", alpha=.88,
               edgecolor="white", linewidth=.8, zorder=2)
    for row in node.itertuples(index=False):
        xy0 = pos[int(row.community)]
        ax.text(xy0[0], xy0[1], f"M{row.community}", ha="center", va="center",
                fontsize=6.5 if row.protein_count < 1000 else 7.5, color="white", weight="bold", zorder=3,
                clip_on=False)
    ax.set_title("PAGA-like topology of 28 Leiden communities", loc="left", fontsize=10, weight="bold")
    ax.text(.0, -0.10, "Node area: protein count. Backbone: union of each community's two strongest normalized cross-community kNN links.\nForce-directed position is topological only; use the MDS panel for centroid distance.", transform=ax.transAxes, fontsize=7, color="#424242", va="top")
    pad = .10 * np.ptp(node_xy, axis=0)
    ax.set_xlim(node_xy[:, 0].min() - pad[0], node_xy[:, 0].max() + pad[0])
    ax.set_ylim(node_xy[:, 1].min() - pad[1], node_xy[:, 1].max() + pad[1])
    ax.set_aspect("equal", adjustable="datalim"); ax.set_axis_off()
    save_all(fig, "paga_like_community_topology")
    plt.close(fig)

    # B. Metric centroid map: a network-free view so spatial distance can be read.
    fig, ax = plt.subplots(figsize=(7.5, 5.6))
    sc = ax.scatter(node.centroid_mds1, node.centroid_mds2, s=node_sizes, c=node.protein_count,
                    cmap="Blues", edgecolor="white", linewidth=.8)
    for row in node.itertuples(index=False):
        ax.annotate(f"M{row.community}", (row.centroid_mds1, row.centroid_mds2), ha="center", va="center", fontsize=6.5, color="#1E2A3A", weight="bold")
    cb = fig.colorbar(sc, ax=ax, shrink=.78, pad=.02); cb.set_label("Proteins per community")
    ax.set_title("Metric MDS of inter-community ESM/topology/family distance", loc="left", fontsize=10, weight="bold")
    ax.text(.0, -0.12, "Each node is a community centroid. Euclidean distances in this panel approximate high-dimensional centroid cosine distances; no graph edges are drawn.", transform=ax.transAxes, fontsize=7, color="#424242", va="top")
    ax.set_xlabel("MDS1"); ax.set_ylabel("MDS2"); ax.set_aspect("equal", adjustable="datalim")
    save_all(fig, "centroid_metric_mds")
    plt.close(fig)

    # C. Quantitative distance heatmap with hierarchical ordering.
    labels = [f"M{c}" for c in community_ids]
    dist_df = pd.DataFrame(cdist, index=labels, columns=labels)
    grid = sns.clustermap(dist_df, row_linkage=linkage_matrix, col_linkage=linkage_matrix,
                          cmap="rocket", linewidths=.12, linecolor="white", figsize=(8.0, 8.0),
                          cbar_kws={"label": "Centroid cosine distance"}, xticklabels=True, yticklabels=True)
    grid.ax_heatmap.tick_params(labelsize=6.5, length=0)
    grid.ax_heatmap.set_xlabel("Leiden community"); grid.ax_heatmap.set_ylabel("Leiden community")
    grid.fig.suptitle("Inter-community distance in the role-blind input space", x=.08, y=.995, ha="left", fontsize=10, weight="bold")
    grid.fig.text(.08, .012, "Average-linkage dendrogram of community-centroid cosine distance. Darker cells represent closer centroids.", fontsize=7, color="#424242")
    for ext in ("png", "pdf", "svg"):
        grid.fig.savefig(DEST / f"community_centroid_distance_heatmap.{ext}", dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(grid.fig)

    manifest = {
        "scope": "all 7800 formal MemPro proteins assigned to 28 pre-existing Leiden communities",
        "no_reclustering": True,
        "graph_rule": "same 25-nearest-neighbour cosine graph and Gaussian distance weighting used for Leiden",
        "topology_edge_rule": "cross-community weighted kNN connectivity divided by a configuration-model expectation; display backbone is union of top two edges per community",
        "metric_position_rule": "MDS of community-centroid cosine distances in the pre-existing 50-D graph input",
        "held_out_annotations": "candidate labels are post-hoc; primary membrane role, disease, and ligand data were not graph inputs",
        "knn_distance_scale": distance_scale,
        "mds_stress": float(mds.stress_),
    }
    (DEST / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (DEST / "README.md").write_text(
        "# Community-topology exploration\n\n"
        "This package gives complementary visual encodings of the same complete 7,800-protein Leiden map. It does not perform a new clustering. "
        "The PAGA-like topology uses observed inter-community kNN links, normalized against a degree-based configuration expectation; it is not a phylogenetic tree. "
        "The MDS panel is the appropriate panel for reading centroid distance. The heatmap displays that distance numerically.\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
