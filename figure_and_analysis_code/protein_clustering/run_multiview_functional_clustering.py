#!/usr/bin/env python3
"""Role-agnostic multi-view clustering of MemPro formal proteins.

This is deliberately distinct from the ESM-only analysis.  It combines five
*pre-existing* evidence blocks, but excludes primary membrane role from model
inputs so the role labels can remain a post-hoc validation overlay.
"""
from __future__ import annotations

import json
from pathlib import Path

import hdbscan
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import umap
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
BASE = ROOT / "results"
OUT = ROOT / "multiview_results"
OUT.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 8,
    "axes.spines.top": False, "axes.spines.right": False,
})


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def sparse_block(rows: list[dict[str, float]], components: int, name: str) -> tuple[np.ndarray, dict]:
    """TF-IDF + SVD for sparse categorical/profile evidence."""
    vec = DictVectorizer(sparse=True)
    x = vec.fit_transform(rows)
    if x.shape[1] < 2:
        return np.zeros((len(rows), 1), dtype=np.float32), {"name": name, "n_features": int(x.shape[1]), "components": 1}
    x = TfidfTransformer(norm="l2", sublinear_tf=True).fit_transform(x)
    n_comp = max(2, min(components, x.shape[0] - 1, x.shape[1] - 1))
    z = TruncatedSVD(n_components=n_comp, algorithm="randomized", n_iter=7, random_state=42).fit_transform(x)
    return z.astype(np.float32), {"name": name, "n_features": int(x.shape[1]), "components": int(n_comp)}


def equal_block_weight(x: np.ndarray) -> np.ndarray:
    """Give each evidence block unit expected total variance after standardising."""
    z = StandardScaler().fit_transform(x)
    return (z / np.sqrt(max(1, z.shape[1]))).astype(np.float32)


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=400, bbox_inches="tight", facecolor="white")


def main() -> None:
    master = read_tsv(INPUT / "protein_master_v72.tsv.gz", usecols=[
        "canonical_uniprot_accession", "sequence_length", "membrane_class_v7",
        "primary_membrane_mode_v7", "secondary_membrane_modes_v7",
    ]).rename(columns={"canonical_uniprot_accession": "uniprot_id"})
    master = master.drop_duplicates("uniprot_id").set_index("uniprot_id")
    ids = read_tsv(BASE / "protein_embedding_ids.tsv")["uniprot_id"].astype(str).tolist()
    assert len(ids) == 7800 and set(ids) == set(master.index), "Formal protein/embedding mismatch"
    master = master.loc[ids].reset_index()

    role = read_tsv(INPUT / "protein_membrane_role_FORMAL.tsv", dtype=str).rename(columns={
        "canonical_uniprot_accession": "uniprot_id", "formal_primary_membrane_role": "posthoc_membrane_role",
    })[["uniprot_id", "posthoc_membrane_role"]]
    master = master.merge(role, on="uniprot_id", how="left")
    master["posthoc_membrane_role"] = master["posthoc_membrane_role"].fillna("Unknown")

    # Block 1: existing label-free ESM sequence features.
    embeddings = np.load(BASE / "protein_embeddings.npy")
    seq = PCA(n_components=50, random_state=42).fit_transform(StandardScaler().fit_transform(embeddings))

    # Block 2: topology, deliberately separate from functional role.
    topo_rows: list[dict[str, float]] = []
    for row in master.itertuples(index=False):
        d = {
            f"class:{row.membrane_class_v7}": 1.0,
            f"mode:{row.primary_membrane_mode_v7}": 1.0,
        }
        for mode in str(row.secondary_membrane_modes_v7 or "").split(";"):
            if mode and mode.lower() != "nan":
                d[f"secondary:{mode}"] = 1.0
        d["length_log1p"] = float(np.log1p(row.sequence_length))
        topo_rows.append(d)
    topo, topo_meta = sparse_block(topo_rows, 12, "topology")

    # Block 3: four of five annotation axes.  Membrane-role labels are excluded
    # from inputs and used only later for interpretation.
    cross = read_tsv(INPUT / "protein_cross_classification_v72.tsv.gz", usecols=[
        "target_uniprot_id", "classification_axis", "classification_label",
    ], dtype=str)
    allowed = {"structural_family", "molecular_function", "biological_process", "specialist_classification"}
    cross = cross[cross["classification_axis"].isin(allowed)]
    annotations: dict[str, dict[str, float]] = {u: {} for u in ids}
    for r in cross.itertuples(index=False):
        if r.target_uniprot_id in annotations and pd.notna(r.classification_label):
            annotations[r.target_uniprot_id][f"{r.classification_axis}:{r.classification_label}"] = 1.0
    annotation, annotation_meta = sparse_block([annotations[u] for u in ids], 30, "role_agnostic_annotation")

    # Block 4: target pharmacology from pair abundance and assay provenance;
    # no membrane-role label is used.
    pairs = read_tsv(INPUT / "protein_compound_pair_v72.tsv.gz", usecols=[
        "target_uniprot_id", "compound_internal_id", "evidence_count", "distinct_database_count",
        "independent_experiment_count", "independent_structure_count", "best_evidence_tier",
    ], dtype={"target_uniprot_id": str, "compound_internal_id": str})
    pairs = pairs[pairs["target_uniprot_id"].isin(ids)].copy()
    chem = read_tsv(INPUT / "compound_chemical_classification_v721.tsv.gz", usecols=[
        "compound_internal_id", "scaffold_representative_butina_cluster", "structural_tags",
    ], dtype=str)
    pairs = pairs.merge(chem, on="compound_internal_id", how="left")
    pharm_rows: dict[str, dict[str, float]] = {u: {} for u in ids}
    for r in pairs.itertuples(index=False):
        d = pharm_rows[r.target_uniprot_id]
        scaf = str(r.scaffold_representative_butina_cluster)
        if scaf and scaf.lower() != "nan":
            d[f"scaffold:{scaf}"] = d.get(f"scaffold:{scaf}", 0.0) + 1.0
        for tag in str(r.structural_tags).split(";"):
            if tag and tag.lower() != "nan":
                d[f"chemical_tag:{tag}"] = d.get(f"chemical_tag:{tag}", 0.0) + 1.0
    pharm, pharm_meta = sparse_block([pharm_rows[u] for u in ids], 30, "compound_profile")
    agg = pairs.groupby("target_uniprot_id").agg(
        pair_count=("compound_internal_id", "nunique"),
        evidence_sum=("evidence_count", "sum"),
        database_mean=("distinct_database_count", "mean"),
        lineage_mean=("independent_experiment_count", "mean"),
        structure_mean=("independent_structure_count", "mean"),
    )
    numeric = np.array([[np.log1p(agg.loc[u, "pair_count"]) if u in agg.index else 0.0,
                         np.log1p(agg.loc[u, "evidence_sum"]) if u in agg.index else 0.0,
                         agg.loc[u, "database_mean"] if u in agg.index else 0.0,
                         agg.loc[u, "lineage_mean"] if u in agg.index else 0.0,
                         agg.loc[u, "structure_mean"] if u in agg.index else 0.0] for u in ids], dtype=np.float32)
    pharmacology = np.hstack([pharm, StandardScaler().fit_transform(numeric)]).astype(np.float32)

    blocks = {
        "ESM_sequence": equal_block_weight(seq),
        "membrane_topology": equal_block_weight(topo),
        "role_agnostic_annotation": equal_block_weight(annotation),
        "compound_profile_and_evidence": equal_block_weight(pharmacology),
    }
    combined = np.hstack(list(blocks.values())).astype(np.float32)
    final_pca = PCA(n_components=50, random_state=42).fit_transform(combined).astype(np.float32)
    np.save(OUT / "multiview_features_pca50.npy", final_pca)

    # Select HDBSCAN strictly from unsupervised quantities; role labels remain unused.
    scans = []
    models = {}
    for mcs in (10, 20, 30, 50, 75, 100):
        for ms in (None, 5, 10, 20):
            model = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms, metric="euclidean", prediction_data=False).fit(final_pca)
            labels = model.labels_
            keep = labels >= 0
            ncl = len(set(labels[keep]))
            sil = float(silhouette_score(final_pca[keep], labels[keep])) if ncl >= 2 and keep.sum() > ncl else float("nan")
            persistence = float(np.mean(model.cluster_persistence_)) if len(model.cluster_persistence_) else 0.0
            # Penalise trivial one-cluster outputs and excessive noise only.
            score = (0 if np.isnan(sil) else sil) + 0.35 * persistence - 0.20 * (1 - keep.mean()) + 0.02 * min(ncl, 10)
            scans.append({"min_cluster_size": mcs, "min_samples": "None" if ms is None else ms,
                          "n_clusters": ncl, "noise_fraction": float(1 - keep.mean()),
                          "silhouette_non_noise": sil, "mean_persistence": persistence, "selection_score": score})
            models[(mcs, ms)] = model
    scan = pd.DataFrame(scans).sort_values("selection_score", ascending=False)
    scan.to_csv(OUT / "hdbscan_multiview_parameter_scan.tsv", sep="\t", index=False)
    best = scan.iloc[0]
    ms = None if best.min_samples == "None" else int(best.min_samples)
    model = models[(int(best.min_cluster_size), ms)]
    master["cluster_id"] = model.labels_

    # Visualisation only: it does not feed into clustering or parameter selection.
    coords = umap.UMAP(n_components=3, n_neighbors=30, min_dist=0.12, metric="cosine", random_state=42).fit_transform(final_pca)
    master[["UMAP1", "UMAP2", "UMAP3"]] = coords
    master.to_csv(OUT / "protein_multiview_umap3d.tsv", sep="\t", index=False)

    known = master["posthoc_membrane_role"].ne("Unknown")
    metrics = {
        "n_proteins": int(len(master)), "n_clusters": int(len(set(master.cluster_id.astype(int)) - {-1})),
        "noise_fraction": float((master.cluster_id < 0).mean()),
        "ARI_posthoc_membrane_role": float(adjusted_rand_score(master.loc[known, "posthoc_membrane_role"], master.loc[known, "cluster_id"])),
        "NMI_posthoc_membrane_role": float(normalized_mutual_info_score(master.loc[known, "posthoc_membrane_role"], master.loc[known, "cluster_id"])),
        "selected_hdbscan": {k: (int(v) if k in {"min_cluster_size", "n_clusters"} else float(v) if k not in {"min_samples"} else str(v)) for k, v in best.items()},
        "feature_blocks": [
            {"name": "ESM_sequence", "dimensions": int(blocks["ESM_sequence"].shape[1])},
            topo_meta, annotation_meta, pharm_meta,
            {"name": "pharmacology_numeric", "dimensions": 5},
        ],
        "interpretation": "Integrative functional/pharmacological map; posthoc membrane role was excluded from inputs.",
        "random_state": 42,
    }
    (OUT / "multiview_clustering_summary.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    master["cluster_display"] = master.cluster_id.map(lambda x: "Noise / unassigned" if x == -1 else f"Cluster {x}")
    fig = px.scatter_3d(master, x="UMAP1", y="UMAP2", z="UMAP3", color="cluster_display",
                        hover_data=["uniprot_id", "posthoc_membrane_role", "cluster_id", "sequence_length"], opacity=0.60,
                        color_discrete_map={"Noise / unassigned": "#B8B8B8"})
    fig.update_layout(template="simple_white", legend_title_text="HDBSCAN cluster")
    fig.write_html(OUT / "protein_multiview_umap3d_interactive.html", include_plotlyjs=True)

    fig2 = plt.figure(figsize=(7.2, 5.7)); ax = fig2.add_subplot(projection="3d")
    colors = {"Noise / unassigned": "#B8B8B8", "Cluster 0": "#4DAF4A", "Cluster 1": "#4BA3D8", "Cluster 2": "#F29E4C", "Cluster 3": "#A67FC5", "Cluster 4": "#E46C7A", "Cluster 5": "#67A583"}
    for label, group in master.groupby("cluster_display", sort=True):
        ax.scatter(group.UMAP1, group.UMAP2, group.UMAP3, s=4, alpha=0.52, edgecolors="none", label=label, c=colors.get(label, "#666666"), rasterized=True)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_zlabel("UMAP3"); ax.view_init(elev=22, azim=42)
    ax.legend(title="HDBSCAN cluster", loc="center left", bbox_to_anchor=(1.04, .5), markerscale=2)
    save(fig2, "protein_multiview_umap3d_cluster"); plt.close(fig2)

    pd.DataFrame([metrics]).to_csv(OUT / "multiview_metrics.tsv", sep="\t", index=False)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
