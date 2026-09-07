from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image


ROOT = Path(r"C:\Users\Administrator\MemPro_narrative_v3")
DATA = ROOT / "data"
FIG = ROOT / "figures"
QA = ROOT / "qa"
REL = Path(r"D:\finale\20_MemPro_V7.1.1_final_20260814\01_release_tables")
ANATOMY = next(Path(r"C:\Users\Administrator\Desktop").rglob("human_anatomy_base_v1.png"))

C = {
    "blue": "#7b95c6", "cyan": "#49c2d9", "lcyan": "#a1d8e8",
    "green": "#67a583", "lgreen": "#a2c986", "yellow": "#fded95",
    "peach": "#ffc1a6", "salmon": "#f59c7c", "coral": "#f47254",
    "red": "#c85e62", "ink": "#27313d", "muted": "#66717e",
    "grid": "#e5e9ee", "missing": "#b9bec6",
}
PAL = [C["blue"], C["cyan"], C["green"], C["lgreen"], C["yellow"], C["peach"], C["salmon"], C["red"]]

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8.4, "axes.titlesize": 9.4,
    "axes.labelsize": 8.3, "xtick.labelsize": 7.3, "ytick.labelsize": 7.3,
    "text.color": C["ink"], "axes.labelcolor": C["ink"],
    "axes.titlecolor": C["ink"], "axes.edgecolor": "#c7cdd4",
    "grid.color": C["grid"], "pdf.fonttype": 42, "svg.fonttype": "none",
    "legend.frameon": False,
})
sns.set_style("whitegrid")


def panel(ax, label: str, title: str) -> None:
    ax.set_title(f"{label}  {title}", loc="left", weight="bold", pad=8)


def clean(ax, grid: str | None = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)
    if grid:
        ax.grid(axis=grid, alpha=.68)


def save(fig, stem: str, source: str) -> None:
    fig.text(.035, .018, "Source: " + source, fontsize=7, color=C["muted"])
    for ext in ("png", "svg", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=450 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def short(value: object, width: int = 18) -> str:
    return textwrap.shorten(str(value).replace("_", " "), width=width, placeholder="…")


COORD = {
    "Brain": (.50, .09), "Heart": (.52, .34), "Lung": (.43, .30),
    "Liver": (.43, .42), "Kidney": (.59, .49), "Stomach": (.56, .43),
    "Large intestine": (.50, .54), "Small intestine": (.50, .57),
    "Spleen": (.61, .43), "Pancreas": (.53, .47), "Thyroid": (.50, .21),
    "Skin": (.68, .36), "Bone marrow": (.61, .73),
    "Reproductive": (.50, .68), "Bladder": (.50, .66),
}


def body_map(ax, table: pd.DataFrame, value: str) -> None:
    image = np.asarray(Image.open(ANATOMY).convert("RGB"))
    h, w = image.shape[:2]
    ax.imshow(image)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    panel(ax, "A", "Ontology-linked disease burden across the body")
    q = table[table.organ.isin(COORD)].nlargest(7, value).copy()
    norm = mpl.colors.Normalize(q[value].min(), q[value].max())
    cmap = mpl.colormaps["YlOrRd"]
    left_y = [.25, .40, .55, .70]
    right_y = [.28, .47, .66]
    li = ri = 0
    for idx, (_, row) in enumerate(q.iterrows()):
        x0, y0 = COORD[row.organ][0] * w, COORD[row.organ][1] * h
        left = idx % 2 == 0
        if left:
            tx, ty = .02 * w, left_y[li] * h
            elbow = .27 * w
            li += 1
            ha = "left"
        else:
            tx, ty = .98 * w, right_y[ri] * h
            elbow = .73 * w
            ri += 1
            ha = "right"
        ax.scatter(x0, y0, s=72, c=[cmap(norm(row[value]))], edgecolor="white", lw=.8, zorder=5)
        ax.plot([x0, elbow, tx], [y0, ty, ty], color=C["muted"], lw=.65)
        ax.text(tx, ty, f"{row.organ}\n{int(row[value]):,} pairs", ha=ha,
                va="center", fontsize=7.1, weight="bold")


def draw_m2() -> None:
    d = pd.read_csv(DATA / "M2_protein_annotation_umap.tsv", sep="\t")
    cc = pd.read_csv(REL / "protein_cross_classification_v71.tsv.gz", sep="\t",
                     usecols=["target_uniprot_id", "classification_axis"], low_memory=False)
    joined = d[["target_uniprot_id", "cluster"]].merge(cc, on="target_uniprot_id", how="left")
    counts = d.cluster.value_counts()
    largest = counts.drop(index=-1, errors="ignore").head(12).index
    tab = pd.crosstab(joined.cluster, joined.classification_axis).reindex(largest).fillna(0)
    preferred = ["structural_family", "molecular_function", "biological_process", "membrane_role", "specialist_classification"]
    tab = tab.reindex(columns=preferred, fill_value=0)
    share = tab.div(tab.sum(axis=1), axis=0)
    axis_codes = ["Family", "Function", "Process", "Role", "Specialist"]

    annotations = pd.read_csv(DATA / "M2_cluster_top_annotations.tsv", sep="\t")
    top_label = annotations.drop_duplicates("cluster").set_index("cluster")["label"].to_dict()

    fig = plt.figure(figsize=(14, 9.4))
    fig.text(.03, .968, "M2 | Multi-axis annotations reveal reproducible membrane-protein neighborhoods",
             fontsize=15, weight="bold", va="top")
    fig.text(.03, .932, "Five independent annotation axes preserve overlapping biology while HDBSCAN leaves diffuse proteins unforced.",
             fontsize=9, color=C["muted"], va="top")
    gs = fig.add_gridspec(2, 3, left=.055, right=.985, top=.87, bottom=.09,
                          hspace=.34, wspace=.33, width_ratios=[1.2, 1.2, 1.0])

    ax = fig.add_subplot(gs[:, :2])
    panel(ax, "A", "Annotation-space UMAP separates overlapping functional neighborhoods")
    for cls, color in zip(["A", "B", "C"], [C["blue"], C["green"], C["peach"]]):
        q = d[d.membrane_class_v7.eq(cls)]
        ax.scatter(q.umap1, q.umap2, s=7, c=color, alpha=.30, label=f"Class {cls}", rasterized=True)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(fontsize=7.2, markerscale=2.2, ncol=3, loc="upper right")
    clean(ax, None)

    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "B", "Dominant communities have distinct annotation identities")
    order = list(reversed(largest))
    y = np.arange(len(order))
    vals = counts.reindex(order).values
    ax.hlines(y, 0, vals, color=C["grid"], lw=3)
    ax.scatter(vals, y, s=45, c=[PAL[int(c) % len(PAL)] for c in order], edgecolor="white", lw=.45)
    labels = [f"C{c} · {short(top_label.get(c, 'unresolved'), 16)}" for c in order]
    ax.set_yticks(y, labels)
    ax.set_xlabel("proteins")
    clean(ax, "x")

    ax = fig.add_subplot(gs[1, 2])
    panel(ax, "C", "Communities differ in which annotation axes support them")
    sns.heatmap(share, ax=ax, cmap="YlGnBu", linewidths=.45, linecolor="white",
                xticklabels=axis_codes, yticklabels=[f"C{x}" for x in largest],
                cbar_kws={"label": "within-community share"})
    ax.set_xlabel("")
    ax.set_ylabel("community")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)

    save(fig, "M2_multi_axis_annotation_space_cleanfinal",
         "MemPro V7.1.1 protein master and five-axis cross-classification; Jaccard UMAP and HDBSCAN, seed=42.")


def draw_m4() -> None:
    burden = pd.read_csv(DATA / "M4_body_disease.tsv", sep="\t")
    concordance = pd.read_csv(DATA / "M4_disease_expression_concordance.tsv", sep="\t")
    role_share = pd.read_csv(DATA / "M4_system_role_share.tsv", sep="\t", index_col=0)

    sys_codes = {
        "nervous system": "Nervous", "musculoskeletal system": "Musculo",
        "digestive system": "Digestive", "alimentary part of gastrointestinal system": "GI tract",
        "sensory organ system": "Sensory", "circulatory system": "Circulatory",
        "integumental system": "Skin", "excretory system": "Renal",
        "reproductive system": "Reproductive", "hematopoietic system": "Heme",
        "respiratory system": "Respiratory", "immune system": "Immune",
    }
    role_codes = {
        "Receptor": "Rec", "Enzyme": "Enz", "Channel": "Chan", "Transporter": "Trans",
        "Role unresolved": "Unres", "Scaffold": "Scaf", "Signaling": "Signal",
        "Adhesion": "Adhes", "Immune": "Immune",
    }
    role_share.index = [sys_codes.get(str(x).lower(), short(x, 12)) for x in role_share.index]
    role_share.columns = [role_codes.get(str(x), short(x, 8)) for x in role_share.columns]

    relation = pd.read_csv(REL / "protein_disease_relation_v71.tsv.gz", sep="\t",
                           usecols=["target_uniprot_id", "canonical_disease_id", "best_evidence_level"], low_memory=False)
    anatomy = pd.read_csv(REL / "disease_anatomical_system_v71.tsv.gz", sep="\t",
                          usecols=["canonical_disease_id", "anatomical_system_name"], low_memory=False).drop_duplicates()
    joined = relation.merge(anatomy, on="canonical_disease_id", how="inner")
    joined["pair"] = joined.target_uniprot_id + "|" + joined.canonical_disease_id
    evidence = pd.crosstab(joined.anatomical_system_name, joined.best_evidence_level,
                           values=joined.pair, aggfunc=pd.Series.nunique).fillna(0)
    evidence = evidence.loc[evidence.sum(axis=1).nlargest(10).index]
    evidence = evidence.reindex(columns=["medium", "high", "very high"], fill_value=0)
    evidence.index = [sys_codes.get(str(x).lower(), short(x, 12)) for x in evidence.index]
    evidence.columns = ["Medium", "High", "Very high"]

    fig = plt.figure(figsize=(14, 9.6))
    fig.text(.03, .968, "M4 | Disease anatomy and normal-tissue expression expose biologically coherent target contexts",
             fontsize=15, weight="bold", va="top")
    fig.text(.03, .932, "Ontology-derived multi-label systems are compared with HPA expression and a 500-permutation null model.",
             fontsize=9, color=C["muted"], va="top")
    gs = fig.add_gridspec(2, 3, left=.035, right=.985, top=.87, bottom=.09,
                          width_ratios=[.9, 1.05, 1.16], hspace=.38, wspace=.34)

    ax = fig.add_subplot(gs[:, 0])
    body_map(ax, burden, "protein_disease_pairs")

    ax = fig.add_subplot(gs[:, 1])
    panel(ax, "B", "Observed disease–expression overlap exceeds matched random expectation")
    q = concordance.sort_values("zscore")
    y = np.arange(len(q))
    ax.hlines(y, 0, q.zscore, color=C["grid"], lw=3)
    ax.scatter(q.zscore, y, c=q.observed_fraction, cmap="YlOrRd", s=72, edgecolor="white", lw=.55)
    ax.axvline(0, color=C["muted"], lw=.8)
    ax.set_yticks(y, q.organ)
    ax.set_xlabel("concordance z score")
    clean(ax, "x")

    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "C", "Disease systems differ in membrane-role composition")
    sns.heatmap(role_share, ax=ax, cmap="YlGnBu", linewidths=.45, linecolor="white",
                cbar_kws={"label": "within-system share"})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)

    ax = fig.add_subplot(gs[1, 2])
    panel(ax, "D", "Evidence strength varies by anatomical system")
    sns.heatmap(np.log1p(evidence), ax=ax, cmap=sns.light_palette(C["red"], as_cmap=True),
                linewidths=.45, linecolor="white", cbar_kws={"label": "log1p(unique pairs)"})
    ax.set_xlabel("best disease-evidence level")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)

    save(fig, "M4_disease_expression_concordance_cleanfinal",
         "Open Targets 26.06, UniProtKB, MONDO, Disease Ontology, Uberon and HPA 25.1; multi-label anatomy; seed=42.")


def draw_m5() -> None:
    ecfp = pd.read_csv(DATA / "M5_ECFP4_chemical_space.tsv", sep="\t")
    pagtn = pd.read_csv(DATA / "M5_PAGTN_pilot_space.tsv", sep="\t")
    qa = json.loads((QA / "M5_embedding_QA.json").read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(14, 8.9))
    fig.text(.03, .966, "M5 | Chemical-space structure depends on both representation and biological status",
             fontsize=15, weight="bold", va="top")
    fig.text(.03, .928, "An interpretable ECFP4 baseline is contrasted with a property-supervised PAGTN pilot; neither 2D map is an identity rule.",
             fontsize=9, color=C["muted"], va="top")
    gs = fig.add_gridspec(1, 3, left=.055, right=.985, top=.84, bottom=.11,
                          width_ratios=[1.25, 1.05, .8], wspace=.32)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A", "ECFP4 maps local chemical-fragment neighborhoods")
    ax.hexbin(ecfp.umap1, ecfp.umap2, gridsize=55, mincnt=1, cmap="Blues", linewidths=0, alpha=.72)
    statuses = [
        ("is_approved_drug", "Approved", C["red"]),
        ("is_clinical_candidate", "Clinical", C["coral"]),
        ("is_endogenous_ligand", "Endogenous", C["green"]),
        ("is_natural_product", "Natural", C["yellow"]),
    ]
    for col, label, color in statuses:
        q = ecfp[pd.to_numeric(ecfp[col], errors="coerce").fillna(0).eq(1)]
        if len(q):
            q = q.sample(min(450, len(q)), random_state=42)
            ax.scatter(q.umap1, q.umap2, s=8, c=color, alpha=.45, label=label, rasterized=True)
    ax.legend(fontsize=7.1, markerscale=1.8, ncol=2, loc="upper right")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    clean(ax, None)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B", "PAGTN pilot emphasizes learned physicochemical gradients")
    top_classes = pagtn.structural_class.fillna("unknown").value_counts().head(6).index
    for idx, cls in enumerate(top_classes):
        q = pagtn[pagtn.structural_class.fillna("unknown").eq(cls)]
        ax.scatter(q.umap1, q.umap2, s=8, color=PAL[idx], alpha=.40,
                   label=short(cls, 18), rasterized=True)
    other = pagtn[~pagtn.structural_class.fillna("unknown").isin(top_classes)]
    ax.scatter(other.umap1, other.umap2, s=5, color=C["missing"], alpha=.22, label="Other", rasterized=True)
    ax.legend(fontsize=6.8, markerscale=1.7, loc="upper right")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    clean(ax, None)

    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "C", "Embedding quality metrics reveal different trade-offs")
    metrics = ["kNN retention", "cluster purity", "assigned fraction"]
    e_vals = [qa["ecfp4_knn_retention"], qa["ecfp4_cluster_purity"], 1 - qa["ecfp4_noise_fraction"]]
    p_vals = [qa["pagtn_knn_retention"], qa["pagtn_cluster_purity"], 1 - qa["pagtn_noise_fraction"]]
    y = np.arange(len(metrics))
    for yy, a, b in zip(y, e_vals, p_vals):
        ax.plot([a, b], [yy, yy], color=C["grid"], lw=3, zorder=1)
    ax.scatter(e_vals, y, s=62, color=C["blue"], label="ECFP4")
    ax.scatter(p_vals, y, s=62, color=C["coral"], marker="D", label="PAGTN pilot")
    ax.set_yticks(y, metrics)
    ax.set_xlim(0, 1.03)
    ax.set_xlabel("score (0–1)")
    ax.legend(fontsize=7.2, loc="lower right")
    clean(ax, "x")
    ax.text(.02, -.12, "PAGTN: property-supervised pilot, n=2,400\nECFP4: deterministic sample, n=15,000",
            transform=ax.transAxes, fontsize=7, color=C["muted"], va="top")

    save(fig, "M5_chemical_space_embedding_comparison_cleanfinal",
         "MemPro V7.1.1 compound master; ECFP4 radius=2/1024 bits; UMAP/HDBSCAN seed=42; PAGTN property-supervised pilot.")


def draw_m6() -> None:
    nodes = pd.read_csv(DATA / "M6_target_network_nodes.tsv", sep="\t")
    edges = pd.read_csv(DATA / "M6_target_network_edges.tsv", sep="\t")
    side = pd.read_csv(DATA / "M6_membrane_side_compound_properties.tsv", sep="\t")
    side_map = {
        "Extra, unresolved": "Extra", "Interface/mixed": "Interface",
        "Intramembrane": "Intramembrane", "Cytoplasmic": "Cytoplasmic",
        "Non-cytoplasmic": "Non-cytoplasmic", "Both/mixed": "Both/mixed",
    }
    side["side"] = side.membrane_side.map(side_map).fillna(side.membrane_side)
    order = [x for x in ["Extra", "Interface", "Intramembrane", "Cytoplasmic", "Non-cytoplasmic"] if x in set(side.side)]

    graph = nx.Graph()
    graph.add_nodes_from(nodes.node.astype(int))
    for row in edges.itertuples(index=False):
        graph.add_edge(int(row.source), int(row.target), weight=float(row.shared_compounds))
    pos = nx.spring_layout(graph, seed=42, weight="weight", k=.38, iterations=220)
    max_w = max((a[2]["weight"] for a in graph.edges(data=True)), default=1)

    fig = plt.figure(figsize=(14, 9.1))
    fig.text(.03, .968, "M6 | Polypharmacology communities connect to membrane-side chemical environments",
             fontsize=15, weight="bold", va="top")
    fig.text(.03, .932, "Shared-ligand target projection is paired with site-level physicochemical distributions; membrane side is not inferred from ligand properties.",
             fontsize=9, color=C["muted"], va="top")
    gs = fig.add_gridspec(2, 3, left=.055, right=.985, top=.86, bottom=.09,
                          width_ratios=[1.2, 1.2, 1.0], hspace=.38, wspace=.34)

    ax = fig.add_subplot(gs[:, :2])
    panel(ax, "A", "Targets sharing ligands form polypharmacology communities")
    nx.draw_networkx_edges(graph, pos, ax=ax,
                           width=[.15 + 1.4 * (data["weight"] / max_w) ** .55 for _, _, data in graph.edges(data=True)],
                           alpha=.18, edge_color=C["muted"])
    node_colors = [PAL[int(nodes.set_index("node").loc[n, "community"]) % len(PAL)] for n in graph.nodes]
    node_sizes = [16 + 5 * np.sqrt(graph.degree(n)) for n in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_color=node_colors, node_size=node_sizes,
                           alpha=.78, edgecolors="white", linewidths=.25)
    ax.axis("off")
    ax.text(.01, .02, "220 highest-degree targets; edge = ≥3 shared compounds; Louvain communities",
            transform=ax.transAxes, fontsize=7.2, color=C["muted"])

    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "B", "Ligand lipophilicity varies by annotated membrane side")
    sns.boxenplot(data=side[side.side.isin(order)], y="side", x="xlogp", order=order,
                  palette=[C["blue"], C["cyan"], C["green"], C["peach"], C["salmon"]][:len(order)],
                  linewidth=.5, k_depth="proportion", ax=ax)
    ax.set_xlabel("XlogP")
    ax.set_ylabel("")
    clean(ax, "x")

    ax = fig.add_subplot(gs[1, 2])
    panel(ax, "C", "Polar surface area shows a complementary site-context pattern")
    sns.boxenplot(data=side[side.side.isin(order)], y="side", x="tpsa", order=order,
                  palette=[C["blue"], C["cyan"], C["green"], C["peach"], C["salmon"]][:len(order)],
                  linewidth=.5, k_depth="proportion", ax=ax)
    ax.set_xlabel("TPSA (Å²)")
    ax.set_ylabel("")
    clean(ax, "x")

    save(fig, "M6_polypharmacology_and_membrane_side_cleanfinal",
         "MemPro V7.1.1 released pairs, interaction sites and compound properties; target projection top 220; Louvain seed=42.")


def draw_upset(bar_ax, matrix_ax, profiles: pd.DataFrame) -> None:
    profiles = profiles.head(7).copy()
    profiles["code"] = [f"P{i+1}" for i in range(len(profiles))]
    x = np.arange(len(profiles))
    bar_ax.bar(x, profiles.evidence_records, color=C["blue"])
    bar_ax.set_yscale("log")
    bar_ax.set_ylabel("records (log)")
    bar_ax.set_xticks(x, [])
    clean(bar_ax, "y")
    modes = ["Quant", "Func", "Struct", "Site", "Identity"]
    matrix_ax.set_xlim(-.5, len(profiles)-.5)
    matrix_ax.set_ylim(-.6, len(modes)-.4)
    matrix_ax.set_yticks(range(len(modes)), modes)
    matrix_ax.set_xticks(x, profiles.code)
    matrix_ax.set_xlabel("evidence-profile code")
    matrix_ax.grid(False)
    matrix_ax.spines[:].set_visible(False)
    for col, row in enumerate(profiles.itertuples(index=False)):
        active = set(str(row.profile).split("+"))
        if row.profile == "identity_only":
            active = {"identity"}
        ys = []
        for yi, mode in enumerate(["quantitative", "functional", "structure", "site", "identity"]):
            on = mode in active
            matrix_ax.scatter(col, yi, s=31, color=C["ink"] if on else "#d9dde2", zorder=3)
            if on:
                ys.append(yi)
        if len(ys) > 1:
            matrix_ax.plot([col, col], [min(ys), max(ys)], color=C["ink"], lw=.9, zorder=2)
    legend = "   ".join(f"{r.code}: {str(r.profile).replace('_', ' ')}" for r in profiles.itertuples(index=False))
    matrix_ax.text(0, -.42, legend, transform=matrix_ax.transAxes, fontsize=6.6,
                   color=C["muted"], va="top", wrap=True)


def draw_m7() -> None:
    profiles = pd.read_csv(DATA / "M7_evidence_profile_intersections.tsv", sep="\t")
    lineage = pd.read_csv(DATA / "M7_source_lineage_matrix.tsv", sep="\t", index_col=0)
    modality = pd.read_csv(DATA / "M7_modality_tier_matrix.tsv", sep="\t", index_col=0)
    lineage.columns = ["DB", "Literature", "PubChem", "Structure"]
    lineage.index = [short(x, 21) for x in lineage.index]
    modality.index = [short(x, 18) for x in modality.index]

    fig = plt.figure(figsize=(14, 9.7))
    fig.text(.03, .968, "M7 | Evidence profiles distinguish record volume from experimental lineage and strength",
             fontsize=15, weight="bold", va="top")
    fig.text(.03, .932, "UpSet combinations are record-level profiles; database count is never presented as independent-experiment count.",
             fontsize=9, color=C["muted"], va="top")
    outer = fig.add_gridspec(2, 2, left=.075, right=.985, top=.87, bottom=.14,
                             height_ratios=[1.12, 1], hspace=.43, wspace=.34)
    upset = outer[0, :].subgridspec(2, 1, height_ratios=[1.25, .92], hspace=.03)
    bar_ax = fig.add_subplot(upset[0])
    panel(bar_ax, "A", "Evidence-profile intersections reveal dominant and multimodal record types")
    matrix_ax = fig.add_subplot(upset[1])
    draw_upset(bar_ax, matrix_ax, profiles)

    ax = fig.add_subplot(outer[1, 0])
    panel(ax, "B", "Sources resolve to different lineage-key types")
    sns.heatmap(np.log1p(lineage), ax=ax, cmap="YlGnBu", linewidths=.45, linecolor="white",
                cbar_kws={"label": "log1p(records)"})
    ax.set_xlabel("lineage-key type")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)

    ax = fig.add_subplot(outer[1, 1])
    panel(ax, "C", "Evidence modality determines the BE-tier composition")
    sns.heatmap(np.log1p(modality), ax=ax, cmap=sns.light_palette(C["red"], as_cmap=True),
                linewidths=.45, linecolor="white", cbar_kws={"label": "log1p(records)"})
    ax.set_xlabel("evidence tier")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)

    save(fig, "M7_evidence_profiles_and_lineage_cleanfinal",
         "MemPro V7.1.1 positive interaction evidence and lineage fields; combinations are record-level, not independent experiments.")


if __name__ == "__main__":
    draw_m2()
    draw_m4()
    draw_m5()
    draw_m6()
    draw_m7()
    (QA / "CLEAN_FINAL_RENDER_COMPLETE.json").write_text(
        json.dumps({"M2": "complete", "M4": "complete", "M5": "complete", "M6": "complete", "M7": "complete"}, indent=2),
        encoding="utf-8",
    )
    print("clean final problem panels complete")
