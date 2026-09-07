from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, PathPatch, Polygon
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr


BASE = Path(r"C:\Users\Administrator\MemPro_narrative_v3")
REL = Path(r"D:\finale\20_MemPro_V7.1.1_final_20260814\01_release_tables")
V3X = Path(r"D:\7.22\v71_publication_candidate_working\work\v3_1_activity_source_crosswalk_v71.tsv.gz")
OUT = Path(r"C:\Users\Administrator\Desktop\gkg sjf summer intern\MemPro最终展示\MemPro_M1-M4_科学修订版_20260816")
DATA = OUT / "01_panel_data"
FIG = OUT / "02_figures"
EDIT = OUT / "03_editable_assets"
QA = OUT / "04_QA"
for p in (DATA, FIG, EDIT, QA):
    p.mkdir(parents=True, exist_ok=True)

C = {
    "blue": "#7b95c6", "cyan": "#49c2d9", "lcyan": "#a1d8e8",
    "green": "#67a583", "lgreen": "#a2c986", "yellow": "#fded95",
    "peach": "#ffc1a6", "salmon": "#f59c7c", "coral": "#f47254",
    "red": "#c85e62", "ink": "#27313d", "muted": "#66717e",
    "grid": "#e5e9ee", "missing": "#b9bec6", "skin": "#f7e8df",
}
PAL = [C["blue"], C["cyan"], C["green"], C["lgreen"], C["yellow"], C["peach"], C["salmon"], C["red"]]
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8.2, "axes.titlesize": 9.2,
    "axes.labelsize": 8.1, "xtick.labelsize": 7.0, "ytick.labelsize": 7.0,
    "text.color": C["ink"], "axes.labelcolor": C["ink"],
    "axes.titlecolor": C["ink"], "axes.edgecolor": "#c7cdd4",
    "grid.color": C["grid"], "pdf.fonttype": 42, "svg.fonttype": "none",
    "legend.frameon": False,
})
sns.set_style("whitegrid")


def panel(ax, letter, title):
    ax.set_title(f"{letter}  {title}", loc="left", weight="bold", pad=8)


def clean(ax, grid="x"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)
    if grid:
        ax.grid(axis=grid, alpha=.7)


def save(fig, stem, source):
    fig.text(.03, .018, "Source: " + source, fontsize=6.8, color=C["muted"])
    for ext in ("png", "svg", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=450 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fmt(v):
    v = float(v)
    return f"{v/1e6:.2f}M" if v >= 1e6 else f"{v/1e3:.1f}k" if v >= 1e3 else f"{int(v):,}"


def bh_fdr(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.empty_like(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked, 1)
    return q


def rebuild_m1_data():
    cols = ["evidence_id", "target_uniprot_id", "compound_internal_id", "source_database", "source_record_id"]
    e = pd.read_csv(REL / "positive_interaction_evidence_v71.tsv.gz", sep="\t", usecols=cols, low_memory=False)
    cross = pd.read_csv(V3X, sep="\t", dtype=str, usecols=["activity_measurement_id", "source_database"])
    amap = cross.drop_duplicates("activity_measurement_id").set_index("activity_measurement_id").source_database
    legacy = e.source_database.eq("v3.1_mixed_sources")
    recovered = e.loc[legacy, "source_record_id"].map(amap)
    e["figure_source"] = e.source_database
    e.loc[legacy, "figure_source"] = recovered.fillna("Legacy v3.1 source unresolved")

    # Expand semicolon-labelled contributing databases into individual provenance labels.
    x = e.assign(figure_source=e.figure_source.fillna("Unspecified").str.split(";"))
    x = x.explode("figure_source")
    x["figure_source"] = x.figure_source.str.strip()
    x["pair"] = x.target_uniprot_id.astype(str) + "|" + x.compound_internal_id.astype(str)
    x = x[["pair", "figure_source"]].drop_duplicates()
    x.to_csv(DATA / "M1_evidence_source_decomposition.tsv.gz", sep="\t", index=False, compression="gzip")
    source_sets = {s: set(q.pair) for s, q in x.groupby("figure_source")}
    sources = sorted(source_sets, key=lambda s: len(source_sets[s]), reverse=True)

    memberships = {}
    pair_sources = x.groupby("pair").figure_source.agg(lambda z: tuple(sorted(set(z))))
    for ss, n in pair_sources.value_counts().items():
        memberships[ss] = int(n)
    inter = pd.DataFrame([
        {"intersection": " & ".join(k), "pair_count": v, "degree": len(k)}
        for k, v in memberships.items()
    ]).sort_values("pair_count", ascending=False)
    inter.to_csv(DATA / "M1_source_intersections_revised.tsv", sep="\t", index=False)

    jac = pd.DataFrame(index=sources, columns=sources, dtype=float)
    for a in sources:
        for b in sources:
            jac.loc[a, b] = len(source_sets[a] & source_sets[b]) / max(1, len(source_sets[a] | source_sets[b]))
    jac.to_csv(DATA / "M1_source_jaccard_revised.tsv", sep="\t")

    remaining = set(pair_sources.index)
    covered = set()
    greedy = []
    while remaining:
        best = max(sources, key=lambda s: len(source_sets[s] - covered) if s not in [r[0] for r in greedy] else -1)
        gain = source_sets[best] - covered
        if not gain or best in [r[0] for r in greedy]:
            break
        covered |= source_sets[best]
        greedy.append((best, len(gain), len(covered)))
        remaining = set(pair_sources.index) - covered
    gd = pd.DataFrame(greedy, columns=["source", "new_pairs", "cumulative_pairs"])
    gd.to_csv(DATA / "M1_greedy_coverage_revised.tsv", sep="\t", index=False)
    uq = pd.DataFrame({
        "source": sources,
        "strictly_unique_pairs": [sum(len(ss) == 1 and ss[0] == s for ss in pair_sources) for s in sources]
    })
    uq.to_csv(DATA / "M1_source_unique_contribution_revised.tsv", sep="\t", index=False)
    report = {
        "legacy_mixed_rows": int(legacy.sum()),
        "legacy_rows_recovered": int(recovered.notna().sum()),
        "legacy_rows_unresolved": int(recovered.isna().sum()),
        "recovered_source_counts": recovered.value_counts(dropna=False).to_dict(),
        "figure_sources": {s: len(source_sets[s]) for s in sources},
    }
    (QA / "M1_SOURCE_DECOMPOSITION_QA.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return inter, jac, gd, uq, sources


def draw_m1():
    inter, jac, greedy, uniq, sources = rebuild_m1_data()
    codes = {s: chr(65 + i) if i < 26 else f"S{i+1}" for i, s in enumerate(sources)}
    arr = inter.head(10).copy()
    arr["active"] = arr.intersection.map(lambda z: str(z).split(" & "))
    fig = plt.figure(figsize=(14, 10))
    fig.text(.03, .972, "M1 | Traceable source decomposition reveals complementary and redundant evidence", fontsize=15, weight="bold", va="top")
    fig.text(.03, .938, "Legacy v3.1 records are restored to their assay-level source before pair-level overlap is calculated.", fontsize=9, color=C["muted"], va="top")
    gs = fig.add_gridspec(2, 2, left=.07, right=.985, top=.87, bottom=.09, hspace=.40, wspace=.34)
    left = gs[0, 0].subgridspec(2, 1, height_ratios=[1.15, .9], hspace=.08)
    ax = fig.add_subplot(left[0]); panel(ax, "A", "Largest pair-level source intersections")
    xx = np.arange(len(arr)); ax.bar(xx, arr.pair_count, color=C["blue"])
    ax.set_yscale("log"); ax.set_ylabel("unique protein–compound pairs (log)"); ax.set_xticks(xx, []); clean(ax, "y")
    am = fig.add_subplot(left[1]); am.set_xlim(-.5, len(arr)-.5); am.set_ylim(-.5, len(sources)-.5); am.axis("off")
    for i, s in enumerate(sources):
        am.text(-.65, len(sources)-1-i, f"{codes[s]}  {s}", ha="right", va="center", fontsize=6.5)
    for j, r in arr.reset_index(drop=True).iterrows():
        ys = []
        for i, s in enumerate(sources):
            y = len(sources)-1-i; on = s in r.active
            am.scatter(j, y, s=22, color=C["ink"] if on else "#d9dde2")
            if on: ys.append(y)
        if len(ys) > 1: am.plot([j, j], [min(ys), max(ys)], color=C["ink"], lw=.9)
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "B", "Pair-set Jaccard similarity quantifies database redundancy")
    dist = np.clip(1-jac.values, 0, 1); order = leaves_list(linkage(squareform(dist, checks=False), method="average"))
    sns.heatmap(jac.values[np.ix_(order, order)], ax=ax, cmap="YlGnBu", vmin=0, vmax=1, square=True,
                xticklabels=[codes[sources[i]] for i in order], yticklabels=[codes[sources[i]] for i in order],
                cbar_kws={"label": "pair-set Jaccard"})
    ax.set_xlabel("source code"); ax.set_ylabel("source code")
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "C", "Marginal gain shows what each added source contributes")
    g = greedy; xx = np.arange(len(g)); ax.bar(xx, g.new_pairs, color=[PAL[i % len(PAL)] for i in xx])
    ax.plot(xx, g.cumulative_pairs, color=C["ink"], marker="o", lw=1.2, label="cumulative")
    ax.set_xticks(xx, [codes.get(s, "?") for s in g.source]); ax.set_ylabel("unique pairs"); ax.set_xlabel("source code in greedy order"); ax.legend(fontsize=7); clean(ax, "y")
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "D", "Strictly unique pairs identify irreplaceable contributions")
    u = uniq.sort_values("strictly_unique_pairs"); yy = np.arange(len(u)); ax.hlines(yy, 0, u.strictly_unique_pairs, color=C["grid"], lw=3)
    ax.scatter(u.strictly_unique_pairs, yy, s=48, color=C["coral"]); ax.set_yticks(yy, [codes.get(s, "?") for s in u.source]); ax.set_xlabel("strictly unique released pairs")
    for y, v in zip(yy, u.strictly_unique_pairs): ax.text(v, y, "  " + fmt(v), va="center", fontsize=6.8)
    clean(ax)
    save(fig, "M1_source_overlap_revised", "V7.1.1 positive evidence; v3.1 ACT→assay crosswalk; contributing databases are not assumed to be independent experiments.")


def draw_m2():
    d = pd.read_csv(BASE / "data" / "M2_protein_annotation_umap.tsv", sep="\t")
    cc = pd.read_csv(REL / "protein_cross_classification_v71.tsv.gz", sep="\t", usecols=["target_uniprot_id", "classification_axis"])
    x = d[["target_uniprot_id", "cluster"]].merge(cc, on="target_uniprot_id", how="left")
    tab = pd.crosstab(x.cluster, x.classification_axis).drop(index=-1, errors="ignore")
    order = d[d.cluster.ne(-1)].cluster.value_counts().sort_values(ascending=False).index
    tab = tab.reindex(order).fillna(0); share = tab.div(tab.sum(1), axis=0)
    counts = d.cluster.value_counts().rename_axis("community").rename("protein_count").reset_index()
    counts.to_csv(DATA / "M2_all_communities.tsv", sep="\t", index=False)
    share.to_csv(DATA / "M2_all_community_axis_share.tsv", sep="\t")

    fig = plt.figure(figsize=(14, 9.6))
    fig.text(.03, .969, "M2 | Five-axis cross-classification preserves membrane-protein multifunctionality", fontsize=15, weight="bold", va="top")
    fig.text(.03, .934, f"All {len(order)} HDBSCAN communities are shown; unclustered proteins remain explicit noise rather than being forced into a label.", fontsize=9, color=C["muted"], va="top")
    gs = fig.add_gridspec(2, 3, left=.055, right=.985, top=.87, bottom=.09, hspace=.35, wspace=.34, width_ratios=[1.15, 1.1, 1.0])
    ax = fig.add_subplot(gs[:, 0:2]); panel(ax, "A", "Annotation-space embedding shows overlapping functional neighborhoods")
    for cl, col in zip(["A", "B", "C"], [C["blue"], C["green"], C["peach"]]):
        q = d[d.membrane_class_v7.eq(cl)]; ax.scatter(q.umap1, q.umap2, s=6, c=col, alpha=.22, label=f"Class {cl}", rasterized=True)
    # Community centroids make the otherwise dense embedding interpretable.
    cent = d[d.cluster.ne(-1)].groupby("cluster")[["umap1", "umap2"]].median()
    for k, r in cent.iterrows(): ax.text(r.umap1, r.umap2, str(k), fontsize=6.5, weight="bold", ha="center", va="center", color=C["ink"], bbox=dict(boxstyle="circle,pad=.18", fc="white", ec=C["grid"], alpha=.9))
    ax.set_xlabel("UMAP 1 (annotation similarity)"); ax.set_ylabel("UMAP 2 (annotation similarity)"); ax.legend(fontsize=7, markerscale=2); clean(ax, None)
    ax = fig.add_subplot(gs[0, 2]); panel(ax, "B", f"All {len(order)} communities plus noise")
    cnt = d.cluster.value_counts(); labels = list(order) + [-1]; vals = [int(cnt.get(k, 0)) for k in labels]
    yy = np.arange(len(labels)); colors = [PAL[i % len(PAL)] for i in range(len(order))] + [C["missing"]]
    ax.barh(yy, vals, color=colors, alpha=.85); ax.set_yticks(yy, [str(k) for k in order] + ["Noise"]); ax.invert_yaxis(); ax.set_xlabel("proteins"); clean(ax, "x")
    ax = fig.add_subplot(gs[1, 2]); panel(ax, "C", "Five-axis composition of every community")
    sns.heatmap(share, ax=ax, cmap="YlGnBu", linewidths=.25, linecolor="white", cbar_kws={"label": "within-community assertion share"})
    ax.set_xlabel(""); ax.set_ylabel("community (ranked by size)"); ax.tick_params(axis="x", rotation=40)
    save(fig, "M2_five_axis_all_communities_revised", "V7.1.1 structural family, molecular function, biological process, membrane role and specialist classification; 23 HDBSCAN communities + noise.")


TISSUE_ORGAN = {
    "brain": "Brain", "cerebral cortex": "Brain", "cerebellum": "Brain", "heart muscle": "Heart",
    "lung": "Lung", "liver": "Liver", "kidney": "Kidney", "stomach": "Stomach",
    "colon": "Large intestine", "rectum": "Large intestine", "small intestine": "Small intestine",
    "duodenum": "Small intestine", "spleen": "Spleen", "pancreas": "Pancreas",
    "thyroid gland": "Thyroid", "skin": "Skin", "bone marrow": "Bone marrow",
    "testis": "Reproductive", "ovary": "Reproductive", "endometrium": "Reproductive",
    "prostate": "Reproductive", "urinary bladder": "Bladder"
}
SYSTEM_ORGAN = {
    "nervous system": "Brain", "circulatory system": "Heart", "cardiovascular system": "Heart",
    "respiratory system": "Lung", "digestive system": "Large intestine",
    "alimentary part of gastrointestinal system": "Large intestine", "excretory system": "Kidney",
    "hematopoietic system": "Bone marrow", "lymphoid system": "Spleen",
    "integumental system": "Skin", "reproductive system": "Reproductive",
    "musculoskeletal system": "Bone marrow", "exocrine system": "Pancreas"
}


def draw_human_vector(ax, values, value_col, title, cmap="YlGnBu", label_suffix="proteins", letter="A"):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.axis("off"); panel(ax, letter, title)
    # Front-facing schematic. Every component remains a separate SVG path/shape.
    skin = dict(fc=C["skin"], ec="#a8afb7", lw=.65, alpha=.78)
    ax.add_patch(Circle((.50, .89), .065, **skin))
    ax.add_patch(FancyBboxPatch((.445, .43), .11, .39, boxstyle="round,pad=.018,rounding_size=.07", **skin))
    ax.add_patch(Polygon([[.45,.78],[.31,.66],[.22,.46],[.255,.445],[.38,.65],[.45,.68]], closed=True, **skin))
    ax.add_patch(Polygon([[.55,.78],[.69,.66],[.78,.46],[.745,.445],[.62,.65],[.55,.68]], closed=True, **skin))
    ax.add_patch(Polygon([[.465,.44],[.405,.12],[.455,.10],[.50,.42]], closed=True, **skin))
    ax.add_patch(Polygon([[.535,.44],[.595,.12],[.545,.10],[.50,.42]], closed=True, **skin))
    ax.add_patch(Ellipse((.50,.89),.10,.085,fc="#e9b4b8",ec="none",alpha=.8))
    ax.add_patch(Ellipse((.50,.795),.035,.025,fc="#d49a62",ec="none"))
    ax.add_patch(Ellipse((.455,.70),.075,.12,fc="#efb4b4",ec="#c98c8c",lw=.4,alpha=.9))
    ax.add_patch(Ellipse((.545,.70),.075,.12,fc="#efb4b4",ec="#c98c8c",lw=.4,alpha=.9))
    ax.add_patch(PathPatch(MplPath([(0.50,.69),(.47,.68),(.48,.63),(.50,.60),(.53,.64),(.53,.68),(0.50,.69)], [MplPath.MOVETO,MplPath.CURVE3,MplPath.CURVE3,MplPath.CURVE3,MplPath.CURVE3,MplPath.CURVE3,MplPath.CLOSEPOLY]), fc="#b95b58", ec="none"))
    ax.add_patch(Polygon([[.41,.61],[.50,.62],[.54,.57],[.49,.54],[.40,.55]], closed=True,fc="#b97859",ec="#8d5745",lw=.4))
    ax.add_patch(Ellipse((.555,.565),.08,.055,angle=-20,fc="#f0b69e",ec="#be8270",lw=.4))
    ax.add_patch(Ellipse((.595,.575),.026,.045,fc="#9d6072",ec="none"))
    ax.add_patch(Ellipse((.46,.505),.026,.055,angle=-12,fc="#a86b63",ec="none")); ax.add_patch(Ellipse((.54,.505),.026,.055,angle=12,fc="#a86b63",ec="none"))
    ax.add_patch(Ellipse((.50,.52),.09,.018,fc="#efad78",ec="none"))
    ax.add_patch(FancyBboxPatch((.44,.39),.12,.10,boxstyle="round,pad=.008,rounding_size=.02",fc="#cf8f82",ec="#9f6c63",lw=.5))
    for yy in [.415,.44,.465]: ax.plot([.46,.54],[yy,yy],color="#efc0ac",lw=2.0,solid_capstyle="round")
    ax.add_patch(Ellipse((.50,.315),.045,.035,fc="#d6a079",ec="none"))
    ax.plot([.42,.42],[.41,.17],color="#c99670",lw=4,alpha=.6); ax.plot([.58,.58],[.41,.17],color="#c99670",lw=4,alpha=.6)

    organ_xy = {"Brain":(.50,.89),"Thyroid":(.50,.795),"Lung":(.455,.70),"Heart":(.50,.65),"Liver":(.44,.575),"Stomach":(.555,.565),"Spleen":(.595,.575),"Pancreas":(.50,.52),"Kidney":(.54,.505),"Large intestine":(.50,.44),"Small intestine":(.50,.445),"Bladder":(.50,.315),"Bone marrow":(.58,.24),"Reproductive":(.50,.275),"Skin":(.31,.63)}
    dd = values.dropna(subset=[value_col]).set_index("organ")[value_col].to_dict()
    finite = [dd[o] for o in organ_xy if o in dd]
    norm = mpl.colors.Normalize(min(finite) if finite else 0, max(finite) if finite else 1); cm = mpl.colormaps[cmap]
    side_y = {"Brain":.92,"Thyroid":.82,"Lung":.74,"Heart":.67,"Liver":.60,"Stomach":.55,"Spleen":.50,"Pancreas":.46,"Kidney":.41,"Large intestine":.36,"Small intestine":.31,"Bladder":.26,"Bone marrow":.20,"Reproductive":.14,"Skin":.08}
    left = {"Brain","Lung","Liver","Pancreas","Large intestine","Bone marrow","Skin"}
    for o, (x,y) in organ_xy.items():
        if o not in dd: continue
        tx = .02 if o in left else .98; elbow = .29 if o in left else .71; ty = side_y[o]
        ax.scatter([x],[y],s=34,c=[cm(norm(dd[o]))],ec="white",lw=.55,zorder=10)
        ax.plot([x,elbow,tx],[y,ty,ty],color="#59636f",lw=.55,zorder=9)
        ax.text(tx,ty,f"{o}  {int(dd[o]):,}",ha="left" if o in left else "right",va="center",fontsize=5.8,weight="bold")
    sm = mpl.cm.ScalarMappable(norm=norm,cmap=cm); cb = plt.colorbar(sm,ax=ax,fraction=.035,pad=.01,shrink=.44); cb.set_label(label_suffix,fontsize=6.2); cb.ax.tick_params(labelsize=5.8)


def build_expression_data():
    R = pd.read_pickle(BASE / "data" / "M3_tissue_rna_matrix.pkl")
    I = pd.read_pickle(BASE / "data" / "M3_ihc_matrix.pkl")
    # Normalize tissue labels for RNA/IHC matching.
    R.columns = [str(x).strip().lower() for x in R.columns]
    I.columns = [str(x).strip().lower() for x in I.columns]
    R = R.groupby(level=0, axis=1).mean(); I = I.groupby(level=0, axis=1).mean()
    pm = pd.read_csv(REL / "canonical_protein_membrane_class_v71.tsv.gz", sep="\t", usecols=["canonical_uniprot_accession","approved_symbol","membrane_class"])
    pm = pm.drop_duplicates("canonical_uniprot_accession").set_index("canonical_uniprot_accession")
    det = R.gt(1.0); breadth = det.sum(axis=1)
    mx = R.max(axis=1); tau = pd.Series(0.0,index=R.index)
    nz = mx.gt(0); tau.loc[nz] = (1 - R.loc[nz].div(mx.loc[nz],axis=0)).sum(axis=1) / max(1,R.shape[1]-1)
    metrics = pd.DataFrame({"target_uniprot_id":R.index,"detected_tissue_count":breadth.values,"tau_specificity":tau.values})
    metrics["membrane_class"] = pm.reindex(R.index).membrane_class.values
    metrics["approved_symbol"] = pm.reindex(R.index).approved_symbol.values
    metrics.to_csv(DATA / "M3_expression_breadth_specificity.tsv",sep="\t",index=False)
    rows=[]
    for t in sorted(set(R.columns)&set(I.columns)):
        q=pd.concat([R[t],I[t]],axis=1,join="inner").dropna();q.columns=["rna","ihc"]
        if len(q)>=100:
            rho,p=spearmanr(np.log1p(q.rna),q.ihc);rows.append((t,len(q),rho,p))
    cor=pd.DataFrame(rows,columns=["tissue","n","spearman_rho","pvalue"]).sort_values("spearman_rho",ascending=False)
    cor["fdr"]=bh_fdr(cor.pvalue.values);cor.to_csv(DATA/"M3_RNA_IHC_concordance_revised.tsv",sep="\t",index=False)
    organ_rows=[]
    for organ in sorted(set(TISSUE_ORGAN.values())):
        ts=[t for t in R.columns if TISSUE_ORGAN.get(t)==organ]
        if ts: organ_rows.append((organ,int(det[ts].any(axis=1).sum())))
    organ=pd.DataFrame(organ_rows,columns=["organ","detected_proteins"]);organ.to_csv(DATA/"M3_body_expression_revised.tsv",sep="\t",index=False)
    return R,I,metrics,cor,organ


def draw_m3():
    R,I,metrics,cor,organ=build_expression_data()
    fig=plt.figure(figsize=(14,9.7))
    fig.text(.03,.969,"M3 | Expression breadth and tissue specificity replace an uninterpretable embedding",fontsize=15,weight="bold",va="top")
    fig.text(.03,.934,"The editable vector anatomy map reports RNA detection; direct axes show how broadly and specifically each protein is expressed.",fontsize=9,color=C["muted"],va="top")
    gs=fig.add_gridspec(2,3,left=.04,right=.985,top=.87,bottom=.09,width_ratios=[1.05,1.25,1.0],hspace=.36,wspace=.32)
    ax=fig.add_subplot(gs[:,0]);draw_human_vector(ax,organ,"detected_proteins","Editable vector map of tissue RNA detection",label_suffix="detected membrane proteins")
    ax=fig.add_subplot(gs[0,1]);panel(ax,"B","Expression breadth vs τ tissue specificity")
    hb=ax.hexbin(metrics.detected_tissue_count,metrics.tau_specificity,gridsize=(28,18),mincnt=1,cmap="YlGnBu",bins="log",alpha=.88)
    cb=fig.colorbar(hb,ax=ax,pad=.02);cb.set_label("log protein density")
    for cl,col in zip(["A","B","C"],[C["blue"],C["green"],C["coral"]]):
        q=metrics[metrics.membrane_class.eq(cl)];ax.scatter(q.detected_tissue_count.median(),q.tau_specificity.median(),s=75,c=col,ec="white",lw=.8,label=f"Class {cl} median",zorder=5)
    ax.set_xlabel("number of HPA tissues with RNA > 1 nTPM");ax.set_ylabel("τ specificity (0 broad → 1 tissue-specific)");ax.legend(fontsize=6.8);clean(ax,None)
    ax=fig.add_subplot(gs[1,1]);panel(ax,"C","Tissues form co-expression blocks without projecting proteins into 2D")
    V=np.log1p(R);var=V.var(axis=1).nlargest(140).index;sel=V.loc[var]
    tissue_var=sel.var(axis=0).nlargest(min(24,sel.shape[1])).index;sel=sel[tissue_var]
    order=leaves_list(linkage(pdist(sel.T,metric="correlation"),method="average"));sel=sel.iloc[:,order]
    z=(sel-sel.mean(axis=0))/sel.std(axis=0).replace(0,np.nan);sns.heatmap(z.T,ax=ax,cmap="vlag",center=0,xticklabels=False,yticklabels=True,cbar_kws={"label":"tissue-wise z score"})
    ax.set_xlabel("140 variable membrane proteins");ax.set_ylabel("")
    ax=fig.add_subplot(gs[:,2]);panel(ax,"D","RNA–IHC concordance is tissue-dependent")
    c=cor.sort_values("spearman_rho").tail(22);yy=np.arange(len(c));ax.hlines(yy,0,c.spearman_rho,color=C["grid"],lw=2);sc=ax.scatter(c.spearman_rho,yy,c=np.log10(c.n),cmap="rocket",s=54,ec="white",lw=.4)
    ax.axvline(0,color=C["muted"],lw=.7);ax.set_yticks(yy,c.tissue);ax.set_xlabel("Spearman ρ: log1p RNA vs IHC detection");clean(ax)
    save(fig,"M3_expression_interpretable_revised","HPA 25.1; mapped + measured normal-tissue records only; RNA detection threshold >1 nTPM; τ is a direct tissue-specificity index.")
    # A fully editable standalone anatomy asset for manual label adjustment.
    f,a=plt.subplots(figsize=(7.2,10));draw_human_vector(a,organ,"detected_proteins","Human tissue RNA detection map",label_suffix="detected membrane proteins",letter="M3A")
    f.savefig(EDIT/"M3A_editable_human_RNA_map.svg",bbox_inches="tight",facecolor="white")
    f.savefig(EDIT/"M3A_editable_human_RNA_map.pdf",bbox_inches="tight",facecolor="white")
    f.savefig(EDIT/"M3A_editable_human_RNA_map.png",dpi=600,bbox_inches="tight",facecolor="white");plt.close(f)
    return R


def draw_m4(R):
    rel=pd.read_csv(REL/"protein_disease_relation_v71.tsv.gz",sep="\t",low_memory=False)
    systems=pd.read_csv(REL/"disease_anatomical_system_v71.tsv.gz",sep="\t",low_memory=False)
    systems["organ"]=systems.anatomical_system_name.str.lower().map(SYSTEM_ORGAN)
    j=rel.merge(systems[["canonical_disease_id","anatomical_system_name","organ"]].drop_duplicates(),on="canonical_disease_id",how="inner")
    j["pair"]=j.target_uniprot_id.astype(str)+"|"+j.canonical_disease_id.astype(str)
    classes=pd.read_csv(REL/"canonical_protein_membrane_class_v71.tsv.gz",sep="\t",usecols=["canonical_uniprot_accession","membrane_class"]).drop_duplicates("canonical_uniprot_accession").set_index("canonical_uniprot_accession")
    R=R.copy();R.columns=[str(x).lower() for x in R.columns];det=R.gt(1.0);breadth=det.sum(axis=1)
    universe=pd.Index(sorted(set(R.index)&set(classes.index)))
    meta=pd.DataFrame(index=universe);meta["class"]=classes.reindex(universe).membrane_class.fillna("U")
    meta["breadth_decile"]=pd.qcut(breadth.reindex(universe).rank(method="first"),10,labels=False,duplicates="drop")
    meta["stratum"]=meta["class"].astype(str)+"|"+meta["breadth_decile"].astype(str)
    strata={s:idx.to_numpy() for s,idx in meta.groupby("stratum").groups.items()}
    organ_sets={}
    for org in set(TISSUE_ORGAN.values()):
        ts=[t for t in R.columns if TISSUE_ORGAN.get(t)==org]
        if ts:organ_sets[org]=set(R.index[det[ts].any(axis=1)])&set(universe)
    rng=np.random.default_rng(42);rows=[];nperm=5000
    for org,q in j.dropna(subset=["organ"]).groupby("organ"):
        dp=set(q.target_uniprot_id)&set(universe);ep=organ_sets.get(org,set())
        if len(dp)<10 or not ep:continue
        obs=len(dp&ep)/len(dp);need=meta.loc[list(dp)].stratum.value_counts();null=np.empty(nperm)
        for i in range(nperm):
            draw=[]
            for s,n in need.items():
                pool=strata[s];draw.extend(rng.choice(pool,size=int(n),replace=len(pool)<int(n)))
            null[i]=len(set(draw)&ep)/len(dp)
        p=(1+int((null>=obs).sum()))/(nperm+1)
        rows.append((org,len(dp),len(ep),obs,float(null.mean()),float(np.quantile(null,.025)),float(np.quantile(null,.975)),obs-float(null.mean()),float((obs-null.mean())/(null.std()+1e-12)),p))
    con=pd.DataFrame(rows,columns=["organ","disease_proteins","expressed_proteins","observed_fraction","matched_null_mean","null_ci_low","null_ci_high","risk_difference","zscore","empirical_p"])
    con["fdr_bh"]=bh_fdr(con.empirical_p.values);con.to_csv(DATA/"M4_matched_disease_expression_concordance.tsv",sep="\t",index=False)
    burden=j.dropna(subset=["organ"]).groupby("organ").pair.nunique().rename("protein_disease_pairs").reset_index();burden.to_csv(DATA/"M4_body_disease_revised.tsv",sep="\t",index=False)
    role=pd.read_csv(REL/"protein_cross_classification_summary_v71.tsv.gz",sep="\t",usecols=["target_uniprot_id","membrane_role_primary_v63"])
    jj=j.merge(role,on="target_uniprot_id",how="left");tab=pd.crosstab(jj.anatomical_system_name,jj.membrane_role_primary_v63)
    tab=tab.loc[tab.sum(1).nlargest(10).index,tab.sum().nlargest(9).index];share=tab.div(tab.sum(1),axis=0);share.to_csv(DATA/"M4_system_role_share_revised.tsv",sep="\t")
    ev=pd.crosstab(j.anatomical_system_name,j.best_evidence_level,values=j.pair,aggfunc=pd.Series.nunique).fillna(0)
    ev=ev.reindex(columns=["medium","high","very_high"],fill_value=0);ev=ev.loc[ev.sum(1).nlargest(11).index];ev.to_csv(DATA/"M4_evidence_levels_revised.tsv",sep="\t")

    fig=plt.figure(figsize=(14,9.7));fig.text(.03,.969,"M4 | Disease anatomy is concordant with normal-tissue expression beyond matched background",fontsize=15,weight="bold",va="top")
    sig=int((con.fdr_bh<.05).sum());fig.text(.03,.934,f"A class- and expression-breadth-matched null controls broad-expression bias; {sig}/{len(con)} mapped organs pass BH-FDR < 0.05.",fontsize=9,color=C["muted"],va="top")
    gs=fig.add_gridspec(2,3,left=.04,right=.985,top=.87,bottom=.09,width_ratios=[1.05,1.18,1.15],hspace=.36,wspace=.34)
    ax=fig.add_subplot(gs[:,0]);draw_human_vector(ax,burden,"protein_disease_pairs","Ontology-mapped disease association burden",cmap="YlOrRd",label_suffix="unique protein–disease pairs")
    ax=fig.add_subplot(gs[:,1]);panel(ax,"B","Observed expression exceeds a matched null in selected organs")
    c=con.sort_values("risk_difference");yy=np.arange(len(c));ax.hlines(yy,c.null_ci_low,c.null_ci_high,color=C["missing"],lw=4,label="matched-null 95% interval")
    ax.scatter(c.matched_null_mean,yy,c="white",ec=C["muted"],s=25,label="null mean",zorder=3)
    colors=np.where(c.fdr_bh<.05,C["red"],C["blue"]);ax.scatter(c.observed_fraction,yy,c=colors,s=50,ec="white",lw=.5,label="observed",zorder=4)
    ax.set_yticks(yy,c.organ);ax.set_xlabel("fraction of organ-disease proteins RNA-detected");ax.legend(fontsize=6.5,loc="lower right");clean(ax,"x")
    ax=fig.add_subplot(gs[0,2]);panel(ax,"C","Disease systems emphasize different membrane roles")
    sns.heatmap(share,ax=ax,cmap="YlGnBu",linewidths=.35,linecolor="white",cbar_kws={"label":"within-system pair share"});ax.set_xlabel("");ax.set_ylabel("");ax.tick_params(axis="x",rotation=40)
    ax=fig.add_subplot(gs[1,2]);panel(ax,"D","Very-high evidence was previously hidden by a label mismatch")
    sns.heatmap(np.log1p(ev),ax=ax,cmap=sns.light_palette(C["red"],as_cmap=True),linewidths=.35,linecolor="white",annot=ev.astype(int),fmt="d",annot_kws={"fontsize":5.7},cbar_kws={"label":"log1p(unique pairs)"})
    ax.set_xticklabels(["Medium","High","Very high"],rotation=0);ax.set_xlabel("best evidence level");ax.set_ylabel("")
    save(fig,"M4_matched_disease_expression_revised","Open Targets 26.06 + UniProtKB disease annotations; MONDO/DO/Uberon multi-label anatomy; HPA 25.1; 5,000 class- and expression-breadth-matched permutations; BH correction.")
    summary={"permutations":nperm,"organs_tested":len(con),"fdr_lt_0_05":sig,"median_risk_difference":float(con.risk_difference.median()),"median_zscore":float(con.zscore.median()),"max_zscore":float(con.zscore.max()),"evidence_level_total_pairs":ev.sum().astype(int).to_dict()}
    (QA/"M4_STATISTICAL_QA.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")


def main():
    draw_m1();draw_m2();R=draw_m3();draw_m4(R)
    readme = """# MemPro M1–M4 scientific revision (2026-08-16)

This revision does not modify the frozen V7.1.1 release tables.

- M1: restores legacy v3.1 `ACT_*` records to assay-level database provenance before overlap analysis.
- M2: displays all 23 HDBSCAN communities plus noise and all five annotation axes.
- M3: replaces the crowded expression UMAP with tissue breadth versus tau specificity; supplies a standalone editable SVG anatomy map.
- M4: replaces the unmatched null with 5,000 membrane-class and expression-breadth matched permutations; fixes `very_high` versus `very high` label mismatch.

All figures are exported as PNG, SVG and PDF. SVG text and vector shapes remain editable.
"""
    (OUT/"README.md").write_text(readme,encoding="utf-8")
    files=[]
    import hashlib
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.name!="SHA256SUMS.tsv":
            h=hashlib.sha256(p.read_bytes()).hexdigest();files.append((str(p.relative_to(OUT)),p.stat().st_size,h))
    pd.DataFrame(files,columns=["file","bytes","sha256"]).to_csv(OUT/"SHA256SUMS.tsv",sep="\t",index=False)
    archive=shutil.make_archive(str(OUT),"zip",root_dir=OUT.parent,base_dir=OUT.name)
    print(json.dumps({"output":str(OUT),"archive":archive,"files":len(files)},ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
