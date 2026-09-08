from __future__ import annotations

import json
import math
import os
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.image as mpimg
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw


ROOT = Path(os.environ["MEMPRO_FIGURE_OUTPUT_ROOT"])
ANALYSIS = ROOT / "02_analysis_data"
SD = ROOT / "03_source_data"
FIG = ROOT / "04_figures"
MAIN = FIG / "main"
SUPP = FIG / "supplementary"
ANATOMY = ROOT / "05_graphics" / "human_anatomy_base_v1.png"
RESULTS = json.loads((ANALYSIS / "analysis_results.json").read_text(encoding="utf-8"))

PALETTE = ["#7b95c6", "#49c2d9", "#a1d8e8", "#67a583", "#a2c986", "#d0e2c0", "#fded95", "#ffc1a6", "#f59c7c", "#f47254", "#c85e62"]
BLUE, CYAN, LTBLUE, GREEN, LTGREEN, PALEGREEN, YELLOW, PEACH, SALMON, ORANGE, RED = PALETTE
DARK = "#26323a"
MID = "#69757d"
LIGHT = "#e8edf0"
ROCKET = sns.color_palette("rocket", as_cmap=True)
VLAAG = sns.color_palette("vlag", as_cmap=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6.4,
    "axes.titlesize": 7.4,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 5.8,
    "ytick.labelsize": 5.8,
    "legend.fontsize": 5.6,
    "figure.titlesize": 11,
    "text.color": DARK,
    "axes.labelcolor": DARK,
    "axes.titlecolor": DARK,
    "axes.edgecolor": "#a8b0b5",
    "axes.linewidth": 0.6,
    "grid.color": "#e2e7ea",
    "grid.linewidth": 0.45,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})
sns.set_theme(style="whitegrid", rc=plt.rcParams)


def read(name: str, directory: Path = SD, **kwargs) -> pd.DataFrame:
    return pd.read_csv(directory / name, sep="\t", compression="infer", **kwargs)


def panel(ax, letter: str, title: str, y: float = 1.0) -> None:
    ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", pad=5, y=y)


def figure_title(fig, title: str, subtitle: str) -> None:
    fig.suptitle(title, x=0.01, y=1.045, ha="left", va="top", fontweight="bold")
    fig.text(0.01, 1.015, subtitle, ha="left", va="top", fontsize=6.8, color=MID)


def source_footer(fig, text: str) -> None:
    fig.text(0.01, -0.022, "Source: " + text, ha="left", va="bottom", fontsize=5.1, color=MID)


def save(fig, stem: str, supplementary: bool = False) -> None:
    base = SUPP if supplementary else MAIN
    png_path = base / "png" / f"{stem}.png"
    svg_path = base / "svg" / f"{stem}.svg"
    pdf_path = base / "pdf" / f"{stem}.pdf"
    for path in [png_path, svg_path, pdf_path]: path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fmt(n: float) -> str:
    if n >= 1_000_000:
        return f"{n/1_000_000:.2f}M"
    if n >= 1_000:
        return f"{n/1_000:.1f}k"
    return f"{int(n)}"


def shorten(value: object, width: int = 27) -> str:
    return textwrap.shorten(str(value).replace("_", " "), width=width, placeholder="…")


def ribbon(ax, x0, y0a, y0b, x1, y1a, y1b, color, alpha=0.35):
    verts = [
        (x0, y0a), (x0 + (x1-x0)*0.45, y0a), (x0 + (x1-x0)*0.55, y1a), (x1, y1a),
        (x1, y1b), (x0 + (x1-x0)*0.55, y1b), (x0 + (x1-x0)*0.45, y0b), (x0, y0b), (x0, y0a),
    ]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
    ax.add_patch(patches.PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha))


def draw_alluvial(ax, stages, edges, colors=None, label_width=20):
    """Draw a compact alluvial from (stage,label) nodes and edge counts."""
    colors = colors or PALETTE
    stage_x = {stage: i / (len(stages)-1) for i, stage in enumerate(stages)}
    totals = defaultdict(float)
    incoming = defaultdict(float)
    outgoing = defaultdict(float)
    for s_stage, s_label, t_stage, t_label, value in edges:
        outgoing[(s_stage, s_label)] += value
        incoming[(t_stage, t_label)] += value
    for node in set(incoming) | set(outgoing):
        totals[node] = max(incoming[node], outgoing[node])
    positions = {}
    for si, stage in enumerate(stages):
        nodes = [(label, value) for (st, label), value in totals.items() if st == stage]
        nodes.sort(key=lambda x: x[1], reverse=True)
        denom = sum(v for _, v in nodes) or 1
        gap = 0.012
        usable = 0.92 - gap * max(len(nodes)-1, 0)
        y = 0.96
        for li, (label, value) in enumerate(nodes):
            h = usable * value / denom
            positions[(stage, label)] = [y-h, y, y-h, y-h]
            color = colors[li % len(colors)] if si == 0 else LIGHT
            ax.add_patch(patches.FancyBboxPatch((stage_x[stage]-0.017, y-h), 0.034, h, boxstyle="round,pad=0.002", facecolor=color, edgecolor="white", linewidth=0.35, zorder=3))
            ha = "right" if si == len(stages)-1 else "left"
            dx = -0.024 if ha == "right" else 0.024
            ax.text(stage_x[stage]+dx, y-h/2, shorten(label, label_width), ha=ha, va="center", fontsize=5.0)
            y -= h + gap
    # Flow offsets are tracked within each node.
    source_offsets = defaultdict(float)
    target_offsets = defaultdict(float)
    stage_totals = {stage: sum(v for (st, _), v in totals.items() if st == stage) or 1 for stage in stages}
    for idx, (ss, sl, ts, tl, value) in enumerate(sorted(edges, key=lambda x: x[4], reverse=True)):
        s0, s1, _, _ = positions[(ss, sl)]
        t0, t1, _, _ = positions[(ts, tl)]
        sh = (s1-s0) * value / max(totals[(ss,sl)], 1)
        th = (t1-t0) * value / max(totals[(ts,tl)], 1)
        sy0 = s0 + source_offsets[(ss,sl)]; sy1 = sy0 + sh
        ty0 = t0 + target_offsets[(ts,tl)]; ty1 = ty0 + th
        source_offsets[(ss,sl)] += sh
        target_offsets[(ts,tl)] += th
        ribbon(ax, stage_x[ss]+0.017, sy0, sy1, stage_x[ts]-0.017, ty0, ty1, colors[idx % len(colors)], alpha=0.22)
    for stage in stages:
        ax.text(stage_x[stage], 1.015, stage.replace("_", " ").title(), ha="center", va="bottom", fontsize=5.6, fontweight="bold")
    ax.set_xlim(-0.12, 1.12); ax.set_ylim(0, 1.06); ax.axis("off")


def draw_upset(ax_bar, ax_matrix, df, set_cols, count_col, top_n=12):
    work = df.sort_values(count_col, ascending=False).head(top_n).reset_index(drop=True)
    x = np.arange(len(work))
    ax_bar.bar(x, work[count_col], color=CYAN, alpha=0.72, width=0.72)
    ax_bar.set_xticks([]); ax_bar.set_ylabel("Count")
    for i, value in enumerate(work[count_col]):
        ax_bar.text(i, value, fmt(value), ha="center", va="bottom", fontsize=5.0, rotation=90, rotation_mode="anchor")
    ax_matrix.set_xlim(-0.5, len(work)-0.5); ax_matrix.set_ylim(-0.5, len(set_cols)-0.5)
    for i, col in enumerate(set_cols):
        y = len(set_cols)-1-i
        ax_matrix.text(-0.75, y, shorten(col, 24), ha="right", va="center", fontsize=5)
        ax_matrix.axhline(y, color="#eef1f3", lw=0.4, zorder=0)
    for j, row in work.iterrows():
        active = []
        for i, col in enumerate(set_cols):
            y = len(set_cols)-1-i
            on = int(row[col]) == 1
            ax_matrix.scatter(j, y, s=14, color=DARK if on else "#d8dde0", zorder=2)
            if on: active.append(y)
        if len(active) > 1:
            ax_matrix.plot([j,j], [min(active), max(active)], color=DARK, lw=0.8, zorder=1)
    ax_matrix.axis("off")


def figure1():
    protein = read("protein_analysis.tsv.gz", ANALYSIS)
    coverage = read("F1D_module_protein_coverage.tsv")
    flow = read("F1C_source_evidence_tier_flow.tsv")
    flow["source_group"] = flow["source_database"].where(flow["source_database"].isin(flow.groupby("source_database")["evidence_record_count"].sum().nlargest(7).index), "Other sources")
    type_map = {
        "direct_quantitative_binding": "Direct quantitative",
        "pubchem_direct_quantitative_binding": "Direct quantitative",
        "pubchem_target_specific_quantitative_activity": "Target-specific activity",
        "functional_or_bioactivity_measurement": "Functional/pharmacology",
        "binding_assay_activity": "Functional/pharmacology",
        "experimental_structure_residue_contact": "Structural contact",
        "compound_specific_structure_context": "Structural context",
        "compound_specific_structure": "Structural context",
    }
    flow["type_group"] = flow["evidence_type"].map(type_map).fillna("Curated/other")
    flow2 = flow.groupby(["source_group", "type_group", "evidence_tier"], as_index=False)["evidence_record_count"].sum()
    edges = []
    for r in flow2.itertuples(index=False):
        edges.append(("source", r.source_group, "evidence_type", r.type_group, r.evidence_record_count))
    for (typ,tier), g in flow2.groupby(["type_group","evidence_tier"]):
        edges.append(("evidence_type", typ, "tier", tier, int(g["evidence_record_count"].sum())))

    fig = plt.figure(figsize=(10.0, 8.8), constrained_layout=True)
    gs = fig.add_gridspec(3, 4, height_ratios=[0.9, 1.8, 1.25])
    figure_title(fig, "Figure 1 | From heterogeneous records to a traceable membrane-target resource", "Identity, entity granularity and evidence provenance are explicit at release time")
    ax = fig.add_subplot(gs[0, :]); panel(ax, "a", "Integration and release contract"); ax.axis("off")
    boxes = [
        (0.02, "Protein · structure\ninteraction · disease\nexpression sources", LTBLUE),
        (0.27, "Canonical identity\nparent/form · isoform\nMONDO exact mapping", PALEGREEN),
        (0.52, "Deduplication · lineage\nE / BE tiers\nmissingness semantics", YELLOW),
        (0.77, "V7.2 frozen entities\nsource IDs · manifest\nreproducible analysis", PEACH),
    ]
    for i, (x, txt, col) in enumerate(boxes):
        ax.add_patch(patches.FancyBboxPatch((x, .18), .20, .62, boxstyle="round,pad=.012,rounding_size=.018", fc=col, ec="white", lw=.8))
        ax.text(x+.10, .49, txt, ha="center", va="center", fontsize=6.2)
        if i < 3:
            ax.annotate("", xy=(boxes[i+1][0]-.01,.49), xytext=(x+.205,.49), arrowprops=dict(arrowstyle="->", color=MID, lw=1))
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax = fig.add_subplot(gs[1, :]); panel(ax, "b", "Source → evidence modality → BE tier"); draw_alluvial(ax, ["source","evidence_type","tier"], edges, label_width=22)
    ax = fig.add_subplot(gs[2, :2]); panel(ax, "c", "Protein-level module coverage")
    cov = coverage[~coverage["module"].eq("formal_protein")].sort_values("coverage_fraction")
    y = np.arange(len(cov)); ax.hlines(y, 0, cov["coverage_fraction"], color="#cbd2d6", lw=1)
    ax.scatter(cov["coverage_fraction"], y, s=38, c=[GREEN, CYAN, BLUE, ORANGE, RED][:len(cov)], alpha=.82)
    ax.set_yticks(y, [shorten(x,25) for x in cov["module"]]); ax.set_xlim(0,1.03); ax.set_xlabel("Fraction of 7,800 formal proteins")
    for yy, frac, n in zip(y, cov["coverage_fraction"], cov["protein_count"]): ax.text(frac+.015, yy, f"{frac:.0%} · n={int(n):,}", va="center", fontsize=5)
    ax = fig.add_subplot(gs[2, 2:]); panel(ax, "d", "Release entities and analysis denominators"); ax.axis("off")
    cards = [("Formal proteins",7800),("Compound registry",646670),("Interaction-linked compounds",240055),("Positive pairs",529168),("Public evidence",942455),("Diseases",3739),("Expression measurements",3046789),("Site assertions",55610)]
    for i,(label,n) in enumerate(cards):
        x=(i%2)*.50; y0=.88-(i//2)*.23
        ax.text(x,y0,fmt(n),fontsize=11,fontweight="bold",color=[BLUE,GREEN,ORANGE,RED][(i//2)%4])
        ax.text(x,y0-.08,label,fontsize=5.7,color=MID)
    source_footer(fig, "MemPro V7.2 frozen release; formal public evidence only. Counts are entities or records, not independent experiments.")
    save(fig, "Figure1_resource_architecture")


def figure2():
    protein = read("protein_analysis.tsv.gz", ANALYSIS)
    rf = read("F2B_role_function_cell_statistics.tsv")
    flow_df = read("F2C_five_axis_alluvial_edges.tsv")
    tau = read("F2D_tissue_specificity_tau.tsv")
    complete = read("F2E_classification_completeness.tsv")

    fig = plt.figure(figsize=(10.0, 11.2))
    fig.subplots_adjust(left=.10, right=.965, top=.91, bottom=.065, hspace=.82, wspace=.95)
    gs = fig.add_gridspec(3, 4, height_ratios=[1.05, 1.72, 1.20])
    figure_title(fig, "Figure 2 | A five-axis framework reveals structured membrane-protein organization", "Primary membrane roles align strongly with molecular function while preserving multi-functional annotations")
    ax = fig.add_subplot(gs[0,:1]); panel(ax,"a","Class × evidence", y=1.05)
    tab = pd.crosstab(protein["membrane_class_v7"], protein["membrane_evidence_level_v7"]).reindex(index=["A","B","C"], columns=["E1","E2","E3"], fill_value=0)
    total=tab.values.sum(); x0=0
    for i,cls in enumerate(tab.index):
        w=tab.loc[cls].sum()/total; y0=0
        for j,ev in enumerate(tab.columns):
            h=tab.loc[cls,ev]/tab.loc[cls].sum() if tab.loc[cls].sum() else 0
            ax.add_patch(patches.Rectangle((x0,y0),w,h,fc=[BLUE,GREEN,ORANGE][i],alpha=[.88,.62,.35][j],ec="white",lw=.7))
            if tab.loc[cls,ev] > 150: ax.text(x0+w/2,y0+h/2,f"{cls}/{ev}\n{tab.loc[cls,ev]:,}",ha="center",va="center",fontsize=5.0)
            y0+=h
        x0+=w
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel("Width: class abundance\nHeight: evidence composition", fontsize=5.2)

    ax = fig.add_subplot(gs[0,1:]); panel(ax,"b","Role × molecular-function enrichment", y=1.05)
    mat = rf.pivot(index="primary_membrane_role",columns="primary_molecular_function",values="standardized_residual")
    mat = mat.loc[mat.abs().max(axis=1).sort_values(ascending=False).index]
    sns.heatmap(mat,ax=ax,cmap=VLAAG,center=0,vmin=-8,vmax=8,linewidths=.25,linecolor="white",cbar_kws={"label":"Standardized residual","shrink":.65})
    ax.set_xlabel("Primary molecular function");ax.set_ylabel("")
    ax.set_xticks(np.arange(len(mat.columns))+.5,[shorten(x.replace('_',' '),15) for x in mat.columns],rotation=38,ha="right",rotation_mode="anchor",fontsize=5.0)
    ax.set_yticks(np.arange(len(mat.index))+.5,[shorten(x.replace('_',' '),17) for x in mat.index],rotation=0,fontsize=5.0)
    v=RESULTS["protein"]["role_function"]["cramers_v_bias_corrected"]
    ax.text(1.0,1.03,f"χ² association · bias-corrected Cramér's V={v:.2f}",transform=ax.transAxes,ha="right",fontsize=5.3,color=MID)

    ax = fig.add_subplot(gs[1,:]); panel(ax,"c","Primary-category flow across the classification system", y=1.03)
    edges=[(r.source_axis,r.source_label,r.target_axis,r.target_label,int(r.protein_count)) for r in flow_df.itertuples(index=False)]
    draw_alluvial(ax,["membrane_class","membrane_role","molecular_function","biological_process"],edges,label_width=18)

    ax = fig.add_subplot(gs[2,:2]); panel(ax,"d","Tissue-RNA specificity by membrane role", y=1.04)
    tau["tau_consensus_tissue_rna"]=pd.to_numeric(tau["tau_consensus_tissue_rna"],errors="coerce")
    counts=tau["primary_membrane_role"].value_counts(); roles=list(counts[counts>=30].index)
    plot=tau[tau["primary_membrane_role"].isin(roles)].copy(); order=list(plot.groupby("primary_membrane_role")["tau_consensus_tissue_rna"].median().sort_values().index)
    sns.violinplot(data=plot,x="tau_consensus_tissue_rna",y="primary_membrane_role",order=order,inner=None,cut=0,linewidth=.4,color=LTBLUE,alpha=.42,ax=ax)
    sns.boxplot(data=plot,x="tau_consensus_tissue_rna",y="primary_membrane_role",order=order,width=.18,showfliers=False,boxprops={"facecolor":"white","alpha":.8},medianprops={"color":RED,"lw":1},ax=ax)
    ax.set_xlim(0,1);ax.set_xlabel("Tau (0 broad expression → 1 tissue-specific)");ax.set_ylabel("")
    ax.tick_params(axis="y",labelsize=5.0)
    ax.text(.99,.02,f"Kruskal–Wallis P={RESULTS['tau']['kruskal_p']:.1e}",transform=ax.transAxes,ha="right",fontsize=5.3,color=MID)

    ax = fig.add_subplot(gs[2,2:]); panel(ax,"e","Axis completeness", y=1.04)
    complete=complete.sort_values("informative")
    y=np.arange(len(complete)); frac=complete["informative"]/complete["denominator"]
    ax.barh(y,frac,color=GREEN,alpha=.65);ax.barh(y,1-frac,left=frac,color=LIGHT)
    ax.set_yticks(y,[shorten(x.replace('_',' '),13) for x in complete["classification_axis"]],fontsize=5.0);ax.set_xlim(0,1);ax.set_xlabel("Fraction")
    for yy,f,n in zip(y,frac,complete["placeholder_or_missing"]):ax.text(min(f+.02,.88),yy,f"{f:.0%} · {int(n):,} unresolved",va="center",fontsize=5.0)
    source_footer(fig,"V7.2 protein master; primary five-axis annotations; HPA consensus tissue RNA. Cell-wise FDR statistics are supplied as source data.")
    save(fig,"Figure2_membrane_functional_architecture")


def figure3():
    scaff=read("F3A_top_scaffolds.tsv")
    rank=read("F3A_scaffold_rank_abundance.tsv")
    props=read("F3B_physicochemical_by_status.tsv")
    entropy=read("F3C_compound_target_entropy.tsv")
    pref=read("F3D_structure_role_cell_statistics.tsv")
    sim=read("F3E_chemical_vs_target_similarity.tsv")
    universe=read("F3_compound_universe_audit.tsv")

    fig=plt.figure(figsize=(10.0,10.0),constrained_layout=True);gs=fig.add_gridspec(3,4,height_ratios=[1.35,1.2,1.25])
    figure_title(fig,"Figure 3 | Chemical diversity is broad, whereas target-family preference is modest","Scaffold diversity, biological status and polypharmacology are separated from the registry-scale compound universe")
    ax=fig.add_subplot(gs[0,:2]);panel(ax,"a","Scaffold rank-abundance")
    rr=rank.head(3000);ax.loglog(rr["rank"],rr["compound_count"],color=BLUE,lw=1.2);ax.fill_between(rr["rank"],rr["compound_count"],1,color=LTBLUE,alpha=.2);ax.set_xlabel("Scaffold rank");ax.set_ylabel("Interaction-linked compounds per scaffold")
    ax.text(.98,.93,f"{len(rank):,} unique scaffold groups",transform=ax.transAxes,ha="right",fontsize=5.4)
    # Representative non-acyclic structures as deterministic inset tiles.
    shown=0
    for _,r in scaff.iterrows():
        smi=r["murcko_scaffold_smiles"]
        if smi=="ACYCLIC_NO_MURCKO_SCAFFOLD" or not isinstance(smi,str):continue
        mol=Chem.MolFromSmiles(smi)
        if mol is None:continue
        x=.05+(shown%4)*.235;y=.05+(shown//4)*.22
        ia=ax.inset_axes([x,y,.20,.18]);img=Draw.MolToImage(mol,size=(180,120));ia.imshow(img);ia.axis("off");ia.text(.5,-.02,f"n={int(r['compound_count']):,}",transform=ia.transAxes,ha="center",va="top",fontsize=5.0)
        shown+=1
        if shown>=8:break

    ax=fig.add_subplot(gs[0,2:]);panel(ax,"b","Physicochemical profiles")
    props["value"]=pd.to_numeric(props["value"],errors="coerce")
    # Use robust within-property scaling so four properties share one compact panel.
    q=props.groupby("property")["value"].quantile([.25,.5,.75]).unstack(); q.columns=["q25","median","q75"]
    props=props.join(q,on="property");props["robust_z"]=(props["value"]-props["median"])/(props["q75"]-props["q25"]).replace(0,np.nan)
    props=props[props["robust_z"].between(-3,3)]
    order=["All interaction-linked","Approved drug","Clinical candidate","Endogenous ligand","Natural product","Chemical probe"]
    sns.violinplot(data=props,x="robust_z",y="biological_status",hue="property",order=order,cut=0,inner=None,linewidth=.3,palette=[BLUE,GREEN,ORANGE,RED],alpha=.28,ax=ax)
    ax.set_xlabel("Robust z-score within property");ax.set_ylabel("");ax.legend(title="Property",loc="lower right",fontsize=5.0,title_fontsize=5)

    ax=fig.add_subplot(gs[1,:2]);panel(ax,"c","Target breadth × role entropy")
    hb=ax.hexbin(np.log10(pd.to_numeric(entropy["target_count"])+1),pd.to_numeric(entropy["target_role_entropy_normalized"]),gridsize=45,mincnt=1,cmap=ROCKET,bins="log",alpha=.82)
    fig.colorbar(hb,ax=ax,label="log10(compound count)",shrink=.75);ax.set_xlabel("log10(target count + 1)");ax.set_ylabel("Normalized membrane-role entropy")
    ax.text(.02,.93,"High right: broad, cross-role polypharmacology",transform=ax.transAxes,fontsize=5.2,color=MID)

    ax=fig.add_subplot(gs[1,2:]);panel(ax,"d","Structure-class preference")
    show=pref[(pref["observed"]>=20)&(pd.to_numeric(pref["fisher_bh_fdr"],errors="coerce")<.05)].copy()
    rows=list(show.groupby("broad_structure_class")["observed"].sum().nlargest(8).index);cols=list(show.groupby("membrane_role")["observed"].sum().nlargest(9).index);show=show[show["broad_structure_class"].isin(rows)&show["membrane_role"].isin(cols)]
    for r in show.itertuples(index=False):
        x=cols.index(r.membrane_role);y=rows.index(r.broad_structure_class);val=np.clip(float(r.log2_odds_ratio),-2,2);size=12+55*min(abs(val)/2,1);ax.scatter(x,y,s=size,c=[val],cmap=VLAAG,vmin=-2,vmax=2,alpha=.76,edgecolors="white",lw=.3)
    ax.set_xticks(range(len(cols)),[shorten(x,15) for x in cols],rotation=45,ha="right");ax.set_yticks(range(len(rows)),[shorten(x,19) for x in rows]);ax.set_xlim(-.6,len(cols)-.4);ax.set_ylim(-.6,len(rows)-.4);ax.invert_yaxis();ax.grid(False)
    v=RESULTS["compound_pair_evidence"]["class_role"]["cramers_v_bias_corrected"]
    ax.text(.99,1.02,f"Overall association is weak: Cramér's V={v:.3f}",transform=ax.transAxes,ha="right",fontsize=5.2,color=MID)

    ax=fig.add_subplot(gs[2,:3]);panel(ax,"e","Chemical similarity × target overlap")
    sim["tanimoto"]=pd.to_numeric(sim["tanimoto"]);sim["target_jaccard"]=pd.to_numeric(sim["target_jaccard"])
    hb=ax.hexbin(sim["tanimoto"],sim["target_jaccard"],gridsize=42,mincnt=1,cmap=ROCKET,bins="log",alpha=.82);fig.colorbar(hb,ax=ax,label="log10(pair count)",shrink=.75);ax.set_xlabel("Morgan/Tanimoto chemical similarity");ax.set_ylabel("Target-set Jaccard")
    ax.axvline(.65,color=MID,ls="--",lw=.7);ax.text(.66,.95,"Butina threshold",transform=ax.get_xaxis_transform(),fontsize=5.0,color=MID,va="top")

    ax=fig.add_subplot(gs[2,3:]);panel(ax,"f","Compound universes");ax.axis("off")
    u=dict(zip(universe["metric"],universe["count"]));vals=[("Registry",u["compound_registry"],BLUE),("Formal-pair linked",u["interaction_linked_compounds"],GREEN),("No formal pair",u["registry_without_formal_pair"],LIGHT)]
    y=.82
    for label,n,c in vals:
        ax.text(.05,y,fmt(n),fontsize=12,fontweight="bold",color=c if c!=LIGHT else MID);ax.text(.05,y-.10,label,fontsize=5.8,color=MID);y-=.27
    source_footer(fig,"V7.2 canonical compound registry and formal pair table. Scaffold/UMAP analyses use deterministic interaction-linked subsets; UMAP is supplementary only.")
    save(fig,"Figure3_chemical_target_landscape")


def figure4():
    upset=read("F4A_source_upset.tsv");jac=read("F4B_source_jaccard.tsv");dbexp=read("F4C_database_vs_putative_experiment.tsv");lor=read("F4D_evidence_lorenz.tsv");ecdf=read("F4E_site_mapping_ecdf.tsv");side=read("F4F_membrane_side.tsv")
    set_cols=[c for c in upset.columns if c not in {"intersection","pair_count"}]
    fig=plt.figure(figsize=(10.0,10.0),constrained_layout=True);gs=fig.add_gridspec(4,4,height_ratios=[1.2,.65,1.25,1.15])
    figure_title(fig,"Figure 4 | Database coverage and evidence independence are distinct dimensions","Pair confidence reflects evidence lineage, modality and structure—not source count alone")
    axb=fig.add_subplot(gs[0,:2]);axm=fig.add_subplot(gs[1,:2]);panel(axb,"a","Major-source intersections");draw_upset(axb,axm,upset,set_cols,"pair_count",top_n=10)
    ax=fig.add_subplot(gs[0:2,2:]);panel(ax,"b","Pair overlap between contributing databases")
    mat=jac.pivot(index="source_a",columns="source_b",values="jaccard").reindex(index=set_cols,columns=set_cols)
    sns.heatmap(mat,ax=ax,cmap=sns.light_palette(BLUE,as_cmap=True),vmin=0,vmax=max(.01,float(np.nanmax(mat.values))),annot=True,fmt=".2f",annot_kws={"fontsize":5.0},square=True,cbar_kws={"label":"Jaccard","shrink":.65});ax.set_xlabel("");ax.set_ylabel("");ax.set_xticks(np.arange(len(mat.columns))+.5,[shorten(x,16) for x in mat.columns],rotation=45,ha="right",rotation_mode="anchor");ax.set_yticks(np.arange(len(mat.index))+.5,[shorten(x,18) for x in mat.index],rotation=0)
    ax=fig.add_subplot(gs[2,:2]);panel(ax,"c","Databases vs putative experiments")
    dbexp[["distinct_database_count","independent_experiment_count","pair_count"]]=dbexp[["distinct_database_count","independent_experiment_count","pair_count"]].apply(pd.to_numeric)
    ax.scatter(dbexp["distinct_database_count"],np.log10(dbexp["independent_experiment_count"]+1),s=8+50*np.log10(dbexp["pair_count"]+1),c=np.log10(dbexp["pair_count"]+1),cmap=ROCKET,alpha=.65,edgecolors="none")
    ax.set_xlabel("Distinct contributing database count");ax.set_ylabel("log10(putative experiment count + 1)");ax.plot([1,9],np.log10(np.array([1,9])+1),ls="--",lw=.6,color=MID)
    ax=fig.add_subplot(gs[2,2:]);panel(ax,"d","Evidence concentration")
    ax.plot(lor["cumulative_pair_fraction"],lor["cumulative_evidence_fraction"],color=RED,lw=1.4,label=f"Observed · Gini={RESULTS['compound_pair_evidence']['pair_evidence_gini']:.2f}");ax.plot([0,1],[0,1],color=MID,ls="--",lw=.7,label="Equality");ax.fill_between(lor["cumulative_pair_fraction"],lor["cumulative_evidence_fraction"],lor["cumulative_pair_fraction"],color=PEACH,alpha=.25);ax.set_xlabel("Cumulative fraction of pairs");ax.set_ylabel("Cumulative fraction of evidence");ax.legend(loc="upper left")
    ax=fig.add_subplot(gs[3,:2]);panel(ax,"e","Residue-coordinate mapping by site tier")
    ecdf["mapping_fraction"]=pd.to_numeric(ecdf["mapping_fraction"]);ecdf["ecdf"]=pd.to_numeric(ecdf["ecdf"])
    for i,(tier,g) in enumerate(ecdf.groupby("site_tier")):
        ax.plot(g["mapping_fraction"],g["ecdf"],lw=1.2,color=PALETTE[i*3%len(PALETTE)],label=shorten(tier,22))
    ax.set_xlabel("Mapped source residues / requested residues");ax.set_ylabel("ECDF");ax.set_xlim(0,1);ax.legend(loc="lower right")
    ax=fig.add_subplot(gs[3,2:]);panel(ax,"f","Membrane-side resolution remains bounded")
    side=side.sort_values("site_count",ascending=False);vals=side["site_count"].to_numpy();labels=side["membrane_side"].tolist();colors=["#cfd5d9" if x=="unknown" else PALETTE[(i+2)%len(PALETTE)] for i,x in enumerate(labels)]
    ax.pie(vals,colors=colors,startangle=90,wedgeprops={"width":.38,"edgecolor":"white"});ax.text(0,0,f"{int(vals.sum()):,}\nsites",ha="center",va="center",fontweight="bold");ax.legend([f"{shorten(l,22)} · {int(v):,}" for l,v in zip(labels,vals)],loc="center left",bbox_to_anchor=(.92,.5),fontsize=5.0)
    source_footer(fig,"V7.2 pair, evidence-lineage and interaction-site tables. 'Independent experiment' is interpreted as a putative lineage key, not a manually audited publication experiment.")
    save(fig,"Figure4_evidence_architecture")


SYSTEM_POS={
    "nervous system":(.50,.86),"entire sense organ system":(.44,.82),"respiratory system":(.48,.67),"cardiovascular system":(.52,.62),"circulatory system":(.56,.61),"vascular system":(.60,.58),"digestive system":(.50,.51),"alimentary part of gastrointestinal system":(.52,.45),"excretory system":(.44,.46),"reproductive system":(.50,.36),"hematopoietic system":(.61,.28),"lymphoid system":(.60,.55),"musculoskeletal system":(.38,.29),"integumental system":(.27,.55),"exocrine system":(.55,.50),"lacrimal apparatus":(.40,.81),"neuroendocrine system":(.50,.72),
}


def figure5():
    anatomy=read("F5A_anatomy_fractional_counts.tsv");heat=read("F5B_therapeutic_anatomy_cell_statistics.tsv");up=read("F5C_disease_channel_upset.tsv");levels=read("F5D_disease_evidence_levels.tsv");edges=read("F5E_disease_shared_target_network.tsv");coverage=read("F5_disease_mapping_coverage.tsv")
    fig=plt.figure(figsize=(10.0,10.0),constrained_layout=True);gs=fig.add_gridspec(3,4,height_ratios=[1.75,1.25,1.2])
    figure_title(fig,"Figure 5 | Ontology-derived disease and anatomical context", "Multi-label anatomy, evidence channels and shared membrane targets are retained without forcing one-disease–one-organ assignments")
    ax=fig.add_subplot(gs[0,:2]);panel(ax,"a","Multi-label disease anatomy (fractional counts)");img=mpimg.imread(ANATOMY);ax.imshow(img);ax.axis("off")
    anatomy=anatomy[anatomy["anatomical_system"].ne("Unmapped")].copy();mx=anatomy["fractional_disease_count"].max()
    for r in anatomy.itertuples(index=False):
        pos=SYSTEM_POS.get(str(r.anatomical_system).lower())
        if not pos:continue
        x=pos[0]*img.shape[1];y=(1-pos[1])*img.shape[0];size=16+170*math.sqrt(r.fractional_disease_count/mx)
        ax.scatter(x,y,s=size,c=[RED],alpha=.35,edgecolors="white",lw=.5)
    ax.text(.01,.01,"Bubble area scales with fractional disease count\nRepresentative anatomical-system position",transform=ax.transAxes,fontsize=5,color=MID)
    ax=fig.add_subplot(gs[0,2:]);panel(ax,"b","Therapeutic area × anatomy enrichment")
    show=heat[(heat["observed"]>=8)&(pd.to_numeric(heat["fisher_bh_fdr"],errors="coerce")<.05)].copy();rows=list(show.groupby("therapeutic_area")["observed"].sum().nlargest(12).index);cols=list(show.groupby("anatomical_system")["observed"].sum().nlargest(11).index);show=show[show["therapeutic_area"].isin(rows)&show["anatomical_system"].isin(cols)];mat=show.pivot(index="therapeutic_area",columns="anatomical_system",values="standardized_residual").reindex(index=rows,columns=cols).fillna(0)
    sns.heatmap(mat,ax=ax,cmap=VLAAG,center=0,vmin=-8,vmax=8,linewidths=.25,linecolor="white",cbar_kws={"label":"Standardized residual","shrink":.58});ax.set_xlabel("");ax.set_ylabel("");ax.set_xticks(np.arange(len(mat.columns))+.5,[shorten(x,15) for x in mat.columns],rotation=45,ha="right");ax.set_yticks(np.arange(len(mat.index))+.5,[shorten(x,24) for x in mat.index],rotation=0)
    ax=fig.add_subplot(gs[1,:2]);panel(ax,"c","Evidence-channel intersections");matrix=ax.inset_axes([0,-.42,1,.38]);set_cols=[c for c in up.columns if c not in {"intersection","relation_count"}];draw_upset(ax,matrix,up,set_cols,"relation_count",top_n=10)
    ax=fig.add_subplot(gs[1,2:]);panel(ax,"d","Evidence and ontology coverage");ax.axis("off")
    colors={"very_high":RED,"high":ORANGE,"medium":BLUE};total=levels["relation_count"].sum();start=0
    for r in levels.itertuples(index=False):
        frac=r.relation_count/total;ax.add_patch(patches.Wedge((.28,.55),.25,start*360,(start+frac)*360,width=.10,facecolor=colors.get(r.evidence_level,MID),edgecolor="white"));start+=frac
    ax.text(.28,.55,f"{int(total):,}\nrelations",ha="center",va="center",fontweight="bold")
    y=.82
    for r in levels.sort_values("relation_count",ascending=False).itertuples(index=False):ax.scatter(.62,y,s=26,color=colors.get(r.evidence_level,MID));ax.text(.68,y,f"{r.evidence_level.replace('_',' ').title()} · {int(r.relation_count):,}",va="center");y-=.14
    cov=dict(zip(coverage["metric"],coverage["count"]));ax.text(.58,.28,f"Anatomy mapped  {int(cov['diseases_with_anatomy']):,} / {int(cov['formal_diseases']):,}\nUnmapped            {int(cov['diseases_without_anatomy']):,}\nTherapeutic area  {int(cov['diseases_with_therapeutic_area']):,}",fontsize=5.4,linespacing=1.5)
    ax=fig.add_subplot(gs[2,:]);panel(ax,"e","Selected shared-target disease communities")
    disease_names=read("disease_relation_analysis.tsv.gz",ANALYSIS).drop_duplicates("canonical_disease_id").set_index("canonical_disease_id")["canonical_disease_name"].to_dict()
    G=nx.Graph();top=edges.head(120)
    for r in top.itertuples(index=False):G.add_edge(r.disease_a,r.disease_b,weight=r.jaccard,shared=r.shared_targets)
    if G.number_of_nodes()>0:
        comm=list(nx.community.greedy_modularity_communities(G,weight="weight"));cmap={n:i for i,c in enumerate(comm) for n in c};pos=nx.spring_layout(G,seed=42,weight="weight",k=.8/math.sqrt(max(G.number_of_nodes(),1)));deg=dict(G.degree());nx.draw_networkx_edges(G,pos,ax=ax,width=[.3+1.4*G[u][v]["weight"] for u,v in G.edges()],alpha=.22,edge_color=MID);nx.draw_networkx_nodes(G,pos,ax=ax,node_size=[10+6*deg[n] for n in G.nodes()],node_color=[PALETTE[cmap[n]%len(PALETTE)] for n in G.nodes()],alpha=.78,linewidths=.2,edgecolors="white")
        label_nodes=sorted(G.nodes(),key=lambda n:deg[n],reverse=True)[:5]
        offsets=[(.04,.05),(.04,-.06),(-.05,.05),(-.05,-.06),(.05,.01)]
        for n,(dx,dy) in zip(label_nodes,offsets):
            ax.annotate(shorten(disease_names.get(n,n),20),xy=pos[n],xytext=(pos[n][0]+dx,pos[n][1]+dy),
                        fontsize=5.0,ha="left" if dx>0 else "right",va="center",
                        arrowprops={"arrowstyle":"-","lw":.35,"color":MID},
                        bbox={"boxstyle":"round,pad=.12","fc":"white","ec":"none","alpha":.78})
    ax.axis("off");ax.text(.99,.02,"Edges: ≥2 shared targets and Jaccard ≥0.20; top 120 shown",transform=ax.transAxes,ha="right",fontsize=5.0,color=MID)
    source_footer(fig,"V7.2 exact canonical disease relations; ontology-derived multi-label anatomy and Open Targets therapeutic areas. Historical out-of-release ontology rows are excluded.")
    save(fig,"Figure5_disease_anatomy_context")


def figure6():
    landscape=read("F6A_integrated_target_landscape.tsv");null=read("F6B_concordance_permutation_null.tsv");candidates=read("F6D_hypothesis_generation_candidates.tsv");role_ta=read("F5B_therapeutic_anatomy_cell_statistics.tsv")
    for col in ["log10_chemical_pair_count_plus1","log10_disease_evidence_index_plus1","mean_anatomical_expression_concordance","coordinate_ready_site_count","hypothesis_generation_candidate"]: landscape[col]=pd.to_numeric(landscape[col],errors="coerce").fillna(0)
    fig=plt.figure(figsize=(10.0,9.6),constrained_layout=True);gs=fig.add_gridspec(3,4,height_ratios=[1.7,1.15,1.05])
    figure_title(fig,"Figure 6 | Integrated evidence exposes chemically underexplored disease-relevant targets","Transparent dimensions—not an opaque score—support hypothesis generation")
    ax=fig.add_subplot(gs[0,:3]);panel(ax,"a","Disease evidence × chemical coverage")
    s=8+55*np.clip(landscape["mean_anatomical_expression_concordance"],0,1);c=np.log10(landscape["coordinate_ready_site_count"]+1);sc=ax.scatter(landscape["log10_chemical_pair_count_plus1"],landscape["log10_disease_evidence_index_plus1"],s=s,c=c,cmap=ROCKET,alpha=.32,edgecolors="none",rasterized=True);fig.colorbar(sc,ax=ax,label="log10(coordinate-ready sites + 1)",shrink=.72)
    cand=landscape[landscape["hypothesis_generation_candidate"].eq(1)];ax.scatter(cand["log10_chemical_pair_count_plus1"],cand["log10_disease_evidence_index_plus1"],s=58,facecolors="none",edgecolors=RED,lw=.9,label="Hypothesis-generation subset")
    label_view=cand.nlargest(4,"transparent_disease_evidence_index")
    label_offsets=[(5,7),(6,-8),(-7,7),(-8,-8)]
    for r,off in zip(label_view.itertuples(index=False),label_offsets):
        ax.annotate(shorten(r.approved_symbol,10),(r.log10_chemical_pair_count_plus1,r.log10_disease_evidence_index_plus1),
                    xytext=off,textcoords="offset points",fontsize=5.0,ha="left" if off[0]>0 else "right")
    ax.set_xlabel("log10(formal protein–compound pairs + 1)");ax.set_ylabel("log10(transparent disease-evidence index + 1)");ax.legend(loc="lower right")
    ax=fig.add_subplot(gs[0,3:]);panel(ax,"b","Visual encoding");ax.axis("off");items=[("X · Chemical coverage","pair count"),("Y · Disease support","3×very high + 2×high + medium"),("SIZE · Anatomical concordance","normal-tissue RNA only"),("COLOUR · Structural readiness","coordinate-ready sites")];y=.88
    for label,desc in items:ax.text(.03,y,label,fontweight="bold",fontsize=5.5,color=BLUE);ax.text(.03,y-.08,desc,fontsize=5.0,color=MID);y-=.22
    ax=fig.add_subplot(gs[1,:2]);panel(ax,"c","Degree-preserving concordance null")
    vals=pd.to_numeric(null["null_mean_concordance"]);sns.histplot(vals,bins=18,color=LTBLUE,edgecolor="white",ax=ax);obs=RESULTS["anatomical_expression_concordance"]["observed_mean"];ax.axvline(obs,color=RED,lw=1.5,label=f"Observed mean={obs:.3f}");ax.axvline(vals.mean(),color=BLUE,ls="--",lw=1,label=f"Null mean={vals.mean():.3f}");ax.set_xlabel("Mean anatomical expression concordance");ax.legend();ax.text(.98,.80,f"z={RESULTS['anatomical_expression_concordance']['z_score']:.1f}\nempirical P={RESULTS['anatomical_expression_concordance']['empirical_p']:.4f}",transform=ax.transAxes,ha="right",fontsize=5.4)
    ax=fig.add_subplot(gs[1,2:]);panel(ax,"d","Bounded candidate subset");ax.axis("off")
    if len(candidates):
        view=candidates.sort_values(["transparent_disease_evidence_index","chemical_pair_count"],ascending=[False,True]).head(9).copy()
        view["approved_symbol"]=view["approved_symbol"].fillna(view["target_uniprot_id"])
        headers=["Target","Disease\nsupport","Formal\npairs","Ready\nsites","Anatomy\nmatch"];xs=[.04,.30,.51,.68,.84]
        for x,h in zip(xs,headers):ax.text(x,.89,h,transform=ax.transAxes,fontsize=5.0,fontweight="bold",color=MID,va="top")
        for i,r in enumerate(view.itertuples(index=False)):
            y=.76-i*.073
            if i%2==0:ax.add_patch(patches.Rectangle((.02,y-.032),.95,.062,transform=ax.transAxes,fc=LIGHT,ec="none",alpha=.55))
            anatomy="NA" if pd.isna(r.mean_anatomical_expression_concordance) else f"{r.mean_anatomical_expression_concordance:.2f}"
            vals=[shorten(r.approved_symbol,10),f"{r.transparent_disease_evidence_index:.0f}",f"{int(r.chemical_pair_count):,}",f"{int(r.coordinate_ready_site_count):,}",anatomy]
            for x,v in zip(xs,vals):ax.text(x,y,v,transform=ax.transAxes,fontsize=5.0,va="center")
        ax.text(.02,.04,"High disease support, low-but-nonzero chemistry and ≥1 coordinate-ready site",transform=ax.transAxes,fontsize=5.0,color=MID)
    else:ax.text(.5,.5,"No protein satisfies the preregistered\nrelative-underexploration rule",ha="center",va="center")
    ax=fig.add_subplot(gs[2,:]);panel(ax,"e","Interpretation and follow-up")
    ax.axis("off");steps=[("Disease support","genetic · clinical · functional",LTBLUE),("Anatomical context","mapped normal-tissue RNA",PALEGREEN),("Chemical gap","relative pair coverage",YELLOW),("Structural route","coordinate-ready site",PEACH),("Follow-up","literature · assay · docking",SALMON)]
    for i,(lab,desc,col) in enumerate(steps):
        x=.015+i*.197;ax.add_patch(patches.FancyBboxPatch((x,.28),.17,.42,boxstyle="round,pad=.012",fc=col,ec="white"));ax.text(x+.085,.56,lab,ha="center",fontweight="bold",fontsize=5.5);ax.text(x+.085,.40,desc,ha="center",fontsize=5.0,color=MID)
        if i<4:ax.annotate("",xy=(x+.193,.49),xytext=(x+.174,.49),arrowprops=dict(arrowstyle="->",lw=.7,color=MID))
    ax.text(.5,.12,"This workflow prioritises hypotheses; it does not infer causality or validate binding.",ha="center",fontsize=5.4,color=RED)
    source_footer(fig,"V7.2 formal protein–compound pairs, protein–disease relations, mapped HPA consensus tissue RNA and coordinate-ready sites. Permutations preserve protein and disease degrees.")
    save(fig,"Figure6_integrated_drug_discovery_insights")


def supplementary_figures():
    emb=read("S3_chemical_umap_deterministic_sample.tsv")
    fig,ax=plt.subplots(figsize=(6.8,5.5),constrained_layout=True);panel(ax,"S1","Deterministic chemical-space sample")
    statuses=list(emb["biological_status_display"].value_counts().index)
    for i,status in enumerate(statuses):
        g=emb[emb["biological_status_display"].eq(status)];ax.scatter(g["UMAP1"],g["UMAP2"],s=4,alpha=.28,color=PALETTE[i%len(PALETTE)],label=f"{status} (n={len(g):,})",rasterized=True)
    ax.set_xlabel("UMAP 1 (no direct chemical meaning)");ax.set_ylabel("UMAP 2 (no direct chemical meaning)");ax.legend(markerscale=2,loc="best");ax.set_title("Supplementary Figure S1 | Morgan-fingerprint chemical atlas",loc="left",fontweight="bold");fig.text(.01,.01,"Deterministic status-enriched sample; local neighbourhoods only. Global UMAP distances are not interpreted as exact Tanimoto similarity.",fontsize=5.2,color=MID)
    save(fig,"Supplementary_Figure_S1_chemical_umap",supplementary=True)


def main():
    figure1();figure2();figure3();figure4();figure5();figure6();supplementary_figures()
    print(json.dumps({"main_figures":6,"supplementary_figures":1,"formats":["png","svg","pdf"]}))


if __name__=="__main__":main()
