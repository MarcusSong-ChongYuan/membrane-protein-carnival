from __future__ import annotations

import json
import math
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.image as mpimg
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CFG = json.loads((ROOT / "00_config" / "nar_config.json").read_text(encoding="utf-8"))
SD = ROOT / "03_source_data"
MAIN = ROOT / "04_main_figures"
SUPP = ROOT / "05_supplementary"
GA = ROOT / "06_graphical_abstract"
ANATOMY = ROOT / "09_assets" / "human_anatomy_highres_v2.png"
MM = 1 / 25.4
P = CFG["palette"]

PROTEIN = P["protein"]; COMPOUND = P["compound"]; DISEASE = P["disease"]
EXPRESSION = P["expression"]; STRUCTURE = P["structure"]; UNKNOWN = P["unknown"]
TEXT = P["text"]; GRID = P["grid"]; NEG = P["negative"]; NEUTRAL = P["neutral"]; POS = P["positive"]
DIVERGING = LinearSegmentedColormap.from_list("mempro_diverging", [NEG, NEUTRAL, POS])

PUBLICATION_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 6.4,
    "axes.titlesize": 7.2,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 5.5,
    "ytick.labelsize": 5.5,
    "legend.fontsize": 5.5,
    "axes.linewidth": .55,
    "axes.edgecolor": "#9AA0A6",
    "axes.labelcolor": TEXT,
    "axes.titlecolor": TEXT,
    "text.color": TEXT,
    "xtick.color": TEXT,
    "ytick.color": TEXT,
    "grid.color": GRID,
    "grid.linewidth": .4,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
}
# Seaborn's context defaults overwrite font sizes even when rcParams are passed
# back into set_theme.  Apply the visual style first, then restore the exact
# journal-scale typography.
sns.set_theme(style="whitegrid")
mpl.rcParams.update(PUBLICATION_RC)


def read(name: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(SD / name, sep="\t", compression="infer", **kwargs)


def wrap(v: object, width: int = 22) -> str:
    return "\n".join(textwrap.wrap(str(v).replace("_", " "), width=width, break_long_words=False, break_on_hyphens=False))


def panel(ax, letter: str, title: str) -> None:
    ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", fontsize=7.2, pad=4)


def source_note(fig, text: str) -> None:
    fig.text(.012, .008, text, ha="left", va="bottom", fontsize=5.0, color="#5F6368")


def save_main(fig, stem: str) -> None:
    for ext in ("pdf", "svg", "png"):
        path = MAIN / ext / f"{stem}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        kwargs = {"dpi": int(CFG["preview_dpi"])} if ext == "png" else {}
        fig.savefig(path, **kwargs)
    plt.close(fig)


def save_supp(fig, stem: str) -> None:
    for ext in ("pdf", "svg", "png"):
        path = SUPP / ext / f"{stem}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        kwargs = {"dpi": int(CFG["preview_dpi"])} if ext == "png" else {}
        fig.savefig(path, **kwargs)
    plt.close(fig)


def fmt(n: float) -> str:
    if n >= 1_000_000: return f"{n/1_000_000:.2f}M"
    if n >= 1_000: return f"{n/1_000:.1f}k"
    return f"{int(n):,}"


def top_other_edges(df: pd.DataFrame, stage_cols: list[str], count_col: str, top_n=6) -> pd.DataFrame:
    work = df.copy()
    for col in stage_cols:
        totals = work.groupby(col)[count_col].sum().sort_values(ascending=False)
        keep = set(totals.head(top_n).index)
        work[col] = work[col].where(work[col].isin(keep), "Other")
    return work.groupby(stage_cols, as_index=False)[count_col].sum()


def draw_alluvial(ax, stages: list[str], edge_frames: list[pd.DataFrame], count_col: str, palette=None,
                  label_fontsize: float = 5.0, label_wrap: int = 16) -> None:
    palette = palette or [PROTEIN, EXPRESSION, COMPOUND, DISEASE, STRUCTURE, "#7B95C6", "#49C2D9"]
    xpos = np.linspace(.06, .94, len(stages))
    node_totals = {s: defaultdict(float) for s in stages}
    for i, e in enumerate(edge_frames):
        a, b = stages[i], stages[i+1]
        for r in e.itertuples(index=False):
            node_totals[a][getattr(r, a)] += getattr(r, count_col)
            node_totals[b][getattr(r, b)] += getattr(r, count_col)
    node_pos = {}
    node_color = {}
    for si, s in enumerate(stages):
        totals = node_totals[s]
        total = sum(totals.values()) or 1
        gap = .020
        scale = (.76 - gap * max(len(totals)-1, 0)) / total
        y = .88
        for j, (label, val) in enumerate(sorted(totals.items(), key=lambda x: x[1], reverse=True)):
            h = max(.012, val * scale)
            node_pos[(s,label)] = (y-h, y, h)
            node_color[(s,label)] = palette[j % len(palette)]
            rect = patches.FancyBboxPatch((xpos[si]-.012, y-h), .024, h, boxstyle="round,pad=.001", fc=node_color[(s,label)], ec="white", lw=.25, transform=ax.transAxes, zorder=3)
            ax.add_patch(rect)
            ha = "right" if si == len(stages)-1 else "left"
            dx = -.018 if si == len(stages)-1 else .018
            ax.text(xpos[si]+dx, y-h/2, wrap(label, label_wrap), ha=ha, va="center", fontsize=label_fontsize,
                    linespacing=1.08, transform=ax.transAxes, clip_on=True)
            y -= h + gap
        ax.text(xpos[si], .955, wrap(s, 18), ha="center", va="bottom", fontsize=5.4,
                fontweight="bold", transform=ax.transAxes, clip_on=True)
    for i, e in enumerate(edge_frames):
        a, b = stages[i], stages[i+1]
        out_offsets = defaultdict(float); in_offsets = defaultdict(float)
        a_totals = node_totals[a]; b_totals = node_totals[b]
        for r in e.sort_values(count_col, ascending=False).itertuples(index=False):
            la, lb, val = getattr(r,a), getattr(r,b), float(getattr(r,count_col))
            ay0, ay1, ah = node_pos[(a,la)]; by0, by1, bh = node_pos[(b,lb)]
            h1 = ah * val / max(a_totals[la],1); h2 = bh * val / max(b_totals[lb],1)
            ya = ay0 + out_offsets[(a,la)] + h1/2; yb = by0 + in_offsets[(b,lb)] + h2/2
            out_offsets[(a,la)] += h1; in_offsets[(b,lb)] += h2
            xs = np.linspace(xpos[i]+.012, xpos[i+1]-.012, 50)
            t = np.linspace(0,1,50); smooth = 3*t*t-2*t*t*t
            center = ya + (yb-ya)*smooth; half = (h1/2)+(h2-h1)/2*smooth
            ax.fill_between(xs, center-half, center+half, color=node_color[(a,la)], alpha=.20, lw=0, transform=ax.transAxes, zorder=1)
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")


def draw_upset(ax_bar, ax_matrix, data: pd.DataFrame, set_cols: list[str], count_col: str, top_n=10) -> None:
    d = data.nlargest(top_n, count_col).reset_index(drop=True)
    x = np.arange(len(d))
    ax_bar.bar(x, d[count_col], color="#49AFC8", width=.72)
    ax_bar.set_xticks([]); ax_bar.set_ylabel("Count")
    for i, v in enumerate(d[count_col]): ax_bar.text(i, v, fmt(v), ha="center", va="bottom", fontsize=5.0)
    ax_matrix.set_xlim(-.5, len(d)-.5); ax_matrix.set_ylim(-.5, len(set_cols)-.5)
    for yi, col in enumerate(set_cols):
        vals = d[col].astype(int).to_numpy()
        ax_matrix.scatter(x, np.full(len(x), yi), s=10, c=np.where(vals==1, TEXT, "#D9DEE2"), zorder=2)
        idx = np.where(vals==1)[0]
        if len(idx):
            for xi in idx:
                active = np.where(d.loc[xi,set_cols].astype(int).to_numpy()==1)[0]
                if len(active)>1: ax_matrix.plot([xi,xi],[active.min(),active.max()], color=TEXT, lw=.6, zorder=1)
    ax_matrix.set_yticks(range(len(set_cols)), [wrap(x,17) for x in set_cols])
    ax_matrix.set_xticks([]); ax_matrix.grid(False)
    for s in ax_matrix.spines.values(): s.set_visible(False)


def figure1() -> None:
    """Render Figure 1 with explicit physical margins.

    This figure deliberately avoids constrained_layout: the combination of an
    alluvial panel, direct labels and a long y-axis label list can otherwise
    collapse axes and inflate text visually.  All labels remain editable and
    are at least 5 pt in the exported PDF.
    """
    fig = plt.figure(figsize=(178*MM, 195*MM))
    gs = fig.add_gridspec(
        3, 4,
        height_ratios=[.66, 1.45, 1.02],
        width_ratios=[1.0, 1.0, 1.0, 1.0],
        left=.055, right=.975, bottom=.115, top=.975,
        hspace=.34, wspace=.68,
    )
    ax = fig.add_subplot(gs[0,:]); panel(ax,"a","Data integration and release framework"); ax.axis("off")
    steps = [
        ("Source records", "Protein · compound\ninteraction · structure\ndisease · expression", PROTEIN),
        ("Entity harmonization", "Canonical protein / isoform\nparent compound / form\nMONDO disease", EXPRESSION),
        ("Evidence harmonization", "Deduplication · lineage\nBE tiers · missingness", COMPOUND),
        ("Frozen V7.2 release", "Tables · source data\nmanifest · SHA-256", DISEASE),
    ]
    for i,(title,body,color) in enumerate(steps):
        x=.02+i*.245
        ax.add_patch(patches.FancyBboxPatch((x,.22),.20,.56,boxstyle="round,pad=.012",fc=color,ec="none",alpha=.18,transform=ax.transAxes))
        ax.text(x+.10,.61,title,ha="center",va="center",fontweight="bold",fontsize=6.3,transform=ax.transAxes)
        ax.text(x+.10,.40,body,ha="center",va="center",fontsize=5.3,linespacing=1.35,transform=ax.transAxes)
        if i<3: ax.annotate("",xy=(x+.24,.50),xytext=(x+.215,.50),xycoords=ax.transAxes,arrowprops={"arrowstyle":"->","lw":.7,"color":"#737B80"})
    flow = read("Figure1B_source_modality_tier_flow.tsv")
    def modality(x):
        s=str(x).lower()
        if "structur" in s or "contact" in s: return "Structural"
        if "direct" in s and "quant" in s: return "Direct quantitative"
        if "functional" in s or "activity" in s or "pharm" in s: return "Functional / pharmacological"
        return "Curated / supporting"
    flow["evidence_modality"] = flow["evidence_type"].map(modality)
    # Four dominant sources remain directly labelled; all smaller contributors
    # are retained quantitatively in an explicit aggregate to avoid collisions.
    source_keep = set(flow.groupby("source_database")["evidence_record_count"].sum().nlargest(4).index)
    flow["source"] = flow["source_database"].where(flow["source_database"].isin(source_keep),"Other sources")
    e1=flow.groupby(["source","evidence_modality"],as_index=False)["evidence_record_count"].sum()
    e2=flow.groupby(["evidence_modality","evidence_tier"],as_index=False)["evidence_record_count"].sum()
    pd.concat([e1.assign(edge="source_to_modality"),e2.assign(edge="modality_to_tier")],ignore_index=True,sort=False).to_csv(SD/"Figure1B_simplified_flow.tsv",sep="\t",index=False)
    ax=fig.add_subplot(gs[1,:]); panel(ax,"b","Source records converge into evidence modalities and BE tiers")
    draw_alluvial(ax,["source","evidence_modality","evidence_tier"],[e1,e2],"evidence_record_count")
    cov=read("Figure1_module_coverage.tsv"); cov=cov[cov["module"].ne("formal_protein")].copy()
    coverage_labels = {
        "expression_mapped_measured": "Expression measured",
        "subcellular_localization": "Subcellular localization",
        "compound_interaction": "Compound interaction",
        "disease": "Disease association",
        "binding_site": "Mapped binding site",
    }
    cov["label"] = cov["module"].map(coverage_labels)
    ax=fig.add_subplot(gs[2,:2]); panel(ax,"c","Protein-level module coverage")
    y=np.arange(len(cov))[::-1]
    ax.hlines(y,0,cov["coverage_fraction"],color="#DDE2E5",lw=1.6)
    ax.scatter(cov["coverage_fraction"],y,s=25,color=[EXPRESSION,EXPRESSION,COMPOUND,DISEASE,STRUCTURE][:len(cov)],edgecolor="white",lw=.4,zorder=3)
    ax.set_yticks(y); ax.set_yticklabels([]); ax.set_xlim(-.72,1.08)
    ax.set_xticks([0,.25,.50,.75,1.0]); ax.tick_params(axis="x",labelsize=5.0)
    ax.set_xlabel("Coverage among 7,800 formal proteins", labelpad=3, fontsize=5.7)
    ax.xaxis.set_label_coords(.70,-.13)
    ax.grid(axis="x",color=GRID,lw=.4); ax.grid(axis="y",visible=False)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y",length=0,pad=3)
    for yi,r in zip(y,cov.itertuples()):
        ax.text(-.69,yi,r.label,ha="left",va="center",fontsize=5.2,clip_on=True)
        ax.text(r.coverage_fraction+.018,yi,f"{r.coverage_fraction:.0%}  n={int(r.protein_count):,}",va="center",fontsize=5.0)
    counts=read("Figure1_resource_counts.tsv")
    ax=fig.add_subplot(gs[2,2:]); panel(ax,"d","Release entities and analysis denominators"); ax.axis("off")
    show=["Formal proteins","Compound registry","Interaction-linked compounds","Formal protein-compound pairs","Public positive evidence","Canonical diseases","Expression measurements","Binding-site assertions"]
    data=counts.set_index("metric").loc[show]
    colors=[PROTEIN,COMPOUND,EXPRESSION,COMPOUND,DISEASE,DISEASE,EXPRESSION,STRUCTURE]
    for i,(idx,r) in enumerate(data.iterrows()):
        x=.02+(i%2)*.50; yy=.86-(i//2)*.23
        ax.text(x,yy,fmt(r["count"]),fontsize=7.8,fontweight="bold",color=colors[i],transform=ax.transAxes)
        ax.text(x,yy-.095,wrap(idx,23),fontsize=5.1,linespacing=1.15,transform=ax.transAxes)
    source_note(fig,"MemPro V7.2 formal release. Counts denote entities or records, not independent experiments.")
    save_main(fig,"Figure1_MemPro_integration_and_coverage")


def figure2() -> None:
    fig=plt.figure(figsize=(178*MM,225*MM)); gs=fig.add_gridspec(
        3,12,height_ratios=[1.00,1.34,1.16],
        left=.17,right=.965,bottom=.082,top=.965,hspace=.50,wspace=.95)
    tab=read("Figure2A_membrane_class_evidence.tsv").set_index("membrane_class_v7")
    ax=fig.add_subplot(gs[0,:5]);panel(ax,"a","Membrane class by evidence tier")
    total=tab.values.sum(); x0=0
    ecolors={"E1":PROTEIN,"E2":"#7B95C6","E3":"#A1D8E8"}
    centers=[]
    for cls,row in tab.iterrows():
        width=row.sum()/total; y0=0
        centers.append((x0+width/2, cls, int(row.sum())))
        for ev,val in row.items():
            h=val/row.sum() if row.sum() else 0
            ax.add_patch(patches.Rectangle((x0,y0),width,h,fc=ecolors[ev],ec="white",lw=.5,alpha=.9))
            if val>=100 and width>.13 and h>.10:
                ax.text(x0+width/2,y0+h/2,f"{ev}\n{int(val):,}",ha="center",va="center",fontsize=5.0)
            y0+=h
        x0+=width
    for xc,cls,n in centers:
        ax.text(xc,-.055,f"{cls}\n{n:,}",ha="center",va="top",fontsize=5.2,fontweight="bold")
    handles=[patches.Patch(facecolor=ecolors[k],edgecolor="none",label=k) for k in ["E1","E2","E3"]]
    ax.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,1.02),ncol=3,frameon=False,
              handlelength=1.0,columnspacing=1.1,fontsize=5.0)
    ax.set_xlim(0,1);ax.set_ylim(-.13,1.05);ax.axis("off")
    cells=read("Figure2B_role_function_enrichment.tsv")
    roles=cells.groupby("primary_membrane_role")["observed"].sum().nlargest(6).index
    funcs=cells.groupby("primary_molecular_function")["observed"].sum().nlargest(6).index
    show=cells[cells["primary_membrane_role"].isin(roles)&cells["primary_molecular_function"].isin(funcs)]
    mat=show.pivot(index="primary_membrane_role",columns="primary_molecular_function",values="log2_odds_ratio").reindex(index=roles,columns=funcs)
    q=show.pivot(index="primary_membrane_role",columns="primary_molecular_function",values="bh_q").reindex(index=roles,columns=funcs)
    ax=fig.add_subplot(gs[0,6:]);panel(ax,"b","Membrane-role by molecular-function enrichment")
    sns.heatmap(mat,ax=ax,cmap=DIVERGING,center=0,vmin=-7,vmax=7,linewidths=.25,linecolor="white",cbar_kws={"label":"log2 odds ratio","shrink":.72})
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(q.iat[i,j]) and q.iat[i,j]<.05: ax.text(j+.5,i+.5,"•",ha="center",va="center",fontsize=6.5,color=TEXT)
    function_display={
        "binding_activity":"binding", "catalytic_activity":"catalytic",
        "signaling_receptor_activity":"receptor",
        "transmembrane_transporter_activity":"transporter",
        "molecular_function_regulator":"function regulator",
        "ion_channel_activity":"ion channel", "unclassified":"unclassified",
    }
    role_display={
        "membrane_associated_enzyme":"membrane-associated enzyme",
        "family_defined_membrane_role_unresolved":"family-defined; role unresolved",
        "membrane_scaffold_or_linker":"scaffold / linker",
    }
    ax.set_xticklabels([wrap(function_display.get(x,x),12) for x in mat.columns],rotation=30,ha="right",rotation_mode="anchor")
    ax.set_yticklabels([wrap(role_display.get(x,x),18) for x in mat.index],rotation=0);ax.set_xlabel("");ax.set_ylabel("")
    ax.text(.015,.985,"Cramer's V = 0.619; dot = BH q < 0.05",ha="left",va="top",
            transform=ax.transAxes,fontsize=5.0,bbox={"fc":"white","ec":"none","alpha":.82,"pad":1.0})
    edges=read("Figure2C_five_axis_alluvial_edges.tsv")
    stage_order=["membrane_class","membrane_role","molecular_function","biological_process"]
    frames=[]
    for a,b in zip(stage_order[:-1],stage_order[1:]):
        e=edges[(edges["source_axis"].eq(a))&(edges["target_axis"].eq(b))][["source_label","target_label","protein_count"]].rename(columns={"source_label":a,"target_label":b})
        frames.append(e)
    all_rows=[]
    simplified=[]
    for i,e in enumerate(frames):
        x=top_other_edges(e,[stage_order[i],stage_order[i+1]],"protein_count",top_n=6);simplified.append(x);all_rows.append(x.assign(edge=f"{stage_order[i]}->{stage_order[i+1]}"))
    pd.concat(all_rows,ignore_index=True,sort=False).to_csv(SD/"Figure2C_simplified_five_axis_flow.tsv",sep="\t",index=False)
    ax=fig.add_subplot(gs[1,:]);panel(ax,"c","Primary-category flow across four classification axes")
    draw_alluvial(ax,stage_order,simplified,"protein_count",label_fontsize=5.0,label_wrap=14)
    tau=read("Figure2D_tissue_specificity_tau.tsv");counts=tau.groupby("primary_membrane_role")["tau_consensus_tissue_rna"].count();keep=counts[counts>=30].index;tau=tau[tau["primary_membrane_role"].isin(keep)].copy();order=tau.groupby("primary_membrane_role")["tau_consensus_tissue_rna"].median().sort_values().index
    ax=fig.add_subplot(gs[2,:]);panel(ax,"d","Normal-tissue RNA specificity by membrane role")
    sns.violinplot(data=tau,x="tau_consensus_tissue_rna",y="primary_membrane_role",order=order,inner=None,cut=0,linewidth=.35,color="#A1D8E8",ax=ax)
    sns.boxplot(data=tau,x="tau_consensus_tissue_rna",y="primary_membrane_role",order=order,width=.16,showfliers=False,boxprops={"facecolor":"white","alpha":.85},medianprops={"color":DISEASE,"lw":.9},ax=ax)
    ax.set_xlim(0,1);ax.set_xlabel("Tau (0 = broad expression; 1 = tissue-specific)");ax.set_ylabel("")
    ax.set_yticks(np.arange(len(order)));ax.set_yticklabels([wrap(x,28) for x in order])
    ax.text(.995,1.012,"Kruskal-Wallis P = 4.34e-92; epsilon-squared = 0.061",ha="right",va="bottom",
            transform=ax.transAxes,fontsize=5.0)
    source_note(fig,"V7.2 canonical membrane proteins and mapped+measured normal HPA tissue RNA. Unresolved annotations remain explicit.")
    save_main(fig,"Figure2_membrane_protein_organization")


def figure3() -> None:
    fig=plt.figure(figsize=(178*MM,210*MM));gs=fig.add_gridspec(
        2,2,height_ratios=[1.0,1.02],left=.065,right=.95,bottom=.08,top=.975,hspace=.34,wspace=.50)
    rank=read("Figure3A_scaffold_rank_abundance.tsv");metrics=read("Figure3A_scaffold_metrics.tsv");top=read("Figure3A_top_scaffolds.tsv")
    outer=gs[0,0].subgridspec(2,1,height_ratios=[3,1]);ax=fig.add_subplot(outer[0]);panel(ax,"a","Bemis–Murcko scaffold diversity")
    ax.loglog(rank["rank"],rank["compound_count"],color=PROTEIN,lw=1.0);ax.fill_between(rank["rank"],rank["compound_count"],1,color="#A1D8E8",alpha=.25);ax.set_xlabel("Scaffold rank");ax.set_ylabel("Interaction-linked compounds per scaffold")
    gini=float(metrics.loc[metrics["metric"].eq("Scaffold Gini"),"value"].iloc[0]);ax.text(.98,.95,f"91,977 groups · Gini={gini:.3f}",ha="right",va="top",transform=ax.transAxes,fontsize=5.2)
    strip=fig.add_subplot(outer[1]);strip.axis("off");molrows=top[top["murcko_scaffold_smiles"].ne("ACYCLIC_NO_MURCKO_SCAFFOLD")].head(5)
    for i,r in enumerate(molrows.itertuples(index=False)):
        mol=Chem.MolFromSmiles(r.murcko_scaffold_smiles)
        if mol is None: continue
        img=np.asarray(Draw.MolToImage(mol,size=(280,160),kekulize=True))
        ia=strip.inset_axes([i*.2,.12,.18,.78]);ia.imshow(img);ia.axis("off");ia.text(.5,-.02,f"n={int(r.compound_count):,}",ha="center",va="top",transform=ia.transAxes,fontsize=5.0)
    enr=read("Figure3B_compound_status_role_enrichment.tsv");roles=enr.groupby("membrane_role")["role_pair_count"].max().nlargest(9).index;show=enr[enr["membrane_role"].isin(roles)];mat=show.pivot(index="compound_status",columns="membrane_role",values="log2_odds_ratio").reindex(columns=roles);q=show.pivot(index="compound_status",columns="membrane_role",values="bh_q").reindex(index=mat.index,columns=mat.columns)
    ax=fig.add_subplot(gs[0,1]);panel(ax,"b","Overlapping compound-status preferences")
    mat=mat.T;q=q.T
    lim=max(1.0,float(np.nanpercentile(np.abs(mat),95)));sns.heatmap(mat,ax=ax,cmap=DIVERGING,center=0,vmin=-lim,vmax=lim,linewidths=.25,linecolor="white",cbar_kws={"label":"log2 odds ratio","shrink":.72})
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(q.iat[i,j]) and q.iat[i,j]<.05: ax.text(j+.5,i+.5,"•",ha="center",va="center",fontsize=6.2)
    ax.set_xticklabels([wrap(x,17) for x in mat.columns],rotation=28,ha="right",rotation_mode="anchor");ax.set_yticklabels([wrap(x,20) for x in mat.index],rotation=0);ax.set_xlabel("");ax.set_ylabel("")
    ax.text(.99,.985,"Binary, non-exclusive labels; dot = BH q < 0.05",ha="right",va="top",transform=ax.transAxes,fontsize=5.0,bbox={"fc":"white","ec":"none","alpha":.82,"pad":1})
    ent=read("Figure3C_polypharmacology_entropy.tsv");ax=fig.add_subplot(gs[1,0]);panel(ax,"c","Target breadth and cross-role polypharmacology")
    hb=ax.hexbin(np.log10(ent["target_count"]+1),ent["target_role_entropy_normalized"],gridsize=42,mincnt=1,bins="log",cmap=sns.color_palette("flare",as_cmap=True));fig.colorbar(hb,ax=ax,label="log10(compound count)",shrink=.72);ax.set_xlabel("log10(formal target count + 1)");ax.set_ylabel("Degree-corrected membrane-role entropy")
    sim=read("Figure3D_chemical_target_similarity.tsv.gz");stats=json.loads((ROOT/"02_analysis_screening"/"chemical_similarity_target_overlap"/"statistics.json").read_text(encoding="utf-8"));ax=fig.add_subplot(gs[1,1]);panel(ax,"d","Chemical similarity and target-set overlap")
    hb=ax.hexbin(sim["tanimoto_similarity"],sim["target_set_jaccard"],gridsize=44,mincnt=1,bins="log",cmap=sns.color_palette("rocket",as_cmap=True));fig.colorbar(hb,ax=ax,label="log10(pair count)",shrink=.72);ax.set_xlabel("Morgan/Tanimoto similarity");ax.set_ylabel("Target-set Jaccard")
    ci=stats["cluster_bootstrap_ci95"]
    ax.text(.99,.985,
            f"rho = {stats['spearman_rho']:.3f} (95% CI {ci[0]:.3f}-{ci[1]:.3f}); ANN recall = {stats['approximate_neighbor_recall_mean']:.3f}",
            ha="right",va="top",transform=ax.transAxes,fontsize=5.0,bbox={"fc":"white","ec":"none","alpha":.82,"pad":1})
    source_note(fig,"V7.2 standard SMILES and formal protein–compound pairs. Compound-status labels overlap; similarity analysis requires ≥2 targets per compound.")
    save_main(fig,"Figure3_chemical_diversity_and_preference")


def figure4() -> None:
    fig=plt.figure(figsize=(178*MM,225*MM));gs=fig.add_gridspec(
        3,4,height_ratios=[1.16,1.0,1.0],left=.145,right=.95,bottom=.075,top=.975,hspace=.52,wspace=.70)
    upset=read("Figure4A_source_upset.tsv")
    source_short={"IUPHAR/BPS Guide to PHARMACOLOGY":"IUPHAR/BPS Guide","PubChem BioAssay":"PubChem BioAssay"}
    upset=upset.rename(columns=source_short)
    outer=gs[0,:2].subgridspec(2,1,height_ratios=[2.2,1]);axb=fig.add_subplot(outer[0]);panel(axb,"a","Major-source intersections");axm=fig.add_subplot(outer[1]);sets=[c for c in upset.columns if c not in {"intersection","pair_count"}];draw_upset(axb,axm,upset,sets,"pair_count",10)
    ov=read("Figure4B_source_overlap_triangular.tsv");sources=ov["source"].tolist();mat=np.full((len(sources),len(sources)),np.nan)
    for i,s in enumerate(sources):
        for j,t in enumerate(sources):
            if i<j: mat[i,j]=ov.loc[i,f"Jaccard::{t}"]
            elif i>j: mat[i,j]=ov.loc[i,f"Overlap::{t}"]
    ax=fig.add_subplot(gs[0,2:]);panel(ax,"b","Source overlap and containment")
    display_sources=[source_short.get(x,x) for x in sources]
    im=ax.imshow(mat,vmin=0,vmax=1,cmap=sns.light_palette(PROTEIN,as_cmap=True));ax.set_xticks(range(len(sources)),[wrap(x,15) for x in display_sources],rotation=38,ha="right",rotation_mode="anchor");ax.set_yticks(range(len(sources)),[wrap(x,18) for x in display_sources]);fig.colorbar(im,ax=ax,label="Upper: Jaccard; lower: overlap coefficient",shrink=.72)
    for i in range(len(sources)):
        for j in range(len(sources)):
            if np.isfinite(mat[i,j]): ax.text(j,i,f"{mat[i,j]:.2f}",ha="center",va="center",fontsize=5.0,color="white" if mat[i,j]>.55 else TEXT)
    agg=read("Figure4C_database_putative_lineage.tsv");x=np.repeat(agg["distinct_database_count"].to_numpy(int),agg["pair_count"].to_numpy(int));y=np.repeat(agg["independent_experiment_count"].to_numpy(int),agg["pair_count"].to_numpy(int))
    ax=fig.add_subplot(gs[1,:2]);panel(ax,"c","Contributing databases versus putative lineages");hb=ax.hexbin(x,np.log10(y+1),gridsize=(18,38),mincnt=1,bins="log",cmap=sns.color_palette("flare",as_cmap=True));cax=ax.inset_axes([.60,.88,.34,.025]);cb=fig.colorbar(hb,cax=cax,orientation="horizontal");cb.set_label("log10(pair count)",fontsize=5.0,labelpad=1);cb.ax.tick_params(labelsize=5.0);xx=np.arange(1,x.max()+1);ax.plot(xx,np.log10(xx+1),ls="--",lw=.7,color="#5F6368",label="One lineage per database");ax.set_xlabel("Distinct contributing databases");ax.set_ylabel("log10(putative lineage count + 1)");ax.legend(loc="upper left")
    contrib=read("Figure4D_source_contribution.tsv").head(6);metrics=["unique_pair_fraction","shared_pair_fraction","quantitative_evidence_fraction","structural_evidence_fraction","resolved_putative_lineage_fraction"]
    ax=fig.add_subplot(gs[1,2:]);panel(ax,"d","Source contribution profile")
    for i,r in contrib.reset_index(drop=True).iterrows():
        for j,m in enumerate(metrics):
            v=float(r[m]);ax.scatter(j,i,s=10+95*v,c=v,cmap=sns.light_palette(STRUCTURE,as_cmap=True),vmin=0,vmax=1,edgecolor="white",lw=.3)
    ax.set_xlim(-.5,len(metrics)-.5);ax.set_ylim(len(contrib)-.5,-.5);ax.set_xticks(range(len(metrics)),[wrap(x.replace("_fraction",""),16) for x in metrics],rotation=38,ha="right",rotation_mode="anchor");ax.set_yticks(range(len(contrib)),[wrap(source_short.get(x,x),18) for x in contrib["source_database"]]);ax.grid(False);sm=mpl.cm.ScalarMappable(norm=Normalize(0,1),cmap=sns.light_palette(STRUCTURE,as_cmap=True));fig.colorbar(sm,ax=ax,label="Fraction",shrink=.72)
    ecdf=read("Figure4E_site_mapping_ecdf.tsv");ax=fig.add_subplot(gs[2,:2]);panel(ax,"e","Residue-coordinate mapping by site tier")
    colors=[PROTEIN,EXPRESSION,COMPOUND]
    tier_display={"S1_EXACT_COCRYSTAL_COORDINATE_COMPLETE":"S1 exact cocrystal",
                  "S2_COORDINATE_COMPLETE_LIGAND_UNRESOLVED":"S2 coordinates complete",
                  "S3_SOURCE_OR_PARTIAL":"S3 source / partial"}
    for (tier,g),c in zip(ecdf.groupby("site_tier"),colors): ax.step(g["mapping_fraction"],g["ecdf"],where="post",label=tier_display.get(tier,tier.replace("_"," ")),color=c,lw=1.0)
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_xlabel("Mapped source residues / requested residues");ax.set_ylabel("ECDF");ax.legend(loc="lower left",fontsize=5.0,frameon=False)
    side=read("Figure4E_membrane_side.tsv").sort_values("site_count");ax=fig.add_subplot(gs[2,2:]);panel(ax,"f","Membrane-side annotation completeness");y=np.arange(len(side));ax.barh(y,side["site_count"],color=[UNKNOWN if str(x).lower()=="unknown" else STRUCTURE for x in side["membrane_side"]],alpha=.82);ax.set_yticks(y,[wrap(x,22) for x in side["membrane_side"]]);ax.set_xlabel("Site assertions")
    side_display={"extramembrane side unresolved":"extramembrane unresolved",
                  "membrane interface or mixed":"interface / mixed",
                  "non cytoplasmic side":"non-cytoplasmic side",
                  "both sides or mixed":"both sides / mixed"}
    ax.set_yticklabels([wrap(side_display.get(str(x).lower(),str(x)),18) for x in side["membrane_side"]])
    ax.set_xlim(0,float(side["site_count"].max())*1.18)
    for yi,r in enumerate(side.itertuples()):ax.text(r.site_count,yi,f"  {int(r.site_count):,} ({r.fraction:.1%})",va="center",fontsize=5.0)
    source_note(fig,"V7.2 formal pairs, public evidence lineages and interaction sites. Putative lineages are not manually audited publication experiments.")
    save_main(fig,"Figure4_evidence_independence_and_sites")


ANCHORS={"nervous system":(.50,.86),"entire sense organ system":(.45,.82),"respiratory system":(.48,.68),"cardiovascular system":(.52,.62),"circulatory system":(.56,.61),"vascular system":(.60,.58),"digestive system":(.50,.51),"alimentary part of gastrointestinal system":(.52,.45),"excretory system":(.44,.46),"reproductive system":(.50,.36),"hematopoietic system":(.60,.28),"lymphoid system":(.60,.55),"musculoskeletal system":(.38,.29),"integumental system":(.28,.55),"exocrine system":(.55,.50),"lacrimal apparatus":(.41,.81),"neuroendocrine system":(.50,.73)}


def anatomy_panel(ax, anatomy: pd.DataFrame) -> None:
    img=mpimg.imread(ANATOMY);ax.imshow(img,extent=(.30,.64,.04,.96),aspect="auto");ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    a=anatomy[~anatomy["anatomical_system"].str.lower().eq("unmapped")].copy()
    a=a[a["anatomical_system"].str.lower().isin(ANCHORS)].nlargest(12,"fractional_disease_count")
    left=[x for x in a["anatomical_system"] if ANCHORS[x.lower()][0]<.5];right=[x for x in a["anatomical_system"] if x not in left]
    left=sorted(left,key=lambda x:ANCHORS[x.lower()][1],reverse=True);right=sorted(right,key=lambda x:ANCHORS[x.lower()][1],reverse=True)
    yleft=np.linspace(.88,.16,len(left));yright=np.linspace(.90,.12,len(right));mx=a["fractional_disease_count"].max()
    lookup=a.set_index("anatomical_system")
    for side,names,ys,xlab in [("left",left,yleft,.06),("right",right,yright,.76)]:
        for name,y in zip(names,ys):
            val=float(lookup.loc[name,"fractional_disease_count"]);anchor=ANCHORS[name.lower()]
            s=28+180*(val/mx)
            ax.scatter([xlab],[y],s=s,facecolor=DISEASE,edgecolor="white",lw=.45,alpha=.48,zorder=4)
            ax.annotate(wrap(name,17),xy=anchor,xytext=(xlab,y),xycoords=ax.transAxes,textcoords=ax.transAxes,ha="right" if side=="left" else "left",va="center",fontsize=5.0,arrowprops={"arrowstyle":"-","lw":.35,"color":"#7A8085"})
    ax.text(.02,.01,"Top 12 mapped systems; bubble area = fractional disease count\nOntology-derived multi-label anatomy; full set in Source Data",transform=ax.transAxes,fontsize=5.0,color="#5F6368")


def figure5() -> None:
    fig=plt.figure(figsize=(178*MM,228*MM));gs=fig.add_gridspec(
        3,4,height_ratios=[1.25,1.0,.72],left=.145,right=.95,bottom=.075,top=.975,hspace=.56,wspace=.72)
    anatomy=read("Figure5A_anatomy_fractional_counts.tsv");ax=fig.add_subplot(gs[0,:2]);panel(ax,"a","Disease anatomy with fractional multi-label counts");anatomy_panel(ax,anatomy)
    rt=read("Figure5B_role_therapeutic_area_permutation.tsv");roles=rt.groupby("membrane_role")["observed"].sum().nlargest(6).index;tas=rt.groupby("therapeutic_area")["observed"].sum().nlargest(6).index;show=rt[rt["membrane_role"].isin(roles)&rt["therapeutic_area"].isin(tas)];mat=show.pivot(index="membrane_role",columns="therapeutic_area",values="permutation_z").reindex(index=roles,columns=tas);q=show.pivot(index="membrane_role",columns="therapeutic_area",values="bh_q").reindex(index=roles,columns=tas)
    ax=fig.add_subplot(gs[0,2:]);panel(ax,"b","Membrane-role × therapeutic-area preference")
    lim=max(2,float(np.nanpercentile(np.abs(mat),95)));sns.heatmap(mat,ax=ax,cmap=DIVERGING,center=0,vmin=-lim,vmax=lim,linewidths=.25,linecolor="white",cbar_kws={"label":"Permutation Z","shrink":.70})
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(q.iat[i,j]) and q.iat[i,j]<.05:ax.text(j+.5,i+.5,"•",ha="center",va="center",fontsize=6.3)
    ta_display={"genetic, familial or congenital disease":"genetic / congenital","cancer or benign tumor":"cancer / benign tumor","nervous system disorder":"nervous system","musculoskeletal or connective tissue disease":"musculoskeletal / connective","nutritional or metabolic disease":"metabolic / nutritional"}
    role_display={"membrane_associated_enzyme":"membrane-associated enzyme","family_defined_membrane_role_unresolved":"family-defined; role unresolved","membrane_scaffold_or_linker":"scaffold / linker","signaling_regulator":"signaling regulator"}
    ax.set_xticklabels([wrap(ta_display.get(x,x),15) for x in mat.columns],rotation=32,ha="right",rotation_mode="anchor");ax.set_yticklabels([wrap(role_display.get(x,x.replace('_',' ')),18) for x in mat.index],rotation=0);ax.set_xlabel("");ax.set_ylabel("")
    ax.text(.99,.985,"10,000 degree-constrained permutations; dot = BH q < 0.05",ha="right",va="top",transform=ax.transAxes,fontsize=5.0,bbox={"fc":"white","ec":"none","alpha":.82,"pad":1})
    upset=read("Figure5C_disease_channel_upset.tsv")
    channel_display={"human_genetic":"human genetics","clinical_genetic":"clinical genetics","somatic_mutation":"somatic mutation","animal_model":"animal model","uniprot_disease":"UniProt disease"}
    upset=upset.rename(columns=channel_display)
    outer=gs[1,:2].subgridspec(2,1,height_ratios=[2.1,1]);axb=fig.add_subplot(outer[0]);panel(axb,"c","Disease-evidence channel combinations");axm=fig.add_subplot(outer[1]);sets=[c for c in upset.columns if c not in {"intersection","relation_count"}];draw_upset(axb,axm,upset,sets,"relation_count",10)
    land=read("Figure5D_disease_chemical_landscape.tsv");land["x"]=np.log10(land["formal_compound_count"]+1);land["y"]=np.log10(land["high_very_high_disease_count"]+1);land["mean_anatomical_expression_concordance"]=pd.to_numeric(land["mean_anatomical_expression_concordance"],errors="coerce")
    ax=fig.add_subplot(gs[1,2:]);panel(ax,"d","Disease support versus chemical exploration")
    no=land[land["coordinate_ready_site_count"].lt(1)];yes=land[land["coordinate_ready_site_count"].ge(1)]
    sc=ax.scatter(no["x"],no["y"],c=no["mean_anatomical_expression_concordance"],cmap=sns.light_palette(EXPRESSION,as_cmap=True),vmin=0,vmax=1,s=9,alpha=.32,marker="o",edgecolors="none",rasterized=True)
    ax.scatter(yes["x"],yes["y"],c=yes["mean_anatomical_expression_concordance"],cmap=sns.light_palette(EXPRESSION,as_cmap=True),vmin=0,vmax=1,s=17,alpha=.55,marker="D",edgecolors="white",lw=.2,rasterized=True)
    cand=land[land["candidate_flag"].eq(1)];ax.scatter(cand["x"],cand["y"],facecolors="none",edgecolors=DISEASE,s=42,lw=.8)
    label_offsets=[(5,10),(18,9),(5,-12)]
    for r,off in zip(cand.sort_values(["high_very_high_disease_count","formal_compound_count"],ascending=[False,True]).head(3).itertuples(),label_offsets):
        ax.annotate(str(r.approved_symbol),(r.x,r.y),xytext=off,textcoords="offset points",fontsize=5.0)
    fig.colorbar(sc,ax=ax,label="Anatomical-expression concordance",shrink=.70);ax.set_xlabel("log10(unique formal compounds + 1)");ax.set_ylabel("log10(High/Very-High diseases + 1)");ax.text(.02,.98,"Diamond: coordinate-ready site; circle: no coordinate-ready site\nCoral outline: transparent candidate rule",ha="left",va="top",transform=ax.transAxes,fontsize=5.0)
    null=read("Figure5E_anatomical_concordance_null.tsv");stats=json.loads((ROOT/"02_analysis_screening"/"anatomical_expression_concordance"/"statistics.json").read_text(encoding="utf-8"));levels=read("Figure5C_disease_evidence_levels.tsv")
    ax=fig.add_subplot(gs[2,:]);panel(ax,"e","Anatomical-expression concordance under a constrained null")
    sns.histplot(null["null_mean_concordance"],bins=38,color="#A1D8E8",edgecolor="white",ax=ax);ax.axvline(stats["observed_mean"],color=DISEASE,lw=1.3,label=f"Observed={stats['observed_mean']:.3f}");ax.axvline(stats["null_mean"],color=PROTEIN,lw=.9,ls="--",label=f"Null mean={stats['null_mean']:.3f}");ax.set_xlabel("Mean anatomical-expression concordance");ax.set_ylabel("Permutations");ax.legend(loc="upper left")
    ax.text(.98,.92,f"z={stats['z_score']:.2f}\nempirical P<1e-4\n10,000 permutations",ha="right",va="top",transform=ax.transAxes,fontsize=5.3)
    lev=" · ".join(f"{r.evidence_level.replace('_',' ').title()} {int(r.relation_count):,}" for r in levels.itertuples())
    ax.text(.98,.10,"Disease evidence: "+lev,ha="right",va="bottom",transform=ax.transAxes,fontsize=5.0,color="#5F6368")
    source_note(fig,"V7.2 canonical diseases, ontology-derived anatomy, Open Targets therapeutic areas and mapped+measured normal HPA tissue RNA. Concordance is not causality.")
    save_main(fig,"Figure5_disease_context_and_hypothesis_generation")


def supplementary() -> None:
    cov=read("Supplementary_role_module_coverage.tsv");mat=cov.pivot(index="membrane_role",columns="module",values="coverage_fraction")
    fig,ax=plt.subplots(figsize=(178*MM,135*MM),layout="constrained");panel(ax,"a","Module coverage varies across membrane roles")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v=mat.iat[i,j];ax.scatter(j,i,s=12+110*v,c=v,cmap=sns.light_palette(PROTEIN,as_cmap=True),vmin=0,vmax=1,edgecolor="white",lw=.25)
    ax.set_xticks(range(len(mat.columns)),[wrap(x,18) for x in mat.columns],rotation=35,ha="right",rotation_mode="anchor");ax.set_yticks(range(len(mat.index)),[wrap(x,26) for x in mat.index]);ax.grid(False);fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0,1),cmap=sns.light_palette(PROTEIN,as_cmap=True)),ax=ax,label="Protein coverage fraction")
    source_note(fig,"Coverage describes MemPro module availability, not biological enrichment.");save_supp(fig,"Supplementary_Figure_S1_role_module_coverage")
    phys=read("Supplementary_physicochemical_by_status.tsv")
    def _robust_z(s):
        med=float(s.median()); mad=float((s-med).abs().median()); scale=1.4826*mad
        return (s-med)/scale if scale>0 else s*0
    phys["robust_z"]=phys.groupby("property")["value"].transform(_robust_z).clip(-6,6)
    fig,ax=plt.subplots(figsize=(178*MM,125*MM),layout="constrained");panel(ax,"a","Physicochemical profiles by compound status")
    sns.violinplot(data=phys,x="robust_z",y="biological_status",hue="property",cut=0,inner=None,linewidth=.25,alpha=.28,ax=ax);ax.set_xlabel("Robust z-score within property");ax.set_ylabel("");ax.legend(title="Property",ncol=2);source_note(fig,"Standardized descriptive properties; not a drug-likeness filter.");save_supp(fig,"Supplementary_Figure_S2_physicochemical_profiles")
    weak=read("Supplementary_structure_role_weak_association.tsv");roles=weak.groupby("membrane_role")["observed"].sum().nlargest(9).index;classes=weak.groupby("broad_structure_class")["observed"].sum().nlargest(8).index;show=weak[weak["membrane_role"].isin(roles)&weak["broad_structure_class"].isin(classes)];mat=show.pivot(index="broad_structure_class",columns="membrane_role",values="log2_odds_ratio").reindex(index=classes,columns=roles)
    fig,ax=plt.subplots(figsize=(178*MM,140*MM),layout="constrained");panel(ax,"a","Broad structure-class association with membrane role is weak")
    sns.heatmap(mat,ax=ax,cmap=DIVERGING,center=0,cbar_kws={"label":"log2 odds ratio"});ax.set_xticklabels([wrap(x,17) for x in mat.columns],rotation=38,ha="right",rotation_mode="anchor");ax.set_yticklabels([wrap(x,22) for x in mat.index],rotation=0);ax.set_xlabel("");ax.set_ylabel("");ax.text(1,1.02,"Overall bias-corrected Cramer's V = 0.046",ha="right",transform=ax.transAxes,fontsize=5.3);save_supp(fig,"Supplementary_Figure_S3_weak_structure_role_association")
    lor=read("Supplementary_evidence_lorenz.tsv");fig,ax=plt.subplots(figsize=(84*MM,75*MM),layout="constrained");panel(ax,"a","Evidence concentration across formal pairs");ax.plot(lor["cumulative_pair_fraction"],lor["cumulative_evidence_fraction"],color=DISEASE,lw=1.1);ax.plot([0,1],[0,1],ls="--",color="#777",lw=.7);ax.fill_between(lor["cumulative_pair_fraction"],lor["cumulative_evidence_fraction"],lor["cumulative_pair_fraction"],color="#FFC1A6",alpha=.25);ax.set_xlabel("Cumulative fraction of pairs");ax.set_ylabel("Cumulative fraction of evidence");ax.text(.05,.90,"Gini=0.361",transform=ax.transAxes);save_supp(fig,"Supplementary_Figure_S4_evidence_lorenz")


def graphical_abstract() -> None:
    fig,ax=plt.subplots(figsize=(150*MM,60*MM));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    # Stage 1: grouped source modalities, shown as scientific symbols rather than database logos.
    ax.text(.105,.92,"SOURCES",ha="center",va="center",fontsize=7.5,fontweight="bold")
    src=[("protein",PROTEIN,"P"),("compound",COMPOUND,"C"),("evidence",DISEASE,"E"),("disease / tissue",EXPRESSION,"D")]
    for i,(lab,col,letter) in enumerate(src):
        yy=.75-i*.16
        ax.add_patch(patches.Circle((.035,yy),.028,fc=col,ec="white",lw=.8,alpha=.88))
        ax.text(.035,yy,letter,ha="center",va="center",fontsize=6.5,fontweight="bold",color="white")
        ax.text(.072,yy,lab,ha="left",va="center",fontsize=6.2)
    ax.annotate("",xy=(.245,.50),xytext=(.205,.50),arrowprops={"arrowstyle":"->","lw":1.4,"color":"#72787C"})
    # Stage 2: deterministic harmonisation funnel.
    ax.text(.35,.92,"HARMONIZATION",ha="center",va="center",fontsize=7.5,fontweight="bold")
    funnel=patches.Polygon([[.255,.78],[.445,.78],[.405,.58],[.385,.40],[.315,.40],[.295,.58]],closed=True,fc="#DCE9E3",ec=EXPRESSION,lw=1.0)
    ax.add_patch(funnel)
    ax.text(.35,.65,"identity",ha="center",va="center",fontsize=6.2)
    ax.text(.35,.54,"deduplicate",ha="center",va="center",fontsize=6.2)
    ax.text(.35,.43,"lineage + ontology",ha="center",va="center",fontsize=6.2)
    ax.text(.35,.27,"protein · form · MONDO",ha="center",va="center",fontsize=5.3,color="#5F666A")
    ax.annotate("",xy=(.505,.50),xytext=(.455,.50),arrowprops={"arrowstyle":"->","lw":1.4,"color":"#72787C"})
    # Stage 3: central MemPro knowledge structure.
    ax.text(.625,.92,"MemPro V7.2",ha="center",va="center",fontsize=7.8,fontweight="bold")
    for x in np.linspace(.51,.74,12):
        ax.add_patch(patches.Circle((x,.64),.011,fc="#9FD3D9",ec="none"));ax.plot([x,x],[.59,.48],color="#9FD3D9",lw=1.1)
        ax.add_patch(patches.Circle((x,.43),.011,fc="#9FD3D9",ec="none"))
    t=np.linspace(0,1,120);xp=.595+.035*np.sin(6*np.pi*t);yp=.40+.30*t
    ax.plot(xp,yp,color=PROTEIN,lw=5.0,solid_capstyle="round",zorder=3)
    hexagon=patches.RegularPolygon((.69,.56),6,radius=.045,orientation=np.pi/6,fc="#F6C38D",ec=COMPOUND,lw=1.2,zorder=4);ax.add_patch(hexagon)
    ax.add_patch(patches.Circle((.56,.27),.052,fc="#E8B0B2",ec=DISEASE,lw=1.1));ax.text(.56,.27,"disease",ha="center",va="center",fontsize=5.5)
    ax.plot([.60,.57],[.43,.31],color="#7B8185",lw=1);ax.plot([.65,.62],[.53,.34],color="#7B8185",lw=1)
    ax.text(.625,.13,"classification · expression · sites · complexes",ha="center",va="center",fontsize=5.5,color="#5F666A")
    ax.annotate("",xy=(.785,.50),xytext=(.755,.50),arrowprops={"arrowstyle":"->","lw":1.4,"color":"#72787C"})
    # Stage 4: interpretable outputs, deliberately separated instead of collapsed to one score.
    ax.text(.89,.92,"OUTPUTS",ha="center",va="center",fontsize=7.5,fontweight="bold")
    outs=[(.82,.66,PROTEIN,"target\npreference"),(.95,.66,COMPOUND,"chemical\npreference"),(.82,.34,DISEASE,"evidence\ncontext"),(.95,.34,EXPRESSION,"underexplored\ntargets")]
    for x,y,c,lab in outs:
        ax.add_patch(patches.Circle((x,y),.048,fc=c,ec="white",lw=.8,alpha=.85))
        ax.text(x,y-.075,lab,ha="center",va="top",fontsize=5.5,linespacing=1.0)
    for ext in ("pdf","svg","png","tif"):
        path=GA/f"graphical_abstract.{ext}";kwargs={}
        if ext in {"png","tif"}:kwargs["dpi"]=int(CFG["graphical_abstract_dpi"])
        fig.savefig(path,**kwargs)
    plt.close(fig)


def main() -> None:
    figure1();figure2();figure3();figure4();figure5();supplementary();graphical_abstract()
    print(json.dumps({"status":"complete","main_figures":5,"supplementary_figures":4,"graphical_abstract":1,"output":str(ROOT)},ensure_ascii=False))


if __name__=="__main__":main()

