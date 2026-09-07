from __future__ import annotations

import json
import math
import os
import re
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr


ROOT = Path(os.environ.get("MEMPRO_REVISION_ROOT", r"D:\finale\07_M1-M8图表修订_20260806"))
V62 = Path(r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730")
V63 = Path(r"D:\finale\02_V6.3_candidate_20260803")
OLD = Path(r"D:\finale\02_展示图表_V6.2")
PNG, SVG, PDF = ROOT / "03_figures/png", ROOT / "03_figures/svg", ROOT / "03_figures/pdf"
DATA, AUDIT, QA = ROOT / "02_data", ROOT / "01_source_audit", ROOT / "09_qa"
for p in (PNG, SVG, PDF, DATA, AUDIT, QA):
    p.mkdir(parents=True, exist_ok=True)

COL = {
    "navy": "#334E60", "blue": "#6F8798", "teal": "#4F8F88",
    "orange": "#C98257", "purple": "#80739E", "gray": "#9AA1A6",
    "light": "#E8ECEF", "text": "#2F3437", "red": "#B76E62",
    "A": "#496A81", "B": "#4F8F88", "C": "#C98257",
    "E1": "#496A81", "E2": "#4F8F88", "E3": "#C98257", "E0": "#A8ADB1",
    "BE1": "#496A81", "BE2": "#4F8F88", "BE3": "#C98257",
}


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Microsoft YaHei", "DejaVu Sans"],
        "font.size": 8.2, "axes.titlesize": 9.6, "axes.titleweight": "bold",
        "axes.labelsize": 8.3, "xtick.labelsize": 7.1, "ytick.labelsize": 7.1,
        "legend.fontsize": 7.0, "figure.titlesize": 13.5,
        "text.color": COL["text"], "axes.labelcolor": COL["text"],
        "xtick.color": "#454A4E", "ytick.color": "#454A4E",
        "axes.edgecolor": "#9AA1A6", "axes.linewidth": 0.6,
        "grid.color": "#D6D9DC", "grid.linewidth": 0.45, "grid.alpha": 0.55,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })


def fig6(title: str):
    setup_style()
    fig, ax = plt.subplots(2, 3, figsize=(14.2, 8.25), constrained_layout=True)
    fig.suptitle(title, x=0.015, y=1.015, ha="left", color=COL["navy"], fontweight="bold")
    return fig, ax.ravel()


def panel(ax, letter: str, title: str, source: str = ""):
    ax.text(-0.08, 1.08, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="top", ha="left", color=COL["navy"])
    ax.set_title(title, loc="left", pad=8)
    if source:
        ax.text(0.0, -0.20, f"Source: {source}", transform=ax.transAxes, fontsize=6.2,
                color="#666D72", va="top", clip_on=False)


def clean(ax, axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)
    ax.grid(True, axis=axis, zorder=0)


def save(fig, stem: str):
    paths = {"png": PNG / f"{stem}.png", "svg": SVG / f"{stem}.svg", "pdf": PDF / f"{stem}.pdf"}
    fig.savefig(paths["png"], dpi=600, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)
    return {k: str(v) for k, v in paths.items()}


def label_barh(ax, vals, fmt=lambda x: f"{int(x):,}"):
    vmax = max(vals) if len(vals) else 1
    for i, v in enumerate(vals):
        ax.text(v + vmax * .012, i, fmt(v), va="center", fontsize=6.4, color=COL["text"])


def truth(s):
    return s.astype(str).str.lower().isin(["1", "true", "yes", "y"])


def hnum(x):
    if x >= 1_000_000: return f"{x/1_000_000:.1f}M"
    if x >= 1_000: return f"{x/1_000:.1f}k"
    return f"{x:.0f}"


def load_protein():
    cols = [
        "target_uniprot_id", "membrane_class_v52", "evidence_level_v52",
        "functional_primary_class_v5", "transmembrane_count_v5", "sequence_length_v53",
        "opm_present_v5", "pdbtm_present_v5", "pdb_ids", "alphafolddb_ids",
        "oligomeric_state_consensus_v62", "subunit_evidence_grade_v62",
        "hpa_mapping_status_v62", "hpa_tissue_detected_count_v62",
        "hpa_cell_type_detected_count_v62", "hpa_ihc_detected_tissue_count_v62",
        "hpa_main_locations_v62", "hpa_additional_locations_v62",
    ]
    p = pd.read_csv(V62 / "human_membrane_protein_master_v6_2.tsv", sep="\t", usecols=cols, low_memory=False)
    p["tm"] = pd.to_numeric(p["transmembrane_count_v5"], errors="coerce").fillna(0)
    p["length"] = pd.to_numeric(p["sequence_length_v53"], errors="coerce")
    p["class"] = p["membrane_class_v52"].fillna("unknown")
    p["grade"] = p["evidence_level_v52"].fillna("E0")
    p["function"] = p["functional_primary_class_v5"].fillna("unclassified")
    membrane_pdb = truth(p["opm_present_v5"]) | truth(p["pdbtm_present_v5"])
    any_pdb = p["pdb_ids"].fillna("").ne("")
    any_af = p["alphafolddb_ids"].fillna("").ne("")
    p["structure"] = np.select([membrane_pdb, any_pdb, any_af],
        ["Membrane PDB", "Other PDB", "AlphaFold only"], default="No structure")
    return p


def source_audit():
    source_map = {
        "M1": ("UniProtKB; HPA; UniTmp HTP/PDBTM; Membranome; PDBe; OPM; PDBbind; ChEMBL; BindingDB; PubChem BioAssay; BRENDA; Open Targets", "records/entities/relations", "integration architecture and release counts"),
        "M2": ("UniProtKB 2026_02; HPA 25.1; UniTmp HTP d.2.2; OPM/PDBTM; PDBe; AlphaFoldDB", "unique canonical proteins", "class, evidence, topology, sequence length, structure and assembly"),
        "M3": ("Human Protein Atlas 25.1", "unique proteins or normalized RNA units", "tissue RNA nTPM; cell-type RNA nCPM; IHC ordinal protein staining; subcellular localization"),
        "M4": ("Open Targets 26.06; UniProtKB; MONDO; DOID; Uberon", "unique canonical protein-disease pairs", "exact-only disease identity, evidence channels, therapeutic area and anatomy multi-label mapping"),
        "M5": ("MemPro compound master v1.3; PubChem; ChEMBL; ChEBI; BindingDB", "unique canonical compounds", "identity, status, RDKit structural rules, descriptors, scaffolds and learned embeddings"),
        "M6": ("PDBe; PDBbind; BioLiP; BindingDB; UniProtKB; OPM; PDBTM; SIFTS/PDBe mapping; UniTmp", "binding evidence or site instances", "BE tiers, site readiness, residue contacts and membrane-side classification"),
        "M7": ("V6.2 positive/negative evidence, review queues and QA", "records or canonical pairs", "negative identity, conflicts, missingness and release validation"),
        "M8": ("V6.3 candidate lineage summary and V6.3.1 docking readiness", "protein-compound candidate pairs", "structure/site/identity/conflict-aware docking tiers"),
    }
    rows = []
    titles = {
        "M1": ["Source architecture","Integration pipeline","Record disposition","Relational model","Source-module coverage","Release snapshot"],
        "M2": ["Membrane class/evidence","Functional hierarchy","TM topology","Length-TM sensitivity","Structure coverage","Assembly state"],
        "M3": ["HPA mapping","Measurement layers","Tissue RNA","Cell-type RNA","Subcellular location","Expression breadth/missingness"],
        "M4": ["Disease sources","Evidence tier","Evidence channels","Disease identity","Anatomical systems","Disease-anatomy matrix"],
        "M5": ["Identity confidence","Biological status","Structural class","Descriptor profile","PAGTN experiment","Scaffold diversity"],
        "M6": ["Source contribution","BE tiers","Site readiness","Membrane side","Contact residue classes","Docking box readiness"],
        "M7": ["Negative identity","Positive-negative conflict","Review queues","Missingness","QA gates","Release accounting"],
        "M8": ["Docking tiers","Structure readiness","Ligand readiness","Site readiness","Score distribution","HPC workload"],
    }
    for m, vals in titles.items():
        src, unit, desc = source_map[m]
        for i, title in enumerate(vals):
            rows.append({"panel": f"{m}{chr(65+i)}", "panel_title": title, "source_databases": src,
                         "source_files": "V6.2 frozen release; V6.3 candidate where explicitly labeled",
                         "version_snapshot": "V6.2 2026-07-30 / V6.3 candidate 2026-08-03",
                         "statistical_unit": unit, "scope_or_filter": desc,
                         "status": "revised_2026-08-06"})
    df = pd.DataFrame(rows)
    df.to_csv(AUDIT / "PANEL_SOURCE_FIELD_VERSION_UNIT.tsv", sep="\t", index=False)
    old_versions = OLD / "data/M0_source_versions_contributions.tsv"
    if old_versions.exists():
        pd.read_csv(old_versions, sep="\t").to_csv(AUDIT / "DATABASE_VERSION_AND_CONTRIBUTION_BASELINE.tsv", sep="\t", index=False)
    return df


def make_m1(protein):
    fig, ax = fig6("M1  Database architecture, integration and traceable release control")
    panel(ax[0], "A", "Four-source architecture", "source registries and frozen snapshots")
    ax[0].axis("off")
    groups = [("Protein", "UniProt · HPA · HTP · Membranome", COL["blue"]),
              ("Structure", "PDBe · OPM · PDBTM · PDBbind", COL["teal"]),
              ("Interaction", "ChEMBL · BindingDB · PubChem · BRENDA", COL["orange"]),
              ("Disease", "Open Targets · UniProtKB · MONDO/DO/Uberon", COL["purple"])]
    for i,(t,s,c) in enumerate(groups):
        y=.83-i*.22
        ax[0].text(.04,y,t,fontweight="bold",fontsize=8,color=COL["text"])
        ax[0].text(.04,y-.075,s,fontsize=6.6,color="#596166")
        ax[0].plot([.04,.92],[y-.105,y-.105],color=c,lw=5,alpha=.48,solid_capstyle="round")
    panel(ax[1], "B", "Integration is a controlled sequence", "MemPro processing policy")
    ax[1].axis("off")
    steps=["Raw\nrecords","ID\nmapping","Parent /\nform","Cross-source\ndedup","Evidence\ntier","Release"]
    xs=np.linspace(.08,.92,len(steps))
    for i,(x,t) in enumerate(zip(xs,steps)):
        ax[1].scatter(x,.52,s=1050,color=sns.color_palette("flare",len(steps),desat=.65)[i],alpha=.72,edgecolor="white")
        ax[1].text(x,.52,t,ha="center",va="center",fontsize=6.5)
        if i<len(steps)-1: ax[1].annotate("",(xs[i+1]-.07,.52),(x+.07,.52),arrowprops={"arrowstyle":"->","lw":.8,"color":COL["navy"]})
    panel(ax[2], "C", "Positive and negative records remain separate", "V6.2 release accounting")
    vals=pd.DataFrame({"Released":[1_575_210,7_630_659],"Review":[1_619_083,2_615_046],"Excluded/duplicate":[96_221,1_243]},index=["Positive input","Negative input"])
    left=np.zeros(2)
    for k,c in zip(vals.columns,[COL["teal"],COL["purple"],COL["gray"]]):
        ax[2].barh(vals.index,vals[k],left=left,color=c,alpha=.86,label=k); left+=vals[k].values
    ax[2].set_xlabel("Records"); ax[2].legend(frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(.5,-.36)); clean(ax[2],"x")
    panel(ax[3], "D", "Entity model prevents forced identity merges", "V6.3 target and disease identity layers")
    ax[3].axis("off")
    nodes=[("Gene",.18,.73),("Canonical\nprotein",.5,.73),("Isoform",.82,.73),("Complex",.18,.28),("Compound\nparent/form",.5,.28),("Disease\nMONDO",.82,.28)]
    for t,x,y in nodes:
        ax[3].scatter(x,y,s=2100,color=COL["light"],edgecolor=COL["navy"],linewidth=.8)
        ax[3].text(x,y,t,ha="center",va="center",fontsize=7)
    for a,b in [(0,1),(1,2),(1,3),(1,4),(1,5)]:
        x1,y1=nodes[a][1:]; x2,y2=nodes[b][1:]; ax[3].plot([x1,x2],[y1,y2],color="#A4AAAE",lw=.8,zorder=0)
    panel(ax[4], "E", "Contributing databases cover distinct modules", "M1 source-module matrix")
    mat=pd.read_csv(OLD / "data/M1_source_module_matrix.tsv",sep="\t").set_index("source")
    ann=mat.map(lambda v:"" if v==0 else hnum(v))
    sns.heatmap(np.log10(mat+1),cmap=sns.light_palette(COL["blue"],as_cmap=True),annot=ann,fmt="",linewidths=.35,cbar_kws={"label":"log10(count+1)"},ax=ax[4])
    ax[4].set_xlabel(""); ax[4].set_ylabel(""); ax[4].tick_params(axis="x",rotation=35)
    panel(ax[5], "F", "Frozen V6.2 scale", "V6.2 manifest and validation report")
    metrics=pd.Series({"Membrane proteins":10_997,"Disease relations":10_264,"Binding sites":95_598,"Positive pairs":1_502_456,"Canonical compounds":2_016_064,"Negative pairs":6_331_307}).sort_values()
    ax[5].barh(metrics.index,metrics.values,color=COL["blue"],alpha=.86); ax[5].set_xscale("log"); ax[5].set_xlabel("Entities / relations (log)"); label_barh(ax[5],metrics.values); clean(ax[5],"x")
    return save(fig,"M1_database_architecture_quality_revised")


def make_m2(protein):
    fig, ax=fig6("M2  Membrane proteome composition, topology and structure coverage")
    panel(ax[0],"A","Membrane class × evidence grade","UniProtKB/HPA/HTP/Membranome/OPM/PDBTM")
    cross=pd.crosstab(protein["class"],protein["grade"]).reindex(index=["A","B","C","unknown"],columns=["E1","E2","E3","E0"],fill_value=0)
    left=np.zeros(len(cross))
    for g in cross.columns:
        ax[0].barh(cross.index,cross[g],left=left,color=COL[g],alpha=.88,label=g); left+=cross[g].values
    ax[0].set_xlabel("Unique canonical proteins"); ax[0].legend(frameon=False,ncol=4,loc="lower center",bbox_to_anchor=(.5,-.31)); clean(ax[0],"x")
    panel(ax[1],"B","Functional hierarchy","UniProtKB and V6.3 cross-classification")
    c=protein["function"].value_counts().head(11).sort_values(); ax[1].barh(c.index.str.replace("_"," "),c.values,color=COL["teal"],alpha=.86); label_barh(ax[1],c.values); ax[1].set_xlabel("Unique proteins"); clean(ax[1],"x")
    panel(ax[2],"C","Transmembrane architecture","UniProtKB/UniTmp topology")
    bins=pd.cut(protein.tm,[-.1,.5,1.5,4.5,8.5,12.5,20.5,np.inf],labels=["0","1","2–4","5–8","9–12","13–20",">20"])
    c=bins.value_counts().sort_index(); ax[2].bar(c.index.astype(str),c.values,color=sns.color_palette("flare",len(c),desat=.72),alpha=.88); ax[2].set_ylabel("Unique proteins"); ax[2].set_xlabel("TM helices"); clean(ax[2],"y")
    panel(ax[3],"D","Length–TM relationship after biological stratification","A-class, E1/E2, canonical multipass proteins")
    base=protein[(protein["class"].eq("A")) & protein["grade"].isin(["E1","E2"]) & protein.tm.between(2,40) & protein.length.between(50,5000)].dropna(subset=["length","tm"])
    hb=ax[3].hexbin(base.tm,base.length,gridsize=(22,35),bins="log",mincnt=1,cmap=sns.light_palette(COL["teal"],as_cmap=True),alpha=.82)
    ax[3].set_yscale("log"); ax[3].set_xlabel("TM helices"); ax[3].set_ylabel("Sequence length (aa, log)")
    rho,p=spearmanr(base.tm,base.length)
    ax[3].text(.97,.96,f"Spearman ρ={rho:.3f}\nn={len(base):,}\n50–5,000 aa; 2–40 TM",transform=ax[3].transAxes,ha="right",va="top",fontsize=7,bbox={"boxstyle":"round,pad=.3","fc":"white","ec":"#D6D9DC"})
    clean(ax[3],"y")
    sens=[]
    for name,q in [("All A/B/C",protein.length.notna()),("A, TM≥1",protein["class"].eq("A")&protein.tm.ge(1)),("A multipass",protein["class"].eq("A")&protein.tm.ge(2)),("Main scope",base.index.to_series().isin(base.index))]:
        d=base if name=="Main scope" else protein.loc[q,["length","tm"]].dropna()
        r=spearmanr(d.length,d.tm).statistic if len(d)>2 else np.nan
        sens.append({"analysis_scope":name,"n":len(d),"spearman_rho":r})
    pd.DataFrame(sens).to_csv(DATA/"M2D_spearman_sensitivity.tsv",sep="\t",index=False)
    panel(ax[4],"E","Structure coverage retains AlphaFold","PDB/OPM/PDBTM and AlphaFoldDB identifiers")
    top=protein.function.value_counts().head(9).index; d=protein[protein.function.isin(top)]
    st=pd.crosstab(d.function,d.structure,normalize="index").reindex(index=top,columns=["Membrane PDB","Other PDB","AlphaFold only","No structure"],fill_value=0)*100
    sns.heatmap(st,cmap=sns.light_palette(COL["blue"],as_cmap=True),annot=True,fmt=".0f",linewidths=.35,cbar_kws={"label":"% of class"},ax=ax[4]); ax[4].set_xlabel(""); ax[4].set_ylabel(""); ax[4].set_xticklabels(["Membrane\nPDB","Other\nPDB","AlphaFold\nonly","No\nstructure"],rotation=0); ax[4].set_yticklabels([x.replace("_"," ") for x in st.index],rotation=0)
    panel(ax[5],"F","Biological assembly evidence","UniProtKB and PDBe biological assemblies")
    c=protein.oligomeric_state_consensus_v62.fillna("unresolved").value_counts().head(10).sort_values(); ax[5].barh(c.index.str.replace("_"," "),c.values,color=COL["purple"],alpha=.82); label_barh(ax[5],c.values); ax[5].set_xlabel("Unique proteins"); clean(ax[5],"x")
    return save(fig,"M2_membrane_proteome_landscape_revised")


def make_m3(protein):
    fig,ax=fig6("M3  HPA expression and localization with explicit measurement layers")
    panel(ax[0],"A","HPA mapping status","HPA 25.1; Ensembl gene → canonical UniProt mapping")
    c=protein.hpa_mapping_status_v62.fillna("missing/unmapped").replace("","missing/unmapped").value_counts().sort_values(); ax[0].barh(c.index,c.values,color=[COL["gray"] if "missing" in x else COL["teal"] for x in c.index],alpha=.86); label_barh(ax[0],c.values); ax[0].set_xlabel("Unique proteins"); clean(ax[0],"x")
    panel(ax[1],"B","RNA and IHC are different measurement layers","HPA 25.1")
    vals=pd.Series({"Tissue RNA\n(nTPM)":protein.hpa_tissue_detected_count_v62.gt(0).sum(),"Cell-type RNA\n(nCPM)":protein.hpa_cell_type_detected_count_v62.gt(0).sum(),"IHC protein\n(ordinal staining)":protein.hpa_ihc_detected_tissue_count_v62.gt(0).sum()})
    ax[1].bar(vals.index,vals.values,color=[COL["blue"],COL["teal"],COL["orange"]],alpha=.84); ax[1].set_ylabel("Proteins with ≥1 detected context"); ax[1].tick_params(axis="x",rotation=12); clean(ax[1],"y")
    tm=pd.read_csv(OLD/"data/M3_tissue_expression_matrix.tsv",sep="\t").set_index("functional_class")
    z=tm.sub(tm.mean(axis=1),axis=0).div(tm.std(axis=1).replace(0,np.nan),axis=0).fillna(0)
    panel(ax[2],"C","Tissue RNA relative pattern","HPA consensus tissue RNA; nTPM; row-wise z-score")
    sns.heatmap(z,cmap="vlag",center=0,vmin=-2.3,vmax=2.3,linewidths=.25,cbar_kws={"label":"within-class z-score"},ax=ax[2]); ax[2].set_xlabel(""); ax[2].set_ylabel(""); ax[2].tick_params(axis="x",rotation=35)
    cm=pd.read_csv(OLD/"data/M3_cell_type_expression_matrix.tsv",sep="\t").set_index("functional_class")
    cz=cm.sub(cm.mean(axis=1),axis=0).div(cm.std(axis=1).replace(0,np.nan),axis=0).fillna(0)
    panel(ax[3],"D","Cell-type RNA relative pattern","HPA single-cell type RNA; nCPM; row-wise z-score")
    sns.heatmap(cz,cmap="vlag",center=0,vmin=-2.3,vmax=2.3,linewidths=.25,cbar_kws={"label":"within-class z-score"},ax=ax[3]); ax[3].set_xlabel(""); ax[3].set_ylabel(""); ax[3].tick_params(axis="x",rotation=35)
    loc=pd.read_csv(V62/"protein_subcellular_localization_v2.tsv.gz",sep="\t",compression="gzip",usecols=["target_uniprot_id","location_term"],low_memory=False).drop_duplicates()
    c=loc.location_term.value_counts().head(12).sort_values()
    panel(ax[4],"E","Subcellular localization is multi-label","HPA subcellular localization")
    ax[4].barh(c.index,c.values,color=COL["blue"],alpha=.84); label_barh(ax[4],c.values); ax[4].set_xlabel("Unique proteins"); clean(ax[4],"x")
    panel(ax[5],"F","Expression breadth separates zero from missing","HPA 25.1 mapping-aware summaries")
    mapped=protein.hpa_mapping_status_v62.eq("mapped")
    dd=pd.DataFrame({"Detected tissues":protein.loc[mapped,"hpa_tissue_detected_count_v62"],"Detected cell types":protein.loc[mapped,"hpa_cell_type_detected_count_v62"]}).melt(var_name="layer",value_name="count")
    sns.violinplot(data=dd,x="layer",y="count",inner=None,color=COL["teal"],alpha=.58,cut=0,ax=ax[5]); sns.boxplot(data=dd,x="layer",y="count",width=.22,showfliers=False,boxprops={"facecolor":"white","alpha":.82},ax=ax[5])
    missing=(~mapped).sum()/len(protein)*100; ax[5].text(.98,.97,f"HPA unmapped/missing: {missing:.1f}%",transform=ax[5].transAxes,ha="right",va="top",fontsize=7)
    ax[5].set_xlabel(""); ax[5].set_ylabel("Detected contexts per mapped protein"); clean(ax[5],"y")
    return save(fig,"M3_expression_localization_atlas_revised")


def make_m4():
    rel=pd.read_csv(V63/"04_disease/protein_disease_relation_v2_1.tsv",sep="\t",low_memory=False)
    rel=rel[rel.default_release_inclusion_v63.eq(1)].copy()
    src=pd.read_csv(V63/"04_disease/disease_source_entity_v0_1.tsv",sep="\t",low_memory=False)
    systems=pd.read_csv(V63/"04_disease/disease_anatomical_system_v0_1.tsv",sep="\t",low_memory=False)
    anat=pd.read_csv(V63/"04_disease/disease_anatomy_v0_1.tsv",sep="\t",low_memory=False)
    fig,ax=fig6("M4  Canonical disease identity, evidence and ontology-derived anatomy")
    panel(ax[0],"A","Disease relation sources","Open Targets 26.06 and UniProtKB")
    c=rel.source_databases.fillna("").str.split(";").explode().replace("",np.nan).dropna().value_counts().sort_values(); ax[0].barh(c.index,c.values,color=COL["purple"],alpha=.84); label_barh(ax[0],c.values); ax[0].set_xlabel("Canonical protein–disease pairs"); clean(ax[0],"x")
    panel(ax[1],"B","Evidence level by therapeutic area","Open Targets 26.06 multi-label therapeutic areas")
    ta=pd.read_csv(V63/"04_disease/disease_therapeutic_area_v0_1.tsv",sep="\t",low_memory=False)
    dm=rel[["canonical_disease_id","best_evidence_level","disease_relation_id_v2"]].merge(ta,on="canonical_disease_id",how="inner")
    name_col="therapeutic_area_name" if "therapeutic_area_name" in dm else [x for x in dm.columns if "area" in x and "name" in x][0]
    tab=pd.crosstab(dm[name_col],dm.best_evidence_level).reindex(columns=["very_high","high","medium"],fill_value=0); tab=tab.loc[tab.sum(axis=1).nlargest(10).index].sort_values(tab.columns.tolist())
    left=np.zeros(len(tab))
    for g,c0 in zip(tab.columns,[COL["navy"],COL["teal"],COL["orange"]]): ax[1].barh(tab.index,tab[g],left=left,label=g.replace("_"," "),color=c0,alpha=.85); left+=tab[g].values
    ax[1].set_xlabel("Unique protein–disease pairs"); ax[1].legend(frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(.5,-.34)); clean(ax[1],"x")
    panel(ax[2],"C","Evidence channels are not interchangeable","Open Targets evidence channels and UniProt disease annotation")
    flags=[("Human genetics","human_genetic_flag"),("Clinical genetics","clinical_genetic_flag"),("Somatic mutation","somatic_mutation_flag"),("Functional","functional_flag"),("Animal model","animal_model_flag"),("Expression","expression_flag"),("Literature","literature_flag"),("UniProt disease","uniprot_disease_flag")]
    cc=pd.Series({n:pd.to_numeric(rel[c],errors="coerce").fillna(0).gt(0).mean()*100 for n,c in flags}).sort_values(); ax[2].barh(cc.index,cc.values,color=COL["teal"],alpha=.84); label_barh(ax[2],cc.values,lambda x:f"{x:.1f}%"); ax[2].set_xlabel("Relations carrying channel (%)"); clean(ax[2],"x")
    panel(ax[3],"D","Disease identity uses official exact/equivalent mappings","MONDO/OMIM/Orphanet/EFO cross-references")
    c=src.canonicalization_status.value_counts().head(8).sort_values(); ax[3].barh(c.index.str.replace("_"," "),c.values,color=[COL["teal"] if "exact" in x or "direct" in x else COL["orange"] for x in c.index],alpha=.84); label_barh(ax[3],c.values); ax[3].set_xlabel("Source disease entities"); clean(ax[3],"x")
    panel(ax[4],"E","Anatomical systems are ontology-derived and multi-label","MONDO/DO → Uberon; no keyword classification")
    pairs=rel[["target_uniprot_id","canonical_disease_id"]].drop_duplicates().merge(systems,on="canonical_disease_id")
    c=pairs.drop_duplicates(["target_uniprot_id","canonical_disease_id","anatomical_system_id"]).anatomical_system_name.value_counts().head(12).sort_values(); ax[4].barh(c.index,c.values,color=COL["blue"],alpha=.84); label_barh(ax[4],c.values); ax[4].set_xlabel("Unique protein–disease memberships"); clean(ax[4],"x")
    panel(ax[5],"F","Disease tissue location uses explicit ontology links","MONDO/DO disease anatomy → Uberon tissue/organ")
    top_sys=c.index[-8:]
    da=rel[["target_uniprot_id","canonical_disease_id"]].drop_duplicates().merge(systems[["canonical_disease_id","anatomical_system_name"]].drop_duplicates(),on="canonical_disease_id").merge(anat[["canonical_disease_id","anatomy_name"]].drop_duplicates(),on="canonical_disease_id")
    da=da[da.anatomical_system_name.isin(top_sys)]; top_an=da.anatomy_name.value_counts().head(10).index
    mat=pd.crosstab(da[da.anatomy_name.isin(top_an)].anatomical_system_name,da[da.anatomy_name.isin(top_an)].anatomy_name).reindex(index=top_sys,columns=top_an,fill_value=0)
    sns.heatmap(np.log1p(mat),cmap="rocket",annot=mat,fmt="d",linewidths=.3,cbar_kws={"label":"log1p(pair memberships)"},ax=ax[5]); ax[5].set_xlabel(""); ax[5].set_ylabel(""); ax[5].tick_params(axis="x",rotation=40)
    mat.to_csv(DATA/"M4F_ontology_disease_anatomy_matrix.tsv",sep="\t")
    return save(fig,"M4_disease_association_landscape_revised")


def structural_class(row):
    formula=str(row.molecular_formula or "")
    rings=float(row.ring_count) if pd.notna(row.ring_count) else 0
    maxring=float(row.max_ring_size) if pd.notna(row.max_ring_size) else 0
    heavy=float(row.heavy_atom_count) if pd.notna(row.heavy_atom_count) else 0
    amide=float(row.amide_bond_count) if pd.notna(row.amide_bond_count) else 0
    if not re.search(r"C",formula): return "inorganic/no carbon"
    if heavy>=25 and amide>=4: return "peptide-like"
    if maxring>=12: return "macrocycle"
    if rings>=4: return "polycyclic (≥4 rings)"
    if rings>=2: return "polycyclic (2–3 rings)"
    if rings==1: return "monocyclic"
    return "acyclic organic"


def scan_compounds():
    path=V62/"small_molecule_master_v1_3.tsv"
    cols=["compound_internal_id","compound_scope_status","record_qc_status","identity_confidence","molecular_formula","molecular_weight","xlogp","tpsa","hbond_donor_count","hbond_acceptor_count","rotatable_bond_count","heavy_atom_count","ring_count","max_ring_size","amide_bond_count","is_approved_drug","is_clinical_candidate","is_endogenous_ligand","is_natural_product","is_chemical_probe"]
    identity=Counter(); classes=Counter(); status=Counter(); desc=[]; n=0
    for ch in pd.read_csv(path,sep="\t",usecols=cols,chunksize=160000,low_memory=False):
        core=ch[ch.compound_scope_status.eq("core") & ch.record_qc_status.eq("ok")].copy(); n+=len(core)
        identity.update(core.identity_confidence.fillna("missing"))
        classes.update(core.apply(structural_class,axis=1))
        for label,col in [("Approved", "is_approved_drug"),("Clinical", "is_clinical_candidate"),("Endogenous", "is_endogenous_ligand"),("Natural", "is_natural_product"),("Probe", "is_chemical_probe")]: status[label]+=int(truth(core[col]).sum())
        if sum(len(x) for x in desc)<50000:
            d=core[["molecular_weight","xlogp","tpsa","rotatable_bond_count"]].dropna().sample(min(5000,len(core.dropna())),random_state=42) if len(core) else core
            desc.append(d)
    return n,pd.Series(identity),pd.Series(classes),pd.Series(status),pd.concat(desc,ignore_index=True).head(50000)


def make_m5():
    n,identity,classes,status,desc=scan_compounds()
    fig,ax=fig6("M5  Compound identity, structural diversity and modeling readiness")
    panel(ax[0],"A","Canonical identity confidence","MemPro compound master v1.3; InChIKey/SMILES identity policy")
    c=identity.sort_values(); ax[0].barh(c.index,c.values,color=COL["blue"],alpha=.84); label_barh(ax[0],c.values); ax[0].set_xlabel("Core QC canonical compounds"); clean(ax[0],"x")
    panel(ax[1],"B","Biological-status annotations overlap","ChEMBL/ChEBI/PubChem and source annotations")
    c=status.sort_values(); ax[1].barh(c.index,c.values,color=sns.color_palette("Paired",len(c),desat=.72),alpha=.84); label_barh(ax[1],c.values); ax[1].set_xlabel("Annotated compounds (non-exclusive)"); clean(ax[1],"x")
    panel(ax[2],"C","One deterministic structural classification","RDKit-compatible descriptor rules; mutually exclusive")
    c=classes.sort_values(); ax[2].barh(c.index,c.values,color=COL["teal"],alpha=.84); label_barh(ax[2],c.values); ax[2].set_xlabel("Core QC canonical compounds"); clean(ax[2],"x")
    classes.rename("compound_count").to_csv(DATA/"M5C_recomputed_structural_classes.tsv",sep="\t")
    panel(ax[3],"D","Descriptor distributions guide docking, not database inclusion","MemPro v1.3 calculated properties")
    q=desc.quantile([.1,.5,.9]).T; q.columns=["p10","median","p90"]
    qn=(q-q.min(axis=1).to_numpy()[:,None])/(q.max(axis=1)-q.min(axis=1)).replace(0,1).to_numpy()[:,None]
    sns.heatmap(qn,cmap="flare",annot=q.round(1),fmt=".1f",linewidths=.4,cbar=False,ax=ax[3]); ax[3].set_xlabel("Quantile"); ax[3].set_ylabel("")
    panel(ax[4],"E","PAGTN clustering is an evaluated method experiment","PAGTN graph encoder; fixed sample; Morgan baseline")
    pagtn=DATA/"M5E_pagtn_cluster_summary.tsv"
    if pagtn.exists():
        d=pd.read_csv(pagtn,sep="\t"); ax[4].bar(d.metric,d.value,color=[COL["blue"],COL["teal"],COL["orange"]][:len(d)],alpha=.84); ax[4].tick_params(axis="x",rotation=25); clean(ax[4],"y")
    else:
        ax[4].axis("off"); ax[4].text(.5,.56,"PAGTN environment/training in progress",ha="center",va="center",fontsize=11,color=COL["purple"]); ax[4].text(.5,.39,"No random-weight embedding is reported",ha="center",va="center",fontsize=8,color="#666D72")
    panel(ax[5],"F","Scaffold diversity complements learned embeddings","Bemis–Murcko scaffold mapping, V6.3 candidate")
    sf=pd.read_csv(V63/"04_disease/scaffold_frequency_v0_1.tsv",sep="\t").head(1000); ax[5].plot(sf.scaffold_rank,sf.compound_count,color=COL["blue"],lw=1.5); ax[5].set_xscale("log"); ax[5].set_yscale("log"); ax[5].set_xlabel("Scaffold rank (log)"); ax[5].set_ylabel("Compounds per scaffold (log)"); clean(ax[5],"both")
    return save(fig,"M5_chemical_identity_space_revised")


def make_m6():
    fig,ax=fig6("M6  Binding evidence, site availability and membrane-side context")
    mat=pd.read_csv(OLD/"data/M6_source_evidence_matrix.tsv",sep="\t").set_index("source")
    panel(ax[0],"A","Sources contribute different evidence tiers","ChEMBL/BindingDB/PubChem/BRENDA/PDBe/PDBbind")
    top=mat.assign(total=mat.sum(axis=1)).nlargest(10,"total").drop(columns="total").sort_values(mat.columns.tolist())
    left=np.zeros(len(top))
    for t in ["BE1","BE2","BE3"]:
        ax[0].barh(top.index,top[t],left=left,label=t,color=COL[t],alpha=.86); left+=top[t].values
    ax[0].set_xscale("log"); ax[0].set_xlabel("Evidence records (log; stacked)"); ax[0].legend(frameon=False,ncol=3); clean(ax[0],"x")
    panel(ax[1],"B","BE1/BE2/BE3 encode evidence type, not potency","V6.2 binding evidence policy")
    totals=mat.sum().reindex(["BE1","BE2","BE3"]); ax[1].bar(totals.index,totals.values,color=[COL[x] for x in totals.index],alpha=.86); ax[1].set_yscale("log"); ax[1].set_ylabel("Evidence records (log)"); clean(ax[1],"y")
    site=pd.read_csv(V62/"binding_site_instances_v6_2.tsv.gz",sep="\t",compression="gzip",low_memory=False)
    panel(ax[2],"C","Binding-site readiness has several levels","PDBe/PDBbind/BioLiP/BindingDB/UniProtKB")
    labels={"experimental_structure_residue_contact":"PDB + residue contacts","experimental_structure_match":"PDB structure match","experimental_complex_pocket":"Experimental pocket","bindingdb_ligand_target_complex":"Complex PDB only","uniprot_curated_binding_site":"UniProt curated site"}
    c=site.site_type.map(labels).fillna("Other").value_counts().sort_values(); ax[2].barh(c.index,c.values,color=COL["teal"],alpha=.84); label_barh(ax[2],c.values); ax[2].set_xlabel("Binding-site instances"); clean(ax[2],"x")
    side_path=ROOT/"06_binding_site_side/binding_site_membrane_side_summary.tsv"
    panel(ax[3],"D","Sites are separated into membrane-relative compartments","PDBe/SIFTS residues + UniProt/UniTmp topology + OPM/PDBTM eligibility")
    if side_path.exists():
        d=pd.read_csv(side_path,sep="\t"); ax[3].barh(d.membrane_side,d.site_count,color=[COL["teal"],COL["blue"],COL["orange"],COL["purple"],COL["gray"]][:len(d)],alpha=.84); label_barh(ax[3],d.site_count.values); ax[3].set_xlabel("Residue-contact site instances"); clean(ax[3],"x")
    else:
        ax[3].axis("off"); ax[3].text(.5,.5,"Membrane-side classifier running",ha="center",fontsize=11,color=COL["purple"])
    panel(ax[4],"E","Contact residues summarize pocket chemistry","PDBe/SIFTS UniProt-coordinate contacts")
    aa=Counter();
    for text in site.residue_or_site_description.fillna(""):
        for a in re.findall(r"\b([A-Z]{3})\d+",text): aa[a]+=1
    groups={"Hydrophobic":set("ALA VAL ILE LEU MET PRO GLY".split()),"Polar":set("SER THR ASN GLN CYS".split()),"Aromatic":set("PHE TYR TRP HIS".split()),"Positive":set("LYS ARG".split()),"Negative":set("ASP GLU".split())}
    gc=pd.Series({g:sum(aa[x] for x in aset) for g,aset in groups.items()}).sort_values(); ax[4].barh(gc.index,gc.values,color=COL["blue"],alpha=.84); label_barh(ax[4],gc.values); ax[4].set_xlabel("Residue-contact occurrences"); clean(ax[4],"x")
    panel(ax[5],"F","Docking-box readiness requires chain, ligand and site context","V6.3.1 docking readiness")
    vals=pd.Series({"PDB + residue contacts":int((site.site_type=="experimental_structure_residue_contact").sum()),"PDB structure match":int((site.site_type=="experimental_structure_match").sum()),"Curated site only":int((site.site_type=="uniprot_curated_binding_site").sum())}).sort_values(); ax[5].barh(vals.index,vals.values,color=[COL["teal"],COL["blue"],COL["orange"]],alpha=.84); label_barh(ax[5],vals.values); ax[5].set_xlabel("Site instances"); clean(ax[5],"x")
    return save(fig,"M6_binding_evidence_sites_revised")


def make_m7():
    fig,ax=fig6("M7  Negative evidence, conflicts, review queues and release QA")
    panels=[
        ("Negative identity disposition",pd.Series({"Mapped release":7_630_659,"Unmapped review":2_615_046,"Duplicates":1_243}),"PubChem BioAssay negative evidence"),
        ("Positive–negative conflicts",pd.Series({"Conflict pairs":68_015,"Non-conflict negative pairs":6_331_307-68_015}),"V6.2 positive/negative pair comparison"),
        ("Review queues",pd.Series({"Positive identity":1_626_364,"Negative identity":2_615_046,"Subunit state":3_348,"Expression/location":496}),"V6.1/V6.2 review layers"),
        ("Key missingness",pd.read_csv(OLD/"data/M7_missingness.tsv",sep="\t").set_index("field").iloc[:,0].nlargest(8)*100,"V6.2 release fields"),
        ("QA gates",pd.Series({"Primary key":1,"Foreign key":1,"Count reconciliation":1,"Hash/manifest":1,"Blocking errors":0}),"V6.2 validation report"),
        ("Release layers",pd.Series({"Default high confidence":1_502_456,"Positive review":1_626_364,"Negative mapped":6_331_307,"Negative review":2_615_046}),"V6.2 frozen release"),
    ]
    for i,(title,s,src) in enumerate(panels):
        panel(ax[i],chr(65+i),title,src); s=s.sort_values(); ax[i].barh(s.index,s.values,color=[COL["red"] if "Conflict" in x or "Blocking" in x else COL["blue"] for x in s.index],alpha=.82); label_barh(ax[i],s.values,lambda x:f"{x:.1f}" if max(s.values)<=100 else f"{int(x):,}"); clean(ax[i],"x")
    return save(fig,"M7_negative_conflicts_quality_revised")


def make_m8():
    fig,ax=fig6("M8  Docking prioritization and HPC execution readiness")
    dock_path=Path(r"D:\finale\03_V6.3.1_readiness_20260803\04_docking_hpc\docking_pilot_manifest_v631.tsv")
    if not dock_path.exists():
        dock_path=V63/"07_docking/docking_priority_v6_3_candidate.tsv.gz"
    d=pd.read_csv(dock_path,sep="\t",compression="gzip" if str(dock_path).endswith(".gz") else None,low_memory=False)
    tier_col="docking_tier" if "docking_tier" in d else [x for x in d.columns if "tier" in x][0]
    score_col="ranking_score_v63" if "ranking_score_v63" in d else [x for x in d.columns if "score" in x][0]
    panels=[]
    panels.append(("Docking tiers",d[tier_col].value_counts()))
    bool_groups={"Candidate receptor PDB":d.get("candidate_receptor_pdb_id_v2",pd.Series(False,index=d.index)).fillna("").astype(str).ne(""),"Pair-specific site":d.get("pair_site_pdb_ids",pd.Series(False,index=d.index)).fillna("").astype(str).ne(""),"Ligand SMILES":d.get("standard_smiles",pd.Series(False,index=d.index)).fillna("").astype(str).ne(""),"Ligand InChIKey":d.get("standard_inchikey",pd.Series(False,index=d.index)).fillna("").astype(str).ne("")}
    panels.append(("Structure/identity readiness",pd.Series({k:int(v.sum()) for k,v in bool_groups.items()})))
    panels.append(("Unique targets and compounds",pd.Series({"Proteins":d.get("target_uniprot_id",pd.Series()).nunique(),"Compounds":d.get("compound_internal_id",pd.Series()).nunique(),"Pairs":len(d)})))
    site_grade=d.get("site_grade",d.get("binding_site_grade",pd.Series("unresolved",index=d.index))).fillna("unresolved").value_counts(); panels.append(("Site readiness grades",site_grade.head(8)))
    score=pd.to_numeric(d[score_col],errors="coerce").dropna(); panels.append(("Ranking score quantiles",score.quantile([.1,.25,.5,.75,.9]).rename(index=lambda x:f"q{x:g}")))
    workload=pd.Series({"Pilot pairs":len(d),"Estimated receptor jobs":d.get("candidate_receptor_pdb_id_v2",pd.Series()).nunique(),"Conflict flagged":truth(d.get("positive_negative_conflict_flag_v62",pd.Series(False,index=d.index))).sum()}); panels.append(("HPC workload",workload))
    for i,(title,s) in enumerate(panels):
        panel(ax[i],chr(65+i),title,"V6.3.1 docking readiness / V6.3 candidate")
        s=s.sort_values(); ax[i].barh(s.index.astype(str),s.values,color=sns.color_palette("flare",max(len(s),3),desat=.72)[:len(s)],alpha=.84); label_barh(ax[i],s.values,lambda x:f"{x:.2f}" if max(s.values)<20 else f"{int(x):,}"); clean(ax[i],"x")
    return save(fig,"M8_docking_prioritization_revised")


def main():
    protein=load_protein(); audit=source_audit()
    outputs={}
    for key,fn in [("M1",lambda:make_m1(protein)),("M2",lambda:make_m2(protein)),("M3",lambda:make_m3(protein)),("M4",make_m4),("M5",make_m5),("M6",make_m6),("M7",make_m7),("M8",make_m8)]:
        try:
            outputs[key]={"status":"PASS","files":fn()}
            print(key,"PASS",flush=True)
        except Exception as exc:
            outputs[key]={"status":"FAIL","error":repr(exc)}
            print(key,"FAIL",repr(exc),flush=True)
    qa={"status":"PASS" if all(x["status"]=="PASS" for x in outputs.values()) else "PARTIAL",
        "release_data_modified":False,"panel_source_rows":len(audit),"outputs":outputs}
    (QA/"FIGURE_REVISION_VALIDATION.json").write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(qa,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
