from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
import build_revision_figures as b

ROOT=Path(os.environ.get("MEMPRO_REVISION_ROOT",r"D:\finale\07_M1-M8图表修订_20260806"))

def main():
    p=b.load_protein(); fig,ax=b.fig6("M2  Membrane proteome composition, topology and structure coverage")
    b.panel(ax[0],"A","Membrane class × evidence grade","UniProtKB/HPA/HTP/Membranome/OPM/PDBTM")
    cross=pd.crosstab(p["class"],p["grade"]).reindex(index=["A","B","C","unknown"],columns=["E1","E2","E3","E0"],fill_value=0); left=np.zeros(len(cross))
    for g in cross.columns: ax[0].barh(cross.index,cross[g],left=left,color=b.COL[g],alpha=.88,label=g); left+=cross[g].values
    ax[0].set_xlabel("Unique canonical proteins"); ax[0].legend(frameon=False,ncol=4,loc="lower center",bbox_to_anchor=(.5,-.31)); b.clean(ax[0],"x")
    b.panel(ax[1],"B","Functional hierarchy","UniProtKB and V6.3 cross-classification")
    c=p.function.value_counts().head(11).sort_values(); ax[1].barh(c.index.str.replace("_"," "),c.values,color=b.COL["teal"],alpha=.86); b.label_barh(ax[1],c.values); ax[1].set_xlabel("Unique proteins"); b.clean(ax[1],"x")
    b.panel(ax[2],"C","Transmembrane architecture","UniProtKB/UniTmp topology")
    bins=pd.cut(p.tm,[-.1,.5,1.5,4.5,8.5,12.5,20.5,np.inf],labels=["0","1","2–4","5–8","9–12","13–20",">20"]); c=bins.value_counts().sort_index(); ax[2].bar(c.index.astype(str),c.values,color=sns.color_palette("flare",len(c),desat=.72),alpha=.88); ax[2].set_ylabel("Unique proteins"); ax[2].set_xlabel("TM helices"); b.clean(ax[2],"y")
    b.panel(ax[3],"D","Length–TM association after biological stratification","A-class E1/E2 canonical multipass; 50–5,000 aa; 2–40 TM")
    d=p[p["class"].eq("A")&p.grade.isin(["E1","E2"])&p.tm.between(2,40)&p.length.between(50,5000)].dropna(subset=["tm","length"]).copy()
    rng=np.random.default_rng(42); x=d.tm.to_numpy()+rng.uniform(-.18,.18,len(d))
    ax[3].scatter(x,d.length,s=8,color=b.COL["teal"],alpha=.26,linewidths=0,rasterized=True)
    med=d.groupby("tm").length.median(); ax[3].plot(med.index,med.values,color=b.COL["navy"],lw=1.5,marker="o",ms=3,label="median by TM count")
    ax[3].set_yscale("log"); ax[3].set_xlim(1.4,40.6); ax[3].set_xlabel("TM helices (jittered)"); ax[3].set_ylabel("Sequence length (aa, log)"); ax[3].legend(frameon=False,loc="lower right")
    rho,pv=spearmanr(d.tm,d.length); ax[3].text(.97,.96,f"Spearman ρ={rho:.3f}\nn={len(d):,}\nbootstrap CI reported in data table",transform=ax[3].transAxes,ha="right",va="top",fontsize=7,bbox={"boxstyle":"round,pad=.3","fc":"white","ec":"#D6D9DC"}); b.clean(ax[3],"y")
    boot=[]
    for _ in range(2000):
        idx=rng.integers(0,len(d),len(d)); boot.append(spearmanr(d.tm.to_numpy()[idx],d.length.to_numpy()[idx]).statistic)
    pd.DataFrame([{"scope":"A; E1/E2; canonical; TM 2–40; length 50–5000 aa","n":len(d),"spearman_rho":rho,"p_value":pv,"bootstrap_ci_low":np.quantile(boot,.025),"bootstrap_ci_high":np.quantile(boot,.975),"seed":42}]).to_csv(ROOT/"02_data/M2D_spearman_main_and_bootstrap.tsv",sep="\t",index=False)
    b.panel(ax[4],"E","Structure coverage retains AlphaFold","PDB/OPM/PDBTM and AlphaFoldDB v6")
    top=p.function.value_counts().head(9).index; q=p[p.function.isin(top)]; st=pd.crosstab(q.function,q.structure,normalize="index").reindex(index=top,columns=["Membrane PDB","Other PDB","AlphaFold only","No structure"],fill_value=0)*100
    sns.heatmap(st,cmap=sns.light_palette(b.COL["blue"],as_cmap=True),annot=True,fmt=".0f",linewidths=.35,cbar_kws={"label":"% of class"},ax=ax[4]); ax[4].set_xlabel(""); ax[4].set_ylabel(""); ax[4].set_xticklabels(["Membrane\nPDB","Other\nPDB","AlphaFold\nonly","No\nstructure"],rotation=0); ax[4].set_yticklabels([x.replace("_"," ") for x in st.index],rotation=0)
    b.panel(ax[5],"F","Biological assembly evidence","UniProtKB and PDBe biological assemblies")
    c=p.oligomeric_state_consensus_v62.fillna("unresolved").value_counts().head(10).sort_values(); ax[5].barh(c.index.str.replace("_"," "),c.values,color=b.COL["purple"],alpha=.82); b.label_barh(ax[5],c.values); ax[5].set_xlabel("Unique proteins"); b.clean(ax[5],"x")
    print(b.save(fig,"M2_membrane_proteome_landscape_revised"))
if __name__=="__main__": main()
