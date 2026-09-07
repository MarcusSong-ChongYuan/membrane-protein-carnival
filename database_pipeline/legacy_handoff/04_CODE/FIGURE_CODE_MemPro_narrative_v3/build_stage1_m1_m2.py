from __future__ import annotations

import json
import math
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

import hdbscan
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import umap
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

ROOT = Path(r"C:\Users\Administrator\MemPro_narrative_v3")
DATA = ROOT / "data"
FIG = ROOT / "figures"
QA = ROOT / "qa"
REL = Path(r"D:\finale\20_MemPro_V7.1.1_final_20260814\01_release_tables")
for p in (DATA, FIG, QA): p.mkdir(parents=True, exist_ok=True)

C = {"blue":"#7b95c6","cyan":"#49c2d9","lcyan":"#a1d8e8","green":"#67a583",
     "lgreen":"#a2c986","pgreen":"#d0e2c0","yellow":"#fded95","peach":"#ffc1a6",
     "salmon":"#f59c7c","coral":"#f47254","red":"#c85e62","ink":"#27313d",
     "muted":"#66717e","grid":"#e5e9ee","missing":"#b9bec6"}
PAL=[C["blue"],C["cyan"],C["green"],C["lgreen"],C["yellow"],C["peach"],C["salmon"],C["red"]]
mpl.rcParams.update({"font.family":"Arial","font.size":8.3,"axes.titlesize":9.4,"axes.labelsize":8.2,
    "xtick.labelsize":7.2,"ytick.labelsize":7.2,"text.color":C["ink"],"axes.labelcolor":C["ink"],
    "axes.titlecolor":C["ink"],"axes.edgecolor":"#c7cdd4","grid.color":C["grid"],
    "pdf.fonttype":42,"svg.fonttype":"none","legend.frameon":False})
sns.set_style("whitegrid")

def fmt(v):
    v=float(v);return f"{v/1e6:.2f}M" if v>=1e6 else f"{v/1e3:.1f}k" if v>=1e3 else f"{int(v):,}"

def panel(ax,l,t): ax.set_title(f"{l}  {t}",loc="left",weight="bold",pad=8)
def clean(ax,grid="x"):
    ax.spines[["top","right"]].set_visible(False);ax.grid(False)
    if grid:ax.grid(axis=grid,alpha=.7)
def save(fig,stem,source):
    fig.text(.035,.018,"Source: "+source,fontsize=7,color=C["muted"])
    for ext in ("png","svg","pdf"):fig.savefig(FIG/f"{stem}.{ext}",dpi=450 if ext=="png" else None,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def source_sets():
    p=REL/"protein_compound_pair_v71.tsv.gz"
    d=pd.read_csv(p,sep="\t",usecols=["pair_id","source_databases"],low_memory=False)
    sets={}
    norm=[]
    for pid,s in zip(d.pair_id,d.source_databases.fillna("")):
        ss=tuple(sorted({x.strip() for x in str(s).replace("|",";").split(";") if x.strip()}))
        norm.append(ss)
        for x in ss:sets.setdefault(x,set()).add(pid)
    top=[x for x,_ in sorted(sets.items(),key=lambda z:len(z[1]),reverse=True)[:10]]
    j=np.eye(len(top))
    for a in range(len(top)):
        for b in range(a+1,len(top)):
            u=len(sets[top[a]]|sets[top[b]]);j[a,b]=j[b,a]=len(sets[top[a]]&sets[top[b]])/u if u else 0
    remaining=set(d.pair_id);covered=set();greedy=[]
    candidates=set(top)
    while candidates:
        best=max(candidates,key=lambda x:len(sets[x]-covered));new=len(sets[best]-covered);covered|=sets[best]
        greedy.append((best,new,len(covered)));candidates.remove(best)
    inter=Counter(tuple(x for x in top if x in ss) for ss in norm)
    inter={" & ".join(k) if k else "none":v for k,v in inter.items() if k}
    unique={x:len(sets[x]-set().union(*(sets[y] for y in sets if y!=x))) for x in top}
    pd.DataFrame(j,index=top,columns=top).to_csv(DATA/"M1_source_jaccard.tsv",sep="\t")
    pd.DataFrame(greedy,columns=["source","new_pairs","cumulative_pairs"]).to_csv(DATA/"M1_greedy_coverage.tsv",sep="\t",index=False)
    pd.DataFrame(sorted(inter.items(),key=lambda z:z[1],reverse=True),columns=["intersection","pair_count"]).to_csv(DATA/"M1_source_intersections.tsv",sep="\t",index=False)
    pd.DataFrame({"source":top,"total_pairs":[len(sets[x]) for x in top],"strictly_unique_pairs":[unique[x] for x in top]}).to_csv(DATA/"M1_source_unique_contribution.tsv",sep="\t",index=False)
    return top,j,greedy,inter,unique,len(d)

def draw_m1():
    top,j,greedy,inter,unique,n=source_sets()
    fig=plt.figure(figsize=(13.6,9.2),facecolor="white")
    fig.text(.035,.968,"M1 | Integration adds unique coverage while exposing source redundancy",fontsize=15,weight="bold",va="top")
    fig.text(.035,.932,"Pair-level overlap distinguishes raw database volume from genuinely new membrane-protein interaction coverage.",fontsize=9,color=C["muted"],va="top")
    gs=fig.add_gridspec(2,2,left=.06,right=.98,top=.87,bottom=.09,hspace=.43,wspace=.33,height_ratios=[1.08,1])
    ax=fig.add_subplot(gs[0,0]);panel(ax,"A","Source intersections reveal shared and source-specific interaction space")
    arr=sorted(inter.items(),key=lambda z:z[1],reverse=True)[:8];x=np.arange(len(arr));ax.bar(x,[v for _,v in arr],color=C["blue"]);ax.set_yscale("log");ax.set_ylabel("unique pairs (log scale)");ax.set_xticks(x,[]);clean(ax,"y")
    ia=ax.inset_axes([0,-.42,1,.34]);ia.set_xlim(-.5,len(arr)-.5);ia.set_ylim(-.5,len(top)-.5);ia.axis("off")
    for i,s in enumerate(top):ia.text(-.65,len(top)-1-i,textwrap.shorten(s,18,placeholder="..."),ha="right",va="center",fontsize=6.8)
    for k,(name,_) in enumerate(arr):
        active=set(name.split(" & "));ys=[]
        for i,s in enumerate(top):
            y=len(top)-1-i;on=s in active;ia.scatter(k,y,s=22,color=C["ink"] if on else "#d9dde2");ys.append(y) if on else None
        if len(ys)>1:ia.plot([k,k],[min(ys),max(ys)],color=C["ink"],lw=.9)
    ax=fig.add_subplot(gs[0,1]);panel(ax,"B","Database similarity clusters sources with overlapping pair coverage")
    dist=np.clip(1-j,0,1);order=leaves_list(linkage(squareform(dist,checks=False),method="average"));jj=j[np.ix_(order,order)];labs=[top[i] for i in order]
    sns.heatmap(jj,ax=ax,cmap="YlGnBu",vmin=0,vmax=1,xticklabels=labs,yticklabels=labs,cbar_kws={"label":"pair-set Jaccard"},square=True)
    ax.tick_params(axis="x",rotation=45);ax.tick_params(axis="y",rotation=0)
    ax=fig.add_subplot(gs[1,0]);panel(ax,"C","Greedy source addition quantifies marginal coverage gain")
    g=pd.DataFrame(greedy,columns=["source","new","cum"]);xx=np.arange(len(g));ax.bar(xx,g.new,color=PAL[:len(g)]);ax.plot(xx,g.cum,color=C["ink"],marker="o",lw=1.4,label="cumulative pairs")
    ax.set_xticks(xx,[textwrap.shorten(x,14,placeholder="...") for x in g.source],rotation=35,ha="right");ax.set_ylabel("pairs");ax.legend(fontsize=7);clean(ax,"y")
    ax=fig.add_subplot(gs[1,1]);panel(ax,"D","Strictly unique pairs identify sources that cannot be replaced by another database")
    uu=sorted(unique.items(),key=lambda z:z[1]);ax.hlines(np.arange(len(uu)),0,[v for _,v in uu],color=C["grid"],lw=3);ax.scatter([v for _,v in uu],np.arange(len(uu)),s=55,color=C["coral"])
    ax.set_yticks(np.arange(len(uu)),[textwrap.shorten(x,24,placeholder="...") for x,_ in uu]);ax.set_xlabel("strictly unique protein-compound pairs")
    for y,(_,v) in enumerate(uu):ax.text(v,y,"  "+fmt(v),va="center",fontsize=7);clean(ax)
    save(fig,"M1_source_overlap_and_marginal_gain","V7.1.1 protein-compound pair table; n="+f"{n:,}"+" released pairs. Source labels represent contributing databases, not independent experiments.")

def protein_space():
    pm=pd.read_csv(REL/"protein_master_v71.tsv.gz",sep="\t",usecols=["canonical_uniprot_accession","membrane_class_v7","membrane_evidence_level_v7"],low_memory=False).rename(columns={"canonical_uniprot_accession":"target_uniprot_id"})
    cc=pd.read_csv(REL/"protein_cross_classification_v71.tsv.gz",sep="\t",usecols=["target_uniprot_id","classification_axis","classification_label"],low_memory=False)
    cc=cc[~cc.classification_label.astype(str).str.contains("unclassified|no_specialist",case=False,regex=True)]
    freq=cc.classification_label.value_counts();keep=set(freq.head(180).index);cc=cc[cc.classification_label.isin(keep)]
    mat=pd.crosstab(cc.target_uniprot_id,cc.classification_label).clip(upper=1)
    mat=mat.reindex(pm.target_uniprot_id,fill_value=0)
    reducer=umap.UMAP(n_neighbors=25,min_dist=.18,n_components=2,metric="jaccard",random_state=42,low_memory=True)
    emb=reducer.fit_transform(mat.values.astype(np.uint8))
    cluster=hdbscan.HDBSCAN(min_cluster_size=90,min_samples=20,cluster_selection_method="eom").fit_predict(emb)
    out=pm.copy();out["umap1"]=emb[:,0];out["umap2"]=emb[:,1];out["cluster"]=cluster
    out.to_csv(DATA/"M2_protein_annotation_umap.tsv",sep="\t",index=False)
    enrich=[]
    for cl in sorted(set(cluster)-{-1}):
        ids=set(out.loc[out.cluster.eq(cl),"target_uniprot_id"]);sub=cc[cc.target_uniprot_id.isin(ids)];cnt=sub.classification_label.value_counts().head(5)
        for lab,v in cnt.items():enrich.append({"cluster":cl,"label":lab,"count":v,"cluster_size":len(ids)})
    pd.DataFrame(enrich).to_csv(DATA/"M2_cluster_top_annotations.tsv",sep="\t",index=False)
    return out,pd.DataFrame(enrich),mat.shape

def draw_m2():
    d,en,shape=protein_space()
    fig=plt.figure(figsize=(13.6,8.8),facecolor="white")
    fig.text(.035,.968,"M2 | Multi-axis annotations organize membrane proteins into reproducible functional neighborhoods",fontsize=15,weight="bold",va="top")
    fig.text(.035,.932,"The embedding uses structural family, molecular function, biological process, membrane role and specialist classifications together.",fontsize=9,color=C["muted"],va="top")
    gs=fig.add_gridspec(2,3,left=.055,right=.98,top=.87,bottom=.10,hspace=.34,wspace=.32,width_ratios=[1.1,1.1,.9])
    ax=fig.add_subplot(gs[:,0:2]);panel(ax,"A","Annotation neighborhoods separate major membrane mechanisms without forcing one label per protein")
    cmap={"A":C["blue"],"B":C["green"],"C":C["peach"]}
    for cl,sub in d.groupby("membrane_class_v7"):
        ax.scatter(sub.umap1,sub.umap2,s=8,c=cmap.get(cl,C["missing"]),alpha=.32,label=f"Class {cl}",rasterized=True)
    ax.set_xlabel("UMAP 1");ax.set_ylabel("UMAP 2");ax.legend(markerscale=2,fontsize=8);clean(ax,None)
    ax=fig.add_subplot(gs[0,2]);panel(ax,"B","HDBSCAN finds dense annotation communities and leaves diffuse proteins unforced")
    counts=d.cluster.value_counts().sort_index();colors=[C["missing"] if x==-1 else PAL[x%len(PAL)] for x in counts.index]
    ax.bar(np.arange(len(counts)),counts.values,color=colors);ax.set_xticks(np.arange(len(counts)),["noise" if x==-1 else str(x) for x in counts.index]);ax.set_xlabel("annotation community");ax.set_ylabel("proteins");clean(ax,"y")
    ax=fig.add_subplot(gs[1,2]);panel(ax,"C","Each community is summarized by enriched biological labels")
    if not en.empty:
        piv=en.pivot_table(index="cluster",columns="label",values="count",fill_value=0);top=en.groupby("label")["count"].sum().nlargest(8).index;piv=piv.reindex(columns=top,fill_value=0)
        row=piv.div(piv.sum(axis=1).replace(0,np.nan),axis=0);sns.heatmap(row,ax=ax,cmap="YlGnBu",cbar_kws={"label":"within-cluster share"},linewidths=.4,linecolor="white")
        ax.set_xlabel("");ax.set_ylabel("community")
    fig.text(.56,.105,f"Embedding input: {shape[0]:,} proteins x {shape[1]:,} recurrent annotation labels; seed=42.",fontsize=7,color=C["muted"])
    save(fig,"M2_multi_axis_annotation_space","V7.1.1 protein master and five-axis cross-classification table; UMAP metric=Jaccard, HDBSCAN clustering.")

if __name__=="__main__":
    draw_m1();draw_m2()
    (QA/"STAGE1_COMPLETE.json").write_text(json.dumps({"M1":"complete","M2":"complete","seed":42},indent=2),encoding="utf-8")
    print("M1/M2 complete")
