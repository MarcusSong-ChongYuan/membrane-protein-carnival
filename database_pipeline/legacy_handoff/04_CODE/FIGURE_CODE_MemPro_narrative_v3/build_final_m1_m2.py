from pathlib import Path
import textwrap
import matplotlib as mpl,matplotlib.pyplot as plt,numpy as np,pandas as pd,seaborn as sns
from scipy.cluster.hierarchy import leaves_list,linkage
from scipy.spatial.distance import squareform

ROOT=Path(r"C:\Users\Administrator\MemPro_narrative_v3");DATA=ROOT/"data";FIG=ROOT/"figures";REL=Path(r"D:\finale\20_MemPro_V7.1.1_final_20260814\01_release_tables")
C={"blue":"#7b95c6","cyan":"#49c2d9","lcyan":"#a1d8e8","green":"#67a583","lgreen":"#a2c986","yellow":"#fded95","peach":"#ffc1a6","salmon":"#f59c7c","coral":"#f47254","red":"#c85e62","ink":"#27313d","muted":"#66717e","grid":"#e5e9ee","missing":"#b9bec6"};PAL=[C["blue"],C["cyan"],C["green"],C["lgreen"],C["yellow"],C["peach"],C["salmon"],C["red"]]
mpl.rcParams.update({"font.family":"Arial","font.size":8.2,"axes.titlesize":9.3,"axes.labelsize":8.2,"xtick.labelsize":7.2,"ytick.labelsize":7.2,"text.color":C["ink"],"axes.labelcolor":C["ink"],"axes.titlecolor":C["ink"],"axes.edgecolor":"#c7cdd4","grid.color":C["grid"],"pdf.fonttype":42,"svg.fonttype":"none","legend.frameon":False});sns.set_style("whitegrid")
def panel(ax,l,t):ax.set_title(f"{l}  {t}",loc="left",weight="bold",pad=8)
def clean(ax,grid="x"):
 ax.spines[["top","right"]].set_visible(False);ax.grid(False)
 if grid:ax.grid(axis=grid,alpha=.7)
def save(fig,stem,src):
 fig.text(.035,.018,"Source: "+src,fontsize=7,color=C["muted"])
 for ext in ("png","svg","pdf"):fig.savefig(FIG/f"{stem}.{ext}",dpi=450 if ext=="png" else None,bbox_inches="tight",facecolor="white")
 plt.close(fig)
def fmt(v):
 v=float(v);return f"{v/1e6:.2f}M" if v>=1e6 else f"{v/1e3:.1f}k" if v>=1e3 else f"{int(v):,}"

def m1():
 inter=pd.read_csv(DATA/"M1_source_intersections.tsv",sep="\t");jac=pd.read_csv(DATA/"M1_source_jaccard.tsv",sep="\t",index_col=0);greedy=pd.read_csv(DATA/"M1_greedy_coverage.tsv",sep="\t");uniq=pd.read_csv(DATA/"M1_source_unique_contribution.tsv",sep="\t")
 sources=list(jac.index);codes={s:chr(65+i) for i,s in enumerate(sources)};arr=inter.head(9).copy();arr["active"]=arr.intersection.map(lambda x:[s for s in sources if s in str(x).split(" & ")])
 fig=plt.figure(figsize=(14,10));fig.text(.03,.972,"M1 | Integration adds unique coverage while exposing source redundancy",fontsize=15,weight="bold",va="top");fig.text(.03,.938,"Pair-level intersections distinguish contributing database volume, shared records and genuinely new coverage.",fontsize=9,color=C["muted"],va="top")
 gs=fig.add_gridspec(2,2,left=.065,right=.985,top=.87,bottom=.09,hspace=.40,wspace=.34)
 left=gs[0,0].subgridspec(2,1,height_ratios=[1.15,.85],hspace=.08);ax=fig.add_subplot(left[0]);panel(ax,"A","A dedicated UpSet view resolves the largest source intersections");x=np.arange(len(arr));ax.bar(x,arr.pair_count,color=C["blue"]);ax.set_yscale("log");ax.set_ylabel("unique pairs (log)");ax.set_xticks(x,[]);clean(ax,"y")
 am=fig.add_subplot(left[1]);am.set_xlim(-.5,len(arr)-.5);am.set_ylim(-.5,len(sources)-.5);am.axis("off")
 for i,s in enumerate(sources):am.text(-.65,len(sources)-1-i,f"{codes[s]}  {textwrap.shorten(s,18,placeholder='...')}",ha="right",va="center",fontsize=6.8)
 for j,r in arr.iterrows():
  ys=[]
  for i,s in enumerate(sources):
   y=len(sources)-1-i;on=s in r.active;am.scatter(j,y,s=25,color=C["ink"] if on else "#d9dde2");ys.append(y) if on else None
  if len(ys)>1:am.plot([j,j],[min(ys),max(ys)],color=C["ink"],lw=.9)
 ax=fig.add_subplot(gs[0,1]);panel(ax,"B","Jaccard similarity groups sources with overlapping pair coverage");dist=np.clip(1-jac.values,0,1);o=leaves_list(linkage(squareform(dist,checks=False),method="average"));sns.heatmap(jac.values[np.ix_(o,o)],ax=ax,cmap="YlGnBu",vmin=0,vmax=1,square=True,xticklabels=[codes[sources[i]] for i in o],yticklabels=[codes[sources[i]] for i in o],cbar_kws={"label":"pair-set Jaccard"});ax.set_xlabel("source code");ax.set_ylabel("source code")
 ax=fig.add_subplot(gs[1,0]);panel(ax,"C","Greedy source addition quantifies marginal coverage gain");g=greedy.head(len(sources));x=np.arange(len(g));ax.bar(x,g.new_pairs,color=PAL[:len(g)]);ax.plot(x,g.cumulative_pairs,color=C["ink"],marker="o",lw=1.2,label="cumulative");ax.set_xticks(x,[codes.get(s,"?") for s in g.source]);ax.set_ylabel("unique pairs");ax.set_xlabel("source code in greedy order");ax.legend(fontsize=7);clean(ax,"y")
 ax=fig.add_subplot(gs[1,1]);panel(ax,"D","Strictly unique pairs identify irreplaceable source contributions");u=uniq.sort_values("strictly_unique_pairs");y=np.arange(len(u));ax.hlines(y,0,u.strictly_unique_pairs,color=C["grid"],lw=3);ax.scatter(u.strictly_unique_pairs,y,s=55,color=C["coral"]);ax.set_yticks(y,[codes.get(s,"?") for s in u.source]);ax.set_xlabel("strictly unique released pairs");
 for yy,v in zip(y,u.strictly_unique_pairs):ax.text(v,yy,"  "+fmt(v),va="center",fontsize=7);clean(ax)
 save(fig,"M1_source_overlap_and_marginal_gain","V7.1.1 released protein-compound pairs; source codes are defined in panel A; contributing databases are not independent experiments.")

def m2():
 d=pd.read_csv(DATA/"M2_protein_annotation_umap.tsv",sep="\t");cc=pd.read_csv(REL/"protein_cross_classification_v71.tsv.gz",sep="\t",usecols=["target_uniprot_id","classification_axis"]);x=d[["target_uniprot_id","cluster"]].merge(cc,on="target_uniprot_id",how="left");tab=pd.crosstab(x.cluster,x.classification_axis);tab=tab.drop(index=-1,errors="ignore");largest=d[d.cluster.ne(-1)].cluster.value_counts().head(12).index;tab=tab.reindex(largest).fillna(0);share=tab.div(tab.sum(1),axis=0)
 fig=plt.figure(figsize=(14,9.4));fig.text(.03,.968,"M2 | Multi-axis annotations organize membrane proteins into reproducible functional neighborhoods",fontsize=15,weight="bold",va="top");fig.text(.03,.932,"The embedding jointly uses structural family, molecular function, biological process, membrane role and specialist classifications.",fontsize=9,color=C["muted"],va="top")
 gs=fig.add_gridspec(2,3,left=.055,right=.985,top=.87,bottom=.09,hspace=.34,wspace=.32,width_ratios=[1.1,1.1,.9])
 ax=fig.add_subplot(gs[:,0:2]);panel(ax,"A","Annotation-space neighborhoods preserve overlapping functions without forcing one label")
 for cl,col in zip(["A","B","C"],[C["blue"],C["green"],C["peach"]]):q=d[d.membrane_class_v7.eq(cl)];ax.scatter(q.umap1,q.umap2,s=7,c=col,alpha=.28,label=f"Class {cl}",rasterized=True)
 ax.set_xlabel("UMAP 1");ax.set_ylabel("UMAP 2");ax.legend(fontsize=7,markerscale=2);clean(ax,None)
 ax=fig.add_subplot(gs[0,2]);panel(ax,"B","HDBSCAN retains diffuse proteins as noise instead of forcing clusters");cnt=d.cluster.value_counts().sort_index();cols=[C["missing"] if i==-1 else PAL[i%len(PAL)] for i in cnt.index];ax.bar(np.arange(len(cnt)),cnt.values,color=cols);ax.set_xticks(np.arange(len(cnt)),["N" if i==-1 else str(i) for i in cnt.index],rotation=90);ax.set_xlabel("community (N=noise)");ax.set_ylabel("proteins");clean(ax,"y")
 ax=fig.add_subplot(gs[1,2]);panel(ax,"C","Large communities differ in the five annotation axes they carry");sns.heatmap(share,ax=ax,cmap="YlGnBu",linewidths=.4,linecolor="white",cbar_kws={"label":"within-community assertion share"});ax.set_xlabel("");ax.set_ylabel("community");ax.tick_params(axis="x",rotation=45)
 save(fig,"M2_multi_axis_annotation_space","V7.1.1 five-axis classification; Jaccard UMAP n_neighbors=25, min_dist=0.18; HDBSCAN min_cluster_size=90; seed=42.")
if __name__=="__main__":m1();m2();print("final M1/M2 complete")
