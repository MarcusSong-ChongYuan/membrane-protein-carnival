from __future__ import annotations

import json
import math
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import umap
from PIL import Image
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

ROOT=Path(r"C:\Users\Administrator\MemPro_narrative_v3");DATA=ROOT/"data";FIG=ROOT/"figures";QA=ROOT/"qa"
REL=Path(r"D:\finale\20_MemPro_V7.1.1_final_20260814\01_release_tables")
ANATOMY=Path(r"C:\Users\Administrator\Desktop\gkg sjf summer intern\MemPro最终展示\M1-M8新叙事终版_20260814\05_assets\human_anatomy_base_v1.png")
for p in (DATA,FIG,QA):p.mkdir(parents=True,exist_ok=True)
C={"blue":"#7b95c6","cyan":"#49c2d9","lcyan":"#a1d8e8","green":"#67a583","lgreen":"#a2c986","yellow":"#fded95","peach":"#ffc1a6","salmon":"#f59c7c","coral":"#f47254","red":"#c85e62","ink":"#27313d","muted":"#66717e","grid":"#e5e9ee","missing":"#b9bec6"}
PAL=[C["blue"],C["cyan"],C["green"],C["lgreen"],C["yellow"],C["peach"],C["salmon"],C["red"]]
mpl.rcParams.update({"font.family":"Arial","font.size":8.2,"axes.titlesize":9.3,"axes.labelsize":8.2,"xtick.labelsize":7.2,"ytick.labelsize":7.2,"text.color":C["ink"],"axes.labelcolor":C["ink"],"axes.titlecolor":C["ink"],"axes.edgecolor":"#c7cdd4","grid.color":C["grid"],"pdf.fonttype":42,"svg.fonttype":"none","legend.frameon":False})
sns.set_style("whitegrid")
def panel(ax,l,t):ax.set_title(f"{l}  {t}",loc="left",weight="bold",pad=8)
def clean(ax,grid="x"):
 ax.spines[["top","right"]].set_visible(False);ax.grid(False)
 if grid:ax.grid(axis=grid,alpha=.7)
def save(fig,stem,src):
 fig.text(.035,.018,"Source: "+src,fontsize=7,color=C["muted"])
 for ext in ("png","svg","pdf"):fig.savefig(FIG/f"{stem}.{ext}",dpi=450 if ext=="png" else None,bbox_inches="tight",facecolor="white")
 plt.close(fig)

TISSUE_ORGAN={"brain":"Brain","cerebral cortex":"Brain","cerebellum":"Brain","heart muscle":"Heart","lung":"Lung","liver":"Liver","kidney":"Kidney","stomach":"Stomach","colon":"Large intestine","rectum":"Large intestine","small intestine":"Small intestine","duodenum":"Small intestine","spleen":"Spleen","pancreas":"Pancreas","thyroid gland":"Thyroid","skin":"Skin","bone marrow":"Bone marrow","testis":"Reproductive","ovary":"Reproductive","endometrium":"Reproductive","prostate":"Reproductive","urinary bladder":"Bladder"}
COORD={"Brain":(.50,.09),"Heart":(.52,.34),"Lung":(.43,.30),"Liver":(.43,.42),"Kidney":(.59,.49),"Stomach":(.56,.43),"Large intestine":(.50,.54),"Small intestine":(.50,.57),"Spleen":(.61,.43),"Pancreas":(.53,.47),"Thyroid":(.50,.21),"Skin":(.68,.36),"Bone marrow":(.61,.73),"Reproductive":(.50,.68),"Bladder":(.50,.66)}
SYSTEM_ORGAN={"nervous system":"Brain","circulatory system":"Heart","cardiovascular system":"Heart","respiratory system":"Lung","digestive system":"Large intestine","alimentary part of gastrointestinal system":"Large intestine","excretory system":"Kidney","hematopoietic system":"Bone marrow","lymphoid system":"Spleen","integumental system":"Skin","reproductive system":"Reproductive","musculoskeletal system":"Bone marrow","exocrine system":"Pancreas"}

def body(ax,table,val,title,cmap,suffix,n=7):
 im=np.asarray(Image.open(ANATOMY).convert("RGB"));h,w=im.shape[:2];ax.imshow(im);ax.set_xlim(0,w);ax.set_ylim(h,0);ax.set_aspect("equal");ax.axis("off");panel(ax,"A",title)
 t=table[table.organ.isin(COORD)].nlargest(n,val);norm=mpl.colors.Normalize(t[val].min(),t[val].max());cm=mpl.colormaps[cmap]
 for i,(_,r) in enumerate(t.iterrows()):
  x0,y0=COORD[r.organ][0]*w,COORD[r.organ][1]*h;left=i%2==0;tx=(.03 if left else .97)*w;ty=(.23+i//2*.15)*h;el=(.28 if left else .72)*w
  ax.scatter(x0,y0,s=70,c=[cm(norm(r[val]))],edgecolor="white",lw=.7,zorder=5);ax.plot([x0,el,tx],[y0,ty,ty],color=C["muted"],lw=.65);ax.text(tx,ty,f"{r.organ}\n{int(r[val]):,} {suffix}",ha="left" if left else "right",va="center",fontsize=7.3,weight="bold")

def expression_matrices():
 p=REL/"expression_measurement_v71.tsv.gz";cols=["target_uniprot_id","expression_layer_v711","tissue","value","measurement_status_v711","detection_status_v711","mapped_measured_denominator_eligible_v71"]
 rnas=[];ihcs=[]
 for c in pd.read_csv(p,sep="\t",usecols=cols,chunksize=220000,low_memory=False):
  ok=pd.to_numeric(c.mapped_measured_denominator_eligible_v71,errors="coerce").fillna(0).eq(1)&c.measurement_status_v711.eq("MEASURED")
  r=c[ok&c.expression_layer_v711.eq("tissue_RNA")].copy();r["value"]=pd.to_numeric(r.value,errors="coerce");rnas.append(r[["target_uniprot_id","tissue","value","detection_status_v711"]])
  i=c[ok&c.expression_layer_v711.eq("normal_tissue_IHC")].copy();ihcs.append(i[["target_uniprot_id","tissue","detection_status_v711"]])
 rna=pd.concat(rnas,ignore_index=True);ihc=pd.concat(ihcs,ignore_index=True)
 R=rna.pivot_table(index="target_uniprot_id",columns="tissue",values="value",aggfunc="mean",fill_value=0)
 I=ihc.assign(det=ihc.detection_status_v711.eq("DETECTED").astype(float)).pivot_table(index="target_uniprot_id",columns="tissue",values="det",aggfunc="mean")
 R.to_pickle(DATA/"M3_tissue_rna_matrix.pkl");I.to_pickle(DATA/"M3_ihc_matrix.pkl")
 return R,I,rna

def draw_m3():
 R,I,rna=expression_matrices();X=np.log1p(R.values);X=(X-X.mean(1,keepdims=True))/(X.std(1,keepdims=True)+1e-6)
 emb=umap.UMAP(n_neighbors=30,min_dist=.12,metric="cosine",random_state=42).fit_transform(X)
 pm=pd.read_csv(REL/"protein_master_v71.tsv.gz",sep="\t",usecols=["canonical_uniprot_accession","membrane_class_v7"]).set_index("canonical_uniprot_accession");meta=pm.reindex(R.index);out=pd.DataFrame({"target_uniprot_id":R.index,"umap1":emb[:,0],"umap2":emb[:,1],"membrane_class":meta.membrane_class_v7.values});out.to_csv(DATA/"M3_expression_umap.tsv",sep="\t",index=False)
 # tissue clustering
 tcorr=np.corrcoef(np.log1p(R.values).T);order=leaves_list(linkage(pdist(np.log1p(R.values).T,metric="correlation"),method="average"));tissues=list(R.columns[order]);top=tissues[:20]
 # RNA/IHC concordance per shared tissue
 rows=[]
 for t in sorted(set(R.columns)&set(I.columns)):
  q=pd.concat([R[t],I[t]],axis=1,join="inner").dropna();q.columns=["rna","ihc"]
  if len(q)>=100:
   rho,p=spearmanr(np.log1p(q.rna),q.ihc);rows.append((t,len(q),rho,p))
 cor=pd.DataFrame(rows,columns=["tissue","n","spearman_rho","pvalue"]).sort_values("spearman_rho",ascending=False);cor.to_csv(DATA/"M3_RNA_IHC_concordance.tsv",sep="\t",index=False)
 organ=rna.assign(organ=rna.tissue.str.lower().map(TISSUE_ORGAN));organ=organ[organ.detection_status_v711.eq("DETECTED")&organ.organ.notna()].groupby("organ").target_uniprot_id.nunique().rename("detected_proteins").reset_index();organ.to_csv(DATA/"M3_body_expression.tsv",sep="\t",index=False)
 fig=plt.figure(figsize=(14,9.5));fig.text(.03,.968,"M3 | Expression patterns connect membrane proteins to tissues and cell contexts",fontsize=15,weight="bold",va="top");fig.text(.03,.932,"Body-wide coverage, expression-space neighborhoods and RNA-IHC concordance are shown as distinct measurement questions.",fontsize=9,color=C["muted"],va="top")
 gs=fig.add_gridspec(2,3,left=.035,right=.985,top=.87,bottom=.09,width_ratios=[.9,1.25,1.0],hspace=.34,wspace=.30)
 ax=fig.add_subplot(gs[:,0]);body(ax,organ,"detected_proteins","Body-wide RNA detection spans major organs","YlGnBu","proteins")
 ax=fig.add_subplot(gs[0,1]);panel(ax,"B","Tissue-expression UMAP reveals membrane-protein neighborhoods")
 for cl,col in zip(["A","B","C"],[C["blue"],C["green"],C["peach"]]):
  q=out[out.membrane_class.eq(cl)];ax.scatter(q.umap1,q.umap2,s=7,c=col,alpha=.28,label=f"Class {cl}",rasterized=True)
 ax.legend(fontsize=7,markerscale=2);ax.set_xlabel("UMAP 1");ax.set_ylabel("UMAP 2");clean(ax,None)
 ax=fig.add_subplot(gs[1,1]);panel(ax,"C","Hierarchical clustering identifies co-expression tissue blocks")
 sel=R[top].copy();sel=np.log1p(sel);sel=(sel-sel.mean())/(sel.std()+1e-6);sample=sel.var(1).nlargest(150).index;sns.heatmap(sel.loc[sample,top].T,ax=ax,cmap="vlag",center=0,xticklabels=False,cbar_kws={"label":"tissue-wise z score"});ax.set_xlabel("150 variable membrane proteins");ax.set_ylabel("")
 ax=fig.add_subplot(gs[:,2]);panel(ax,"D","RNA-IHC concordance varies across tissues")
 c=cor.sort_values("spearman_rho").tail(18);y=np.arange(len(c));ax.hlines(y,0,c.spearman_rho,color=C["grid"],lw=2);ax.scatter(c.spearman_rho,y,c=c.n,cmap="YlGnBu",s=55);ax.axvline(0,color=C["muted"],lw=.7);ax.set_yticks(y,c.tissue);ax.set_xlabel("Spearman rho: log1p RNA vs IHC detection");clean(ax)
 save(fig,"M3_expression_space_and_RNA_IHC","HPA 25.1 V7.1.1 corrected measurement layers; only mapped and measured records enter denominators; UMAP seed=42.")

def draw_m4():
 rel=pd.read_csv(REL/"protein_disease_relation_v71.tsv.gz",sep="\t",low_memory=False);sys=pd.read_csv(REL/"disease_anatomical_system_v71.tsv.gz",sep="\t",low_memory=False);sys["organ"]=sys.anatomical_system_name.str.lower().map(SYSTEM_ORGAN)
 j=rel.merge(sys[["canonical_disease_id","anatomical_system_name","organ"]].drop_duplicates(),on="canonical_disease_id",how="inner");j["pair"]=j.target_uniprot_id+"|"+j.canonical_disease_id
 R=pd.read_pickle(DATA/"M3_tissue_rna_matrix.pkl");det=R.gt(1.0);organ_sets={}
 for org in set(TISSUE_ORGAN.values()):
  ts=[t for t in R.columns if TISSUE_ORGAN.get(str(t).lower())==org]
  if ts:organ_sets[org]=set(R.index[det[ts].any(axis=1)])
 rng=np.random.default_rng(42);universe=np.array(sorted(set(R.index)&set(rel.target_uniprot_id)));rows=[]
 for org,q in j.dropna(subset=["organ"]).groupby("organ"):
  dp=set(q.target_uniprot_id)&set(universe);ep=organ_sets.get(org,set())&set(universe)
  if not dp or not ep:continue
  obs=len(dp&ep)/len(dp);null=np.empty(500)
  for i in range(500):null[i]=len(set(rng.choice(universe,size=len(dp),replace=False))&ep)/len(dp)
  z=(obs-null.mean())/(null.std()+1e-9);p=(1+(null>=obs).sum())/501;rows.append((org,len(dp),len(ep),obs,null.mean(),z,p))
 con=pd.DataFrame(rows,columns=["organ","disease_proteins","expressed_proteins","observed_fraction","null_mean","zscore","empirical_p"]);con["fdr"]=np.minimum(1,con.empirical_p.rank(method="max")/len(con)*con.empirical_p.sort_values().values[-1] if len(con) else 1);con.to_csv(DATA/"M4_disease_expression_concordance.tsv",sep="\t",index=False)
 burden=j.dropna(subset=["organ"]).groupby("organ").pair.nunique().rename("protein_disease_pairs").reset_index();burden.to_csv(DATA/"M4_body_disease.tsv",sep="\t",index=False)
 # system x membrane-role enrichment (relative share)
 role=pd.read_csv(REL/"protein_cross_classification_summary_v71.tsv.gz",sep="\t",usecols=["target_uniprot_id","membrane_role_primary_v63"]);jj=j.merge(role,on="target_uniprot_id",how="left");tab=pd.crosstab(jj.anatomical_system_name,jj.membrane_role_primary_v63);toprow=tab.sum(1).nlargest(10).index;topcol=tab.sum().nlargest(9).index;tab=tab.loc[toprow,topcol];share=tab.div(tab.sum(1),axis=0);share.to_csv(DATA/"M4_system_role_share.tsv",sep="\t")
 fig=plt.figure(figsize=(14,9.4));fig.text(.03,.968,"M4 | Disease anatomy aligns unevenly with normal-tissue membrane-protein expression",fontsize=15,weight="bold",va="top");fig.text(.03,.932,"Ontology-derived multi-label anatomy is integrated with HPA expression and evaluated against a permutation null model.",fontsize=9,color=C["muted"],va="top")
 gs=fig.add_gridspec(2,3,left=.035,right=.985,top=.87,bottom=.09,width_ratios=[.9,1.1,1.15],hspace=.35,wspace=.30)
 ax=fig.add_subplot(gs[:,0]);body(ax,burden,"protein_disease_pairs","Disease-associated pairs map across organ systems","YlOrRd","pairs")
 ax=fig.add_subplot(gs[:,1]);panel(ax,"B","Permutation testing identifies disease-expression concordance above chance")
 c=con.sort_values("zscore");y=np.arange(len(c));ax.hlines(y,0,c.zscore,color=C["grid"],lw=2);ax.scatter(c.zscore,y,c=-np.log10(c.empirical_p),cmap="YlOrRd",s=70,edgecolor="white",lw=.5);ax.axvline(0,color=C["muted"],lw=.8);ax.set_yticks(y,c.organ);ax.set_xlabel("concordance z score (500 permutations)");clean(ax)
 ax=fig.add_subplot(gs[0,2]);panel(ax,"C","Disease systems emphasize different membrane roles")
 sns.heatmap(share,ax=ax,cmap="YlGnBu",linewidths=.4,linecolor="white",cbar_kws={"label":"within-system pair share"});ax.set_xlabel("");ax.set_ylabel("");ax.tick_params(axis="x",rotation=45)
 ax=fig.add_subplot(gs[1,2]);panel(ax,"D","Evidence strength is not uniform across anatomical systems")
 ev=pd.crosstab(j.anatomical_system_name,j.best_evidence_level,values=j.pair,aggfunc=pd.Series.nunique).fillna(0);ev=ev.loc[ev.sum(1).nlargest(10).index];sns.heatmap(np.log1p(ev),ax=ax,cmap=sns.light_palette(C["red"],as_cmap=True),linewidths=.4,linecolor="white",cbar_kws={"label":"log1p(unique pairs)"});ax.set_xlabel("best evidence level");ax.set_ylabel("")
 save(fig,"M4_disease_expression_concordance","Open Targets 26.06, UniProtKB, MONDO, Disease Ontology, Uberon and HPA 25.1; ontology multi-label mapping; permutation seed=42.")

if __name__=="__main__":draw_m3();draw_m4();(QA/"STAGE2_COMPLETE.json").write_text(json.dumps({"M3":"complete","M4":"complete","seed":42,"permutations":500},indent=2),encoding="utf-8");print("M3/M4 complete")
