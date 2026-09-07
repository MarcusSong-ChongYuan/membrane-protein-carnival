from __future__ import annotations
import json,textwrap
from collections import Counter
from pathlib import Path
import matplotlib as mpl,matplotlib.pyplot as plt,numpy as np,pandas as pd,seaborn as sns

ROOT=Path(r"C:\Users\Administrator\MemPro_narrative_v3");DATA=ROOT/"data";FIG=ROOT/"figures";QA=ROOT/"qa";REL=Path(r"D:\finale\20_MemPro_V7.1.1_final_20260814\01_release_tables")
for p in (DATA,FIG,QA):p.mkdir(parents=True,exist_ok=True)
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

def evidence_profiles():
 cols=["evidence_id","target_uniprot_id","compound_internal_id","source_database","evidence_tier_recomputed_v7","evidence_modality_v7","standard_value_nM","pdb_ids","binding_site_residues","lineage_resolution_v7"]
 d=pd.read_csv(REL/"positive_interaction_evidence_v71.tsv.gz",sep="\t",usecols=cols,low_memory=False)
 d["quantitative"]=d.standard_value_nM.notna();d["structure"]=d.pdb_ids.fillna("").str.strip().ne("");d["site"]=d.binding_site_residues.fillna("").str.strip().ne("");d["functional"]=d.evidence_modality_v7.fillna("").str.contains("functional|pharmacology|enzyme",case=False,regex=True);d["lineage_keyed"]=d.lineage_resolution_v7.fillna("").str.contains("unique|detected|resolved",case=False,regex=True)
 flags=["quantitative","functional","structure","site","lineage_keyed"];combo=Counter("+".join(x for x in flags if bool(r[x])) or "identity_only" for _,r in d.iterrows());pd.DataFrame(combo.most_common(),columns=["profile","evidence_records"]).to_csv(DATA/"M7_evidence_profile_intersections.tsv",sep="\t",index=False)
 sl=pd.crosstab(d.source_database,d.lineage_resolution_v7).fillna(0);sl=sl.loc[sl.sum(1).nlargest(10).index];sl.to_csv(DATA/"M7_source_lineage_matrix.tsv",sep="\t")
 mt=pd.crosstab(d.evidence_modality_v7,d.evidence_tier_recomputed_v7).fillna(0);mt=mt.loc[mt.sum(1).nlargest(10).index];mt.to_csv(DATA/"M7_modality_tier_matrix.tsv",sep="\t")
 return d,combo,sl,mt,flags
def draw_upset(ax,combo,flags):
 arr=combo.most_common(11);x=np.arange(len(arr));vals=[v for _,v in arr];ax.bar(x,vals,color=C["blue"]);ax.set_yscale("log");ax.set_ylabel("evidence records (log scale)");ax.set_xticks(x,[]);clean(ax,"y")
 ia=ax.inset_axes([0,-.43,1,.34]);ia.set_xlim(-.5,len(arr)-.5);ia.set_ylim(-.5,len(flags)-.5);ia.axis("off")
 for i,s in enumerate(flags):ia.text(-.6,len(flags)-1-i,s.replace("_"," "),ha="right",va="center",fontsize=7)
 for j,(key,_) in enumerate(arr):
  active=set(key.split("+"));ys=[]
  for i,s in enumerate(flags):
   y=len(flags)-1-i;on=s in active;ia.scatter(j,y,s=28,color=C["ink"] if on else "#d9dde2");ys.append(y) if on else None
  if len(ys)>1:ia.plot([j,j],[min(ys),max(ys)],color=C["ink"],lw=.9)
def draw_m7():
 d,combo,sl,mt,flags=evidence_profiles();fig=plt.figure(figsize=(14,10));fig.text(.03,.972,"M7 | Evidence profiles and lineage expose depth, redundancy and release reliability",fontsize=15,weight="bold",va="top");fig.text(.03,.938,"Evidence combinations are preserved explicitly; contributing databases are not mislabeled as independent experiments.",fontsize=9,color=C["muted"],va="top")
 gs=fig.add_gridspec(2,2,left=.075,right=.985,top=.87,bottom=.09,hspace=.58,wspace=.34,height_ratios=[1.05,1])
 ax=fig.add_subplot(gs[0,:]);panel(ax,"A","Evidence records occupy recurring quantitative, functional, structural and site profiles");draw_upset(ax,combo,flags)
 ax=fig.add_subplot(gs[1,0]);panel(ax,"B","Lineage status differs across contributing databases");sns.heatmap(np.log1p(sl),ax=ax,cmap="YlGnBu",linewidths=.4,linecolor="white",cbar_kws={"label":"log1p(records)"});ax.set_xlabel("lineage resolution");ax.set_ylabel("");ax.tick_params(axis="x",rotation=35)
 ax=fig.add_subplot(gs[1,1]);panel(ax,"C","Evidence modalities map unevenly to recomputed BE tiers");sns.heatmap(np.log1p(mt),ax=ax,cmap=sns.light_palette(C["red"],as_cmap=True),linewidths=.4,linecolor="white",cbar_kws={"label":"log1p(records)"});ax.set_xlabel("recomputed BE tier");ax.set_ylabel("")
 save(fig,"M7_evidence_profiles_and_lineage","V7.1.1 positive evidence and experiment-lineage fields; combinations are record-level, non-exclusive evidence profiles.")

def pct(s):return s.rank(pct=True,method="average").fillna(0)
def priority_space():
 pm=pd.read_csv(REL/"protein_master_v71.tsv.gz",sep="\t",usecols=["canonical_uniprot_accession","approved_symbol","membrane_class_v7"]).rename(columns={"canonical_uniprot_accession":"target_uniprot_id"})
 pair=pd.read_csv(REL/"protein_compound_pair_v71.tsv.gz",sep="\t",usecols=["target_uniprot_id","compound_internal_id","independent_structure_count","independent_evidence_modality_count"],low_memory=False)
 pg=pair.groupby("target_uniprot_id").agg(pair_count=("compound_internal_id","nunique"),structure_count=("independent_structure_count","sum"),modality_count=("independent_evidence_modality_count","max"))
 dis=pd.read_csv(REL/"protein_disease_relation_v71.tsv.gz",sep="\t",usecols=["target_uniprot_id","canonical_disease_id"]).groupby("target_uniprot_id").canonical_disease_id.nunique().rename("disease_count")
 R=pd.read_pickle(DATA/"M3_tissue_rna_matrix.pkl");breadth=R.gt(1).sum(1).rename("expression_breadth")
 site=pd.read_csv(REL/"interaction_site_v71.tsv.gz",sep="\t",usecols=["target_uniprot_id","compound_internal_id","coordinate_docking_eligible"],low_memory=False);sg=site.groupby("target_uniprot_id").agg(site_count=("compound_internal_id","nunique"),coordinate_pair_count=("coordinate_docking_eligible","sum"))
 d=pm.set_index("target_uniprot_id").join([pg,dis,breadth,sg]).fillna(0).reset_index();d["biology_score"]=(pct(np.log1p(d.pair_count))+pct(np.log1p(d.disease_count))+pct(d.modality_count))/3;d["readiness_score"]=(pct(np.log1p(d.structure_count))+pct(np.log1p(d.site_count))+pct(np.log1p(d.coordinate_pair_count)))/3
 d["quadrant"]=np.select([(d.biology_score>=.7)&(d.readiness_score>=.7),(d.biology_score>=.7)&(d.readiness_score<.3),(d.biology_score<.3)&(d.readiness_score>=.7)], ["actionable high-value","high-value knowledge gap","structure-rich lower-priority"],default="intermediate")
 d.to_csv(DATA/"M8_protein_priority_space.tsv",sep="\t",index=False)
 # pair-level G tiers
 pairkey=pair.target_uniprot_id+"|"+pair.compound_internal_id;g1=set(site.loc[pd.to_numeric(site.coordinate_docking_eligible,errors="coerce").fillna(0).eq(1),"target_uniprot_id"]+"|"+site.loc[pd.to_numeric(site.coordinate_docking_eligible,errors="coerce").fillna(0).eq(1),"compound_internal_id"]);g2=set(site.target_uniprot_id+"|"+site.compound_internal_id)-g1
 ev=pd.read_csv(REL/"positive_interaction_evidence_v71.tsv.gz",sep="\t",usecols=["target_uniprot_id","compound_internal_id","pdb_ids"],low_memory=False);pdb=set(ev.loc[ev.pdb_ids.fillna("").str.strip().ne(""),"target_uniprot_id"]+"|"+ev.loc[ev.pdb_ids.fillna("").str.strip().ne(""),"compound_internal_id"]);allp=set(pairkey);g3=(pdb-g1-g2)&allp;g4=allp-g1-g2-g3;tiers={"G1 coordinate-ready":len(g1&allp),"G2 site/no complete coordinates":len(g2&allp),"G3 structure/no site":len(g3),"G4 interaction-only":len(g4)};pd.DataFrame(tiers.items(),columns=["tier","pairs"]).to_csv(DATA/"M8_pair_readiness_tiers.tsv",sep="\t",index=False)
 return d,tiers
def pareto(d):
 q=d.sort_values(["biology_score","readiness_score"],ascending=False);front=[];best=-1
 for _,r in q.iterrows():
  if r.readiness_score>best:front.append(r);best=r.readiness_score
 return pd.DataFrame(front)
def draw_m8():
 d,tiers=priority_space();front=pareto(d);fig=plt.figure(figsize=(14,9.4));fig.text(.03,.968,"M8 | Biological value and structural readiness expose actionable targets and knowledge gaps",fontsize=15,weight="bold",va="top");fig.text(.03,.932,"Transparent component percentiles define a two-dimensional priority space; no single opaque score controls release inclusion.",fontsize=9,color=C["muted"],va="top")
 gs=fig.add_gridspec(2,3,left=.055,right=.985,top=.87,bottom=.09,hspace=.36,wspace=.32)
 ax=fig.add_subplot(gs[:,0:2]);panel(ax,"A","The priority landscape separates immediately actionable proteins from high-value coverage deserts")
 hb=ax.hexbin(d.biology_score,d.readiness_score,gridsize=40,cmap="YlGnBu",bins="log",mincnt=1);fig.colorbar(hb,ax=ax,fraction=.025,pad=.015,label="log protein density");ax.axvline(.7,color=C["muted"],ls="--",lw=.8);ax.axhline(.7,color=C["muted"],ls="--",lw=.8);ax.axvline(.3,color=C["muted"],ls=":",lw=.7);ax.axhline(.3,color=C["muted"],ls=":",lw=.7);ax.plot(front.biology_score,front.readiness_score,color=C["red"],lw=1.2,label="Pareto frontier")
 gaps=d[d.quadrant.eq("high-value knowledge gap")].nlargest(10,"biology_score");ax.scatter(gaps.biology_score,gaps.readiness_score,s=35,c=C["coral"],edgecolor="white",lw=.5); 
 for _,r in gaps.iterrows():ax.text(r.biology_score+.008,r.readiness_score,str(r.approved_symbol),fontsize=6.8)
 ax.set_xlabel("biological-value percentile score");ax.set_ylabel("structural-readiness percentile score");ax.legend(fontsize=7);clean(ax,None)
 ax=fig.add_subplot(gs[0,2]);panel(ax,"B","Four protein-level regions define distinct follow-up strategies");qc=d.quadrant.value_counts();ax.pie(qc.values,labels=[x.replace(" ","\n") for x in qc.index],colors=PAL[:len(qc)],startangle=90,autopct=lambda p:f"{p:.1f}%" if p>5 else "",textprops={"fontsize":7},wedgeprops={"edgecolor":"white"})
 ax=fig.add_subplot(gs[1,2]);panel(ax,"C","Pair-level readiness tiers preserve interactions that are not yet box-ready");td=pd.Series(tiers).sort_values();ax.hlines(np.arange(len(td)),0,td.values,color=C["grid"],lw=3);ax.scatter(td.values,np.arange(len(td)),s=65,color=[C["green"],C["cyan"],C["yellow"],C["peach"]]);ax.set_yticks(np.arange(len(td)),td.index);ax.set_xscale("log");ax.set_xlabel("unique protein-compound pairs (log scale)");clean(ax)
 save(fig,"M8_value_readiness_and_knowledge_gaps","V7.1.1 proteins, pairs, diseases, HPA expression, structures and interaction sites; component scores are equal-weight percentile summaries shown transparently.")
if __name__=="__main__":draw_m7();draw_m8();(QA/"STAGE4_COMPLETE.json").write_text(json.dumps({"M7":"complete","M8":"complete"},indent=2),encoding="utf-8");print("M7/M8 complete")
