from __future__ import annotations
import json,math,textwrap
from collections import Counter,defaultdict
from pathlib import Path
import hdbscan,matplotlib as mpl,matplotlib.pyplot as plt,networkx as nx,numpy as np,pandas as pd,seaborn as sns,umap
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy import sparse
from scipy.stats import kruskal
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors

ROOT=Path(r"C:\Users\Administrator\MemPro_narrative_v3");DATA=ROOT/"data";FIG=ROOT/"figures";QA=ROOT/"qa";REL=Path(r"D:\finale\20_MemPro_V7.1.1_final_20260814\01_release_tables")
PAG=Path(r"D:\finale\07_M1-M8图表修订_20260806\04_pagtn")
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
def knn_retention(X,Y,k=15,n=1000,seed=42,metric="jaccard"):
 rng=np.random.default_rng(seed);idx=rng.choice(len(X),min(n,len(X)),replace=False);xx=X[idx];yy=Y[idx]
 hi=NearestNeighbors(n_neighbors=k+1,metric=metric).fit(xx).kneighbors(return_distance=False)[:,1:];lo=NearestNeighbors(n_neighbors=k+1).fit(yy).kneighbors(return_distance=False)[:,1:]
 return float(np.mean([len(set(a)&set(b))/k for a,b in zip(hi,lo)]))
def scaffold_purity(labels,classes):
 vals=[]
 for c in set(labels)-{-1}:
  q=pd.Series(classes[np.array(labels)==c]);vals.append(q.value_counts(normalize=True).iloc[0] if len(q) else np.nan)
 return float(np.nanmean(vals)) if vals else np.nan

def chemical_space():
 cols=["compound_internal_id","standard_smiles","computed_structural_class","is_approved_drug","is_clinical_candidate","is_endogenous_ligand","is_natural_product","is_chemical_probe","molecular_weight","xlogp","tpsa","rotatable_bond_count"]
 d=pd.read_csv(REL/"compound_master_v71.tsv.gz",sep="\t",usecols=cols,low_memory=False);d=d[d.standard_smiles.notna()].copy();flags=["is_approved_drug","is_clinical_candidate","is_endogenous_ligand","is_natural_product","is_chemical_probe"]
 special=d[np.any(d[flags].fillna(0).astype(bool),axis=1)];rest=d.drop(index=special.index);nmax=15000;take=max(0,nmax-len(special));sample=pd.concat([special,rest.sample(min(take,len(rest)),random_state=42)]).drop_duplicates("compound_internal_id").head(nmax).reset_index(drop=True)
 fps=[];valid=[]
 for i,s in enumerate(sample.standard_smiles):
  m=Chem.MolFromSmiles(str(s))
  if m is None:continue
  fp=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=1024);arr=np.zeros(1024,dtype=np.uint8);Chem.DataStructs.ConvertToNumpyArray(fp,arr);fps.append(arr);valid.append(i)
 sample=sample.iloc[valid].reset_index(drop=True);X=np.asarray(fps,dtype=np.uint8)
 emb=umap.UMAP(n_neighbors=25,min_dist=.08,metric="jaccard",random_state=42,low_memory=True).fit_transform(X);cl=hdbscan.HDBSCAN(min_cluster_size=100,min_samples=20).fit_predict(emb)
 sample["umap1"]=emb[:,0];sample["umap2"]=emb[:,1];sample["cluster"]=cl;sample.to_csv(DATA/"M5_ECFP4_chemical_space.tsv",sep="\t",index=False)
 ret=knn_retention(X,emb,k=15,n=1200,metric="jaccard");pur=scaffold_purity(cl,sample.computed_structural_class.fillna("unknown").values)
 # PAGTN pilot already trained on four descriptors; analyze but retain boundary
 ps=pd.read_csv(PAG/"pagtn_sample.tsv",sep="\t");pe=np.load(PAG/"pagtn_embeddings.npy");pem=umap.UMAP(n_neighbors=25,min_dist=.08,random_state=42).fit_transform(pe);pcl=hdbscan.HDBSCAN(min_cluster_size=35,min_samples=10).fit_predict(pem);ps["umap1"]=pem[:,0];ps["umap2"]=pem[:,1];ps["cluster"]=pcl;ps.to_csv(DATA/"M5_PAGTN_pilot_space.tsv",sep="\t",index=False);pret=knn_retention(pe,pem,k=15,n=1200,metric="euclidean");ppur=scaffold_purity(pcl,ps.structural_class.fillna("unknown").values)
 qa={"ecfp4_n":len(sample),"ecfp4_knn_retention":ret,"ecfp4_cluster_purity":pur,"ecfp4_noise_fraction":float(np.mean(cl==-1)),"pagtn_n":len(ps),"pagtn_knn_retention":pret,"pagtn_cluster_purity":ppur,"pagtn_noise_fraction":float(np.mean(pcl==-1)),"pagtn_boundary":"property-supervised pilot, not universal pretrained embedding"};(QA/"M5_embedding_QA.json").write_text(json.dumps(qa,indent=2),encoding="utf-8")
 return sample,ps,qa
def draw_m5():
 d,p,qa=chemical_space();fig=plt.figure(figsize=(14,9.4));fig.text(.03,.968,"M5 | Chemical space is structured by scaffold, biological status and target context",fontsize=15,weight="bold",va="top");fig.text(.03,.932,"ECFP4 provides the interpretable baseline; the property-supervised PAGTN pilot is retained only as a quality-controlled comparison.",fontsize=9,color=C["muted"],va="top")
 gs=fig.add_gridspec(2,3,left=.05,right=.985,top=.87,bottom=.09,hspace=.34,wspace=.30)
 ax=fig.add_subplot(gs[:,0:2]);panel(ax,"A","ECFP4-UMAP reveals dense chemical neighborhoods without plotting the full 646k compounds")
 hb=ax.hexbin(d.umap1,d.umap2,gridsize=75,cmap="YlGnBu",bins="log",mincnt=1);cb=fig.colorbar(hb,ax=ax,fraction=.025,pad=.015);cb.set_label("log density")
 status=[("Approved",d.is_approved_drug,C["red"]),("Clinical",d.is_clinical_candidate,C["coral"]),("Endogenous",d.is_endogenous_ligand,C["green"])]
 for lab,flag,col in status:
  q=d[pd.to_numeric(flag,errors="coerce").fillna(0).eq(1)];ax.scatter(q.umap1,q.umap2,s=9,c=col,alpha=.38,label=lab,rasterized=True)
 ax.legend(fontsize=7,markerscale=2);ax.set_xlabel("UMAP 1");ax.set_ylabel("UMAP 2");clean(ax,None)
 ax=fig.add_subplot(gs[0,2]);panel(ax,"B","PAGTN pilot clusters a property-supervised representation")
 ax.scatter(p.umap1,p.umap2,s=8,c=[C["missing"] if x==-1 else PAL[x%len(PAL)] for x in p.cluster],alpha=.5,rasterized=True);ax.set_xlabel("UMAP 1");ax.set_ylabel("UMAP 2");clean(ax,None)
 ax=fig.add_subplot(gs[1,2]);panel(ax,"C","Embedding quality determines which map enters the main interpretation")
 met=pd.DataFrame({"metric":["kNN retention","cluster purity","assigned fraction"]*2,"value":[qa["ecfp4_knn_retention"],qa["ecfp4_cluster_purity"],1-qa["ecfp4_noise_fraction"],qa["pagtn_knn_retention"],qa["pagtn_cluster_purity"],1-qa["pagtn_noise_fraction"]],"method":["ECFP4"]*3+["PAGTN pilot"]*3});sns.barplot(met,x="metric",y="value",hue="method",palette=[C["blue"],C["coral"]],ax=ax);ax.set_ylim(0,1);ax.set_ylabel("score");ax.tick_params(axis="x",rotation=25);clean(ax,"y")
 save(fig,"M5_chemical_space_embedding_comparison","V7.1.1 compound master; ECFP4 radius=2/1024 bits; deterministic stratified sample; UMAP and HDBSCAN seed=42. PAGTN pilot trained on four physicochemical descriptors.")

def interaction_network():
 pair=pd.read_csv(REL/"protein_compound_pair_v71.tsv.gz",sep="\t",usecols=["pair_id","target_uniprot_id","compound_internal_id","best_evidence_tier","independent_structure_count"],low_memory=False);deg=pair.target_uniprot_id.value_counts();top=list(deg.head(220).index);q=pair[pair.target_uniprot_id.isin(top)];pi={p:i for i,p in enumerate(top)};cs={c:i for i,c in enumerate(q.compound_internal_id.unique())};rows=q.target_uniprot_id.map(pi).values;cols=q.compound_internal_id.map(cs).values;M=sparse.csr_matrix((np.ones(len(q)),(rows,cols)),shape=(len(top),len(cs)));shared=(M@M.T).toarray();np.fill_diagonal(shared,0);G=nx.Graph();G.add_nodes_from(range(len(top)))
 for i in range(len(top)):
  ids=np.argsort(shared[i])[-5:]
  for j in ids:
   if shared[i,j]>=3:G.add_edge(i,int(j),weight=float(shared[i,j]))
 comm=nx.algorithms.community.louvain_communities(G,weight="weight",seed=42);cid={n:k for k,c in enumerate(comm) for n in c};nd=pd.DataFrame({"node":range(len(top)),"target_uniprot_id":top,"degree":[G.degree(i) for i in range(len(top))],"community":[cid.get(i,-1) for i in range(len(top))]});nd.to_csv(DATA/"M6_target_network_nodes.tsv",sep="\t",index=False);pd.DataFrame([(a,b,x["weight"]) for a,b,x in G.edges(data=True)],columns=["source","target","shared_compounds"]).to_csv(DATA/"M6_target_network_edges.tsv",sep="\t",index=False)
 return G,nd,pair
def draw_m6():
 G,nd,pair=interaction_network();site=pd.read_csv(REL/"interaction_site_v71.tsv.gz",sep="\t",usecols=["compound_internal_id","membrane_side","coordinate_docking_eligible","site_residue_status"],low_memory=False);cmp=pd.read_csv(REL/"compound_master_v71.tsv.gz",sep="\t",usecols=["compound_internal_id","xlogp","tpsa"],low_memory=False);sc=site.merge(cmp,on="compound_internal_id",how="left");sc=sc[~sc.membrane_side.fillna("unknown").str.lower().eq("unknown")];sc.to_csv(DATA/"M6_membrane_side_compound_properties.tsv",sep="\t",index=False)
 groups=[x.dropna().values for _,x in sc.groupby("membrane_side").xlogp if len(x.dropna())>=10];kw=kruskal(*groups) if len(groups)>=2 else None;(QA/"M6_membrane_side_stats.json").write_text(json.dumps({"groups":sc.membrane_side.value_counts().to_dict(),"kruskal_xlogp_H":None if kw is None else kw.statistic,"kruskal_xlogp_p":None if kw is None else kw.pvalue},indent=2),encoding="utf-8")
 fig=plt.figure(figsize=(14,9.4));fig.text(.03,.968,"M6 | Protein-ligand communities and membrane-side pockets reveal interaction organization",fontsize=15,weight="bold",va="top");fig.text(.03,.932,"Shared-ligand communities summarize polypharmacology, while site-side comparisons retain coordinate and assay limitations.",fontsize=9,color=C["muted"],va="top")
 gs=fig.add_gridspec(2,3,left=.05,right=.985,top=.87,bottom=.09,hspace=.36,wspace=.32)
 ax=fig.add_subplot(gs[:,0:2]);panel(ax,"A","Shared ligands organize high-degree targets into polypharmacology communities")
 pos=nx.spring_layout(G,seed=42,weight="weight",iterations=100);sizes=35+np.sqrt(nd.degree.values)*20;cols=[PAL[x%len(PAL)] for x in nd.community];nx.draw_networkx_edges(G,pos,ax=ax,width=.25,alpha=.18,edge_color=C["muted"]);nx.draw_networkx_nodes(G,pos,ax=ax,node_size=sizes,node_color=cols,alpha=.75,linewidths=.3,edgecolors="white");ax.axis("off")
 ax=fig.add_subplot(gs[0,2]);panel(ax,"B","Ligand lipophilicity differs across annotated membrane-side contexts")
 if len(sc):sns.violinplot(sc,x="membrane_side",y="xlogp",ax=ax,palette="Set2",inner="quart",cut=0);ax.tick_params(axis="x",rotation=35);ax.set_xlabel("");ax.set_ylabel("XlogP (1st-99th pct shown)");
 ax=fig.add_subplot(gs[1,2]);panel(ax,"C","Polar surface area provides an independent pocket-context contrast")
 if len(sc):sns.boxenplot(sc,x="membrane_side",y="tpsa",ax=ax,palette="Set2",showfliers=False);ax.tick_params(axis="x",rotation=35);ax.set_xlabel("");ax.set_ylabel("TPSA (A²)")
 save(fig,"M6_polypharmacology_and_membrane_side","V7.1.1 released pairs, interaction sites and compound properties; target projection restricted to 220 highest-degree proteins; shared-ligand edges >=3; Louvain seed=42.")
if __name__=="__main__":draw_m5();draw_m6();(QA/"STAGE3_COMPLETE.json").write_text(json.dumps({"M5":"complete","M6":"complete","seed":42},indent=2),encoding="utf-8");print("M5/M6 complete")
