from pathlib import Path
import numpy as np,pandas as pd,matplotlib as mpl,matplotlib.pyplot as plt,seaborn as sns,plotly.express as px
from pipeline_lib import load_config
mpl.rcParams.update({"font.family":"Arial","font.size":8,"svg.fonttype":"none","pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
def save(fig,path):
 for e in ("png","pdf","svg"):fig.savefig(path.with_suffix('.'+e),dpi=600,bbox_inches="tight",facecolor="white")
def main():
 cfg=load_config();out=cfg["output_dir"];fd=out/"figures";d3=pd.read_csv(out/"protein_umap3d.tsv",sep="\t");d2=pd.read_csv(out/"protein_umap2d.tsv",sep="\t");
 # class colours: target labels only color the post-hoc display.
 classes=sorted(d3.target_class.fillna("Unknown").unique()); pal=dict(zip(classes,sns.color_palette("tab20",n_colors=len(classes))))
 def p3(color_col,name,iscluster=False):
  fig=plt.figure(figsize=(6.8,5.5));ax=fig.add_subplot(projection="3d");vals=d3[color_col]
  cats=sorted(vals.unique())
  cp=dict(zip(cats,sns.color_palette("husl",len(cats)))) if iscluster else pal
  if iscluster and -1 in cats:cp[-1]="#b8b8b8"
  for c in cats:
   z=d3[vals==c];ax.scatter(z.UMAP1,z.UMAP2,z.UMAP3,s=4,alpha=.55,c=[cp[c]],label=("Noise" if c==-1 else str(c)),edgecolors="none",rasterized=True)
  ax.set(xlabel="UMAP1",ylabel="UMAP2",zlabel="UMAP3");ax.view_init(elev=22,azim=42);ax.legend(title=("HDBSCAN cluster" if iscluster else "Target class"),loc="center left",bbox_to_anchor=(1.04,.5),markerscale=2,fontsize=6);save(fig,fd/name);plt.close(fig)
 p3("target_class","protein_umap3d_target_class");p3("cluster_id","protein_umap3d_cluster",True)
 def p2(col,name,cl=False):
  fig,ax=plt.subplots(figsize=(5.4,4.4));vals=d2[col];cats=sorted(vals.unique());cp=dict(zip(cats,sns.color_palette("husl",len(cats)))) if cl else pal
  if cl and -1 in cats:cp[-1]="#b8b8b8"
  for c in cats:
   z=d2[vals==c];ax.scatter(z.UMAP1,z.UMAP2,s=5,alpha=.5,c=[cp[c]],label=("Noise" if c==-1 else str(c)),edgecolors="none",rasterized=True)
  ax.set(xlabel="UMAP1",ylabel="UMAP2");ax.legend(loc="center left",bbox_to_anchor=(1.02,.5),fontsize=6,frameon=False);save(fig,fd/name);plt.close(fig)
 p2("target_class","protein_umap2d_target_class");p2("cluster_id","protein_umap2d_cluster",True)
 comp=pd.read_csv(out/"cluster_target_composition.tsv",sep="\t");piv=comp.pivot(index="cluster_id",columns="target_class",values="proportion").fillna(0).sort_index();fig,ax=plt.subplots(figsize=(7,3.8));piv.plot(kind="bar",stacked=True,ax=ax,colormap="tab20",width=.9);ax.set(xlabel="HDBSCAN cluster",ylabel="Target-class proportion",ylim=(0,1));ax.legend(loc="center left",bbox_to_anchor=(1.02,.5),fontsize=6,title="Target class");save(fig,fd/"cluster_target_composition");plt.close(fig)
 html=px.scatter_3d(d3,x="UMAP1",y="UMAP2",z="UMAP3",color="target_class",hover_data=["uniprot_id","target_class","cluster_id","sequence_length"],opacity=.65);html.update_layout(template="simple_white",legend_title_text="Target class");html.write_html(fd/"protein_umap3d_interactive.html",include_plotlyjs=True)
if __name__=="__main__":main()
