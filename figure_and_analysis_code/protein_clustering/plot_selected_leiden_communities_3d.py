"""3D PCA community-shape view for selected, independently interpretable modules.

This is a curated community gallery, not the full-proteome map.  It avoids the
false visual density caused by plotting all 7,800 proteins in one 3D panel.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"leiden_structural_results"
DEST=OUT/"selected_community_3d"
DEST.mkdir(exist_ok=True)
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"svg.fonttype":"none","pdf.fonttype":42,"font.size":8})

# Communities pass the independent enrichment criteria and have a specific
# label covering >=50% of their members.  The choice is frozen here rather than
# made visually after plotting.
SELECTED=[9,12,14,18,20,21,22,23]
PALETTE=["#7B95C6","#49C2D9","#67A583","#A2C986","#F59C7C","#C85E62","#A67FC5","#FFC1A6"]

def main():
    d=pd.read_csv(OUT/"protein_leiden_structural_umap3d.tsv",sep="\t")
    labels=pd.read_csv(OUT/"leiden_macro_module_candidate_labels.tsv",sep="\t")
    selected=labels[labels.community.isin(SELECTED)].copy().set_index("community")
    x=np.load(OUT/"leiden_input_pca50.npy")
    pca=PCA(n_components=3,random_state=42).fit_transform(x)
    d[["PC1","PC2","PC3"]]=pca
    d=d[d.macro_community.isin(SELECTED)].copy()
    names={c: f"M{c}: {selected.loc[c,'candidate_label']}" for c in SELECTED}
    d["module_label"]=d.macro_community.map(names)
    d["color"]=d.macro_community.map(dict(zip(SELECTED,PALETTE)))
    d.to_csv(DEST/"selected_community_pca3d_source_data.tsv",sep="\t",index=False)
    evr=PCA(n_components=3,random_state=42).fit(x).explained_variance_ratio_*100

    fig=plt.figure(figsize=(7.2,5.5));ax=fig.add_subplot(projection="3d")
    for community,color in zip(SELECTED,PALETTE):
        g=d[d.macro_community.eq(community)]; pts=g[["PC1","PC2","PC3"]].to_numpy()
        # Render the complete module; selected modules are all small enough that
        # no point downsampling is needed.
        ax.scatter(pts[:,0],pts[:,1],pts[:,2],s=16,alpha=.78,c=color,edgecolors="white",linewidths=.18,label=f"{names[community]} (n={len(g)})",depthshade=True)
        if len(pts)>=4:
            try:
                hull=ConvexHull(pts)
                faces=[pts[s] for s in hull.simplices]
                mesh=Poly3DCollection(faces,facecolors=color,edgecolors=color,linewidths=.28,alpha=.055)
                ax.add_collection3d(mesh)
            except Exception:
                pass
    ax.set_xlabel(f"PC1 ({evr[0]:.1f}%)");ax.set_ylabel(f"PC2 ({evr[1]:.1f}%)");ax.set_zlabel(f"PC3 ({evr[2]:.1f}%)")
    ax.view_init(elev=22,azim=42);ax.legend(title="Selected Leiden communities",loc="center left",bbox_to_anchor=(1.03,.5),frameon=False,fontsize=7,markerscale=1.3)
    for ext in ("png","pdf","svg"):
        fig.savefig(DEST/f"selected_leiden_communities_pca3d.{ext}",dpi=600,bbox_inches="tight",facecolor="white")
    plt.close(fig)

    interactive=go.Figure()
    for community,color in zip(SELECTED,PALETTE):
        g=d[d.macro_community.eq(community)]; pts=g[["PC1","PC2","PC3"]].to_numpy()
        interactive.add_trace(go.Scatter3d(x=pts[:,0],y=pts[:,1],z=pts[:,2],mode="markers",name=f"{names[community]} (n={len(g)})",marker=dict(size=4,color=color,opacity=.82,line=dict(color="white",width=.3)),customdata=g[["uniprot_id","posthoc_membrane_role","fine_community"]],hovertemplate="%{customdata[0]}<br>role: %{customdata[1]}<br>fine community: %{customdata[2]}<extra></extra>"))
        if len(pts)>=4:
            try:
                hull=ConvexHull(pts)
                interactive.add_trace(go.Mesh3d(x=pts[:,0],y=pts[:,1],z=pts[:,2],i=hull.simplices[:,0],j=hull.simplices[:,1],k=hull.simplices[:,2],name=f"{names[community]} hull",color=color,opacity=.09,hoverinfo="skip",showlegend=False))
            except Exception:
                pass
    interactive.update_layout(template="simple_white",scene=dict(xaxis_title=f"PC1 ({evr[0]:.1f}%)",yaxis_title=f"PC2 ({evr[1]:.1f}%)",zaxis_title=f"PC3 ({evr[2]:.1f}%)"),legend_title_text="Selected community")
    interactive.write_html(DEST/"selected_leiden_communities_pca3d_interactive.html",include_plotlyjs=True)
    (DEST/"selection_manifest.json").write_text(json.dumps({"selection_rule":"independent post-hoc candidate label with specific term covering >=50% of module","selected_communities":SELECTED,"n_points":len(d),"pca_explained_variance_percent":evr.tolist()},indent=2),encoding="utf-8")

if __name__=="__main__":main()
