#!/usr/bin/env python3
"""All-protein Leiden map from sequence, topology and structural-family data.

Primary membrane role, GO molecular function/BP, pathway and ligand evidence are
not used as input. They are retained for post-hoc community interpretation.
"""
from __future__ import annotations

import json
from pathlib import Path

import igraph as ig
import leidenalg
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import umap
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.feature_extraction import DictVectorizer
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent
INPUT, BASE, OUT = ROOT / "inputs", ROOT / "results", ROOT / "leiden_structural_results"
OUT.mkdir(exist_ok=True)
mpl.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","DejaVu Sans"], "svg.fonttype":"none", "pdf.fonttype":42, "font.size":8})


def read(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def block_sparse(rows: list[dict[str,float]], n_components: int) -> tuple[np.ndarray, int, int]:
    x = DictVectorizer(sparse=True).fit_transform(rows)
    if x.shape[1] < 2:
        return np.zeros((len(rows),1),dtype=np.float32), int(x.shape[1]), 1
    # Structural families appearing once are kept in raw metadata but not used to
    # define global graph geometry: they cannot support reproducible communities.
    support = np.asarray((x > 0).sum(axis=0)).ravel()
    x = x[:, support >= 5]
    n = max(2, min(n_components, x.shape[0]-1, x.shape[1]-1))
    z = TruncatedSVD(n_components=n, random_state=42, n_iter=7).fit_transform(x)
    return z.astype(np.float32), int(x.shape[1]), int(n)


def equal_weight(z: np.ndarray) -> np.ndarray:
    z = StandardScaler().fit_transform(z)
    return (z / np.sqrt(max(1,z.shape[1]))).astype(np.float32)


def graph_from_knn(z: np.ndarray, k: int=25) -> ig.Graph:
    nn = NearestNeighbors(n_neighbors=k+1, metric="cosine", n_jobs=-1).fit(z)
    d, ind = nn.kneighbors(z)
    nonzero = d[:,1:][d[:,1:] > 0]
    scale = float(np.median(nonzero)) if len(nonzero) else 1.0
    edge_w: dict[tuple[int,int],float] = {}
    for i in range(len(z)):
        for j, dist in zip(ind[i,1:], d[i,1:]):
            a,b = sorted((int(i),int(j)))
            weight = float(np.exp(-((float(dist)/scale)**2)))
            edge_w[(a,b)] = max(edge_w.get((a,b),0.0),weight)
    g = ig.Graph(n=len(z), edges=list(edge_w))
    g.es["weight"] = list(edge_w.values())
    return g


def consensus_partition(g: ig.Graph, resolution: float, seeds=range(5)) -> tuple[leidenalg.VertexPartition, float]:
    parts=[]
    for seed in seeds:
        parts.append(leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
            weights="weight", resolution_parameter=resolution, seed=int(seed)))
    # Choose the partition with median quality, then report stability separately.
    quality=np.array([p.quality() for p in parts]); best=parts[int(np.argsort(quality)[len(parts)//2])]
    labels=[np.array(p.membership) for p in parts]
    from sklearn.metrics import adjusted_rand_score
    ari=np.mean([adjusted_rand_score(labels[0],x) for x in labels[1:]]) if len(labels)>1 else 1.0
    return best, float(ari)


def choose_resolution(g: ig.Graph, target_low: int, target_high: int, name: str) -> tuple[leidenalg.VertexPartition, pd.DataFrame, float]:
    rows=[]; saved={}
    # Include genuinely coarse resolutions for the macro map.  The previous
    # scan began at 0.15, which was already too fine for this graph.
    for r in (0.005,0.01,0.02,0.03,0.05,0.075,0.10,0.15,0.25,0.35,0.5,0.7,0.9,1.2,1.6,2.0,2.5):
        p, stability = consensus_partition(g,r)
        n=len(p)
        # Quality/stability decide inside the requested descriptive scale; roles
        # and other held-out biology are never used.
        in_range=target_low<=n<=target_high
        distance=0 if in_range else min(abs(n-target_low),abs(n-target_high))
        rows.append({"map_level":name,"resolution":r,"n_communities":n,"quality":p.quality(),"seed_stability_ARI":stability,"in_requested_scale":in_range,"distance_to_scale":distance})
        saved[r]=p
    scan=pd.DataFrame(rows)
    viable=scan[scan.in_requested_scale]
    pick=(viable if len(viable) else scan).sort_values(["seed_stability_ARI","quality","distance_to_scale"],ascending=[False,False,True]).iloc[0]
    return saved[float(pick.resolution)],scan,float(pick.resolution)


def enrich_label(meta: pd.DataFrame, annotation: pd.DataFrame, cluster_col: str, axis: str) -> pd.DataFrame:
    """Held-out label enrichment: observed vs background fraction, no auto naming."""
    a=annotation[annotation.classification_axis.eq(axis)][["target_uniprot_id","classification_label"]].drop_duplicates()
    x=meta[["uniprot_id",cluster_col]].merge(a,left_on="uniprot_id",right_on="target_uniprot_id",how="inner")
    if x.empty: return pd.DataFrame()
    bg=x.classification_label.value_counts(); rows=[]
    for (c,label),n in x.groupby([cluster_col,"classification_label"]).size().items():
        cs=(meta[cluster_col]==c).sum(); total=len(meta); b=bg[label]
        if n < 5 or b < 5: continue
        expected=cs*b/total
        rows.append({"community":c,"held_out_axis":axis,"label":label,"observed":int(n),"expected":expected,"fold_enrichment":n/expected if expected else np.nan})
    return pd.DataFrame(rows)


def main() -> None:
    master=read(INPUT/"protein_master_v72.tsv.gz",usecols=["canonical_uniprot_accession","sequence_length","membrane_class_v7","primary_membrane_mode_v7","secondary_membrane_modes_v7"]).rename(columns={"canonical_uniprot_accession":"uniprot_id"})
    ids=read(BASE/"protein_embedding_ids.tsv")["uniprot_id"].astype(str).tolist()
    master=master.drop_duplicates("uniprot_id").set_index("uniprot_id").loc[ids].reset_index()
    role=read(INPUT/"protein_membrane_role_FORMAL.tsv",dtype=str).rename(columns={"canonical_uniprot_accession":"uniprot_id","formal_primary_membrane_role":"posthoc_membrane_role"})[["uniprot_id","posthoc_membrane_role"]]
    master=master.merge(role,on="uniprot_id",how="left").fillna({"posthoc_membrane_role":"Unknown"})
    xseq=np.load(BASE/"protein_embeddings.npy")
    sequence=PCA(n_components=50,random_state=42).fit_transform(StandardScaler().fit_transform(xseq))

    topo=[]
    for r in master.itertuples(index=False):
        d={f"class:{r.membrane_class_v7}":1.0,f"mode:{r.primary_membrane_mode_v7}":1.0,"length_log1p":float(np.log1p(r.sequence_length))}
        for m in str(r.secondary_membrane_modes_v7 or "").split(";"):
            if m and m.lower()!="nan": d[f"secondary:{m}"]=1.0
        topo.append(d)
    topology,topo_features,topo_dims=block_sparse(topo,12)

    ann=read(INPUT/"protein_cross_classification_v72.tsv.gz",usecols=["target_uniprot_id","classification_axis","classification_label"],dtype=str)
    structural=ann[(ann.classification_axis=="structural_family") & ann.target_uniprot_id.isin(ids)].copy()
    structural=structural[~structural.classification_label.str.lower().isin(["unclassified","unknown","nan"])]
    sd={u:{} for u in ids}
    for r in structural.itertuples(index=False): sd[r.target_uniprot_id][f"family:{r.classification_label}"]=1.0
    family,family_features,family_dims=block_sparse([sd[u] for u in ids],30)

    feature=np.hstack([equal_weight(sequence),equal_weight(topology),equal_weight(family)]).astype(np.float32)
    z=PCA(n_components=50,random_state=42).fit_transform(feature).astype(np.float32)
    np.save(OUT/"leiden_input_pca50.npy",z)
    g=graph_from_knn(z,k=25)
    macro,macro_scan,macro_r=choose_resolution(g,8,14,"macro")
    fine,fine_scan,fine_r=choose_resolution(g,30,60,"fine")
    scan=pd.concat([macro_scan,fine_scan],ignore_index=True); scan.to_csv(OUT/"leiden_resolution_scan.tsv",sep="\t",index=False)
    master["macro_community"]=macro.membership; master["fine_community"]=fine.membership

    coords=umap.UMAP(n_components=3,n_neighbors=30,min_dist=.12,metric="cosine",random_state=42).fit_transform(z)
    master[["UMAP1","UMAP2","UMAP3"]]=coords
    master.to_csv(OUT/"protein_leiden_structural_umap3d.tsv",sep="\t",index=False)
    enrich=pd.concat([enrich_label(master,ann,"macro_community",axis) for axis in ["molecular_function","biological_process","specialist_classification"]],ignore_index=True)
    enrich.sort_values(["community","fold_enrichment","observed"],ascending=[True,False,False]).to_csv(OUT/"macro_community_heldout_enrichment.tsv",sep="\t",index=False)

    master["macro_display"]=master.macro_community.map(lambda x:f"Module {x}")
    fig=px.scatter_3d(master,x="UMAP1",y="UMAP2",z="UMAP3",color="macro_display",hover_data=["uniprot_id","posthoc_membrane_role","macro_community","fine_community","sequence_length"],opacity=.62)
    fig.update_layout(template="simple_white",legend_title_text="Leiden macro-module")
    fig.write_html(OUT/"protein_leiden_structural_interactive.html",include_plotlyjs=True)
    # Static companion: largest 10 macro modules are readable; all labels stay in TSV/HTML.
    size=master.macro_display.value_counts(); top=size.head(10).index.tolist()
    master["static_group"]=master.macro_display.where(master.macro_display.isin(top),"Other macro-modules")
    colors=dict(zip(top,["#4DAF4A","#4BA3D8","#F29E4C","#A67FC5","#E46C7A","#67A583","#7B95C6","#C85E62","#A2C986","#FDed95"])); colors["Other macro-modules"]="#C5C5C5"
    f=plt.figure(figsize=(7,5.3)); ax=f.add_subplot(projection="3d")
    for lab in ["Other macro-modules"]+top:
        q=master[master.static_group==lab]
        ax.scatter(q.UMAP1,q.UMAP2,q.UMAP3,s=4 if lab=="Other macro-modules" else 6,alpha=.22 if lab=="Other macro-modules" else .62,c=colors[lab],edgecolors="none",label=f"{lab} (n={len(q):,})",rasterized=True)
    ax.set_xlabel("UMAP1");ax.set_ylabel("UMAP2");ax.set_zlabel("UMAP3");ax.view_init(elev=22,azim=42);ax.legend(title="Leiden macro-module",loc="center left",bbox_to_anchor=(1.02,.5),frameon=False,fontsize=7,markerscale=2)
    for ext in ("png","pdf","svg"):f.savefig(OUT/f"protein_leiden_structural_top10.{ext}",dpi=500,bbox_inches="tight",facecolor="white")
    plt.close(f)
    summary={"n_proteins":len(master),"graph":"25-nearest-neighbour cosine graph","input_blocks":{"ESM_sequence_PCA":50,"membrane_topology":topo_dims,"structural_family_SVD":family_dims},"structural_features_with_support_ge5":family_features,"primary_role_used_as_input":False,"macro_resolution":macro_r,"macro_communities":len(macro),"fine_resolution":fine_r,"fine_communities":len(fine),"note":"All proteins receive a Leiden community. Held-out functional/role annotations are for interpretation, not inputs.","random_state":42}
    (OUT/"leiden_structural_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2))

if __name__=="__main__": main()
