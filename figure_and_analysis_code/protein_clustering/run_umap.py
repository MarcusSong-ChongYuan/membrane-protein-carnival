import numpy as np,pandas as pd,umap
from sklearn.manifold import trustworthiness
from pipeline_lib import load_config,write_tsv,save_json
def main() -> None:
 cfg=load_config();out=cfg["output_dir"];x=np.load(out/"pca_embeddings.npy");meta=pd.read_csv(out/"protein_clusters.tsv",sep="\t",dtype=str);meta["cluster_id"]=meta.cluster_id.astype(int); scans=[]
 for nn in cfg["umap"]["n_neighbors_sensitivity"]:
  u=umap.UMAP(n_components=3,n_neighbors=nn,min_dist=cfg["umap"]["min_dist"],metric=cfg["umap"]["metric"],random_state=cfg["random_state"]).fit_transform(x);np.save(out/f"umap3d_neighbors_{nn}.npy",u.astype(np.float32));scans.append({"n_neighbors":nn,"n_components":3,"random_state":cfg["random_state"]})
 final=np.load(out/f"umap3d_neighbors_{cfg['umap']['final_n_neighbors']}.npy"); o=meta.copy();o[["UMAP1","UMAP2","UMAP3"]]=final;write_tsv(o,out/"protein_umap3d.tsv")
 u2=umap.UMAP(n_components=2,n_neighbors=cfg["umap"]["final_n_neighbors"],min_dist=cfg["umap"]["min_dist"],metric=cfg["umap"]["metric"],random_state=cfg["random_state"]).fit_transform(x);o2=meta.copy();o2[["UMAP1","UMAP2"]]=u2;write_tsv(o2,out/"protein_umap2d.tsv");save_json({"parameter_scan":scans,"final_n_neighbors":cfg["umap"]["final_n_neighbors"],"note":"UMAP is visualization-only; HDBSCAN ran on PCA embeddings."},out/"umap_metadata.json")
if __name__=="__main__":main()
