import numpy as np,pandas as pd,hdbscan
from sklearn.metrics import silhouette_score
from pipeline_lib import load_config,write_tsv,save_json
def main() -> None:
 cfg=load_config();out=cfg["output_dir"];x=np.load(out/"pca_embeddings.npy");records=[];models={}
 for mcs in cfg["clustering"]["min_cluster_size"]:
  for ms in cfg["clustering"]["min_samples"]:
   model=hdbscan.HDBSCAN(min_cluster_size=mcs,min_samples=ms,metric="euclidean",prediction_data=True).fit(x);lab=model.labels_; valid=lab>=0;ncl=len(set(lab[valid]));noise=float((~valid).mean());sil=np.nan
   if ncl>=2 and valid.sum()>ncl: sil=float(silhouette_score(x[valid],lab[valid]))
   rec={"min_cluster_size":mcs,"min_samples":"None" if ms is None else ms,"n_clusters":ncl,"n_noise":int((~valid).sum()),"noise_fraction":noise,"silhouette_non_noise":sil,"cluster_persistence_mean":float(np.mean(model.cluster_persistence_)) if len(model.cluster_persistence_) else np.nan}
   records.append(rec);models[(mcs,ms)]=model
 scan=pd.DataFrame(records);write_tsv(scan,out/"hdbscan_parameter_scan.tsv")
 # Selection uses unsupervised diagnostics only: avoid degenerate one-cluster/noise outcomes, then optimize persistence and silhouette with modest noise penalty.
 viable=scan[(scan.n_clusters>=2)&(scan.noise_fraction<=0.70)].copy(); viable["rank_score"]=viable.cluster_persistence_mean.fillna(0)+viable.silhouette_non_noise.fillna(-1)-0.15*viable.noise_fraction
 best=(viable.sort_values(["rank_score","min_cluster_size"],ascending=[False,True]).iloc[0] if len(viable) else scan.sort_values(["noise_fraction","n_clusters"],ascending=[True,False]).iloc[0]); ms=None if str(best.min_samples)=="None" else int(best.min_samples); model=models[(int(best.min_cluster_size),ms)]
 ids=pd.read_csv(out/"protein_unique_table.tsv",sep="\t",dtype=str); ids["cluster_id"]=model.labels_;write_tsv(ids[["uniprot_id","cluster_id","target_class","sequence_length"]],out/"protein_clusters.tsv"); save_json({"selection_rule":"unsupervised persistence + silhouette - noise penalty; target class not used","selected":best.to_dict()},out/"selected_hdbscan_model.json")
if __name__=="__main__":main()
