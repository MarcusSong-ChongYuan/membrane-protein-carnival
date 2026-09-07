import numpy as np,pandas as pd
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from pipeline_lib import load_config,write_tsv
def main() -> None:
 cfg=load_config();out=cfg["output_dir"];d=pd.read_csv(out/"protein_clusters.tsv",sep="\t",dtype=str);d.cluster_id=d.cluster_id.astype(int);known=d.target_class.ne("Unknown");ari=adjusted_rand_score(d.loc[known,"target_class"],d.loc[known,"cluster_id"]);nmi=normalized_mutual_info_score(d.loc[known,"target_class"],d.loc[known,"cluster_id"])
 rows=[]; comp=[]
 for c,g in d[d.cluster_id>=0].groupby("cluster_id"):
  x=g.target_class.value_counts();n=len(g);dominant=x.index[0];purity=x.iloc[0]/n;rows.append({"cluster_id":c,"cluster_size":n,"dominant_target_class":dominant,"dominant_class_fraction":purity,"cluster_purity":purity})
  for klass,count in x.items():comp.append({"cluster_id":c,"target_class":klass,"count":count,"proportion":count/n})
 write_tsv(pd.DataFrame(comp),out/"cluster_target_composition.tsv");write_tsv(pd.DataFrame(rows),out/"cluster_summary.tsv");write_tsv(pd.DataFrame([{"metric":"ARI_known_target_classes","value":ari},{"metric":"NMI_known_target_classes","value":nmi},{"metric":"n_known_target_classes","value":int(known.sum())}]),out/"clustering_metrics.tsv")
if __name__=="__main__":main()
