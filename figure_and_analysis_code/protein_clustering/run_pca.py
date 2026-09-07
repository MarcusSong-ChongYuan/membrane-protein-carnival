import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib as mpl; import matplotlib.pyplot as plt
from pipeline_lib import load_config, write_tsv

def main() -> None:
    cfg=load_config();out=cfg["output_dir"]; x=np.load(out/"protein_embeddings.npy"); n=min(cfg["pca"]["max_components"],x.shape[0]-1,x.shape[1])
    z=StandardScaler().fit_transform(x); p=PCA(n_components=n,random_state=cfg["random_state"]); pc=p.fit_transform(z)
    np.save(out/"pca_embeddings.npy",pc.astype(np.float32)); write_tsv(pd.DataFrame({"PC":range(1,n+1),"individual_explained_variance":p.explained_variance_ratio_,"cumulative_explained_variance":np.cumsum(p.explained_variance_ratio_)}),out/"pca_explained_variance.tsv")
    mpl.rcParams.update({"font.family":"Arial","font.size":8,"svg.fonttype":"none","pdf.fonttype":42}); fig,ax=plt.subplots(figsize=(6.2,3.5)); k=np.arange(1,n+1);ax.bar(k,p.explained_variance_ratio_*100,color="#7b95c6",label="Individual"); ax.set(xlabel="Principal component",ylabel="Explained variance (%)"); ax2=ax.twinx();ax2.plot(k,np.cumsum(p.explained_variance_ratio_)*100,color="#c85e62",lw=1.5,label="Cumulative");ax2.set_ylabel("Cumulative variance (%)");fig.tight_layout()
    for ext in ("png","pdf","svg"):fig.savefig(out/f"pca_explained_variance.{ext}",dpi=600,bbox_inches="tight")
if __name__=="__main__": main()
