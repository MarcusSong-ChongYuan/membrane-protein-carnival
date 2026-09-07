from pathlib import Path
import numpy as np, pandas as pd
from pipeline_lib import load_config, write_tsv

def main() -> None:
    cfg=load_config(); out=cfg["output_dir"]; x=np.load(out/"protein_embeddings.npy"); ids=pd.read_csv(out/"protein_embedding_ids.tsv",sep="\t")
    duplicate=[]
    _, first, counts=np.unique(x,axis=0,return_index=True,return_counts=True)
    for i,c in zip(first,counts):
        if c>1: duplicate.append(ids.iloc[i].uniprot_id)
    bad=np.where(~np.isfinite(x).all(axis=1))[0]
    problems=pd.DataFrame({"uniprot_id":ids.iloc[bad].uniprot_id,"problem":"non_finite_embedding"})
    if duplicate: problems=pd.concat([problems,pd.DataFrame({"uniprot_id":duplicate,"problem":"duplicate_embedding_vector"})])
    write_tsv(problems,out/"embedding_problematic_proteins.tsv")
    text=(f"N proteins\t{x.shape[0]}\nraw embedding dimensions\t{x.shape[1]}\nNaN present\t{bool(np.isnan(x).any())}\n"
          f"Inf present\t{bool(np.isinf(x).any())}\nduplicate vectors\t{int(sum(counts>1))}\nzero-variance dimensions\t{int((x.var(axis=0)==0).sum())}\n")
    (out/"embedding_qc.txt").write_text(text); print(text)
if __name__=="__main__": main()
