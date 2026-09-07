"""Generate reproducible protein-level ESM-2 embeddings with restart-safe caching."""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import pandas as pd
from pipeline_lib import load_config, save_json, write_tsv

def chunks(seq: str, maximum: int):
    # Non-overlapping chunks avoid silently truncating long membrane proteins.
    return [seq[i:i+maximum] for i in range(0, len(seq), maximum)]

def main() -> None:
    cfg = load_config(); out = cfg["output_dir"]
    table = pd.read_csv(out / "protein_unique_table.tsv", sep="\t", dtype=str)
    import torch, esm
    torch.manual_seed(cfg["random_state"]); np.random.seed(cfg["random_state"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, alphabet = esm.pretrained.load_model_and_alphabet(cfg["embedding"]["model_name"])
    model.eval().to(device)
    if device.type == "cuda": model.half()
    converter = alphabet.get_batch_converter(); layer = int(cfg["embedding"]["representation_layer"])
    cache = out / "embedding_chunks"; cache.mkdir(exist_ok=True)
    vectors=[]; issues=[]; maximum=int(cfg["embedding"]["max_residues_per_chunk"])
    for idx, row in table.iterrows():
        uid, seq = row.uniprot_id, row.sequence
        cp = cache / f"{uid}.npy"
        if cp.exists():
            vec=np.load(cp)
        else:
            segment_vectors=[]; weights=[]
            for part in chunks(seq, maximum):
                _, _, tokens = converter([(uid, part)])
                tokens=tokens.to(device)
                with torch.no_grad():
                    rep=model(tokens, repr_layers=[layer], return_contacts=False)["representations"][layer][0, 1:len(part)+1]
                v=rep.float().mean(dim=0).cpu().numpy()
                segment_vectors.append(v); weights.append(len(part))
            vec=np.average(np.vstack(segment_vectors), axis=0, weights=np.array(weights))
            np.save(cp, vec.astype(np.float32))
        if not np.isfinite(vec).all(): issues.append({"uniprot_id":uid,"problem":"non_finite_embedding"})
        vectors.append(vec)
        if (idx+1) % 25 == 0: print(f"embedded {idx+1}/{len(table)}", flush=True)
    arr=np.vstack(vectors).astype(np.float32)
    np.save(out / "protein_embeddings.npy", arr)
    write_tsv(table[["uniprot_id"]], out / "protein_embedding_ids.tsv")
    write_tsv(pd.DataFrame(issues), out / "embedding_problematic_proteins.tsv")
    save_json({"model_name":cfg["embedding"]["model_name"],"embedding_dimension":int(arr.shape[1]),
               "representation_layer":layer,"pooling_method":cfg["embedding"]["pooling_method"],
               "long_protein_handling":"length-weighted mean of non-overlapping <=1022-residue ESM segments",
               "device":str(device),"n_proteins":int(arr.shape[0])}, out / "method_metadata.json")
    print(arr.shape)

if __name__ == "__main__": main()
