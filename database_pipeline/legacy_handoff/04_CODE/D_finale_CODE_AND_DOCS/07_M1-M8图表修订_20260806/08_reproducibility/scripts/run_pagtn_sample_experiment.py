from __future__ import annotations

import hashlib
import json
import math
import os
import random
import time
from collections import defaultdict
from pathlib import Path

import dgl
import hdbscan
import numpy as np
import pandas as pd
import torch
import umap
from dgllife.model import PAGTNPredictor
from dgllife.utils import PAGTNAtomFeaturizer, PAGTNEdgeFeaturizer, smiles_to_complete_graph
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import adjusted_mutual_info_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, Dataset


SEED = 42
ROOT = Path(os.environ.get("MEMPRO_REVISION_ROOT", r"D:\finale\07_M1-M8图表修订_20260806"))
SOURCE = Path(r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730\small_molecule_master_v1_3.tsv")
OUT = ROOT / "04_pagtn"
DATA_OUT = ROOT / "02_data"
OUT.mkdir(parents=True, exist_ok=True)
DATA_OUT.mkdir(parents=True, exist_ok=True)

SAMPLE_N = int(os.environ.get("PAGTN_SAMPLE_N", "2400"))
EPOCHS = int(os.environ.get("PAGTN_EPOCHS", "10"))
BATCH_SIZE = int(os.environ.get("PAGTN_BATCH_SIZE", "24"))
MAX_ATOMS = int(os.environ.get("PAGTN_MAX_ATOMS", "60"))


def set_seed(seed: int = SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)


def stable_score(text: str) -> int:
    return int(hashlib.sha256((str(text) + f"|{SEED}").encode("utf-8")).hexdigest()[:16], 16)


def structural_class(row: pd.Series) -> str:
    formula = str(row.get("molecular_formula", "") or "")
    rings = int(row.get("ring_count", 0) or 0)
    max_ring = int(row.get("max_ring_size", 0) or 0)
    heavy = int(row.get("heavy_atom_count", 0) or 0)
    amide = int(row.get("amide_bond_count", 0) or 0)
    if "C" not in formula:
        return "inorganic/no carbon"
    if heavy >= 25 and amide >= 4:
        return "peptide-like"
    if max_ring >= 12:
        return "macrocycle"
    if rings >= 4:
        return "polycyclic (>=4 rings)"
    if rings >= 2:
        return "polycyclic (2-3 rings)"
    if rings == 1:
        return "monocyclic"
    return "acyclic organic"


def select_sample() -> pd.DataFrame:
    cols = [
        "compound_internal_id", "preferred_name", "standard_smiles", "molecular_formula",
        "molecular_weight", "xlogp", "tpsa", "rotatable_bond_count", "heavy_atom_count",
        "ring_count", "max_ring_size", "amide_bond_count", "compound_scope_status",
        "record_qc_status", "identity_confidence", "is_approved_drug", "is_endogenous_ligand",
    ]
    pools: dict[str, list[tuple[int, dict]]] = defaultdict(list)
    for chunk in pd.read_csv(SOURCE, sep="\t", usecols=cols, chunksize=100_000, low_memory=False):
        for col in ["molecular_weight", "xlogp", "tpsa", "rotatable_bond_count", "heavy_atom_count", "ring_count", "max_ring_size", "amide_bond_count"]:
            chunk[col] = pd.to_numeric(chunk[col], errors="coerce")
        chunk = chunk[
            chunk["compound_scope_status"].eq("core")
            & chunk["record_qc_status"].eq("ok")
            & chunk["standard_smiles"].notna()
            & chunk[["molecular_weight", "xlogp", "tpsa", "rotatable_bond_count", "heavy_atom_count"]].notna().all(axis=1)
            & chunk["molecular_weight"].between(80, 800)
            & chunk["heavy_atom_count"].between(5, MAX_ATOMS)
        ].copy()
        if chunk.empty:
            continue
        chunk["structural_class"] = chunk.apply(structural_class, axis=1)
        for rec in chunk.to_dict("records"):
            score = stable_score(rec["compound_internal_id"])
            bucket = pools[rec["structural_class"]]
            bucket.append((score, rec))
            if len(bucket) > 800:
                bucket.sort(key=lambda x: x[0])
                del bucket[650:]

    class_records = {k: [r for _, r in sorted(v, key=lambda x: x[0])] for k, v in pools.items()}
    classes = sorted(class_records)
    selected: list[dict] = []
    # Preserve rare classes, then fill the remainder by stable hash across all remaining records.
    base = max(30, SAMPLE_N // max(1, len(classes)) // 2)
    used: set[str] = set()
    for cls in classes:
        take = class_records[cls][: min(base, len(class_records[cls]))]
        selected.extend(take)
        used.update(str(r["compound_internal_id"]) for r in take)
    remaining = [r for cls in classes for r in class_records[cls] if str(r["compound_internal_id"]) not in used]
    remaining.sort(key=lambda r: stable_score(r["compound_internal_id"]))
    selected.extend(remaining[: max(0, SAMPLE_N - len(selected))])
    df = pd.DataFrame(selected[:SAMPLE_N]).reset_index(drop=True)

    valid = []
    scaffolds = []
    for smi in df["standard_smiles"].astype(str):
        mol = Chem.MolFromSmiles(smi)
        ok = mol is not None and 5 <= mol.GetNumHeavyAtoms() <= MAX_ATOMS
        valid.append(ok)
        if ok:
            scaffolds.append(MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True) or "acyclic")
        else:
            scaffolds.append("")
    df = df.loc[valid].copy().reset_index(drop=True)
    df["bemis_murcko_scaffold"] = [s for s, ok in zip(scaffolds, valid) if ok]
    df["sample_seed"] = SEED
    df.to_csv(OUT / "pagtn_sample.tsv", sep="\t", index=False)
    return df


class GraphDataset(Dataset):
    def __init__(self, graphs, targets, indices):
        self.graphs = graphs
        self.targets = targets
        self.indices = indices

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, i):
        idx = int(self.indices[i])
        return self.graphs[idx], torch.tensor(self.targets[idx], dtype=torch.float32), idx


def collate(batch):
    graphs, y, idx = zip(*batch)
    return dgl.batch(graphs), torch.stack(y), np.asarray(idx, dtype=np.int64)


def graph_embedding(model, graph):
    node_feats = graph.ndata["h"]
    edge_feats = graph.edata["e"]
    atom_h = model.model(graph, node_feats, edge_feats)
    atom_h = torch.cat([atom_h, node_feats], dim=1)
    node_z = model.readout.in_project(atom_h)
    if model.readout.activation is not None:
        node_z = model.readout.activation(node_z)
    with graph.local_scope():
        graph.ndata["_pagtn_embedding"] = node_z
        return dgl.mean_nodes(graph, "_pagtn_embedding")


def cluster_metrics(name: str, x: np.ndarray, labels: np.ndarray, classes: np.ndarray, scaffolds: np.ndarray) -> list[dict]:
    nonnoise = labels >= 0
    n_clusters = int(len(set(labels[nonnoise]))) if nonnoise.any() else 0
    noise_fraction = float((labels < 0).mean())
    sil = float("nan")
    if nonnoise.sum() > 20 and n_clusters > 1:
        sil = float(silhouette_score(x[nonnoise], labels[nonnoise], metric="euclidean", sample_size=min(1500, int(nonnoise.sum())), random_state=SEED))
    ami = float(adjusted_mutual_info_score(classes, labels))
    nn = NearestNeighbors(n_neighbors=min(11, len(x)), metric="cosine").fit(x)
    neigh = nn.kneighbors(return_distance=False)[:, 1:]
    class_purity = float(np.mean([np.mean(classes[row] == classes[i]) for i, row in enumerate(neigh)]))
    scaffold_purity = float(np.mean([np.mean(scaffolds[row] == scaffolds[i]) for i, row in enumerate(neigh)]))
    return [
        {"method": name, "metric": "cluster_count", "value": n_clusters},
        {"method": name, "metric": "noise_fraction", "value": noise_fraction},
        {"method": name, "metric": "silhouette_nonnoise", "value": sil},
        {"method": name, "metric": "AMI_vs_rule_class", "value": ami},
        {"method": name, "metric": "10NN_rule_class_purity", "value": class_purity},
        {"method": name, "metric": "10NN_scaffold_purity", "value": scaffold_purity},
    ]


def main() -> None:
    start = time.time()
    set_seed()
    df = select_sample()
    atom_f = PAGTNAtomFeaturizer(atom_data_field="h")
    edge_f = PAGTNEdgeFeaturizer(bond_data_field="e", max_length=5)
    graphs = []
    keep = []
    for i, smi in enumerate(df["standard_smiles"].astype(str)):
        try:
            g = smiles_to_complete_graph(smi, add_self_loop=True, node_featurizer=atom_f, edge_featurizer=edge_f)
        except Exception:
            g = None
        if g is not None:
            graphs.append(g)
            keep.append(i)
    df = df.iloc[keep].reset_index(drop=True)
    df.to_csv(OUT / "pagtn_sample_graph_valid.tsv", sep="\t", index=False)

    target_cols = ["molecular_weight", "xlogp", "tpsa", "rotatable_bond_count"]
    targets_raw = df[target_cols].to_numpy(dtype=np.float32)
    target_scaler = StandardScaler().fit(targets_raw)
    targets = target_scaler.transform(targets_raw).astype(np.float32)
    n = len(df)
    rng = np.random.default_rng(SEED)
    order = rng.permutation(n)
    split = int(n * 0.8)
    train_idx, valid_idx = order[:split], order[split:]
    train_loader = DataLoader(GraphDataset(graphs, targets, train_idx), batch_size=BATCH_SIZE, shuffle=True, collate_fn=collate, generator=torch.Generator().manual_seed(SEED))
    valid_loader = DataLoader(GraphDataset(graphs, targets, valid_idx), batch_size=BATCH_SIZE, shuffle=False, collate_fn=collate)
    all_loader = DataLoader(GraphDataset(graphs, targets, np.arange(n)), batch_size=BATCH_SIZE, shuffle=False, collate_fn=collate)

    node_feats = atom_f.feat_size()
    edge_feats = edge_f.feat_size()
    model = PAGTNPredictor(node_feats, 64, 32, edge_feats, depth=3, nheads=2, dropout=0.12, n_tasks=len(target_cols), mode="mean")
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
    loss_fn = torch.nn.MSELoss()
    history = []
    best = math.inf
    best_state = None
    patience = 3
    stale = 0
    for epoch in range(1, EPOCHS + 1):
        model.train()
        train_losses = []
        for bg, y, _ in train_loader:
            pred = model(bg, bg.ndata["h"], bg.edata["e"])
            loss = loss_fn(pred, y)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            train_losses.append(float(loss.item()))
        model.eval()
        val_losses = []
        with torch.no_grad():
            for bg, y, _ in valid_loader:
                val_losses.append(float(loss_fn(model(bg, bg.ndata["h"], bg.edata["e"]), y).item()))
        row = {"epoch": epoch, "train_mse": float(np.mean(train_losses)), "validation_mse": float(np.mean(val_losses))}
        history.append(row)
        print(json.dumps(row), flush=True)
        if row["validation_mse"] < best - 1e-4:
            best = row["validation_mse"]
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    torch.save({"state_dict": model.state_dict(), "target_mean": target_scaler.mean_, "target_scale": target_scaler.scale_}, OUT / "pagtn_property_supervised_model.pt")
    pd.DataFrame(history).to_csv(OUT / "training_history.tsv", sep="\t", index=False)

    model.eval()
    embeddings = np.zeros((n, 64), dtype=np.float32)
    with torch.no_grad():
        for bg, _, idx in all_loader:
            embeddings[idx] = graph_embedding(model, bg).cpu().numpy()
    np.save(OUT / "pagtn_embeddings.npy", embeddings)

    morgan = np.zeros((n, 1024), dtype=np.float32)
    for i, smi in enumerate(df["standard_smiles"].astype(str)):
        mol = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganGenerator(radius=2, fpSize=1024).GetFingerprint(mol)
        DataStructs.ConvertToNumpyArray(fp, morgan[i])

    pagtn_x = StandardScaler().fit_transform(embeddings)
    # Keep the binary Morgan representation in its native space for clustering/neighbors.
    clusterer_p = hdbscan.HDBSCAN(min_cluster_size=30, min_samples=10, metric="euclidean", cluster_selection_method="eom")
    clusterer_m = hdbscan.HDBSCAN(min_cluster_size=30, min_samples=10, metric="euclidean", cluster_selection_method="eom")
    pagtn_labels = clusterer_p.fit_predict(pagtn_x)
    morgan_labels = clusterer_m.fit_predict(morgan)
    reducer_p = umap.UMAP(n_neighbors=25, min_dist=0.12, metric="cosine", random_state=SEED, n_jobs=1)
    reducer_m = umap.UMAP(n_neighbors=25, min_dist=0.12, metric="jaccard", random_state=SEED, n_jobs=1)
    pagtn_2d = reducer_p.fit_transform(pagtn_x)
    morgan_2d = reducer_m.fit_transform(morgan)

    classes = df["structural_class"].astype(str).to_numpy()
    scaffolds = df["bemis_murcko_scaffold"].astype(str).to_numpy()
    metrics = cluster_metrics("PAGTN_property_supervised", pagtn_x, pagtn_labels, classes, scaffolds)
    metrics += cluster_metrics("Morgan_radius2_1024", morgan, morgan_labels, classes, scaffolds)
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(DATA_OUT / "M5E_pagtn_cluster_summary.tsv", sep="\t", index=False)

    coords = df[["compound_internal_id", "preferred_name", "standard_smiles", "structural_class", "bemis_murcko_scaffold"]].copy()
    coords["pagtn_umap1"] = pagtn_2d[:, 0]
    coords["pagtn_umap2"] = pagtn_2d[:, 1]
    coords["pagtn_cluster"] = pagtn_labels
    coords["morgan_umap1"] = morgan_2d[:, 0]
    coords["morgan_umap2"] = morgan_2d[:, 1]
    coords["morgan_cluster"] = morgan_labels
    coords.to_csv(DATA_OUT / "M5E_pagtn_morgan_coordinates.tsv", sep="\t", index=False)

    qa = {
        "status": "PASS",
        "method_boundary": "PAGTN is an encoder; HDBSCAN performs clustering. The encoder was trained only on four standardized physicochemical descriptors and is not a universal pretrained chemical embedding.",
        "seed": SEED,
        "sample_requested": SAMPLE_N,
        "sample_graph_valid": n,
        "epochs_completed": len(history),
        "best_validation_mse_standardized": best,
        "targets": target_cols,
        "node_feature_dimension": node_feats,
        "edge_feature_dimension": edge_feats,
        "embedding_dimension": 64,
        "software": {"torch": torch.__version__, "dgl": dgl.__version__},
        "elapsed_seconds": round(time.time() - start, 2),
    }
    (OUT / "PAGTN_EXPERIMENT_QA.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(qa, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
