from __future__ import annotations

import hashlib
import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import umap
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina


ROOT = Path(os.environ["MEMPRO_FIGURE_OUTPUT_ROOT"])
ANALYSIS = ROOT / "02_analysis_data"
SD = ROOT / "03_source_data"
QA = ROOT / "08_QA"
SEED = 42
RNG = np.random.default_rng(SEED)


def stable_hash(value: str) -> int:
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:16], 16)


def fingerprint(mol: Chem.Mol):
    return AllChem.GetMorganGenerator(radius=2, fpSize=1024).GetFingerprint(mol)


def bit_array(fp) -> np.ndarray:
    arr = np.zeros((1024,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def main() -> None:
    chunks = []
    for chunk in pd.read_csv(ANALYSIS / "compound_features.tsv.gz", sep="\t", compression="gzip", dtype=str, chunksize=100_000):
        take = chunk[chunk["interaction_linked"].eq("1")].copy()
        if not take.empty:
            chunks.append(take)
    compounds = pd.concat(chunks, ignore_index=True)

    scaffold_rows = []
    parse_failures = 0
    for row in compounds[["compound_internal_id", "standard_smiles"]].itertuples(index=False):
        mol = Chem.MolFromSmiles(row.standard_smiles)
        if mol is None:
            parse_failures += 1
            scaffold = ""
        else:
            scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True) if scaffold_mol.GetNumAtoms() else "ACYCLIC_NO_MURCKO_SCAFFOLD"
        scaffold_rows.append({"compound_internal_id": row.compound_internal_id, "murcko_scaffold_smiles": scaffold})
    scaffold_df = pd.DataFrame(scaffold_rows)
    scaffold_df["scaffold_id"] = scaffold_df["murcko_scaffold_smiles"].map(lambda x: "SCF-" + hashlib.sha1(x.encode()).hexdigest()[:12] if x else "PARSE_FAILED")
    scaffold_df.to_csv(ANALYSIS / "compound_scaffold_mapping.tsv.gz", sep="\t", index=False, compression="gzip")

    scaffold_counts = scaffold_df.groupby(["scaffold_id", "murcko_scaffold_smiles"]).size().reset_index(name="compound_count").sort_values("compound_count", ascending=False)
    scaffold_counts["rank"] = np.arange(1, len(scaffold_counts) + 1)
    scaffold_counts.to_csv(SD / "F3A_scaffold_rank_abundance.tsv", sep="\t", index=False)
    scaffold_counts.head(30).to_csv(SD / "F3A_top_scaffolds.tsv", sep="\t", index=False)

    # Deterministic, status-enriched sample for the supplementary chemical atlas.
    priority_mask = (
        compounds["is_approved_drug"].eq("1")
        | compounds["is_clinical_candidate"].eq("1")
        | compounds["is_endogenous_ligand"].eq("1")
        | compounds["is_natural_product"].eq("1")
        | compounds["is_chemical_probe"].eq("1")
    )
    priority = compounds[priority_mask]
    remaining = compounds[~priority_mask].copy()
    remaining["_hash"] = remaining["compound_internal_id"].map(stable_hash)
    remaining = remaining.nsmallest(max(0, 12_000 - len(priority)), "_hash")
    sample = pd.concat([priority, remaining], ignore_index=True).drop_duplicates("compound_internal_id").head(12_000)

    valid_records = []
    fps = []
    for row in sample.itertuples(index=False):
        mol = Chem.MolFromSmiles(row.standard_smiles)
        if mol is None:
            continue
        fp = fingerprint(mol)
        fps.append(fp)
        status = "Other interaction-linked"
        for label, field in [
            ("Approved drug", "is_approved_drug"),
            ("Clinical candidate", "is_clinical_candidate"),
            ("Endogenous ligand", "is_endogenous_ligand"),
            ("Natural product", "is_natural_product"),
            ("Chemical probe", "is_chemical_probe"),
        ]:
            if getattr(row, field) == "1":
                status = label
                break
        valid_records.append({
            "compound_internal_id": row.compound_internal_id,
            "biological_status_display": status,
            "broad_structure_class": row.broad_structure_class_recomputed,
        })
    matrix = np.vstack([bit_array(fp) for fp in fps])
    reducer = umap.UMAP(n_neighbors=30, min_dist=0.18, n_components=2, metric="jaccard", random_state=SEED, low_memory=True)
    coords = reducer.fit_transform(matrix)
    embedding = pd.DataFrame(valid_records)
    embedding["UMAP1"] = coords[:, 0]
    embedding["UMAP2"] = coords[:, 1]
    embedding["sampling_note"] = "all priority-status compounds plus deterministic hash sample of other interaction-linked compounds"
    embedding.to_csv(SD / "S3_chemical_umap_deterministic_sample.tsv", sep="\t", index=False)

    # Butina is run on a bounded deterministic subset because its distance matrix
    # is quadratic. This benchmark tests neighbourhood structure; it is not used
    # as a release classification.
    butina_n = min(4_000, len(fps))
    bidx = np.argsort([stable_hash(x["compound_internal_id"]) for x in valid_records])[:butina_n]
    bfps = [fps[i] for i in bidx]
    distances = []
    for i in range(1, len(bfps)):
        sims = DataStructs.BulkTanimotoSimilarity(bfps[i], bfps[:i])
        distances.extend(1 - x for x in sims)
    clusters = Butina.ClusterData(distances, len(bfps), 0.35, isDistData=True)
    cluster_label = {}
    for cid, members in enumerate(clusters, 1):
        for local_idx in members:
            cluster_label[int(bidx[local_idx])] = cid
    butina_rows = []
    for original_idx, cid in cluster_label.items():
        rec = valid_records[original_idx].copy()
        rec["butina_cluster"] = cid
        rec["threshold_tanimoto"] = 0.65
        butina_rows.append(rec)
    pd.DataFrame(butina_rows).to_csv(SD / "S3_butina_cluster_sample.tsv", sep="\t", index=False)

    # Chemical similarity versus target-set similarity on a degree-bounded subset.
    pairs = pd.read_csv(ANALYSIS / "pair_analysis.tsv.gz", sep="\t", compression="gzip", dtype=str, usecols=["compound_internal_id", "target_uniprot_id"])
    target_sets = pairs.groupby("compound_internal_id")["target_uniprot_id"].agg(set).to_dict()
    eligible = [i for i, rec in enumerate(valid_records) if 2 <= len(target_sets.get(rec["compound_internal_id"], set())) <= 50]
    eligible = sorted(eligible, key=lambda i: stable_hash(valid_records[i]["compound_internal_id"]))[:3_000]
    efps = [fps[i] for i in eligible]
    comparison_rows = []
    for local_i, global_i in enumerate(eligible):
        sims = DataStructs.BulkTanimotoSimilarity(efps[local_i], efps)
        order = np.argsort(sims)[::-1]
        neighbours = [j for j in order if j != local_i][:3]
        for local_j in neighbours:
            global_j = eligible[local_j]
            a = valid_records[global_i]["compound_internal_id"]
            b = valid_records[global_j]["compound_internal_id"]
            ta, tb = target_sets[a], target_sets[b]
            comparison_rows.append({"compound_a": a, "compound_b": b, "pair_type": "nearest_chemical_neighbour", "tanimoto": sims[local_j], "target_jaccard": len(ta & tb) / len(ta | tb)})
    if len(eligible) >= 2:
        for _ in range(20_000):
            li, lj = RNG.choice(len(eligible), size=2, replace=False)
            gi, gj = eligible[li], eligible[lj]
            a = valid_records[gi]["compound_internal_id"]
            b = valid_records[gj]["compound_internal_id"]
            ta, tb = target_sets[a], target_sets[b]
            comparison_rows.append({"compound_a": a, "compound_b": b, "pair_type": "random_pair", "tanimoto": DataStructs.TanimotoSimilarity(efps[li], efps[lj]), "target_jaccard": len(ta & tb) / len(ta | tb)})
    pd.DataFrame(comparison_rows).to_csv(SD / "F3E_chemical_vs_target_similarity.tsv", sep="\t", index=False)

    report = {
        "interaction_linked_compounds": len(compounds),
        "smiles_parse_failures": parse_failures,
        "unique_murcko_scaffolds_including_acyclic": len(scaffold_counts),
        "chemical_umap_sample": len(embedding),
        "butina_sample": butina_n,
        "butina_clusters": len(clusters),
        "similarity_target_pairs": len(comparison_rows),
        "random_seed": SEED,
    }
    (QA / "chemical_space_build.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
