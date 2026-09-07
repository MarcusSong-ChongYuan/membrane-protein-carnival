"""Audit whether 3D community envelopes are driven by systematic outliers."""
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.covariance import MinCovDet
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "leiden_structural_results"
DATA = OUT / "full_proteome_3d_atlas" / "full_proteome_pca3d_source_data.tsv"

df = pd.read_csv(DATA, sep="\t")
rows, all_dist = [], []
for community, group in df.groupby("macro_community"):
    xyz = group[["PC1", "PC2", "PC3"]].to_numpy()
    if len(group) < 8:
        continue
    mcd = MinCovDet(random_state=42).fit(xyz)
    distance = mcd.mahalanobis(xyz)
    idx = int(np.argmax(distance))
    row = group.iloc[idx]
    rows.append({
        "community": int(community), "community_n": int(len(group)),
        "outlier_uniprot_id": row.uniprot_id,
        "outlier_sequence_length": int(row.sequence_length),
        "outlier_role": row.posthoc_membrane_role,
        "outlier_membrane_mode": row.primary_membrane_mode_v7,
        "outlier_robust_distance_squared": float(distance[idx]),
        "community_95pct_robust_distance_squared": float(np.quantile(distance, .95)),
        "outlier_to_p95_ratio": float(distance[idx] / np.quantile(distance, .95)),
    })
    all_dist.append(pd.DataFrame({"sequence_length": group.sequence_length.to_numpy(), "robust_distance_squared": distance}))

out = pd.DataFrame(rows).sort_values("outlier_to_p95_ratio", ascending=False)
all_dist = pd.concat(all_dist, ignore_index=True)
rho, p = spearmanr(all_dist.sequence_length, all_dist.robust_distance_squared)
summary = {
    "n_proteins": int(len(df)), "n_communities": int(df.macro_community.nunique()),
    "sequence_length_vs_robust_3d_distance_spearman_rho": float(rho),
    "sequence_length_vs_robust_3d_distance_p": float(p),
    "median_length_of_one_top_outlier_per_community": float(out.outlier_sequence_length.median()),
    "median_length_all_proteins": float(df.sequence_length.median()),
    "communities_with_top_outlier_over_5x_community_p95": int((out.outlier_to_p95_ratio > 5).sum()),
}
out.to_csv(OUT / "full_proteome_3d_atlas" / "community_outlier_audit.tsv", sep="\t", index=False)
pd.Series(summary).to_json(OUT / "full_proteome_3d_atlas" / "community_outlier_audit.json", indent=2)
print(pd.Series(summary).to_string())
print("\nTop outlier per community:")
print(out.to_string(index=False))
