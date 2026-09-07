"""Readable static companion for the full interactive multi-view map.

All HDBSCAN labels remain in the TSV/HTML.  The raster preview only highlights
the eight largest resolved communities; otherwise a 30-item legend would hide
the data cloud rather than explain it.
"""
from pathlib import Path
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "multiview_results"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 8})

d = pd.read_csv(OUT / "protein_multiview_umap3d.tsv", sep="\t")
sizes = d.loc[d.cluster_id >= 0, "cluster_id"].value_counts().sort_values(ascending=False)
top = sizes.head(8).index.tolist()
d["preview_group"] = d.cluster_id.map(lambda x: f"Cluster {x}" if x in top else ("Noise / unassigned" if x == -1 else "Other resolved communities"))
palette = dict(zip([f"Cluster {x}" for x in top], ["#4DAF4A", "#4BA3D8", "#F29E4C", "#A67FC5", "#E46C7A", "#67A583", "#7B95C6", "#C85E62"]))
palette.update({"Other resolved communities": "#737373", "Noise / unassigned": "#D1D1D1"})
order = ["Noise / unassigned", "Other resolved communities"] + [f"Cluster {x}" for x in top]
fig = plt.figure(figsize=(7.0, 5.3)); ax = fig.add_subplot(projection="3d")
for label in order:
    g = d[d.preview_group == label]
    if g.empty: continue
    ax.scatter(g.UMAP1, g.UMAP2, g.UMAP3, s=5 if label.startswith("Cluster") else 3,
               alpha=.56 if label.startswith("Cluster") else .25, edgecolors="none", c=palette[label], label=f"{label} (n={len(g):,})", rasterized=True)
ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_zlabel("UMAP3"); ax.view_init(elev=22, azim=42)
ax.legend(title="HDBSCAN display", loc="center left", bbox_to_anchor=(1.03, .5), frameon=False, markerscale=2, fontsize=7)
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUT / f"protein_multiview_umap3d_cluster_top8.{ext}", dpi=500, bbox_inches="tight", facecolor="white")
plt.close(fig)
sizes.rename_axis("cluster_id").reset_index(name="protein_count").to_csv(OUT / "multiview_cluster_sizes.tsv", sep="\t", index=False)
