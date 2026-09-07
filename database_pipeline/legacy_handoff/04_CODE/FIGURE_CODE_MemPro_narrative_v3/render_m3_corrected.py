from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr


ROOT = Path(r"C:\Users\Administrator\MemPro_narrative_v3")
DATA = ROOT / "data"
FIG = ROOT / "figures"
ANATOMY = next(Path(r"C:\Users\Administrator\Desktop").rglob("human_anatomy_base_v1.png"))

C = {
    "blue": "#7b95c6", "cyan": "#49c2d9", "green": "#67a583",
    "lgreen": "#a2c986", "peach": "#ffc1a6", "red": "#c85e62",
    "ink": "#27313d", "muted": "#66717e", "grid": "#e5e9ee",
}
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8.4, "axes.titlesize": 9.4,
    "axes.labelsize": 8.3, "xtick.labelsize": 7.3, "ytick.labelsize": 7.3,
    "text.color": C["ink"], "axes.labelcolor": C["ink"],
    "axes.titlecolor": C["ink"], "axes.edgecolor": "#c7cdd4",
    "grid.color": C["grid"], "pdf.fonttype": 42, "svg.fonttype": "none",
    "legend.frameon": False,
})
sns.set_style("whitegrid")

COORD = {
    "Brain": (.50, .09), "Heart": (.52, .34), "Lung": (.43, .30),
    "Liver": (.43, .42), "Kidney": (.59, .49), "Stomach": (.56, .43),
    "Large intestine": (.50, .54), "Small intestine": (.50, .57),
    "Spleen": (.61, .43), "Pancreas": (.53, .47), "Thyroid": (.50, .21),
    "Skin": (.68, .36), "Bone marrow": (.61, .73),
    "Reproductive": (.50, .68), "Bladder": (.50, .66),
}


def panel(ax, label: str, title: str) -> None:
    ax.set_title(f"{label}  {title}", loc="left", weight="bold", pad=8)


def clean(ax, grid: str | None = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)
    if grid:
        ax.grid(axis=grid, alpha=.68)


def body_expression(ax, table: pd.DataFrame) -> None:
    image = np.asarray(Image.open(ANATOMY).convert("RGB"))
    h, w = image.shape[:2]
    ax.imshow(image)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    panel(ax, "A", "RNA-detected membrane proteins span major organs")
    q = table[table.organ.isin(COORD)].nlargest(7, "detected_proteins")
    norm = mpl.colors.Normalize(q.detected_proteins.min(), q.detected_proteins.max())
    cmap = mpl.colormaps["YlGnBu"]
    left_y = [.25, .40, .55, .70]
    right_y = [.28, .47, .66]
    li = ri = 0
    for idx, (_, row) in enumerate(q.iterrows()):
        x0, y0 = COORD[row.organ][0] * w, COORD[row.organ][1] * h
        left = idx % 2 == 0
        if left:
            tx, ty, elbow, ha = .02 * w, left_y[li] * h, .27 * w, "left"
            li += 1
        else:
            tx, ty, elbow, ha = .98 * w, right_y[ri] * h, .73 * w, "right"
            ri += 1
        ax.scatter(x0, y0, s=72, c=[cmap(norm(row.detected_proteins))], edgecolor="white", lw=.8, zorder=5)
        ax.plot([x0, elbow, tx], [y0, ty, ty], color=C["muted"], lw=.65)
        ax.text(tx, ty, f"{row.organ}\n{int(row.detected_proteins):,} proteins", ha=ha,
                va="center", fontsize=7.1, weight="bold")


def main() -> None:
    rna = pd.read_pickle(DATA / "M3_tissue_rna_matrix.pkl")
    ihc = pd.read_pickle(DATA / "M3_ihc_matrix.pkl")
    rna.columns = rna.columns.astype(str).str.strip().str.lower()
    ihc.columns = ihc.columns.astype(str).str.strip().str.lower()
    umap = pd.read_csv(DATA / "M3_expression_umap.tsv", sep="\t")
    body = pd.read_csv(DATA / "M3_body_expression.tsv", sep="\t")

    rows = []
    for tissue in sorted(set(rna.columns) & set(ihc.columns)):
        q = pd.concat([rna[tissue], ihc[tissue]], axis=1, join="inner").dropna()
        q.columns = ["rna", "ihc"]
        if len(q) >= 100 and q.ihc.nunique() > 1:
            rho, pvalue = spearmanr(np.log1p(q.rna), q.ihc)
            rows.append((tissue, len(q), float(rho), float(pvalue)))
    concordance = pd.DataFrame(rows, columns=["tissue", "n", "spearman_rho", "pvalue"]).sort_values("spearman_rho")
    concordance.to_csv(DATA / "M3_RNA_IHC_concordance.tsv", sep="\t", index=False)

    tissue_order = leaves_list(linkage(pdist(np.log1p(rna.values).T, metric="correlation"), method="average"))
    tissue_names = list(rna.columns[tissue_order][:20])
    heat = np.log1p(rna[tissue_names])
    heat = (heat - heat.mean()) / (heat.std() + 1e-6)
    proteins = heat.var(axis=1).nlargest(150).index

    fig = plt.figure(figsize=(14, 9.5))
    fig.text(.03, .968, "M3 | Tissue RNA and protein staining define complementary membrane-protein contexts",
             fontsize=15, weight="bold", va="top")
    fig.text(.03, .932, "Only HPA mapped-and-measured records enter denominators; RNA and IHC remain distinct measurement layers.",
             fontsize=9, color=C["muted"], va="top")
    gs = fig.add_gridspec(2, 3, left=.035, right=.985, top=.87, bottom=.09,
                          width_ratios=[.9, 1.25, 1.0], hspace=.36, wspace=.31)

    ax = fig.add_subplot(gs[:, 0])
    body_expression(ax, body)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B", "Tissue-expression UMAP reveals protein neighborhoods")
    for cls, color in zip(["A", "B", "C"], [C["blue"], C["green"], C["peach"]]):
        q = umap[umap.membrane_class.eq(cls)]
        ax.scatter(q.umap1, q.umap2, s=7, c=color, alpha=.30, label=f"Class {cls}", rasterized=True)
    ax.legend(fontsize=7.1, markerscale=2, ncol=3)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    clean(ax, None)

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "C", "Co-expression clustering identifies tissue blocks")
    sns.heatmap(heat.loc[proteins, tissue_names].T, ax=ax, cmap="vlag", center=0,
                xticklabels=False, yticklabels=[x.title() for x in tissue_names],
                cbar_kws={"label": "tissue-wise z score"})
    ax.set_xlabel("150 variable membrane proteins")
    ax.set_ylabel("")
    ax.tick_params(axis="y", rotation=0)

    ax = fig.add_subplot(gs[:, 2])
    panel(ax, "D", "RNA–IHC concordance is tissue-dependent")
    q = concordance.tail(20)
    y = np.arange(len(q))
    ax.hlines(y, 0, q.spearman_rho, color=C["grid"], lw=3)
    marks = ax.scatter(q.spearman_rho, y, c=q.n, cmap="YlGnBu", s=62, edgecolor="white", lw=.5)
    ax.axvline(0, color=C["muted"], lw=.8)
    ax.set_yticks(y, [x.title() for x in q.tissue])
    ax.set_xlabel("Spearman ρ: log1p RNA vs IHC detection")
    cbar = fig.colorbar(marks, ax=ax, fraction=.04, pad=.02)
    cbar.set_label("matched proteins")
    clean(ax, "x")

    fig.text(.035, .018,
             "Source: HPA 25.1 through MemPro V7.1.1 expression layers; mapped + measured only; case-normalized tissue matching; UMAP seed=42.",
             fontsize=7, color=C["muted"])
    for ext in ("png", "svg", "pdf"):
        fig.savefig(FIG / f"M3_expression_space_and_RNA_IHC_cleanfinal.{ext}",
                    dpi=450 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"M3 corrected: {len(concordance)} matched tissues")


if __name__ == "__main__":
    main()
