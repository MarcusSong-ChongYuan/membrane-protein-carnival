"""Visual-only external-source redraw of the existing C1 source data.

The input TSV is immutable source data. This script only excludes the
MemPro-derived row and modules without direct external-source attribution.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import FixedFormatter, FixedLocator
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compound_figure_pool"
DATA = OUT / "C1_source_module_compound_coverage.tsv"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

PALETTE = ["#EAF2F6", "#B9D4DF", "#79B4C8", "#3B86AA", "#174A73"]
BLANK, DARK_TEXT, SLATE = "#F5F6F7", "#23364D", "#5E7182"
EXTERNAL_MODULES = [
    "Identity / chemical structure",
    "Protein interaction",
    "Quantitative activity",
    "Binding evidence",
    "Structural binding context",
]
DISPLAY_SOURCE = {"IUPHAR/BPS Guide to PHARMACOLOGY": "GtoPdb (IUPHAR/BPS)"}


def fmt(value):
    value = int(round(value))
    return f"{value / 1_000_000:.1f}M" if value >= 1_000_000 else (f"{value / 1_000:.1f}k" if value >= 1_000 else str(value))


def relative_luminance(rgba):
    channel = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in rgba[:3]]
    return 0.2126 * channel[0] + 0.7152 * channel[1] + 0.0722 * channel[2]


def best_text_color(background):
    white, dark = mpl.colors.to_rgba("#FFFFFF"), mpl.colors.to_rgba(DARK_TEXT)
    lum_bg = relative_luminance(background)
    white_ratio = (max(relative_luminance(white), lum_bg) + .05) / (min(relative_luminance(white), lum_bg) + .05)
    dark_ratio = (max(relative_luminance(dark), lum_bg) + .05) / (min(relative_luminance(dark), lum_bg) + .05)
    return "#FFFFFF" if white_ratio >= dark_ratio else DARK_TEXT


def main():
    table = pd.read_csv(DATA, sep="\t")
    # The existing C1 table already contains the frozen values. "MemPro-derived"
    # is the only non-external row and is deliberately excluded for this view.
    external_rows = table.loc[table["data_source"].ne("MemPro-derived")].copy()
    sources = external_rows["data_source"].drop_duplicates().tolist()
    matrix = external_rows.pivot(index="data_source", columns="compound_module", values="unique_canonical_compounds")
    matrix = matrix.reindex(index=sources, columns=EXTERNAL_MODULES)
    values = matrix.stack().astype(float)
    transformed = np.log1p(matrix.astype(float))
    normalizer = Normalize(vmin=float(np.log1p(values.min())), vmax=float(np.log1p(values.max())), clip=True)
    positions = 0.08 + 0.92 * normalizer(transformed.to_numpy())
    plot_values = np.ma.masked_where(np.isnan(matrix.to_numpy(dtype=float)), positions)
    cmap = LinearSegmentedColormap.from_list("external_source_teal", PALETTE)
    cmap.set_bad(BLANK)

    fig = plt.figure(figsize=(1800 / 150, 920 / 150), dpi=600)
    ax = fig.add_axes([.29, .19, .57, .63])
    image = ax.imshow(plot_values, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            value = matrix.iat[y, x]
            if pd.notna(value):
                ax.text(x, y, fmt(value), ha="center", va="center", fontsize=8.8, weight="bold", color=best_text_color(cmap(positions[y, x])))

    labels = [
        "Identity / chemical structure", "Protein\ninteraction", "Quantitative\nactivity",
        "Binding\nevidence", "Structural binding\ncontext",
    ]
    ax.set_xticks(range(len(EXTERNAL_MODULES)), labels)
    ax.set_yticks(range(len(sources)), [DISPLAY_SOURCE.get(source, source) for source in sources])
    ax.tick_params(axis="x", labelrotation=30, labelsize=9, length=0, pad=9)
    ax.tick_params(axis="y", labelsize=9.5, length=0, pad=8)
    plt.setp(ax.get_xticklabels(), ha="right", rotation_mode="anchor")
    ax.set_xticks(np.arange(-.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="#FFFFFF", linewidth=1.35)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    references = np.array([1_000, 10_000, 100_000, 500_000], dtype=float)
    references = references[(references >= values.min()) & (references <= values.max())]
    cbar = fig.colorbar(image, ax=ax, pad=.028, fraction=.048)
    cbar.locator = FixedLocator(0.08 + 0.92 * normalizer(np.log1p(references)))
    cbar.formatter = FixedFormatter([fmt(value) for value in references])
    cbar.update_ticks()
    cbar.set_label("Canonical compound coverage", fontsize=9, color=DARK_TEXT)
    cbar.ax.tick_params(labelsize=8, colors=SLATE)
    cbar.outline.set_visible(False)

    fig.text(.29, .955, "External compound-source coverage by data module", fontsize=15, weight="bold", color=DARK_TEXT, ha="left")
    fig.text(.29, .910, "Direct external-source contributions only; internally derived annotations are excluded.", fontsize=9, color=SLATE, ha="left")
    fig.text(.29, .875, "Cell = direct source-specific coverage; blank = no direct source-level coverage counted.", fontsize=8.3, color=SLATE, ha="left")
    for suffix, kwargs in {".svg": {}, ".pdf": {}, ".png": {"dpi": 600}}.items():
        fig.savefig(OUT / f"C1_external_source_module_compound_coverage{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
