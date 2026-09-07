"""Visual-only v2 redraw of C1 using the existing source-data TSV.

No release table is opened and no module/source statistic is recomputed.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import FixedLocator, FixedFormatter
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
BLANK = "#F5F6F7"
DARK_TEXT = "#23364D"
SLATE = "#5E7182"
DISPLAY_SOURCE = {"IUPHAR/BPS Guide to PHARMACOLOGY": "GtoPdb (IUPHAR/BPS)"}


def fmt(value):
    value = int(round(value))
    return f"{value / 1_000_000:.1f}M" if value >= 1_000_000 else (f"{value / 1_000:.1f}k" if value >= 1_000 else str(value))


def relative_luminance(rgba):
    channels = []
    for channel in rgba[:3]:
        channels.append(channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4)
    return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2]


def contrast_ratio(foreground, background):
    return (max(relative_luminance(foreground), relative_luminance(background)) + 0.05) / (min(relative_luminance(foreground), relative_luminance(background)) + 0.05)


def main():
    table = pd.read_csv(DATA, sep="\t")
    # Preserve the exact ordering saved by the existing C1 output table.
    sources = table["data_source"].drop_duplicates().tolist()
    # This is the pre-existing C1 processing order, retained verbatim from the
    # current heatmap rather than inferred from TSV row order.
    modules = [
        "Identity / chemical structure",
        "Physicochemical properties",
        "Scaffold / structural annotation",
        "Status annotation",
        "Protein interaction",
        "Quantitative activity",
        "Binding evidence",
        "Structural binding context",
    ]
    matrix = table.pivot(index="data_source", columns="compound_module", values="unique_canonical_compounds")
    matrix = matrix.reindex(index=sources, columns=modules)
    values = matrix.stack().astype(float)
    transformed = np.log1p(matrix.astype(float))
    vmin, vmax = float(np.log1p(values.min())), float(np.log1p(values.max()))
    normalizer = Normalize(vmin=vmin, vmax=vmax, clip=True)
    cmap = LinearSegmentedColormap.from_list("mempro_soft_blue", PALETTE)
    cmap.set_bad(BLANK)

    # Map non-zero coverage to the 0.08–1.00 interval of the palette so the
    # lowest observed non-zero cell remains visibly blue rather than white.
    normalized = normalizer(transformed.to_numpy())
    color_positions = 0.08 + 0.92 * normalized
    plot_values = np.ma.masked_where(np.isnan(matrix.to_numpy(dtype=float)), color_positions)

    fig = plt.figure(figsize=(1800 / 150, 1050 / 150), dpi=600)
    ax = fig.add_axes([0.27, 0.18, 0.60, 0.65])
    image = ax.imshow(plot_values, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")

    dark_rgba = mpl.colors.to_rgba(DARK_TEXT)
    white_rgba = mpl.colors.to_rgba("#FFFFFF")
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            value = matrix.iat[y, x]
            if pd.isna(value):
                continue
            background = cmap(color_positions[y, x])
            text_color = "#FFFFFF" if contrast_ratio(white_rgba, background) >= contrast_ratio(dark_rgba, background) else DARK_TEXT
            ax.text(x, y, fmt(value), ha="center", va="center", fontsize=8.2, color=text_color, weight="bold")

    display_modules = [module.replace(" properties", "\nproperties").replace(" structural annotation", " structural\nannotation").replace(" annotation", "\nannotation").replace(" interaction", "\ninteraction").replace(" activity", "\nactivity").replace(" evidence", "\nevidence").replace(" context", "\ncontext") for module in modules]
    ax.set_xticks(range(len(modules)), display_modules)
    ax.set_yticks(range(len(sources)), [DISPLAY_SOURCE.get(source, source) for source in sources])
    ax.tick_params(axis="x", labelrotation=30, labelsize=8.5, length=0, pad=8)
    ax.tick_params(axis="y", labelsize=9, length=0, pad=8)
    plt.setp(ax.get_xticklabels(), ha="right", rotation_mode="anchor")
    ax.set_xticks(np.arange(-.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="#FFFFFF", linewidth=1.35)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # The colorbar shares exactly the same compressed log1p mapping as cells.
    reference_counts = np.array([1_000, 10_000, 100_000, 500_000], dtype=float)
    reference_counts = reference_counts[(reference_counts >= values.min()) & (reference_counts <= values.max())]
    reference_positions = 0.08 + 0.92 * normalizer(np.log1p(reference_counts))
    cbar = fig.colorbar(image, ax=ax, pad=0.025, fraction=0.045)
    cbar.locator = FixedLocator(reference_positions)
    cbar.formatter = FixedFormatter([fmt(value) for value in reference_counts])
    cbar.update_ticks()
    cbar.set_label("Canonical compound coverage", fontsize=8.7, color=DARK_TEXT)
    cbar.ax.tick_params(labelsize=7.5, colors=SLATE)
    cbar.outline.set_visible(False)

    analysis_n = int(table["analysis_universe_n"].iloc[0])
    fig.text(0.27, 0.955, "Source-by-module compound coverage", fontsize=15, weight="bold", color=DARK_TEXT, ha="left")
    fig.text(0.27, 0.910, f"All formal canonical compounds (Analysis N = {analysis_n:,}); colour intensity is log1p-scaled.", fontsize=9, color=SLATE, ha="left")
    fig.text(0.27, 0.875, "Cell = direct source-specific coverage; blank = no direct source-level coverage counted.", fontsize=8.3, color=SLATE, ha="left")

    for suffix, kwargs in {".svg": {}, ".pdf": {}, ".png": {"dpi": 600}}.items():
        fig.savefig(OUT / f"C1_source_module_compound_coverage_v2{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
