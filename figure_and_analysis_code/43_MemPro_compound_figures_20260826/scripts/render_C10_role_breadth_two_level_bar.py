"""Visual-only two-level composition redraw from the existing C10 statistics."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compound_figure_pool"
DATA = OUT / "C10_compound_role_breadth_statistics.tsv"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

INK, SLATE = "#24364B", "#5E7182"
NAVY, TEAL = "#173B6C", "#168C8C"
TWO_ROLE, THREE_ROLE, FOUR_ROLE, FIVE_ROLE = "#3B86AA", "#7B95C6", "#887AA9", "#D6A23A"


def main():
    # This TSV contains the already-computed C10 counts. No pair, target or role
    # table is read, and no role assignment is recalculated.
    stat = pd.read_csv(DATA, sep="\t").set_index("role_breadth")
    count = {key: int(stat.loc[key, "canonical_compound_count"]) for key in ["1", "2", "3", "4", "5+"]}
    total = sum(count.values())
    cross_total = total - count["1"]
    assert total == 240_055 and cross_total == 21_513, "Unexpected frozen C10 source statistics"

    top = [("1 role", count["1"], NAVY), ("2+ roles", cross_total, TEAL)]
    bottom = [("2 roles", count["2"], TWO_ROLE), ("3 roles", count["3"], THREE_ROLE), ("4 roles", count["4"], FOUR_ROLE), ("5+ roles", count["5+"], FIVE_ROLE)]
    fig = plt.figure(figsize=(178 / 25.4, 108 / 25.4), dpi=600)
    ax = fig.add_axes([.16, .18, .73, .65])
    top_y, bottom_y, height = 2.12, .52, .52

    left = 0.0
    for label, n, color in top:
        width = 100 * n / total
        ax.barh(top_y, width, left=left, height=height, color=color, edgecolor="white", linewidth=1.3)
        if label == "1 role":
            ax.text(left + width / 2, top_y, f"1 role\n{n:,}\n{width:.1f}%", ha="center", va="center", color="white", fontsize=8.5, weight="bold", linespacing=1.16)
        else:
            center = left + width / 2
            ax.annotate(f"2+ roles\n{n:,}\n{width:.1f}%", xy=(center, top_y + height / 2), xytext=(105.4, 2.87), ha="left", va="bottom", color=INK, fontsize=8.1, weight="bold", arrowprops={"arrowstyle":"-", "color":TEAL, "lw":.85, "connectionstyle":"angle,angleA=90,angleB=0,rad=0"})
        left += width

    # A thin bracket/connector signals that the second composition expands the
    # 2+ role segment, without using a heavy arrow or an inset panel.
    ax.plot([91.04, 102.4, 102.4, 100.0], [top_y - height / 2, top_y - height / 2, 1.18, bottom_y + height / 2], color=TEAL, lw=.9, solid_capstyle="round")

    left = 0.0
    centers = {}
    for label, n, color in bottom:
        width = 100 * n / cross_total
        ax.barh(bottom_y, width, left=left, height=height, color=color, edgecolor="white", linewidth=1.3)
        centers[label] = left + width / 2
        if label == "2 roles":
            ax.text(centers[label], bottom_y, f"2 roles\n{n:,}\n{width:.1f}%", ha="center", va="center", color="white", fontsize=8.4, weight="bold", linespacing=1.16)
        left += width

    # Narrow segments receive external labels with staggered thin leader lines.
    external = [
        ("3 roles", count["3"], 100 * count["3"] / cross_total, 84.5, 1.45, "right"),
        ("4 roles", count["4"], 100 * count["4"] / cross_total, 96.0, 1.77, "center"),
        ("5+ roles", count["5+"], 100 * count["5+"] / cross_total, 106.2, 1.36, "left"),
    ]
    for label, n, pct, text_x, text_y, ha in external:
        target_x = centers[label]
        ax.annotate(f"{label}\n{n:,} · {pct:.1f}%", xy=(target_x, bottom_y + height / 2), xytext=(text_x, text_y), ha=ha, va="bottom", fontsize=7.1, color=INK, weight="semibold", arrowprops={"arrowstyle":"-", "color":SLATE, "lw":.75, "connectionstyle":"angle,angleA=90,angleB=0,rad=0"})

    ax.text(-1.7, top_y, f"All interaction-linked\ncanonical compounds\nN = {total:,}", ha="right", va="center", fontsize=8.2, color=INK, weight="semibold", linespacing=1.08)
    ax.text(-1.7, bottom_y, f"Cross-role compounds only\nN = {cross_total:,}", ha="right", va="center", fontsize=8.2, color=INK, weight="semibold", linespacing=1.08)
    ax.set_xlim(-2, 115.5)
    ax.set_ylim(-.02, 3.34)
    ax.set_xticks([0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%"])
    ax.tick_params(axis="x", labelsize=7.5, colors=SLATE, length=3, pad=3)
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.spines["bottom"].set_color("#B8C4CC")
    ax.grid(False)

    fig.text(.16, .955, "Membrane-protein role breadth of interaction-linked compounds", fontsize=12.8, color=INK, weight="bold", ha="left")
    fig.text(.16, .906, "Most compounds are linked to a single membrane-protein role; the cross-role subset is dominated by compounds spanning two roles.", fontsize=7.9, color=SLATE, ha="left")
    fig.text(.16, .070, "Role breadth describes database linkage breadth across membrane-protein roles and should not be interpreted as pharmacological selectivity.", fontsize=6.8, color=SLATE, ha="left")
    for suffix, kwargs in {".svg": {}, ".pdf": {}, ".png": {"dpi": 600}}.items():
        fig.savefig(OUT / f"C10_role_breadth_two_level_bar{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
