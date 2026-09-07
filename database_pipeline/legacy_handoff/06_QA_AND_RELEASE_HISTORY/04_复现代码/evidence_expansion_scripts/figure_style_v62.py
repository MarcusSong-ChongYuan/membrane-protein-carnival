from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Circle, FancyBboxPatch, Polygon
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working\figures\v62_main_20260730")
SVG_DIR = ROOT / "svg"
PDF_DIR = ROOT / "pdf"
PNG_DIR = ROOT / "png"
DATA_DIR = ROOT / "data"
CAPTION_DIR = ROOT / "captions"
QA_DIR = ROOT / "qa"
for directory in (SVG_DIR, PDF_DIR, PNG_DIR, DATA_DIR, CAPTION_DIR, QA_DIR):
    directory.mkdir(parents=True, exist_ok=True)

COLORS = {
    "A": "#0072B2",
    "B": "#009E73",
    "C": "#E69F00",
    "unknown": "#B7B7B7",
    "E1": "#2F5597",
    "E2": "#70AD47",
    "E3": "#ED7D31",
    "E0": "#A6A6A6",
    "BE1": "#2F5597",
    "BE2": "#70AD47",
    "BE3": "#ED7D31",
    "positive": "#1B9E77",
    "negative": "#7F7F7F",
    "conflict": "#D55E00",
    "review": "#CC79A7",
    "excluded": "#B7B7B7",
    "duplicate": "#4D4D4D",
    "R0": "#6A3D9A",
    "D1": "#1B9E77",
    "D2": "#D9A400",
    "D3": "#5B7FA3",
    "navy": "#17365D",
    "blue": "#5B9BD5",
    "light_blue": "#D9EAF7",
    "grid": "#D9D9D9",
    "text": "#202124",
}


def configure_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.titlesize": 9.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 8.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 7.0,
            "figure.titlesize": 13.0,
            "figure.titleweight": "bold",
            "axes.edgecolor": "#A6A6A6",
            "axes.linewidth": 0.6,
            "grid.color": COLORS["grid"],
            "grid.linewidth": 0.45,
            "grid.alpha": 0.65,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def make_figure(title: str):
    configure_style()
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.25), constrained_layout=True)
    fig.suptitle(title, x=0.015, y=1.015, ha="left", color=COLORS["navy"])
    return fig, axes.ravel()


def panel(ax, letter: str, title: str) -> None:
    ax.text(
        -0.08,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
        color=COLORS["navy"],
    )
    ax.set_title(title, loc="left", pad=8, color=COLORS["text"])


def clean(ax, grid_axis: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(grid_axis == "x", axis="x")
    ax.grid(grid_axis == "y", axis="y")


def label_bars(ax, fmt="{:,.0f}", fontsize=6.5, padding=2) -> None:
    for container in ax.containers:
        try:
            labels = [fmt.format(v) if np.isfinite(v) else "" for v in container.datavalues]
            ax.bar_label(container, labels=labels, padding=padding, fontsize=fontsize)
        except Exception:
            continue


def save(fig, stem: str) -> dict[str, str]:
    paths = {
        "svg": SVG_DIR / f"{stem}.svg",
        "pdf": PDF_DIR / f"{stem}.pdf",
        "png": PNG_DIR / f"{stem}.png",
    }
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    fig.savefig(paths["png"], bbox_inches="tight", dpi=600)
    plt.close(fig)
    return {key: str(value) for key, value in paths.items()}


def write_table(frame: pd.DataFrame, name: str) -> str:
    path = DATA_DIR / name
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def write_json(data: dict, name: str) -> str:
    path = QA_DIR / name
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    return str(path)


def clustered_order(frame: pd.DataFrame, axis: int = 0) -> list[int]:
    values = frame.to_numpy(dtype=float)
    if axis == 1:
        values = values.T
    if values.shape[0] < 3:
        return list(range(values.shape[0]))
    values = np.nan_to_num(values, nan=np.nanmedian(values))
    distance = pdist(values, metric="correlation")
    if not np.isfinite(distance).all() or np.allclose(distance, 0):
        return list(range(values.shape[0]))
    return leaves_list(linkage(distance, method="average")).tolist()


def draw_flow_boxes(
    ax,
    labels: Iterable[str],
    counts: Iterable[int | float | None] | None = None,
    colors: Iterable[str] | None = None,
) -> None:
    labels = list(labels)
    counts = list(counts) if counts is not None else [None] * len(labels)
    colors = list(colors) if colors is not None else [COLORS["light_blue"]] * len(labels)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    xs = np.linspace(0.07, 0.93, len(labels))
    width = min(0.2, 0.75 / len(labels))
    for index, (x, label, count, color) in enumerate(zip(xs, labels, counts, colors)):
        box = FancyBboxPatch(
            (x - width / 2, 0.38),
            width,
            0.24,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            linewidth=0.9,
            edgecolor=COLORS["navy"],
            facecolor=color,
        )
        ax.add_patch(box)
        text = label if count is None else f"{label}\n{count:,.0f}"
        ax.text(x, 0.5, text, ha="center", va="center", fontsize=7, color=COLORS["text"])
        if index < len(labels) - 1:
            ax.annotate(
                "",
                xy=(xs[index + 1] - width / 2 - 0.008, 0.5),
                xytext=(x + width / 2 + 0.008, 0.5),
                arrowprops=dict(arrowstyle="->", color=COLORS["navy"], lw=1.0),
            )


def draw_human_navigation(
    ax,
    system_values: dict[str, float],
    color: str = "#5B9BD5",
    max_systems: int = 7,
) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    # Deliberately schematic: a navigation silhouette, not an anatomical measurement.
    ax.add_patch(Circle((0.22, 0.85), 0.07, facecolor="#E8EEF5", edgecolor=COLORS["navy"], lw=0.8))
    torso = Polygon(
        [(0.15, 0.76), (0.29, 0.76), (0.33, 0.43), (0.27, 0.38), (0.17, 0.38), (0.11, 0.43)],
        closed=True,
        facecolor="#E8EEF5",
        edgecolor=COLORS["navy"],
        lw=0.8,
    )
    ax.add_patch(torso)
    for x1, x2 in [(0.15, 0.05), (0.29, 0.39)]:
        ax.plot([x1, x2], [0.7, 0.42], color=COLORS["navy"], lw=5, solid_capstyle="round")
    for x1, x2 in [(0.19, 0.15), (0.25, 0.29)]:
        ax.plot([x1, x2], [0.39, 0.08], color=COLORS["navy"], lw=6, solid_capstyle="round")
    ordered = sorted(system_values.items(), key=lambda item: item[1], reverse=True)[:max_systems]
    if not ordered:
        return
    vmax = max(value for _, value in ordered) or 1
    positions = {
        "Nervous": (0.22, 0.85),
        "Cardiovascular": (0.22, 0.64),
        "Respiratory": (0.19, 0.67),
        "Digestive/hepatic": (0.22, 0.52),
        "Renal/urinary": (0.25, 0.48),
        "Immune/hematologic": (0.18, 0.46),
        "Endocrine/metabolic": (0.22, 0.57),
        "Reproductive": (0.22, 0.40),
        "Musculoskeletal": (0.30, 0.31),
        "Skin": (0.12, 0.58),
        "Other/multisystem": (0.32, 0.78),
    }
    for name, value in ordered:
        x, y = positions.get(name, (0.32, 0.78))
        radius = 0.012 + 0.025 * math.sqrt(value / vmax)
        ax.add_patch(Circle((x, y), radius, facecolor=color, edgecolor="white", lw=0.6, alpha=0.9))
    bar_ax = ax.inset_axes([0.49, 0.08, 0.49, 0.84])
    names = [name for name, _ in ordered][::-1]
    values = [value for _, value in ordered][::-1]
    positions_y = np.arange(len(names))
    bars = bar_ax.barh(positions_y, values, color=color, alpha=0.88)
    bar_ax.set_yticks(positions_y, names)
    bar_ax.bar_label(bars, labels=[f"{int(v):,}" for v in values], padding=2, fontsize=6.2)
    bar_ax.set_xlabel("Unique proteins / pairs")
    clean(bar_ax, "x")


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.array([]), np.array([])
    x = np.sort(values)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y
