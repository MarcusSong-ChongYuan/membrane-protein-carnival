from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


ROOT = Path(r"D:\finale")
V63 = ROOT / "02_V6.3_candidate_20260803"
V62 = V63 / "01_baseline_v6_2"
OLD_DATA = ROOT / "02_展示图表_V6.2" / "data"
READY = ROOT / "03_V6.3.1_readiness_20260803"
OUT = ROOT / "05_V6.3.1_双页模块PPT重制_20260803"

PNG = OUT / "figures" / "png"
PDF = OUT / "figures" / "pdf"
SVG = OUT / "figures" / "svg"
DATA = OUT / "figures" / "data"
for d in (PNG, PDF, SVG, DATA):
    d.mkdir(parents=True, exist_ok=True)


C = {
    "navy": "#173B57",
    "blue": "#3E78A8",
    "teal": "#2A9D8F",
    "green": "#69A65B",
    "orange": "#F0A35E",
    "coral": "#D96C63",
    "purple": "#7A5AA6",
    "gray": "#A9B3BC",
    "dark": "#44515C",
    "light": "#E8EEF2",
    "pale": "#F5F7F8",
}

mpl.rcParams.update(
    {
        "font.family": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.unicode_minus": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)
sns.set_theme(style="whitegrid", rc=mpl.rcParams)


def wrap(value: object, width: int = 24) -> str:
    return "\n".join(textwrap.wrap(str(value).replace("_", " "), width=width))


def panel(ax, letter: str, title: str) -> None:
    ax.set_title(f"{letter}  {title}", loc="left", pad=11, color=C["navy"], fontweight="bold")


def finish_ax(ax, axis="x") -> None:
    ax.grid(False)
    if axis in {"x", "both"}:
        ax.grid(axis="x", color=C["light"], lw=0.8, zorder=0)
    if axis in {"y", "both"}:
        ax.grid(axis="y", color=C["light"], lw=0.8, zorder=0)


def save(fig, stem: str) -> dict[str, str]:
    paths = {}
    for directory, ext in ((PNG, "png"), (PDF, "pdf"), (SVG, "svg")):
        p = directory / f"{stem}.{ext}"
        fig.savefig(p, dpi=220 if ext == "png" else None, bbox_inches=None, facecolor="white")
        paths[ext] = str(p)
    plt.close(fig)
    return paths


def annotate_barh(ax, values, formatter=lambda x: f"{int(x):,}", pad=0.01):
    mx = max(values) if len(values) else 1
    for i, v in enumerate(values):
        ax.text(v + mx * pad, i, formatter(v), va="center", fontsize=9, color=C["dark"])


def load_core():
    protein_cols = [
        "target_uniprot_id", "membrane_class_v52", "evidence_level_v52",
        "website_default_v52", "membrane_role_primary_v63",
        "structural_family_primary_v63", "molecular_function_primary_v63",
        "biological_process_primary_v63", "specialist_classification_primary_v63",
        "reactome_top_level_v63", "hpa_mapping_status_v62",
        "hpa_tissue_detected_count_v62", "hpa_cell_type_detected_count_v62",
    ]
    proteins = pd.read_csv(
        V63 / "05_protein_annotation" / "human_membrane_protein_master_v6_3_candidate.tsv.gz",
        sep="\t", compression="gzip", dtype=str, keep_default_na=False, usecols=protein_cols,
    )
    expr = pd.read_csv(
        V62 / "expression_location_summary_v2.tsv.gz",
        sep="\t", compression="gzip", dtype=str, keep_default_na=False,
    )
    return proteins, expr


def make_protein(proteins: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(14, 6.1), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.02, 1.22, 1.0]})

    ax = axes[0]
    panel(ax, "A", "膜结合类别 × 证据等级")
    ct = pd.crosstab(proteins["membrane_class_v52"], proteins["evidence_level_v52"])
    ct = ct.reindex(index=["A", "B", "C", "unknown"], columns=["E1", "E2", "E3", "E0"], fill_value=0)
    y = np.arange(len(ct))
    left = np.zeros(len(ct))
    colors = [C["teal"], C["blue"], C["orange"], C["gray"]]
    for level, color in zip(ct.columns, colors):
        vals = ct[level].to_numpy()
        ax.barh(y, vals, left=left, color=color, height=0.65, label=level)
        for i, (lft, value) in enumerate(zip(left, vals)):
            if value >= 180:
                ax.text(lft + value / 2, i, f"{value:,}", ha="center", va="center",
                        fontsize=8.5, color="white" if level in {"E1", "E2"} else C["dark"], fontweight="bold")
        left += vals
    ax.set_yticks(y, ["A  整合/跨膜", "B  直接嵌入", "C  外周膜相关", "E0/未知"])
    ax.invert_yaxis()
    ax.set_xlabel("蛋白数量")
    ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.20))
    finish_ax(ax, "x")

    ax = axes[1]
    panel(ax, "B", "五轴分类中的主要膜功能角色")
    role_map = {
        "membrane_associated_enzyme": "膜相关酶",
        "family_defined_membrane_role_unresolved": "家族明确、角色未决",
        "receptor": "受体",
        "transporter": "转运体",
        "membrane_scaffold_or_linker": "膜支架/连接",
        "signaling_regulator": "信号调节",
        "ion_channel": "离子通道",
        "adhesion_recognition": "黏附/识别",
        "immune_or_cell_recognition": "免疫/细胞识别",
        "membrane_trafficking": "膜运输",
        "membrane_organizer": "膜组织",
        "junction_or_adhesion": "连接/黏附",
        "unclassified": "未分类",
    }
    roles = proteins["membrane_role_primary_v63"].replace("", "unclassified").value_counts().head(11).sort_values()
    labels = [role_map.get(x, wrap(x, 20)) for x in roles.index]
    colors = [C["orange"] if x in {"family_defined_membrane_role_unresolved", "unclassified"} else C["blue"] for x in roles.index]
    ax.barh(np.arange(len(roles)), roles.values, color=colors, zorder=2)
    ax.set_yticks(np.arange(len(roles)), labels)
    ax.set_xlabel("蛋白数量")
    annotate_barh(ax, roles.values)
    finish_ax(ax, "x")

    ax = axes[2]
    panel(ax, "C", "独立分类轴的注释覆盖率")
    checks = {
        "结构家族": ~proteins["structural_family_primary_v63"].isin(["", "unclassified"]),
        "分子功能": ~proteins["molecular_function_primary_v63"].isin(["", "unclassified"]),
        "生物过程": ~proteins["biological_process_primary_v63"].isin(["", "unclassified"]),
        "具体膜角色": ~proteins["membrane_role_primary_v63"].isin(["", "unclassified", "family_defined_membrane_role_unresolved"]),
        "专科层级": ~proteins["specialist_classification_primary_v63"].isin(["", "no_specialist_classification"]),
        "Reactome通路": proteins["reactome_top_level_v63"].ne(""),
    }
    cov = pd.Series({k: v.mean() * 100 for k, v in checks.items()}).sort_values()
    ax.barh(np.arange(len(cov)), cov.values, color=[C["purple"] if x == "专科层级" else C["teal"] for x in cov.index])
    ax.set_yticks(np.arange(len(cov)), cov.index)
    ax.set_xlim(0, 105)
    ax.set_xlabel("具有注释的蛋白比例 (%)")
    annotate_barh(ax, cov.values, formatter=lambda x: f"{x:.1f}%", pad=0.008)
    finish_ax(ax, "x")
    paths = save(fig, "M2_protein_classification_rebuilt")
    ct.to_csv(DATA / "M2_class_evidence_counts.tsv", sep="\t")
    cov.rename("coverage_percent").to_csv(DATA / "M2_axis_coverage.tsv", sep="\t")
    return paths, {"rows": len(proteins), "class_evidence": ct.to_dict(), "coverage": cov.to_dict()}


def make_expression(proteins: pd.DataFrame, expr: pd.DataFrame):
    fig = plt.figure(figsize=(14, 6.1), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[0.90, 1.45], height_ratios=[0.82, 1.18])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])]

    mapped = expr["hpa_mapping_status_v62"].eq("mapped")
    tissue_count = pd.to_numeric(expr["hpa_tissue_detected_count_v62"], errors="coerce")
    categories = pd.Series({
        "已映射且检测到RNA": int((mapped & tissue_count.gt(0)).sum()),
        "已映射但检测组织数=0": int((mapped & tissue_count.eq(0)).sum()),
        "未映射/缺失": int((~mapped).sum()),
    })
    ax = axes[0]
    panel(ax, "A", "HPA状态：0值与缺失分开")
    vals = categories.values
    left = 0
    for v, color, label in zip(vals, [C["teal"], C["orange"], C["gray"]], categories.index):
        ax.barh([0], [v], left=[left], color=color, height=0.5, label=label)
        if v > 100:
            ax.text(left + v / 2, 0, f"{v:,}", ha="center", va="center", fontsize=10,
                    color="white" if color == C["teal"] else C["dark"], fontweight="bold")
        left += v
    ax.set_xlim(0, len(expr))
    ax.set_yticks([])
    ax.set_xlabel("蛋白数量")
    ax.legend(frameon=False, fontsize=8.5, loc="lower center", bbox_to_anchor=(0.5, -0.55))
    finish_ax(ax, "x")

    ax = axes[1]
    panel(ax, "B", "已映射蛋白的组织表达广度")
    joined = proteins[["target_uniprot_id", "membrane_role_primary_v63"]].merge(
        expr[["target_uniprot_id", "hpa_mapping_status_v62", "hpa_tissue_detected_count_v62"]],
        on="target_uniprot_id", suffixes=("", "_expr"), how="left"
    )
    joined = joined[joined["hpa_mapping_status_v62_expr"].eq("mapped")].copy()
    joined["detected_tissues"] = pd.to_numeric(joined["hpa_tissue_detected_count_v62_expr"], errors="coerce")
    wanted = ["receptor", "transporter", "ion_channel", "membrane_associated_enzyme"]
    role_cn = {"receptor": "受体", "transporter": "转运体", "ion_channel": "离子通道", "membrane_associated_enzyme": "膜相关酶"}
    data = [joined.loc[joined["membrane_role_primary_v63"].eq(r), "detected_tissues"].dropna().to_numpy() for r in wanted]
    bp = ax.boxplot(data, vert=False, labels=[role_cn[r] for r in wanted], patch_artist=True,
                    showfliers=False, widths=0.6, medianprops={"color": C["navy"], "linewidth": 2})
    for patch, color in zip(bp["boxes"], [C["teal"], C["blue"], C["purple"], C["green"]]):
        patch.set_facecolor(color); patch.set_alpha(0.75)
    ax.set_xlabel("每个蛋白检测到的组织数量（nTPM ≥ 1）")
    finish_ax(ax, "x")

    matrix = pd.read_csv(OLD_DATA / "M3_tissue_expression_matrix.tsv", sep="\t")
    idx_col = matrix.columns[0]
    matrix = matrix.set_index(idx_col)
    means = matrix.mean(axis=1)
    stds = matrix.std(axis=1).replace(0, np.nan)
    z = matrix.sub(means, axis=0).div(stds, axis=0).fillna(0)
    ax = axes[2]
    panel(ax, "C", "组织表达热图（类内相对模式）")
    sns.heatmap(z, cmap="vlag", center=0, vmin=-2.3, vmax=2.3, linewidths=0.35,
                cbar_kws={"label": "每个蛋白类别内 z-score", "shrink": 0.80}, ax=ax)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_xticklabels([wrap(x, 13) for x in z.columns], rotation=32, ha="right")
    ax.set_yticklabels([wrap(x, 19) for x in z.index], rotation=0)

    loc = pd.read_csv(V62 / "protein_subcellular_localization_v2.tsv.gz", sep="\t", compression="gzip", dtype=str, keep_default_na=False)
    loc_counts = loc.drop_duplicates(["target_uniprot_id", "location_term"])["location_term"].value_counts().head(10).sort_values()
    ax = axes[3]
    panel(ax, "D", "亚细胞定位（非互斥）")
    ax.barh(np.arange(len(loc_counts)), loc_counts.values, color=C["blue"])
    ax.set_yticks(np.arange(len(loc_counts)), [wrap(x, 26) for x in loc_counts.index])
    ax.set_xlabel("唯一蛋白数量")
    annotate_barh(ax, loc_counts.values)
    finish_ax(ax, "x")
    paths = save(fig, "M3_expression_localization_rebuilt")
    categories.rename("protein_count").to_csv(DATA / "M3_hpa_mapping_zero_missing.tsv", sep="\t")
    z.to_csv(DATA / "M3_tissue_expression_zscore.tsv", sep="\t")
    return paths, {"mapping": categories.to_dict(), "localization_top": loc_counts.to_dict()}


def make_disease():
    fig, axes = plt.subplots(2, 2, figsize=(14, 6.1), constrained_layout=True)
    axes = axes.ravel()
    source = pd.read_csv(V63 / "04_disease" / "disease_source_entity_v0_1.tsv", sep="\t", dtype=str, keep_default_na=False)
    rel = pd.read_csv(V63 / "04_disease" / "protein_disease_relation_v2_1.tsv", sep="\t", dtype=str, keep_default_na=False)
    ta_counts = pd.read_csv(V63 / "08_figures" / "disease_v63_working" / "M4_therapeutic_area_counts.tsv", sep="\t")
    an_counts = pd.read_csv(V63 / "08_figures" / "disease_v63_working" / "M4_anatomical_system_counts.tsv", sep="\t")
    channels = pd.read_csv(OLD_DATA / "S4_disease_channel_prevalence.tsv", sep="\t")

    ax = axes[0]
    panel(ax, "A", "疾病身份统一结果")
    status = source["canonicalization_status"].value_counts()
    labels = {
        "unique_official_exact_mapping": "官方1:1 exact/equivalent",
        "direct_mondo_id": "直接MONDO ID",
        "non_disease_mondo_entity_review": "非疾病MONDO实体复核",
        "unmapped_source_only": "仅保留来源ID",
        "unique_mapping_to_obsolete_mondo_review": "obsolete MONDO复核",
    }
    order = status.sort_values()
    colors = [C["teal"] if x in {"unique_official_exact_mapping", "direct_mondo_id"} else C["orange"] for x in order.index]
    ax.barh(np.arange(len(order)), order.values, color=colors)
    ax.set_yticks(np.arange(len(order)), [wrap(labels.get(x, x), 23) for x in order.index])
    ax.set_xlabel("来源疾病实体数量")
    annotate_barh(ax, order.values)
    finish_ax(ax, "x")

    ax = axes[1]
    panel(ax, "B", "治疗领域（Open Targets 26.06，多标签）")
    ta = ta_counts.sort_values(ta_counts.columns[-1]).tail(10)
    ax.barh(np.arange(len(ta)), ta.iloc[:, -1].values, color=C["purple"])
    ax.set_yticks(np.arange(len(ta)), [wrap(x, 24) for x in ta.iloc[:, 0]])
    ax.set_xlabel("唯一蛋白–疾病关系")
    annotate_barh(ax, ta.iloc[:, -1].values)
    finish_ax(ax, "x")

    ax = axes[2]
    panel(ax, "C", "解剖系统（DO/MONDO → Uberon，多标签）")
    an = an_counts.sort_values(an_counts.columns[-1]).tail(10)
    ax.barh(np.arange(len(an)), an.iloc[:, -1].values, color=C["blue"])
    ax.set_yticks(np.arange(len(an)), [wrap(x, 24) for x in an.iloc[:, 0]])
    ax.set_xlabel("唯一蛋白–疾病关系")
    annotate_barh(ax, an.iloc[:, -1].values)
    finish_ax(ax, "x")

    ax = axes[3]
    panel(ax, "D", "默认疾病关系中的证据渠道")
    ch = channels.sort_values("relation_percent")
    ax.barh(np.arange(len(ch)), ch["relation_percent"], color=[C["green"] if x != "Literature" else C["orange"] for x in ch["evidence_channel"]])
    ax.set_yticks(np.arange(len(ch)), [wrap(x, 23) for x in ch["evidence_channel"]])
    ax.set_xlim(0, 102)
    ax.set_xlabel("关系覆盖比例 (%)")
    annotate_barh(ax, ch["relation_percent"].values, formatter=lambda x: f"{x:.1f}%", pad=0.006)
    finish_ax(ax, "x")
    paths = save(fig, "M4_disease_ontology_rebuilt")
    return paths, {"source_entities": len(source), "relations": len(rel), "status": status.to_dict()}


def scan_compounds():
    path = V62 / "small_molecule_master_v1_3.tsv"
    usecols = [
        "compound_scope_status", "record_qc_status", "computed_structural_class",
        "molecular_weight", "xlogp", "hbond_donor_count", "hbond_acceptor_count",
        "is_approved_drug", "is_clinical_candidate", "is_endogenous_ligand",
        "is_natural_product", "is_chemical_probe",
    ]
    combos = Counter()
    criteria_pass = Counter()
    complete = 0
    core_ok = 0
    names = ["Approved", "Clinical", "Probe", "Endogenous", "Natural"]
    flag_cols = ["is_approved_drug", "is_clinical_candidate", "is_chemical_probe", "is_endogenous_ligand", "is_natural_product"]
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=180000, low_memory=False):
        flags = chunk[flag_cols].fillna(0).astype(str).isin(["1", "True", "true"])
        for bits, count in flags.value_counts().items():
            if any(bits):
                combos[tuple(bool(x) for x in bits)] += int(count)
        core = chunk[chunk["compound_scope_status"].eq("core") & chunk["record_qc_status"].eq("ok")].copy()
        core_ok += len(core)
        for col in ["molecular_weight", "xlogp", "hbond_donor_count", "hbond_acceptor_count"]:
            core[col] = pd.to_numeric(core[col], errors="coerce")
        valid = core.dropna(subset=["molecular_weight", "xlogp", "hbond_donor_count", "hbond_acceptor_count"])
        complete += len(valid)
        conds = {
            "MW ≤ 500 Da": valid["molecular_weight"].le(500),
            "XlogP ≤ 5": valid["xlogp"].le(5),
            "HBD ≤ 5": valid["hbond_donor_count"].le(5),
            "HBA ≤ 10": valid["hbond_acceptor_count"].le(10),
        }
        violations = np.zeros(len(valid), dtype=int)
        for label, cond in conds.items():
            criteria_pass[label] += int(cond.sum())
            violations += (~cond.to_numpy()).astype(int)
        criteria_pass["≤1项违反（Lipinski口径）"] += int((violations <= 1).sum())
    return names, combos, criteria_pass, complete, core_ok


def draw_upset(ax, names, combos):
    top = combos.most_common(9)
    x = np.arange(len(top))
    values = [count for _, count in top]
    inset = ax.inset_axes([0.10, 0.42, 0.86, 0.53])
    inset.bar(x, values, color=C["blue"], width=0.68)
    inset.set_yscale("log")
    inset.set_ylabel("化合物数量（log）")
    inset.set_xticks([])
    inset.grid(axis="y", color=C["light"], lw=0.7)
    inset.spines[["top", "right"]].set_visible(False)
    for i, v in enumerate(values):
        inset.text(i, v * 1.10, f"{v:,}", ha="center", fontsize=7.5)
    ax.set_xlim(-0.7, len(top) - 0.3)
    ax.set_ylim(-0.6, len(names) - 0.4)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(len(names)), names)
    ax.set_xticks(x, [str(i + 1) for i in x])
    ax.set_xlabel("前9个生物学状态组合")
    ax.grid(False)
    for i, (bits, _) in enumerate(top):
        active = [j for j, bit in enumerate(bits) if bit]
        if active:
            ax.plot([i, i], [min(active), max(active)], color=C["dark"], lw=1.2, zorder=1)
        for j in range(len(names)):
            ax.scatter(i, j, s=34 if bits[j] else 20, color=C["navy"] if bits[j] else C["light"], zorder=2)


def make_compound():
    names, combos, criteria_pass, complete, core_ok = scan_compounds()
    fig, axes = plt.subplots(2, 2, figsize=(14, 6.1), constrained_layout=True)
    axes = axes.ravel()

    ax = axes[0]
    panel(ax, "A", "有生物学状态注释的化合物交集")
    draw_upset(ax, names, combos)

    ax = axes[1]
    panel(ax, "B", "Lipinski性质满足度（描述性，不作收录门槛）")
    crit = pd.Series(criteria_pass) / max(complete, 1) * 100
    crit = crit.sort_values()
    ax.barh(np.arange(len(crit)), crit.values, color=[C["purple"] if "Lipinski" in x else C["teal"] for x in crit.index])
    ax.set_yticks(np.arange(len(crit)), crit.index)
    ax.set_xlim(0, 103)
    ax.set_xlabel(f"完整性质记录中的比例 (%)，n={complete:,}")
    annotate_barh(ax, crit.values, formatter=lambda x: f"{x:.1f}%", pad=0.005)
    finish_ax(ax, "x")

    ax = axes[2]
    panel(ax, "C", "Bemis–Murcko骨架频率与累积覆盖")
    freq = pd.read_csv(V63 / "04_disease" / "scaffold_frequency_v0_1.tsv", sep="\t").head(1000)
    x = freq["scaffold_rank"].to_numpy()
    y = freq["compound_count"].to_numpy()
    cumulative = np.cumsum(y) / 896240 * 100
    ax.plot(x, y, color=C["blue"], lw=1.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("骨架排名（log）")
    ax.set_ylabel("每个骨架的化合物数（log）", color=C["blue"])
    ax.tick_params(axis="y", labelcolor=C["blue"])
    ax2 = ax.twinx()
    ax2.plot(x, cumulative, color=C["orange"], lw=1.8)
    ax2.set_ylabel("前N个骨架的累积覆盖率 (%)", color=C["orange"])
    ax2.tick_params(axis="y", labelcolor=C["orange"])
    ax2.spines["top"].set_visible(False)
    finish_ax(ax, "both")

    ax = axes[3]
    panel(ax, "D", "核心QC化合物的结构类别")
    cls = pd.read_csv(OLD_DATA / "M5_structural_class_counts.tsv", sep="\t")
    value_col = [c for c in cls.columns if c != "structural_class"][0]
    cls = cls.sort_values(value_col).tail(10)
    ax.barh(np.arange(len(cls)), cls[value_col], color=C["blue"])
    ax.set_yticks(np.arange(len(cls)), [wrap(x, 24) for x in cls["structural_class"]])
    ax.set_xlabel("核心QC canonical化合物")
    annotate_barh(ax, cls[value_col].values)
    finish_ax(ax, "x")
    paths = save(fig, "M5_compound_identity_properties_rebuilt")
    combo_rows = []
    for bits, count in combos.most_common():
        combo_rows.append({"combination": "+".join(n for n, b in zip(names, bits) if b), "count": count})
    pd.DataFrame(combo_rows).to_csv(DATA / "M5_biological_status_intersections.tsv", sep="\t", index=False)
    crit.rename("percent").to_csv(DATA / "M5_lipinski_descriptor_pass.tsv", sep="\t")
    return paths, {"core_ok": core_ok, "complete_lipinski": complete, "annotated_combinations": sum(combos.values())}


def lineage_coverage():
    path = V63 / "06_evidence_lineage" / "protein_compound_summary_v6_3_candidate.tsv.gz"
    cols = ["independent_structure_count_v63", "independent_pubchem_assay_count_v63",
            "independent_literature_experiment_count_v63", "independent_evidence_modality_count_v63"]
    total = 0
    counts = Counter()
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", usecols=cols, chunksize=250000, low_memory=False):
        total += len(chunk)
        num = {c: pd.to_numeric(chunk[c], errors="coerce").fillna(0) for c in cols}
        counts["结构谱系"] += int(num[cols[0]].gt(0).sum())
        counts["PubChem assay谱系"] += int(num[cols[1]].gt(0).sum())
        counts["文献实验代理谱系"] += int(num[cols[2]].gt(0).sum())
        counts["≥2种证据模态"] += int(num[cols[3]].ge(2).sum())
    return total, counts


def make_evidence():
    fig, axes = plt.subplots(2, 2, figsize=(14, 6.1), constrained_layout=True)
    axes = axes.ravel()
    source_tier = pd.read_csv(OLD_DATA / "M6_source_evidence_matrix.tsv", sep="\t").set_index("source")
    source_tier["total"] = source_tier.sum(axis=1)
    top = source_tier.sort_values("total").tail(10)

    ax = axes[0]
    panel(ax, "A", "主要来源的BE1/BE2/BE3贡献")
    y = np.arange(len(top))
    offsets = [-0.22, 0, 0.22]
    for tier, color, off in zip(["BE1", "BE2", "BE3"], [C["blue"], C["green"], C["orange"]], offsets):
        vals = top[tier].replace(0, np.nan)
        ax.scatter(vals, y + off, label=tier, color=color, s=42, zorder=3)
    ax.set_xscale("log")
    ax.set_yticks(y, [wrap(x, 24) for x in top.index])
    ax.set_xlabel("证据记录数（log；0不绘制）")
    ax.legend(frameon=False, ncol=3, loc="lower right")
    finish_ax(ax, "x")

    total, cov_counts = lineage_coverage()
    ax = axes[1]
    panel(ax, "B", "蛋白–化合物关系的谱系覆盖")
    cov = pd.Series({k: v / total * 100 for k, v in cov_counts.items()}).sort_values()
    ax.barh(np.arange(len(cov)), cov.values, color=[C["teal"], C["purple"], C["green"], C["blue"]])
    ax.set_yticks(np.arange(len(cov)), [wrap(x, 26) for x in cov.index])
    ax.set_xlim(0, max(cov.max() * 1.22, 5))
    ax.set_xlabel(f"阳性蛋白–化合物关系比例 (%)，n={total:,}")
    annotate_barh(ax, cov.values, formatter=lambda x: f"{x:.1f}%", pad=0.008)
    finish_ax(ax, "x")

    ax = axes[2]
    panel(ax, "C", "结合位点实例的可用类型")
    site = pd.read_csv(V62 / "binding_site_instances_v6_2.tsv.gz", sep="\t", compression="gzip", dtype=str, keep_default_na=False,
                       usecols=["site_type", "record_qc_status"])
    site_map = {
        "experimental_structure_residue_contact": "PDB + 残基接触",
        "experimental_structure_match": "PDB结构匹配",
        "experimental_complex_pocket": "实验复合物口袋",
        "bindingdb_ligand_target_complex": "仅复合物PDB",
        "uniprot_curated_binding_site": "UniProt人工位点",
    }
    readiness = site["site_type"].map(site_map).fillna("其他").value_counts().sort_values()
    ax.barh(np.arange(len(readiness)), readiness.values, color=C["green"])
    ax.set_yticks(np.arange(len(readiness)), [wrap(x, 24) for x in readiness.index])
    ax.set_xlabel("结合位点实例数")
    annotate_barh(ax, readiness.values)
    finish_ax(ax, "x")

    ax = axes[3]
    panel(ax, "D", "负证据身份解析与冲突边界")
    stages = pd.Series({"原始负证据": 10_246_948, "成功映射到正式层": 7_630_659, "未映射复核层": 2_615_046, "重复记录移除": 1_243})
    y = np.arange(len(stages))
    ax.barh(y, stages.values, color=[C["gray"], C["teal"], C["orange"], C["dark"]])
    ax.set_yticks(y, stages.index)
    ax.invert_yaxis(); ax.set_xscale("log")
    ax.set_xlabel("记录数（log）")
    for i, v in enumerate(stages.values):
        ax.text(v * 1.06, i, f"{v:,}", va="center", fontsize=9)
    ax.text(0.98, 0.04, "正式负关系 6,331,307\n正–负冲突关系 68,015",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10, color=C["coral"], fontweight="bold")
    finish_ax(ax, "x")
    paths = save(fig, "M6_binding_evidence_lineage_rebuilt")
    cov.rename("coverage_percent").to_csv(DATA / "M6_lineage_coverage.tsv", sep="\t")
    return paths, {"positive_pairs": total, "lineage_counts": dict(cov_counts), "site_rows": len(site)}


def make_docking():
    fig, axes = plt.subplots(2, 2, figsize=(14, 6.1), constrained_layout=True)
    axes = axes.ravel()
    dock = pd.read_csv(READY / "04_docking_hpc" / "docking_pilot_manifest_v631.tsv", sep="\t", dtype=str, keep_default_na=False)
    order = [x for x in ["R0", "D1", "D2", "D3"] if x in set(dock["docking_tier"])]

    ax = axes[0]
    panel(ax, "A", "HPC Pilot分层")
    tiers = dock["docking_tier"].value_counts().reindex(order).fillna(0)
    ax.bar(np.arange(len(tiers)), tiers.values, color=[C["navy"], C["teal"], C["blue"], C["orange"]][:len(tiers)])
    ax.set_xticks(np.arange(len(tiers)), tiers.index)
    ax.set_ylabel("蛋白–小分子对")
    for i, v in enumerate(tiers.values):
        ax.text(i, v + max(tiers.values) * 0.025, f"{int(v):,}", ha="center", fontsize=10)
    finish_ax(ax, "y")

    ax = axes[1]
    panel(ax, "B", "进入实际Docking前的可用信息")
    readiness = pd.Series({
        "结构化合物身份\n(SMILES + InChIKey)": (dock["standard_smiles"].ne("") & dock["standard_inchikey"].ne("")).sum(),
        "候选受体PDB": dock["candidate_receptor_pdb_id_v2"].ne("").sum(),
        "成对实验结合位点": dock["pair_site_pdb_ids"].ne("").sum(),
        "无正–负冲突": ~dock["positive_negative_conflict_flag_v62"].isin(["1", "true", "True"]),
    })
    readiness.iloc[-1] = int(readiness.iloc[-1].sum()) if hasattr(readiness.iloc[-1], "sum") else int(readiness.iloc[-1])
    vals = readiness.astype(int)
    ax.barh(np.arange(len(vals)), vals.values, color=[C["teal"], C["blue"], C["green"], C["purple"]])
    ax.set_yticks(np.arange(len(vals)), readiness.index)
    ax.set_xlim(0, len(dock) * 1.08)
    ax.set_xlabel("Pilot pairs")
    annotate_barh(ax, vals.values)
    finish_ax(ax, "x")

    ax = axes[2]
    panel(ax, "C", "谱系加分与冲突惩罚后的分数变化")
    old = pd.to_numeric(dock["ranking_score_v2"], errors="coerce").fillna(0)
    new = pd.to_numeric(dock["ranking_score_v63"], errors="coerce").fillna(0)
    delta = new - old
    bins = np.arange(math.floor(delta.min()) - 0.5, math.ceil(delta.max()) + 1.5, 1)
    ax.hist(delta, bins=bins, color=C["teal"], edgecolor="white", lw=0.5)
    ax.axvline(0, color=C["dark"], lw=1.2)
    ax.set_yscale("log")
    ax.set_xlabel("V6.3分数 − V6.2分数")
    ax.set_ylabel("候选对数量（log）")
    finish_ax(ax, "both")

    ax = axes[3]
    panel(ax, "D", "各层候选的V6.3排序分数")
    score = pd.to_numeric(dock["ranking_score_v63"], errors="coerce")
    values = [score[dock["docking_tier"].eq(t)].dropna().to_numpy() for t in order]
    bp = ax.boxplot(values, labels=order, patch_artist=True, showfliers=False,
                    medianprops={"color": C["navy"], "linewidth": 2})
    for p, color in zip(bp["boxes"], [C["navy"], C["teal"], C["blue"], C["orange"]]):
        p.set_facecolor(color); p.set_alpha(0.72)
    ax.set_ylabel("ranking_score_v63")
    finish_ax(ax, "y")
    paths = save(fig, "M8_docking_readiness_rebuilt")
    pd.DataFrame({"tier": tiers.index, "pairs": tiers.values}).to_csv(DATA / "M8_pilot_tiers.tsv", sep="\t", index=False)
    return paths, {"pilot_pairs": len(dock), "targets": dock["target_uniprot_id"].nunique(), "tiers": tiers.to_dict(), "delta_median": float(delta.median())}


def main():
    proteins, expr = load_core()
    metrics = {"status": "PASS", "figures": {}, "counts": {}}
    for key, maker in [
        ("M2", lambda: make_protein(proteins)),
        ("M3", lambda: make_expression(proteins, expr)),
        ("M4", make_disease),
        ("M5", make_compound),
        ("M6", make_evidence),
        ("M8", make_docking),
    ]:
        paths, counts = maker()
        metrics["figures"][key] = paths
        metrics["counts"][key] = counts
        print(key, paths["png"])
    (OUT / "FIGURE_REBUILD_VALIDATION.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
