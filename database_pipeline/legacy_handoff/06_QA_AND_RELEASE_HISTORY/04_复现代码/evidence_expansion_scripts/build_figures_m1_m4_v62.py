from __future__ import annotations

import json
import math
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
from scipy.stats import spearmanr

from figure_style_v62 import (
    CAPTION_DIR,
    COLORS,
    DATA_DIR,
    ROOT,
    clean,
    clustered_order,
    draw_flow_boxes,
    draw_human_navigation,
    ecdf,
    make_figure,
    panel,
    save,
    write_json,
    write_table,
)


REL = Path(
    r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730"
)
RUN_V62 = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\v62_completion_20260729\qa"
)
V61 = Path(
    r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_1_20260729"
)


def truth(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"1", "true", "yes", "y"})


def human_number(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return f"{value:.0f}"


SYSTEM_KEYWORDS = {
    "Nervous": [
        "brain",
        "cerebr",
        "neuro",
        "alzheimer",
        "parkinson",
        "schizoph",
        "epilep",
        "autism",
        "dementia",
        "glioma",
        "astrocyt",
        "spinal",
    ],
    "Cardiovascular": [
        "cardio",
        "heart",
        "vascular",
        "arter",
        "hypertension",
        "arrhythm",
        "aortic",
    ],
    "Respiratory": ["lung", "pulmonary", "asthma", "bronch", "respirat"],
    "Digestive/hepatic": [
        "liver",
        "hepato",
        "colon",
        "colorectal",
        "gastric",
        "stomach",
        "pancrea",
        "crohn",
        "colitis",
        "bowel",
        "intestinal",
        "esoph",
    ],
    "Renal/urinary": ["renal", "kidney", "nephro", "urinary", "bladder"],
    "Endocrine/metabolic": [
        "diabetes",
        "thyroid",
        "adrenal",
        "obesity",
        "metabolic",
        "pituitary",
        "insulin",
    ],
    "Immune/hematologic": [
        "leukemia",
        "lymphoma",
        "myeloma",
        "immune",
        "immuno",
        "autoimmune",
        "blood",
        "anemia",
        "hemat",
        "thymoma",
    ],
    "Musculoskeletal": ["muscle", "skeletal", "bone", "osteo", "muscular"],
    "Reproductive": [
        "breast",
        "ovarian",
        "uterine",
        "prostate",
        "testicular",
        "cervical",
        "endometr",
    ],
    "Skin": ["melanoma", "skin", "epiderm", "dermat"],
}


def map_system(value: str) -> str:
    text = str(value).lower()
    for system, keywords in SYSTEM_KEYWORDS.items():
        if any(keyword in text for keyword in keywords):
            return system
    return "Other/multisystem"


def structure_grade(frame: pd.DataFrame) -> pd.Series:
    membrane_structure = (
        truth(frame["opm_present_v5"])
        | truth(frame["pdbtm_present_v5"])
        | frame["opm_pdb_ids_v5"].fillna("").ne("")
        | frame["pdbtm_pdb_ids_v5"].fillna("").ne("")
    )
    any_pdb = frame["pdb_ids"].fillna("").ne("")
    any_af = frame["alphafolddb_ids"].fillna("").ne("")
    return pd.Series(
        np.select(
            [membrane_structure, any_pdb, any_af],
            ["S2 membrane PDB", "S3 other PDB", "S4 AlphaFold only"],
            default="S0 no structure",
        ),
        index=frame.index,
    )


protein_cols = [
    "target_uniprot_id",
    "approved_symbol",
    "protein_name",
    "membrane_class_v52",
    "evidence_level_v52",
    "functional_primary_class_v5",
    "transmembrane_count_v5",
    "sequence_length_v53",
    "opm_present_v5",
    "opm_pdb_ids_v5",
    "pdbtm_present_v5",
    "pdbtm_pdb_ids_v5",
    "pdb_ids",
    "alphafolddb_ids",
    "hpa_predicted_membrane_v5",
    "htp_present_v5",
    "membranome_present_v5",
    "oligomeric_state_consensus_v62",
    "subunit_evidence_grade_v62",
    "subunit_conflict_flag_v62",
    "hpa_mapping_status_v62",
    "hpa_tissue_detected_count_v62",
    "hpa_cell_type_detected_count_v62",
    "hpa_ihc_detected_tissue_count_v62",
    "hpa_main_locations_v62",
    "hpa_additional_locations_v62",
]
protein = pd.read_csv(
    REL / "human_membrane_protein_master_v6_2.tsv",
    sep="\t",
    usecols=protein_cols,
    low_memory=False,
)
protein["functional_class"] = (
    protein["functional_primary_class_v5"].fillna("unclassified").astype(str)
)
protein["membrane_class"] = protein["membrane_class_v52"].fillna("unknown")
protein["evidence_level"] = protein["evidence_level_v52"].fillna("E0")
protein["tm"] = pd.to_numeric(protein["transmembrane_count_v5"], errors="coerce").fillna(0)
protein["length"] = pd.to_numeric(protein["sequence_length_v53"], errors="coerce")
protein["structure_grade"] = structure_grade(protein)
protein_index = protein.set_index("target_uniprot_id")
class_map = protein_index["functional_class"]

disease = pd.read_csv(
    REL / "protein_gene_disease_relations_v6_2.tsv",
    sep="\t",
    low_memory=False,
)
disease["organ_system"] = disease["disease_name"].map(map_system)

pair = pd.read_csv(
    REL / "protein_compound_summary_v6_2.tsv.gz",
    sep="\t",
    compression="gzip",
    usecols=[
        "target_uniprot_id",
        "compound_internal_id",
        "best_binding_evidence_level",
        "independent_source_count",
        "best_standard_value_nM",
    ],
    low_memory=False,
)

site = pd.read_csv(
    REL / "binding_site_instances_v6_2.tsv.gz",
    sep="\t",
    compression="gzip",
    usecols=["target_uniprot_id", "source_database", "site_type", "pdb_ids"],
    low_memory=False,
)

assembly = pd.read_csv(
    REL / "protein_biological_assembly_v2.tsv.gz",
    sep="\t",
    compression="gzip",
    usecols=["target_uniprot_id", "pdb_id", "biological_assembly_id"],
    low_memory=False,
)

localization = pd.read_csv(
    REL / "protein_subcellular_localization_v2.tsv.gz",
    sep="\t",
    compression="gzip",
    usecols=["target_uniprot_id", "location_term", "location_role"],
    low_memory=False,
)

ihc = pd.read_csv(
    REL / "protein_tissue_cell_ihc_expression_v2.tsv.gz",
    sep="\t",
    compression="gzip",
    usecols=["target_uniprot_id", "tissue", "detection_status"],
    low_memory=False,
)


def evidence_source_counts() -> pd.Series:
    counts = Counter()
    for chunk in pd.read_csv(
        REL / "binding_evidence_master_v6_2.tsv.gz",
        sep="\t",
        compression="gzip",
        usecols=["source_database"],
        chunksize=300_000,
        low_memory=False,
    ):
        counts.update(chunk["source_database"].fillna("Unknown").astype(str))
    return pd.Series(dict(counts), dtype=int).sort_values(ascending=False)


evidence_sources = evidence_source_counts()


def load_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    tissue = pd.read_csv(
        REL / "protein_tissue_expression_v2.tsv.gz",
        sep="\t",
        compression="gzip",
        usecols=["target_uniprot_id", "source_dataset", "tissue", "value"],
        low_memory=False,
    )
    tissue = tissue.loc[tissue["source_dataset"].eq("HPA_consensus_tissue_RNA")].copy()
    tissue["value"] = pd.to_numeric(tissue["value"], errors="coerce").fillna(0).astype("float32")
    tissue["functional_class"] = tissue["target_uniprot_id"].map(class_map).fillna("unclassified")

    cell = pd.read_csv(
        REL / "protein_cell_type_expression_v2.tsv.gz",
        sep="\t",
        compression="gzip",
        usecols=["target_uniprot_id", "source_dataset", "cell_type", "value"],
        low_memory=False,
    )
    cell = cell.loc[cell["source_dataset"].eq("HPA_single_cell_type_RNA")].copy()
    cell["value"] = pd.to_numeric(cell["value"], errors="coerce").fillna(0).astype("float32")
    cell["functional_class"] = cell["target_uniprot_id"].map(class_map).fillna("unclassified")
    return tissue, cell


tissue, cell = load_expression()


def m1_database_architecture() -> dict:
    fig, axes = make_figure("M1  Database architecture, integration and quality control")

    ax = axes[0]
    panel(ax, "A", "Source architecture")
    ax.axis("off")
    groups = [
        ("Protein", ["UniProt", "HPA", "HTP", "Membranome"], "#D9EAF7"),
        ("Structure", ["PDBe", "OPM", "PDBTM", "PDBbind"], "#E2F0D9"),
        ("Interaction", ["ChEMBL", "BindingDB", "PubChem", "BRENDA"], "#FCE4D6"),
        ("Disease", ["Open Targets", "UniProtKB"], "#E4DFEC"),
    ]
    y_positions = [0.78, 0.57, 0.36, 0.15]
    for (title, names, color), y in zip(groups, y_positions):
        ax.add_patch(
            FancyBboxPatch(
                (0.01, y - 0.075),
                0.28,
                0.13,
                boxstyle="round,pad=0.01",
                facecolor=color,
                edgecolor=COLORS["navy"],
                lw=0.7,
            )
        )
        ax.text(0.03, y + 0.035, title, fontweight="bold", fontsize=7)
        ax.text(0.03, y - 0.025, " · ".join(names), fontsize=6.4)
        ax.annotate("", xy=(0.42, 0.5), xytext=(0.30, y), arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.add_patch(
        FancyBboxPatch(
            (0.42, 0.34),
            0.22,
            0.32,
            boxstyle="round,pad=0.015",
            facecolor="#FFF2CC",
            edgecolor=COLORS["navy"],
            lw=0.9,
        )
    )
    ax.text(0.53, 0.57, "MemPro integration", ha="center", fontweight="bold", fontsize=8)
    ax.text(
        0.53,
        0.47,
        "Identity mapping\nCanonicalization\nEvidence grading\nConflict control",
        ha="center",
        va="center",
        fontsize=6.8,
    )
    outputs = ["Protein", "Disease", "Binding & sites", "Compound"]
    for index, output in enumerate(outputs):
        y = 0.78 - index * 0.19
        ax.annotate("", xy=(0.74, y), xytext=(0.65, 0.5), arrowprops=dict(arrowstyle="->", lw=0.7))
        ax.add_patch(
            FancyBboxPatch(
                (0.74, y - 0.045),
                0.24,
                0.09,
                boxstyle="round,pad=0.01",
                facecolor="#EAF2F8",
                edgecolor=COLORS["navy"],
                lw=0.7,
            )
        )
        ax.text(0.86, y, output, ha="center", va="center", fontsize=7)

    ax = axes[1]
    panel(ax, "B", "Standardization and evidence pipeline")
    draw_flow_boxes(
        ax,
        ["Raw\nrecords", "ID\nmapping", "Parent /\nform", "Cross-source\ndedup", "Evidence\ntier", "Release"],
        colors=["#F2F2F2", "#D9EAF7", "#D9EAD3", "#FFF2CC", "#FCE4D6", "#DDEBF7"],
    )

    ax = axes[2]
    panel(ax, "C", "Record disposition and identity resolution")
    positive = {
        "Positive": 1_575_210,
        "Review": 1_619_083,
        "Excluded": 96_205,
        "Duplicate": 16,
    }
    negative = {
        "Mapped": 7_630_659,
        "Review": 2_615_046,
        "Duplicate": 1_243,
    }
    y = [1, 0]
    left = np.zeros(2)
    colors = {
        "Positive": COLORS["positive"],
        "Mapped": COLORS["positive"],
        "Review": COLORS["review"],
        "Excluded": COLORS["excluded"],
        "Duplicate": COLORS["duplicate"],
    }
    all_categories = ["Positive", "Mapped", "Review", "Excluded", "Duplicate"]
    for category in all_categories:
        values = [positive.get(category, 0), negative.get(category, 0)]
        ax.barh(y, values, left=left, color=colors[category], label=category, height=0.58)
        left += values
    ax.set_yticks(y, ["Incremental positive\n3,290,514 records", "Negative identity\n10,246,948 records"])
    ax.set_xlabel("Input records")
    ax.xaxis.set_major_formatter(lambda value, _: human_number(value))
    ax.legend(ncol=3, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.35))
    clean(ax, "x")

    ax = axes[3]
    panel(ax, "D", "Relational data model")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    nodes = {
        "Membrane\nprotein": (0.5, 0.55, "#D9EAF7"),
        "Gene–disease": (0.18, 0.78, "#E4DFEC"),
        "Compound\nparent / form": (0.82, 0.78, "#FCE4D6"),
        "Binding\nevidence": (0.82, 0.3, "#FFF2CC"),
        "Binding site": (0.5, 0.12, "#E2F0D9"),
        "Expression &\nlocalization": (0.18, 0.3, "#DDEBF7"),
    }
    for label, (x, y0, color) in nodes.items():
        ax.add_patch(
            FancyBboxPatch(
                (x - 0.12, y0 - 0.065),
                0.24,
                0.13,
                boxstyle="round,pad=0.012",
                facecolor=color,
                edgecolor=COLORS["navy"],
                lw=0.8,
            )
        )
        ax.text(x, y0, label, ha="center", va="center", fontsize=7)
    for label in nodes:
        if label == "Membrane\nprotein":
            continue
        x, y0, _ = nodes[label]
        ax.plot([0.5, x], [0.55, y0], color="#7F8C8D", lw=0.8, zorder=0)
    ax.text(0.5, 0.95, "Primary key: UniProt accession", ha="center", fontsize=6.5, color="#666666")

    ax = axes[4]
    panel(ax, "E", "Source-by-module coverage")
    sources = [
        "UniProt",
        "HPA",
        "HTP",
        "Membranome",
        "OPM/PDBTM",
        "PDBe",
        "PDBbind",
        "BindingDB",
        "ChEMBL",
        "PubChem BioAssay",
        "BRENDA",
        "Open Targets",
    ]
    modules = ["Protein", "Structure", "Expression", "Disease", "Interaction", "Site"]
    matrix = pd.DataFrame(0.0, index=sources, columns=modules)
    matrix.loc["UniProt", ["Protein", "Disease", "Site"]] = [
        len(protein),
        int((disease["source_database"] == "UniProtKB").sum()),
        int((site["source_database"] == "UniProt").sum()),
    ]
    matrix.loc["HPA", ["Protein", "Expression"]] = [
        int((protein["hpa_mapping_status_v62"] == "mapped").sum()),
        int(tissue["target_uniprot_id"].nunique()),
    ]
    matrix.loc["HTP", "Protein"] = int(truth(protein["htp_present_v5"]).sum())
    matrix.loc["Membranome", "Protein"] = int(truth(protein["membranome_present_v5"]).sum())
    matrix.loc["OPM/PDBTM", "Structure"] = int(
        (truth(protein["opm_present_v5"]) | truth(protein["pdbtm_present_v5"])).sum()
    )
    matrix.loc["PDBe", ["Structure", "Site"]] = [
        assembly["target_uniprot_id"].nunique(),
        int((site["source_database"] == "PDBe").sum()),
    ]
    for source in ["PDBbind", "BindingDB"]:
        matrix.loc[source, "Interaction"] = int(evidence_sources.get(source, 0))
        matrix.loc[source, "Site"] = int((site["source_database"] == source).sum())
    for source in ["ChEMBL", "PubChem BioAssay", "BRENDA"]:
        matrix.loc[source, "Interaction"] = int(evidence_sources.get(source, 0))
    matrix.loc["Open Targets", "Disease"] = int((disease["source_database"] == "Open Targets").sum())
    log_matrix = np.log10(matrix + 1)
    annotations = matrix.map(lambda value: "" if value == 0 else human_number(value))
    sns.heatmap(
        log_matrix,
        cmap=sns.light_palette(COLORS["blue"], as_cmap=True),
        annot=annotations,
        fmt="",
        cbar_kws={"label": "log10(count + 1)"},
        linewidths=0.4,
        linecolor="white",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=35)

    ax = axes[5]
    panel(ax, "F", "Frozen V6.2 release snapshot")
    metrics = pd.Series(
        {
            "Membrane proteins": 10_997,
            "Disease relations": 10_264,
            "Binding sites": 95_598,
            "Positive pairs": 1_502_456,
            "Canonical compounds": 2_016_064,
            "Positive evidence": 3_003_306,
            "Negative pairs": 6_331_307,
            "Negative evidence": 7_784_223,
        }
    ).sort_values()
    bars = ax.barh(metrics.index, metrics.values, color=COLORS["blue"])
    ax.set_xscale("log")
    ax.set_xlabel("Entities / relations (log scale)")
    ax.bar_label(bars, labels=[f"{value:,.0f}" for value in metrics.values], padding=3, fontsize=6.5)
    clean(ax, "x")
    fig_paths = save(fig, "M1_database_architecture_quality")
    write_table(matrix.reset_index(names="source"), "M1_source_module_matrix.tsv")
    return fig_paths


def m2_membrane_proteome() -> dict:
    fig, axes = make_figure("M2  Landscape of the human membrane proteome")
    ax = axes[0]
    panel(ax, "A", "Membrane class and evidence grade")
    class_order = ["A", "B", "C", "unknown"]
    evidence_order = ["E1", "E2", "E3", "E0"]
    cross = (
        pd.crosstab(protein["membrane_class"], protein["evidence_level"])
        .reindex(index=class_order, columns=evidence_order, fill_value=0)
    )
    left = np.zeros(len(cross))
    for evidence in evidence_order:
        values = cross[evidence].to_numpy()
        ax.barh(cross.index, values, left=left, label=evidence, color=COLORS[evidence])
        left += values
    ax.set_xlabel("Unique proteins")
    ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.28))
    clean(ax, "x")

    ax = axes[1]
    panel(ax, "B", "Functional hierarchy")
    counts = protein["functional_class"].value_counts().head(10).sort_values()
    bars = ax.barh(counts.index.str.replace("_", " "), counts.values, color=COLORS["blue"])
    ax.bar_label(bars, labels=[f"{value:,}" for value in counts.values], padding=2, fontsize=6.3)
    ax.set_xlabel("Unique proteins")
    clean(ax, "x")

    ax = axes[2]
    panel(ax, "C", "Transmembrane topology")
    bins = pd.cut(
        protein["tm"],
        bins=[-0.1, 0.5, 1.5, 6.5, 7.5, 12.5, np.inf],
        labels=["0", "1", "2–6", "7", "8–12", ">12"],
    )
    tm_counts = bins.value_counts().sort_index()
    colors = [COLORS["C"], COLORS["B"], "#69B3A2", COLORS["A"], "#4C78A8", "#2F5597"]
    bars = ax.bar(tm_counts.index.astype(str), tm_counts.values, color=colors)
    ax.bar_label(bars, labels=[f"{value:,}" for value in tm_counts.values], fontsize=6.3)
    ax.set_ylabel("Unique proteins")
    ax.set_xlabel("Annotated transmembrane helices")
    clean(ax, "y")

    ax = axes[3]
    panel(ax, "D", "Sequence length versus TM count")
    valid = protein[["length", "tm"]].dropna()
    ax.hexbin(
        valid["length"],
        valid["tm"],
        gridsize=42,
        mincnt=1,
        bins="log",
        cmap=sns.light_palette(COLORS["A"], as_cmap=True),
    )
    rho, p_value = spearmanr(valid["length"], valid["tm"])
    ax.text(
        0.98,
        0.97,
        f"Spearman ρ = {rho:.2f}\nn = {len(valid):,}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#D9D9D9"),
    )
    ax.set_xlabel("Sequence length (aa)")
    ax.set_ylabel("TM helices")
    clean(ax, "y")

    ax = axes[4]
    panel(ax, "E", "Structure coverage by functional class")
    top_classes = protein["functional_class"].value_counts().head(9).index
    subset = protein[protein["functional_class"].isin(top_classes)]
    struct = (
        pd.crosstab(subset["functional_class"], subset["structure_grade"], normalize="index")
        .reindex(top_classes)
        .fillna(0)
    )
    struct = struct.reindex(
        columns=["S2 membrane PDB", "S3 other PDB", "S4 AlphaFold only", "S0 no structure"],
        fill_value=0,
    )
    sns.heatmap(
        struct * 100,
        cmap=sns.light_palette(COLORS["blue"], as_cmap=True),
        annot=True,
        fmt=".0f",
        linewidths=0.4,
        cbar_kws={"label": "% of class"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticklabels(["Membrane\nPDB", "Other\nPDB", "AlphaFold\nonly", "No\nstructure"], rotation=0)
    ax.set_yticklabels([x.replace("_", " ") for x in struct.index], rotation=0)

    ax = axes[5]
    panel(ax, "F", "Oligomeric and assembly states")
    def oligomer_group(value: str) -> str:
        value = str(value)
        if value == "unknown":
            return "Unknown"
        if "multiple" in value:
            return "State-dependent"
        if "monomer" in value:
            return "Monomer"
        if "dimer" in value:
            return "Dimer"
        if "trimer" in value:
            return "Trimer"
        if "tetramer" in value:
            return "Tetramer"
        return "Larger / unresolved"

    protein["oligomer_group"] = protein["oligomeric_state_consensus_v62"].fillna("unknown").map(oligomer_group)
    oligo = pd.crosstab(protein["oligomer_group"], protein["subunit_evidence_grade_v62"].fillna("SA0_unresolved"))
    order = ["Monomer", "Dimer", "Trimer", "Tetramer", "Larger / unresolved", "State-dependent", "Unknown"]
    oligo = oligo.reindex(order, fill_value=0)
    oligo_colors = ["#2F5597", "#5B9BD5", "#70AD47", "#ED7D31", "#B7B7B7"]
    left = np.zeros(len(oligo))
    for color, grade in zip(oligo_colors, oligo.columns):
        ax.barh(oligo.index, oligo[grade], left=left, color=color, label=grade.replace("_", " "))
        left += oligo[grade].to_numpy()
    ax.set_xlabel("Unique proteins")
    ax.legend(frameon=False, fontsize=5.5, loc="lower center", bbox_to_anchor=(0.5, -0.34), ncol=2)
    clean(ax, "x")
    paths = save(fig, "M2_membrane_proteome_landscape")
    write_table(cross.reset_index(), "M2_class_evidence_counts.tsv")
    write_table(struct.reset_index(), "M2_structure_coverage.tsv")
    return paths


def expression_matrices() -> tuple[pd.DataFrame, pd.DataFrame]:
    top_classes = protein["functional_class"].value_counts().head(8).index
    tissue_work = tissue[tissue["functional_class"].isin(top_classes)].copy()
    tissue_work["log_value"] = np.log1p(tissue_work["value"])
    tissue_matrix = tissue_work.pivot_table(
        index="functional_class", columns="tissue", values="log_value", aggfunc="median"
    ).reindex(top_classes)
    tissue_select = tissue_matrix.var(axis=0).sort_values(ascending=False).head(12).index
    tissue_matrix = tissue_matrix[tissue_select]

    cell_work = cell[cell["functional_class"].isin(top_classes)].copy()
    cell_work["log_value"] = np.log1p(cell_work["value"])
    cell_matrix = cell_work.pivot_table(
        index="functional_class", columns="cell_type", values="log_value", aggfunc="median"
    ).reindex(top_classes)
    cell_select = cell_matrix.var(axis=0).sort_values(ascending=False).head(12).index
    cell_matrix = cell_matrix[cell_select]
    return tissue_matrix, cell_matrix


tissue_matrix, cell_matrix = expression_matrices()


def row_zscore(frame: pd.DataFrame) -> pd.DataFrame:
    means = frame.mean(axis=1)
    stds = frame.std(axis=1).replace(0, np.nan)
    return frame.sub(means, axis=0).div(stds, axis=0).fillna(0)


def m3_expression_localization() -> dict:
    fig, axes = make_figure("M3  Tissue, cell-type and subcellular localization atlas")
    ax = axes[0]
    panel(ax, "A", "Anatomical expression coverage")
    tissue["organ_system"] = tissue["tissue"].map(map_system)
    detected = tissue.loc[tissue["value"] >= 1]
    system_counts = detected.groupby("organ_system")["target_uniprot_id"].nunique().to_dict()
    draw_human_navigation(ax, system_counts, color=COLORS["blue"])

    ax = axes[1]
    panel(ax, "B", "Tissue expression patterns")
    z = row_zscore(tissue_matrix)
    col_order = clustered_order(z, axis=1)
    z = z.iloc[:, col_order]
    sns.heatmap(
        z,
        cmap="vlag",
        center=0,
        cbar_kws={"label": "Within-class z-score"},
        linewidths=0.35,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_yticklabels([x.replace("_", " ") for x in z.index], rotation=0)

    ax = axes[2]
    panel(ax, "C", "Cell-type expression patterns")
    zc = row_zscore(cell_matrix)
    col_order = clustered_order(zc, axis=1)
    zc = zc.iloc[:, col_order]
    sns.heatmap(
        zc,
        cmap="vlag",
        center=0,
        cbar_kws={"label": "Within-class z-score"},
        linewidths=0.35,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_yticklabels([x.replace("_", " ") for x in zc.index], rotation=0)

    ax = axes[3]
    panel(ax, "D", "Expression breadth")
    top_classes = ["enzyme", "receptor", "transporter", "ion_channel"]
    for name, color in zip(top_classes, [COLORS["A"], COLORS["B"], COLORS["C"], "#6A3D9A"]):
        values = pd.to_numeric(
            protein.loc[protein["functional_class"].eq(name), "hpa_tissue_detected_count_v62"],
            errors="coerce",
        )
        x, y = ecdf(values.to_numpy())
        ax.plot(x, y, label=name.replace("_", " "), color=color, lw=1.6)
    ax.set_xlabel("Detected tissues per protein")
    ax.set_ylabel("Cumulative fraction")
    ax.legend(frameon=False)
    clean(ax, "both")

    ax = axes[4]
    panel(ax, "E", "Subcellular localization")
    loc = localization.drop_duplicates(["target_uniprot_id", "location_term"])
    loc_counts = loc["location_term"].value_counts().head(12).sort_values()
    bars = ax.barh(loc_counts.index, loc_counts.values, color=COLORS["blue"])
    ax.bar_label(bars, labels=[f"{v:,}" for v in loc_counts.values], padding=2, fontsize=6.2)
    ax.set_xlabel("Unique proteins (non-exclusive)")
    clean(ax, "x")

    ax = axes[5]
    panel(ax, "F", "RNA–IHC concordance at protein–tissue level")
    rna_binary = tissue[["target_uniprot_id", "tissue", "value"]].copy()
    rna_binary["tissue_key"] = rna_binary["tissue"].str.lower().str.strip()
    rna_binary["RNA detected"] = rna_binary["value"].ge(1)
    rna_binary = (
        rna_binary.groupby(["target_uniprot_id", "tissue_key"], as_index=False)["RNA detected"]
        .max()
    )
    ihc_binary = ihc.copy()
    ihc_binary["tissue_key"] = ihc_binary["tissue"].str.lower().str.strip()
    ihc_binary["IHC detected"] = ~ihc_binary["detection_status"].fillna("").eq("not_detected")
    ihc_binary = (
        ihc_binary.groupby(["target_uniprot_id", "tissue_key"], as_index=False)["IHC detected"]
        .max()
    )
    concordance = rna_binary.merge(
        ihc_binary, on=["target_uniprot_id", "tissue_key"], how="inner"
    )
    confusion = pd.crosstab(
        concordance["RNA detected"],
        concordance["IHC detected"],
    ).reindex(index=[False, True], columns=[False, True], fill_value=0)
    sns.heatmap(
        confusion,
        cmap=sns.light_palette(COLORS["B"], as_cmap=True),
        annot=confusion.map(lambda value: f"{value:,.0f}"),
        fmt="",
        linewidths=1,
        linecolor="white",
        cbar_kws={"label": "Protein–tissue observations"},
        ax=ax,
    )
    ax.set_xticklabels(["Not detected", "Detected"], rotation=0)
    ax.set_yticklabels(["Not detected", "Detected"], rotation=0)
    ax.set_xlabel("IHC")
    ax.set_ylabel("RNA (nTPM ≥ 1)")
    observed_agreement = (
        confusion.loc[False, False] + confusion.loc[True, True]
    ) / confusion.to_numpy().sum()
    ax.text(
        0.98,
        -0.22,
        f"Observed agreement = {observed_agreement:.1%}; n = {len(concordance):,}",
        transform=ax.transAxes,
        ha="right",
        fontsize=6.8,
    )
    paths = save(fig, "M3_expression_localization_atlas")
    write_table(tissue_matrix.reset_index(), "M3_tissue_expression_matrix.tsv")
    write_table(cell_matrix.reset_index(), "M3_cell_type_expression_matrix.tsv")
    return paths


def m4_disease_landscape() -> dict:
    fig, axes = make_figure("M4  Disease-association landscape")
    ax = axes[0]
    panel(ax, "A", "Disease associations by organ system")
    system_pair_counts = (
        disease.drop_duplicates(["target_uniprot_id", "disease_id"])
        .groupby("organ_system")
        .size()
        .to_dict()
    )
    draw_human_navigation(ax, system_pair_counts, color="#9C7DBA")

    ax = axes[1]
    panel(ax, "B", "Disease systems and evidence grade")
    disease_unique = disease.drop_duplicates(["target_uniprot_id", "disease_id"])
    cross = pd.crosstab(disease_unique["organ_system"], disease_unique["disease_evidence_level"])
    systems = cross.sum(axis=1).sort_values(ascending=False).head(10).index
    cross = cross.reindex(systems)
    order = ["very_high", "high", "medium"]
    left = np.zeros(len(cross))
    colors = ["#2F5597", "#70AD47", "#ED7D31"]
    for evidence, color in zip(order, colors):
        values = cross.get(evidence, pd.Series(0, index=cross.index)).to_numpy()
        ax.barh(cross.index, values, left=left, label=evidence.replace("_", " "), color=color)
        left += values
    ax.invert_yaxis()
    ax.set_xlabel("Unique protein–disease pairs")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.28))
    clean(ax, "x")

    ax = axes[2]
    panel(ax, "C", "Evidence architecture by disease system")
    flags = [
        "human_genetic_flag",
        "clinical_genetic_flag",
        "somatic_mutation_flag",
        "functional_flag",
        "animal_model_flag",
        "expression_flag",
        "literature_flag",
        "uniprot_disease_flag",
    ]
    evidence_matrix = disease_unique.groupby("organ_system")[flags].mean().reindex(systems)
    evidence_matrix.columns = [
        "Human\ngenetic",
        "Clinical\ngenetic",
        "Somatic",
        "Functional",
        "Animal",
        "Expression",
        "Literature",
        "UniProt",
    ]
    sns.heatmap(
        evidence_matrix * 100,
        cmap=sns.light_palette("#6A3D9A", as_cmap=True),
        annot=True,
        fmt=".0f",
        cbar_kws={"label": "% of pairs"},
        linewidths=0.35,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    ax = axes[3]
    panel(ax, "D", "High-confidence protein–disease network")
    high = disease_unique[disease_unique["disease_evidence_level"].isin(["very_high", "high"])]
    top_diseases = high["disease_name"].value_counts().head(9).index
    top_proteins = high["approved_symbol"].value_counts().head(9).index
    edges = high[
        high["disease_name"].isin(top_diseases) & high["approved_symbol"].isin(top_proteins)
    ][["approved_symbol", "disease_name"]].drop_duplicates()
    graph = nx.Graph()
    graph.add_nodes_from(top_proteins, bipartite=0)
    graph.add_nodes_from(top_diseases, bipartite=1)
    graph.add_edges_from(edges.itertuples(index=False, name=None))
    positions = {}
    for idx, node in enumerate(top_proteins):
        positions[node] = (0.1, 0.9 - idx * 0.1)
    for idx, node in enumerate(top_diseases):
        positions[node] = (0.9, 0.9 - idx * 0.1)
    nx.draw_networkx_edges(graph, positions, ax=ax, alpha=0.24, width=0.7, edge_color="#7F8C8D")
    nx.draw_networkx_nodes(graph, positions, nodelist=list(top_proteins), node_color=COLORS["blue"], node_size=90, ax=ax)
    nx.draw_networkx_nodes(graph, positions, nodelist=list(top_diseases), node_color="#B39DDB", node_size=90, ax=ax)
    for node in top_proteins:
        ax.text(0.13, positions[node][1], node, va="center", fontsize=5.8)
    for node in top_diseases:
        ax.text(0.87, positions[node][1], str(node)[:24], va="center", ha="right", fontsize=5.4)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax = axes[4]
    panel(ax, "E", "Disease burden per protein")
    degree = disease_unique.groupby("target_uniprot_id").size()
    protein_degree = protein[["target_uniprot_id", "functional_class"]].copy()
    protein_degree["disease_count"] = protein_degree["target_uniprot_id"].map(degree).fillna(0)
    for name, color in zip(
        ["receptor", "enzyme", "transporter", "ion_channel"],
        [COLORS["A"], COLORS["B"], COLORS["C"], "#6A3D9A"],
    ):
        x, y = ecdf(protein_degree.loc[protein_degree["functional_class"].eq(name), "disease_count"])
        ax.step(x + 1, 1 - y + 1 / max(len(y), 1), where="post", label=name, color=color, lw=1.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Disease relations per protein + 1")
    ax.set_ylabel("Complementary cumulative fraction")
    ax.legend(frameon=False)
    clean(ax, "both")

    ax = axes[5]
    panel(ax, "F", "Expression–disease organ concordance")
    expr = tissue[["target_uniprot_id", "tissue", "value"]].copy()
    expr["organ_system"] = expr["tissue"].map(map_system)
    expr_profile = (
        expr.groupby(["target_uniprot_id", "organ_system"], as_index=False)["value"]
        .median()
    )
    expr_profile = expr_profile.loc[expr_profile["value"] >= 1].copy()
    expr_profile["max_value"] = expr_profile.groupby("target_uniprot_id")["value"].transform("max")
    expr = expr_profile.loc[expr_profile["value"].eq(expr_profile["max_value"])]
    expr_sets = {
        system: set(group["target_uniprot_id"])
        for system, group in expr.groupby("organ_system")
    }
    disease_sets = {
        system: set(group["target_uniprot_id"])
        for system, group in disease_unique.groupby("organ_system")
    }
    common_systems = [
        system
        for system in SYSTEM_KEYWORDS
        if system in expr_sets and system in disease_sets
    ]
    jaccard = pd.DataFrame(0.0, index=common_systems, columns=common_systems)
    for row_system in common_systems:
        for col_system in common_systems:
            left_set = expr_sets[row_system]
            right_set = disease_sets[col_system]
            union = left_set | right_set
            jaccard.loc[row_system, col_system] = len(left_set & right_set) / len(union) if union else 0
    sns.heatmap(
        jaccard,
        cmap=sns.light_palette("#6A3D9A", as_cmap=True),
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Jaccard index"},
        linewidths=0.35,
        ax=ax,
    )
    ax.set_xlabel("Disease organ system")
    ax.set_ylabel("Expression organ system")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    paths = save(fig, "M4_disease_association_landscape")
    write_table(cross.reset_index(), "M4_disease_system_evidence_counts.tsv")
    write_table(evidence_matrix.reset_index(), "M4_disease_evidence_matrix.tsv")
    write_table(jaccard.reset_index(), "M4_expression_disease_jaccard.tsv")
    return paths


outputs = {
    "M1": m1_database_architecture(),
    "M2": m2_membrane_proteome(),
    "M3": m3_expression_localization(),
    "M4": m4_disease_landscape(),
}

caption = f"""# Main figure captions — M1–M4

## Figure M1. Database architecture, integration and quality control

Sources are grouped by the data module to which they contribute. The pipeline
distinguishes record-level identity resolution from pair-level aggregation.
The incremental positive review input comprised 3,290,514 records; the negative
identity input comprised 10,246,948 records. Counts in the release snapshot are
from the frozen MemPro V6.2 release dated 2026-07-30. Heatmap values are
log10(count + 1), with raw counts printed in cells.

## Figure M2. Landscape of the human membrane proteome

All panels use one row per UniProt accession (n = {len(protein):,}). Membrane
classes A, B and C are biological scope categories; E1–E3 are evidence grades
and should not be interpreted as the same variable. Structure grade separates
membrane-curated experimental structures, other experimental PDB entries,
AlphaFold-only models and proteins without an assigned structure. Oligomeric
states marked unresolved or state-dependent remain explicit.

## Figure M3. Tissue, cell-type and subcellular localization atlas

HPA tissue and cell-type values were log1p transformed and summarized by
functional class. Heatmap colors are within-class z-scores and therefore show
relative patterns rather than absolute expression. The anatomical silhouette is
a navigation schematic; quantitative counts are shown in the adjacent bars.
RNA–IHC concordance is evaluated as detection breadth per protein, not as a
claim of assay equivalence.

## Figure M4. Disease-association landscape

Disease analyses use unique protein–disease pairs. Organ systems are a
transparent keyword-derived grouping of disease names for visualization; they
are not prevalence estimates or a replacement for the original MONDO/OMIM
identifiers. The network is filtered to high-confidence, high-degree nodes for
readability. Expression–disease concordance is descriptive Jaccard overlap and
does not imply tissue causality.
"""
(CAPTION_DIR / "CAPTIONS_M1_M4.md").write_text(caption, encoding="utf-8")
write_json(
    {
        "status": "PASS",
        "figures": outputs,
        "protein_rows": int(len(protein)),
        "disease_rows": int(len(disease)),
        "positive_pair_rows": int(len(pair)),
        "site_rows": int(len(site)),
        "tissue_rows_used": int(len(tissue)),
        "cell_type_rows_used": int(len(cell)),
    },
    "FIGURES_M1_M4_VALIDATION.json",
)
print(json.dumps(outputs, ensure_ascii=False, indent=2))
