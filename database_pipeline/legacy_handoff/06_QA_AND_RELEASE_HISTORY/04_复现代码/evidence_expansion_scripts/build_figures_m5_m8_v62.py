from __future__ import annotations

import gzip
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA

from figure_style_v62 import (
    CAPTION_DIR,
    COLORS,
    DATA_DIR,
    clean,
    draw_flow_boxes,
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
DOCK = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\docking_selection_v2_v62_20260730"
)
V61_REPORT = Path(
    r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_1_20260729\V61_VALIDATION_REPORT.json"
)


def truth(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"1", "true", "yes", "y"})


def human_number(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return f"{value:.0f}"


protein = pd.read_csv(
    REL / "human_membrane_protein_master_v6_2.tsv",
    sep="\t",
    usecols=[
        "target_uniprot_id",
        "functional_primary_class_v5",
        "membrane_class_v52",
        "evidence_level_v52",
        "pdb_ids",
        "alphafolddb_ids",
        "opm_present_v5",
        "pdbtm_present_v5",
    ],
    low_memory=False,
)
protein["functional_class"] = (
    protein["functional_primary_class_v5"].fillna("unclassified").astype(str)
)
class_map = protein.set_index("target_uniprot_id")["functional_class"]

pair = pd.read_csv(
    REL / "protein_compound_summary_v6_2.tsv.gz",
    sep="\t",
    compression="gzip",
    usecols=[
        "target_uniprot_id",
        "compound_internal_id",
        "best_binding_evidence_level",
        "BE1_evidence_count",
        "BE2_evidence_count",
        "BE3_evidence_count",
        "binding_evidence_count",
        "independent_source_count",
        "best_standard_value_nM",
    ],
    low_memory=False,
)
pair["functional_class"] = pair["target_uniprot_id"].map(class_map).fillna("unclassified")
pair_compounds = set(pair["compound_internal_id"].dropna().astype(str))

site = pd.read_csv(
    REL / "binding_site_instances_v6_2.tsv.gz",
    sep="\t",
    compression="gzip",
    usecols=[
        "target_uniprot_id",
        "source_database",
        "site_type",
        "residue_or_site_description",
        "pdb_ids",
    ],
    low_memory=False,
)


compound_columns = [
    "compound_internal_id",
    "compound_scope_status",
    "record_qc_status",
    "computed_structural_class",
    "standard_smiles",
    "molecular_weight",
    "xlogp",
    "tpsa",
    "hbond_donor_count",
    "hbond_acceptor_count",
    "rotatable_bond_count",
    "formal_charge",
    "form_count",
    "source_database_count",
    "protein_target_count",
    "best_binding_evidence_level",
    "is_approved_drug",
    "is_clinical_candidate",
    "is_endogenous_ligand",
    "is_natural_product",
    "is_chemical_probe",
]


def scan_compounds():
    class_counts = Counter()
    status_combos = Counter()
    chem_class_map: dict[str, str] = {}
    sampled = []
    target_hist = {
        "Approved": Counter(),
        "Clinical": Counter(),
        "Probe": Counter(),
        "Endogenous": Counter(),
        "Other": Counter(),
    }
    missing = Counter()
    total = 0
    core_ok_total = 0
    for chunk in pd.read_csv(
        REL / "small_molecule_master_v1_3.tsv",
        sep="\t",
        usecols=compound_columns,
        chunksize=150_000,
        low_memory=False,
    ):
        total += len(chunk)
        core_ok = chunk["compound_scope_status"].eq("core") & chunk["record_qc_status"].eq("ok")
        filtered = chunk.loc[core_ok].copy()
        core_ok_total += len(filtered)
        structural_class = filtered["computed_structural_class"].fillna("unclassified").astype(str)
        class_counts.update(structural_class)
        ids = chunk["compound_internal_id"].astype(str)
        needed = ids.isin(pair_compounds)
        chem_class_map.update(
            zip(
                ids[needed],
                chunk.loc[needed, "computed_structural_class"]
                .fillna("unclassified")
                .astype(str),
            )
        )
        approved = truth(filtered["is_approved_drug"])
        clinical = truth(filtered["is_clinical_candidate"])
        endogenous = truth(filtered["is_endogenous_ligand"])
        natural = truth(filtered["is_natural_product"])
        probe = truth(filtered["is_chemical_probe"])
        combo = (
            approved.astype(int).astype(str)
            + clinical.astype(int).astype(str)
            + probe.astype(int).astype(str)
            + endogenous.astype(int).astype(str)
            + natural.astype(int).astype(str)
        )
        status_combos.update(combo)
        targets = pd.to_numeric(filtered["protein_target_count"], errors="coerce").fillna(0).astype(int)
        category_masks = {
            "Approved": approved,
            "Clinical": clinical & ~approved,
            "Probe": probe & ~approved & ~clinical,
            "Endogenous": endogenous & ~approved & ~clinical & ~probe,
            "Other": ~(approved | clinical | probe | endogenous),
        }
        for category, mask in category_masks.items():
            target_hist[category].update(targets[mask])
        for field in [
            "standard_smiles",
            "molecular_weight",
            "xlogp",
            "tpsa",
            "rotatable_bond_count",
            "formal_charge",
            "computed_structural_class",
        ]:
            missing[field] += int(filtered[field].isna().sum() + filtered[field].astype(str).eq("").sum())
        hash_values = pd.util.hash_pandas_object(
            filtered["compound_internal_id"].astype(str), index=False
        )
        sample_mask = hash_values.mod(160).eq(0) & filtered["standard_smiles"].notna()
        sampled.append(filtered.loc[sample_mask])
    sample = pd.concat(sampled, ignore_index=True)
    if len(sample) > 14_000:
        sample = sample.sample(14_000, random_state=42)
    return {
        "class_counts": class_counts,
        "status_combos": status_combos,
        "chem_class_map": chem_class_map,
        "sample": sample,
        "target_hist": target_hist,
        "missing": missing,
        "total": total,
        "core_ok_total": core_ok_total,
    }


compound_scan = scan_compounds()
compound_sample = compound_scan["sample"]
pair["chemical_class"] = pair["compound_internal_id"].map(
    compound_scan["chem_class_map"]
).fillna("unclassified")


def chemical_embedding(sample: pd.DataFrame) -> pd.DataFrame:
    cache = DATA_DIR / "M5_morgan_umap_coordinates.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    import numpy as np
    sys.path.append(
        r"D:\7.22\evidence_expansion_v2_working\tools\python_v61"
    )
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    records = []
    fingerprints = []
    selected = sample.sample(min(len(sample), 7_500), random_state=42)
    for row in selected.itertuples(index=False):
        smiles = getattr(row, "standard_smiles")
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            continue
        fingerprint = generator.GetFingerprint(mol)
        array = np.zeros((1024,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fingerprint, array)
        fingerprints.append(array)
        records.append(
            {
                "compound_internal_id": getattr(row, "compound_internal_id"),
                "computed_structural_class": getattr(row, "computed_structural_class")
                if pd.notna(getattr(row, "computed_structural_class"))
                else "unclassified",
            }
        )
    matrix = np.asarray(fingerprints, dtype=np.uint8)
    try:
        import umap

        embedding = umap.UMAP(
            n_neighbors=25,
            min_dist=0.12,
            n_components=2,
            metric="jaccard",
            random_state=42,
            n_jobs=1,
            low_memory=True,
        ).fit_transform(matrix)
        method = "Morgan UMAP"
    except Exception:
        embedding = PCA(n_components=2, random_state=42).fit_transform(matrix)
        method = "Morgan PCA fallback"
    result = pd.DataFrame(records)
    result["x"] = embedding[:, 0]
    result["y"] = embedding[:, 1]
    result["method"] = method
    result.to_csv(cache, sep="\t", index=False)
    return result


embedding = chemical_embedding(compound_sample)


def evidence_aggregates():
    source_tier = Counter()
    source_counts = Counter()
    for chunk in pd.read_csv(
        REL / "binding_evidence_master_v6_2.tsv.gz",
        sep="\t",
        compression="gzip",
        usecols=["source_database", "evidence_tier"],
        chunksize=300_000,
        low_memory=False,
    ):
        sources = chunk["source_database"].fillna("Unknown").astype(str)
        tiers = chunk["evidence_tier"].fillna("Unknown").astype(str)
        source_counts.update(sources)
        source_tier.update(zip(sources, tiers))
    matrix = pd.Series(source_tier).rename("count").reset_index()
    matrix.columns = ["source", "tier", "count"]
    pivot = matrix.pivot_table(index="source", columns="tier", values="count", fill_value=0)
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).head(12).index]
    return pivot, pd.Series(source_counts).sort_values(ascending=False)


source_tier, evidence_source_counts = evidence_aggregates()


AA3 = {
    "ALA": "Hydrophobic",
    "VAL": "Hydrophobic",
    "LEU": "Hydrophobic",
    "ILE": "Hydrophobic",
    "MET": "Hydrophobic",
    "PRO": "Hydrophobic",
    "PHE": "Aromatic",
    "TYR": "Aromatic",
    "TRP": "Aromatic",
    "SER": "Polar",
    "THR": "Polar",
    "ASN": "Polar",
    "GLN": "Polar",
    "CYS": "Polar",
    "GLY": "Special",
    "LYS": "Positive",
    "ARG": "Positive",
    "HIS": "Positive",
    "ASP": "Negative",
    "GLU": "Negative",
}
AA1 = {
    "A": "Hydrophobic",
    "V": "Hydrophobic",
    "L": "Hydrophobic",
    "I": "Hydrophobic",
    "M": "Hydrophobic",
    "P": "Hydrophobic",
    "F": "Aromatic",
    "Y": "Aromatic",
    "W": "Aromatic",
    "S": "Polar",
    "T": "Polar",
    "N": "Polar",
    "Q": "Polar",
    "C": "Polar",
    "G": "Special",
    "K": "Positive",
    "R": "Positive",
    "H": "Positive",
    "D": "Negative",
    "E": "Negative",
}


def residue_counts() -> Counter:
    counts = Counter()
    for description in site["residue_or_site_description"].dropna().astype(str):
        for residue in re.findall(r"\b([A-Z]{3})\d+\b", description.upper()):
            if residue in AA3:
                counts[AA3[residue]] += 1
        for residue in re.findall(r"\b([ACDEFGHIKLMNPQRSTVWY])\d+\b", description.upper()):
            counts[AA1[residue]] += 1
    return counts


residues = residue_counts()


def negative_aggregates():
    class_counts = Counter()
    conflict_class_counts = Counter()
    conflicts = []
    for chunk in pd.read_csv(
        REL / "protein_compound_negative_summary_v6_2.tsv.gz",
        sep="\t",
        compression="gzip",
        chunksize=500_000,
        low_memory=False,
    ):
        classes = chunk["target_uniprot_id"].map(class_map).fillna("unclassified")
        class_counts.update(classes)
        flag = pd.to_numeric(chunk["positive_negative_conflict_flag"], errors="coerce").fillna(0).astype(bool)
        conflict_class_counts.update(classes[flag])
        if flag.any():
            subset = chunk.loc[
                flag,
                [
                    "target_uniprot_id",
                    "compound_internal_id",
                    "negative_evidence_count",
                ],
            ].copy()
            conflicts.append(subset)
    return class_counts, conflict_class_counts, pd.concat(conflicts, ignore_index=True)


negative_class_counts, conflict_class_counts, conflict_pairs = negative_aggregates()
conflict_context = conflict_pairs.merge(
    pair[
        [
            "target_uniprot_id",
            "compound_internal_id",
            "best_binding_evidence_level",
            "functional_class",
        ]
    ],
    on=["target_uniprot_id", "compound_internal_id"],
    how="left",
)


def load_docking():
    full = pd.read_csv(
        DOCK / "staging" / "docking_full_nonconflict_v2_v62.tsv.gz",
        sep="\t",
        compression="gzip",
        low_memory=False,
    )
    standard = pd.read_csv(
        DOCK / "staging" / "docking_hpc_standard_v2_v62.tsv.gz",
        sep="\t",
        compression="gzip",
        low_memory=False,
    )
    pilot = pd.read_csv(
        DOCK / "staging" / "docking_hpc_pilot_v2_v62.tsv",
        sep="\t",
        low_memory=False,
    )
    return full, standard, pilot


docking_full, docking_standard, docking_pilot = load_docking()
docking_standard["functional_class"] = docking_standard["target_uniprot_id"].map(
    class_map
).fillna("unclassified")


def m5_chemical_space() -> dict:
    fig, axes = make_figure("M5  Chemical identity and small-molecule space")

    ax = axes[0]
    panel(ax, "A", "Canonical parent and form hierarchy")
    draw_flow_boxes(
        ax,
        ["Source\nrecords", "Exact\nstructure", "Canonical\nparent", "Parent / form\nhierarchy", "QC / review"],
        counts=[10_246_948, 1_621_892, 2_016_064, 2_020_898, 444_309],
        colors=["#F2F2F2", "#D9EAF7", "#D9EAD3", "#FFF2CC", "#F4CCCC"],
    )

    ax = axes[1]
    panel(ax, "B", "Biological-status intersections")
    names = ["Approved", "Clinical", "Probe", "Endogenous", "Natural"]
    combo_labels = []
    combo_values = []
    for bits, count in compound_scan["status_combos"].most_common(9):
        active = [name for name, bit in zip(names, bits) if bit == "1"]
        combo_labels.append(" + ".join(active) if active else "Other")
        combo_values.append(count)
    order = np.arange(len(combo_labels))[::-1]
    bars = ax.barh(order, combo_values[::-1], color=COLORS["blue"])
    ax.set_yticks(order, combo_labels[::-1])
    ax.bar_label(bars, labels=[human_number(v) for v in combo_values[::-1]], padding=2, fontsize=6.2)
    ax.set_xlabel("Canonical compounds")
    ax.set_xscale("log")
    clean(ax, "x")

    ax = axes[2]
    panel(ax, "C", "Computed structural classes")
    class_counts = pd.Series(compound_scan["class_counts"]).sort_values(ascending=False).head(12).sort_values()
    bars = ax.barh(class_counts.index.str.replace("_", " "), class_counts.values, color=COLORS["blue"])
    ax.bar_label(bars, labels=[human_number(v) for v in class_counts.values], padding=2, fontsize=6.2)
    ax.set_xlabel("Core QC-passed canonical compounds")
    clean(ax, "x")

    ax = axes[3]
    panel(ax, "D", "Physicochemical landscape")
    ax.axis("off")
    properties = [
        ("molecular_weight", "MW (Da)", (0, 900)),
        ("xlogp", "XlogP", (-6, 12)),
        ("tpsa", "TPSA (Å²)", (0, 300)),
        ("rotatable_bond_count", "Rotatable bonds", (0, 35)),
    ]
    for index, (field, label, limits) in enumerate(properties):
        inset = ax.inset_axes([0.04 + (index % 2) * 0.49, 0.55 - (index // 2) * 0.5, 0.44, 0.38])
        values = pd.to_numeric(compound_sample[field], errors="coerce")
        values = values[values.between(*limits)]
        x, y = ecdf(values.to_numpy())
        inset.plot(x, y, color=COLORS["A"], lw=1.4)
        inset.set_xlabel(label, fontsize=6.4)
        inset.set_ylabel("ECDF", fontsize=6.2)
        clean(inset, "both")

    ax = axes[4]
    panel(ax, "E", "Morgan-fingerprint chemical space")
    top_classes = (
        embedding["computed_structural_class"].value_counts().head(6).index
    )
    plot_data = embedding.copy()
    plot_data["class_display"] = np.where(
        plot_data["computed_structural_class"].isin(top_classes),
        plot_data["computed_structural_class"],
        "Other",
    )
    palette = sns.color_palette("colorblind", n_colors=len(top_classes))
    palette_map = dict(zip(top_classes, palette))
    palette_map["Other"] = "#C7C7C7"
    for class_name, group in plot_data.groupby("class_display"):
        ax.scatter(
            group["x"],
            group["y"],
            s=3,
            alpha=0.45 if class_name != "Other" else 0.18,
            color=palette_map[class_name],
            label=str(class_name).replace("_", " "),
            rasterized=True,
        )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(frameon=False, markerscale=2, fontsize=5.2, ncol=2)
    clean(ax, "both")

    ax = axes[5]
    panel(ax, "F", "Target promiscuity")
    status_colors = {
        "Approved": COLORS["A"],
        "Clinical": COLORS["B"],
        "Probe": COLORS["C"],
        "Endogenous": "#6A3D9A",
        "Other": "#7F7F7F",
    }
    for category, histogram in compound_scan["target_hist"].items():
        if not histogram:
            continue
        counts = pd.Series(histogram, dtype=float).sort_index()
        x = counts.index.to_numpy(dtype=float) + 1
        y = counts.iloc[::-1].cumsum().iloc[::-1].to_numpy() / counts.sum()
        ax.step(
            x,
            y,
            where="post",
            label=category,
            color=status_colors[category],
            lw=1.4,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Supported protein targets + 1")
    ax.set_ylabel("Complementary cumulative fraction")
    ax.legend(frameon=False)
    clean(ax, "both")
    paths = save(fig, "M5_chemical_identity_space")
    write_table(
        pd.Series(compound_scan["class_counts"], name="count")
        .rename_axis("structural_class")
        .reset_index(),
        "M5_structural_class_counts.tsv",
    )
    write_table(embedding, "M5_morgan_embedding.tsv")
    return paths


def m6_binding_evidence() -> dict:
    fig, axes = make_figure("M6  Binding evidence and binding-site landscape")
    ax = axes[0]
    panel(ax, "A", "Evidence source and tier")
    tier_order = ["BE1", "BE2", "BE3"]
    source_tier_plot = source_tier.reindex(columns=tier_order, fill_value=0)
    left = np.zeros(len(source_tier_plot))
    for tier in tier_order:
        values = source_tier_plot[tier].to_numpy()
        ax.barh(
            source_tier_plot.index,
            values,
            left=left,
            color=COLORS[tier],
            label=tier,
        )
        left += values
    ax.invert_yaxis()
    ax.set_xlabel("Evidence records")
    ax.set_xscale("log")
    ax.legend(frameon=False, ncol=3)
    clean(ax, "x")

    ax = axes[1]
    panel(ax, "B", "Independent source support")
    source_dist = pair["independent_source_count"].value_counts().sort_index()
    bars = ax.bar(source_dist.index.astype(str), source_dist.values, color=COLORS["blue"])
    ax.set_yscale("log")
    ax.set_xlabel("Independent sources per pair")
    ax.set_ylabel("Unique positive pairs (log)")
    for bar, value in zip(bars, source_dist.values):
        if value > 1000:
            ax.text(bar.get_x() + bar.get_width() / 2, value * 1.1, human_number(value), ha="center", fontsize=6)
    clean(ax, "y")

    ax = axes[2]
    panel(ax, "C", "Quantitative potency by evidence tier")
    potency = pair.copy()
    potency["best_standard_value_nM"] = pd.to_numeric(
        potency["best_standard_value_nM"], errors="coerce"
    )
    potency = potency[
        potency["best_standard_value_nM"].between(1e-3, 1e8)
    ]
    for tier in ["BE1", "BE2", "BE3"]:
        values = potency.loc[
            potency["best_binding_evidence_level"].eq(tier),
            "best_standard_value_nM",
        ]
        if values.empty:
            continue
        x, y = ecdf(values.to_numpy())
        ax.plot(x, y, label=f"{tier} (n={len(values):,})", color=COLORS[tier], lw=1.5)
    ax.set_xscale("log")
    ax.set_xlabel("Best standardized value (nM; log scale)")
    ax.set_ylabel("Cumulative fraction")
    ax.legend(frameon=False)
    clean(ax, "both")

    ax = axes[3]
    panel(ax, "D", "Binding-site readiness")
    readiness_map = {
        "experimental_structure_residue_contact": "PDB + residue contacts",
        "experimental_structure_match": "PDB structure match",
        "experimental_complex_pocket": "Experimental pocket",
        "bindingdb_ligand_target_complex": "Complex PDB only",
        "uniprot_curated_binding_site": "Curated site annotation",
    }
    readiness = site["site_type"].map(readiness_map).fillna("Other").value_counts().sort_values()
    bars = ax.barh(readiness.index, readiness.values, color="#70AD47")
    ax.bar_label(bars, labels=[f"{v:,}" for v in readiness.values], padding=2, fontsize=6.1)
    ax.set_xlabel("Binding-site instances")
    clean(ax, "x")

    ax = axes[4]
    panel(ax, "E", "Binding-residue composition")
    residue_series = pd.Series(residues).sort_values()
    residue_colors = {
        "Hydrophobic": "#4C78A8",
        "Aromatic": "#6A3D9A",
        "Polar": "#72B7B2",
        "Positive": "#E45756",
        "Negative": "#F2CF5B",
        "Special": "#B7B7B7",
    }
    bars = ax.barh(
        residue_series.index,
        residue_series.values,
        color=[residue_colors.get(x, "#B7B7B7") for x in residue_series.index],
    )
    ax.bar_label(bars, labels=[f"{v:,}" for v in residue_series.values], padding=2, fontsize=6.2)
    ax.set_xlabel("Parsed residue mentions")
    clean(ax, "x")

    ax = axes[5]
    panel(ax, "F", "Protein class × chemical class enrichment")
    top_protein = pair["functional_class"].value_counts().head(8).index
    top_chemical = pair["chemical_class"].value_counts().head(9).index
    subset = pair[
        pair["functional_class"].isin(top_protein)
        & pair["chemical_class"].isin(top_chemical)
    ]
    observed = pd.crosstab(subset["functional_class"], subset["chemical_class"]).reindex(
        index=top_protein, columns=top_chemical, fill_value=0
    )
    expected = np.outer(observed.sum(axis=1), observed.sum(axis=0)) / observed.to_numpy().sum()
    residual = (observed.to_numpy() - expected) / np.sqrt(expected + 1e-9)
    residual = pd.DataFrame(
        np.clip(residual, -5, 5),
        index=observed.index,
        columns=observed.columns,
    )
    sns.heatmap(
        residual,
        cmap="vlag",
        center=0,
        vmin=-5,
        vmax=5,
        cbar_kws={"label": "Pearson residual (clipped)"},
        linewidths=0.35,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticklabels([x.replace("_", " ") for x in residual.columns], rotation=45, ha="right")
    ax.set_yticklabels([x.replace("_", " ") for x in residual.index], rotation=0)
    paths = save(fig, "M6_binding_evidence_sites")
    write_table(source_tier_plot.reset_index(names="source"), "M6_source_evidence_matrix.tsv")
    write_table(observed.reset_index(), "M6_protein_chemical_observed.tsv")
    write_table(residual.reset_index(), "M6_protein_chemical_residual.tsv")
    return paths


def m7_negative_quality() -> dict:
    fig, axes = make_figure("M7  Negative evidence, conflicts and data reliability")
    top_classes = pair["functional_class"].value_counts().head(9).index

    ax = axes[0]
    panel(ax, "A", "Positive and negative pair coverage")
    positive_counts = pair["functional_class"].value_counts().reindex(top_classes).fillna(0)
    negative_counts = pd.Series(negative_class_counts).reindex(top_classes).fillna(0)
    y = np.arange(len(top_classes))
    height = 0.36
    ax.barh(y - height / 2, positive_counts.values, height=height, color=COLORS["positive"], label="Positive")
    ax.barh(y + height / 2, negative_counts.values, height=height, color=COLORS["negative"], label="Negative")
    ax.set_yticks(y, [x.replace("_", " ") for x in top_classes])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Unique pairs (log scale)")
    ax.legend(frameon=False)
    clean(ax, "x")

    ax = axes[1]
    panel(ax, "B", "Positive–negative conflict rate")
    conflict = pd.Series(conflict_class_counts).reindex(top_classes).fillna(0)
    rate = conflict / positive_counts.replace(0, np.nan) * 100
    bars = ax.barh(
        [x.replace("_", " ") for x in top_classes[::-1]],
        rate[::-1],
        color=COLORS["conflict"],
    )
    ax.bar_label(bars, labels=[f"{v:.1f}%" for v in rate[::-1]], padding=2, fontsize=6.2)
    ax.set_xlabel("Conflict pairs / positive pairs (%)")
    clean(ax, "x")

    ax = axes[2]
    panel(ax, "C", "Conflict context by evidence tier")
    conflict_matrix = pd.crosstab(
        conflict_context["functional_class"],
        conflict_context["best_binding_evidence_level"],
    ).reindex(index=top_classes, columns=["BE1", "BE2", "BE3"], fill_value=0)
    sns.heatmap(
        np.log10(conflict_matrix + 1),
        cmap=sns.light_palette(COLORS["conflict"], as_cmap=True),
        annot=conflict_matrix.map(lambda value: f"{value:,.0f}" if value else ""),
        fmt="",
        cbar_kws={"label": "log10(count + 1)"},
        linewidths=0.35,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels([x.replace("_", " ") for x in conflict_matrix.index], rotation=0)

    ax = axes[3]
    panel(ax, "D", "Negative identity resolution")
    draw_flow_boxes(
        ax,
        ["Source\nrecords", "Mapped\nrelease", "Review\nqueue", "Duplicates\nremoved"],
        counts=[10_246_948, 7_630_659, 2_615_046, 1_243],
        colors=["#F2F2F2", "#D9EAD3", "#F4CCCC", "#D9D9D9"],
    )

    ax = axes[4]
    panel(ax, "E", "Core-field missingness")
    compound_missing = {
        field: compound_scan["missing"][field] / compound_scan["core_ok_total"] * 100
        for field in compound_scan["missing"]
    }
    protein_missing = {
        "Protein: class": protein["functional_class"].isna().mean() * 100,
        "Protein: PDB IDs": protein["pdb_ids"].isna().mean() * 100,
        "Protein: AlphaFold": protein["alphafolddb_ids"].isna().mean() * 100,
    }
    missing_values = {
        **protein_missing,
        **{f"Compound: {field}": value for field, value in compound_missing.items()},
        "Pair: quantitative value": pair["best_standard_value_nM"].isna().mean() * 100,
        "Site: PDB ID": site["pdb_ids"].isna().mean() * 100,
        "Site: residue description": site["residue_or_site_description"].isna().mean() * 100,
    }
    missing_series = pd.Series(missing_values).sort_values().tail(12)
    bars = ax.barh(missing_series.index, missing_series.values, color="#A5A5A5")
    ax.bar_label(bars, labels=[f"{v:.1f}%" for v in missing_series.values], padding=2, fontsize=6.1)
    ax.set_xlabel("Missing / not available (%)")
    clean(ax, "x")

    ax = axes[5]
    panel(ax, "F", "Review queues and validation gates")
    queues = pd.Series(
        {
            "V6.1 positive identity review": 1_626_364,
            "V6.2 negative identity review": 2_615_046,
            "Positive–negative conflict pairs": 68_015,
            "Subunit review": 3_348,
            "Expression/location review": 496,
            "Blocking QA errors": 0,
        }
    ).sort_values()
    bars = ax.barh(
        queues.index,
        queues.values + 1,
        color=[COLORS["review"] if value else COLORS["positive"] for value in queues.values],
    )
    ax.set_xscale("log")
    ax.bar_label(bars, labels=[f"{v:,.0f}" for v in queues.values], padding=2, fontsize=6.1)
    ax.set_xlabel("Records / pairs + 1 (log scale)")
    clean(ax, "x")
    paths = save(fig, "M7_negative_conflicts_quality")
    write_table(conflict_matrix.reset_index(), "M7_conflict_evidence_matrix.tsv")
    write_table(pd.DataFrame({"field": missing_series.index, "missing_percent": missing_series.values}), "M7_missingness.tsv")
    return paths


def m8_docking_prioritization() -> dict:
    fig, axes = make_figure("M8  Structure-aware docking prioritization")

    ax = axes[0]
    panel(ax, "A", "Docking selection funnel")
    funnel = pd.Series(
        {
            "All positive pairs": 1_502_456,
            "R0–D3 eligible": 192_592,
            "Non-conflict": 186_625,
            "Standard": 14_713,
            "Pilot": 6_019,
        }
    )
    y = np.arange(len(funnel))
    bars = ax.barh(y, funnel.values, color=["#A5A5A5", "#7F8C8D", "#5B9BD5", "#1B9E77", "#6A3D9A"])
    ax.set_yticks(y, funnel.index)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.bar_label(bars, labels=[f"{v:,.0f}" for v in funnel.values], padding=3, fontsize=6.4)
    ax.set_xlabel("Unique protein–compound pairs (log)")
    clean(ax, "x")

    ax = axes[1]
    panel(ax, "B", "Standard-set tier architecture")
    tier_order = ["R0", "D1", "D2", "D3"]
    tier_counts = docking_standard["docking_tier"].value_counts().reindex(tier_order)
    bars = ax.bar(
        tier_order,
        tier_counts.values,
        color=[COLORS[tier] for tier in tier_order],
    )
    ax.bar_label(bars, labels=[f"{v:,}" for v in tier_counts.values], fontsize=6.4)
    ax.set_ylabel("Standard pairs")
    ax.set_xlabel("R0 validation; D1–D3 prospective / mechanistic")
    clean(ax, "y")

    ax = axes[2]
    panel(ax, "C", "Target class × docking tier")
    top_classes = docking_standard["functional_class"].value_counts().head(9).index
    target_tier = pd.crosstab(
        docking_standard["functional_class"],
        docking_standard["docking_tier"],
    ).reindex(index=top_classes, columns=tier_order, fill_value=0)
    sns.heatmap(
        np.log10(target_tier + 1),
        cmap=sns.light_palette(COLORS["blue"], as_cmap=True),
        annot=target_tier.map(lambda value: f"{value:,.0f}" if value else ""),
        fmt="",
        cbar_kws={"label": "log10(pair count + 1)"},
        linewidths=0.35,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels([x.replace("_", " ") for x in target_tier.index], rotation=0)

    ax = axes[3]
    panel(ax, "D", "Receptor structure readiness")
    structure_order = [
        "S1_pair_specific_experimental_site",
        "S2_membrane_experimental_OPM_PDBTM",
        "S3_other_experimental_PDB",
    ]
    readiness = pd.crosstab(
        docking_standard["docking_tier"],
        docking_standard["structure_grade"],
        normalize="index",
    ).reindex(index=tier_order, columns=structure_order, fill_value=0)
    bottom = np.zeros(len(readiness))
    structure_colors = ["#6A3D9A", "#1B9E77", "#5B9BD5"]
    for column, color in zip(readiness.columns, structure_colors):
        values = readiness[column].to_numpy() * 100
        ax.bar(readiness.index, values, bottom=bottom, label=column.split("_")[0], color=color)
        bottom += values
    ax.set_ylabel("Pairs within tier (%)")
    ax.legend(frameon=False, ncol=3)
    clean(ax, "y")

    ax = axes[4]
    panel(ax, "E", "Ligand tractability across tiers")
    ax.axis("off")
    fields = [
        ("molecular_weight", "MW (Da)", (80, 750)),
        ("xlogp", "XlogP", (-5, 10)),
        ("tpsa", "TPSA (Å²)", (0, 220)),
        ("rotatable_bond_count", "Rotatable bonds", (0, 25)),
    ]
    for index, (field, label, limits) in enumerate(fields):
        inset = ax.inset_axes([0.04 + (index % 2) * 0.49, 0.55 - (index // 2) * 0.5, 0.44, 0.38])
        for tier in tier_order:
            values = pd.to_numeric(
                docking_standard.loc[docking_standard["docking_tier"].eq(tier), field],
                errors="coerce",
            )
            values = values[values.between(*limits)]
            x, y = ecdf(values.to_numpy())
            inset.plot(x, y, color=COLORS[tier], label=tier, lw=1.1)
        inset.set_xlabel(label, fontsize=6.2)
        inset.set_ylabel("ECDF", fontsize=6)
        clean(inset, "both")
    axes[4].legend(
        [Line2D([0], [0], color=COLORS[tier]) for tier in tier_order],
        tier_order,
        frameon=False,
        ncol=4,
        loc="lower center",
    )

    ax = axes[5]
    panel(ax, "F", "HPC workload choices")
    sets = {
        "Pilot": docking_pilot,
        "Standard": docking_standard,
        "Full": docking_full,
    }
    workload = pd.DataFrame(
        {
            name: {
                "Pairs": len(frame),
                "Targets": frame["target_uniprot_id"].nunique(),
                "Compounds": frame["compound_internal_id"].nunique(),
            }
            for name, frame in sets.items()
        }
    ).T
    x = np.arange(len(workload))
    width = 0.24
    for index, metric in enumerate(["Pairs", "Targets", "Compounds"]):
        ax.bar(
            x + (index - 1) * width,
            workload[metric],
            width,
            label=metric,
            color=[COLORS["A"], COLORS["B"], COLORS["C"]][index],
        )
    ax.set_xticks(x, workload.index)
    ax.set_yscale("log")
    ax.set_ylabel("Count (log scale)")
    ax.legend(frameon=False, ncol=3)
    clean(ax, "y")
    paths = save(fig, "M8_docking_prioritization")
    write_table(target_tier.reset_index(), "M8_target_tier_matrix.tsv")
    write_table(workload.reset_index(names="set"), "M8_workload.tsv")
    return paths


outputs = {
    "M5": m5_chemical_space(),
    "M6": m6_binding_evidence(),
    "M7": m7_negative_quality(),
    "M8": m8_docking_prioritization(),
}

caption = f"""# Main figure captions — M5–M8

## Figure M5. Chemical identity and small-molecule space

Canonical parents and exact forms are separate entities. The chemical-space
embedding uses Morgan fingerprints (radius 2, 1024 bits) and a fixed random seed
on a deterministic sample of QC-passed core compounds. UMAP is descriptive:
distances should not be interpreted as an absolute chemical-similarity scale.
Physicochemical ECDFs retain outliers within the displayed ranges.

## Figure M6. Binding evidence and binding-site landscape

Evidence records are grouped by source and mutually exclusive highest evidence
tier. Independent-source counts and quantitative potency use unique
protein–canonical compound pairs. Standardized nM values are displayed on a log
scale but activity types are not assumed to be thermodynamically equivalent.
Residue composition is calculated from parsable residue mentions in
coordinate-supported site descriptions.

## Figure M7. Negative evidence, conflicts and data reliability

Positive and negative coverage is reported at the unique-pair level. Conflict
rate is the number of positive pairs with at least one negative record divided
by positive pairs within each functional class. Conflicts may reflect assay or
context differences and are flagged rather than automatically adjudicated.
Missing, not-applicable and review states remain explicit.

## Figure M8. Structure-aware docking prioritization

The docking funnel begins with {len(pair):,} positive pairs. R0 is a protocol
validation/redocking tier and is not a prospective discovery tier. D1–D3 are
ranked prospective or mechanistic sets. Candidate PDB identifiers indicate
triage readiness only; chain, state, biological assembly, binding box and ligand
microstates still require preparation before HPC submission. Docking scores
must not be interpreted as binding affinities without protocol validation.
"""
(CAPTION_DIR / "CAPTIONS_M5_M8.md").write_text(caption, encoding="utf-8")
write_json(
    {
        "status": "PASS",
        "figures": outputs,
        "compound_rows_scanned": int(compound_scan["total"]),
        "compound_core_qc_rows": int(compound_scan["core_ok_total"]),
        "compound_embedding_rows": int(len(embedding)),
        "positive_pair_rows": int(len(pair)),
        "negative_pair_rows": int(sum(negative_class_counts.values())),
        "conflict_pair_rows": int(len(conflict_pairs)),
        "site_rows": int(len(site)),
        "docking_standard_rows": int(len(docking_standard)),
        "docking_pilot_rows": int(len(docking_pilot)),
        "docking_full_rows": int(len(docking_full)),
    },
    "FIGURES_M5_M8_VALIDATION.json",
)
print(json.dumps(outputs, ensure_ascii=False, indent=2))
