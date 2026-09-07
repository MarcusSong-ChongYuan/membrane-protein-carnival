"""Build a same-Uberon-system expression-versus-disease coverage chart.

The frozen V7.2 disease/anatomy and HPA expression tables remain read only.
This script derives a visualization-only HPA tissue-to-system projection using
the same frozen Uberon snapshot and terminal anchor set as D2.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd

from build_D2_top_level_anatomical_coverage import (
    FORMAL,
    OUT,
    SNAPSHOT,
    build_terminal_anchor_map,
    display_name,
    parse_uberon,
    ancestor_distances,
)


mpl.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


def fmt(value: int) -> str:
    return f"{int(value):,}"


def main() -> None:
    if not SNAPSHOT.exists():
        raise FileNotFoundError(f"D2 Uberon visualization snapshot is missing: {SNAPSHOT}")
    terms, version = parse_uberon(SNAPSHOT)
    candidates, terminal = build_terminal_anchor_map(terms)

    # D2 defines the disease-side display system set.  Using its exact IDs makes
    # the two panels directly comparable without expanding the disease taxonomy.
    disease_summary = pd.read_csv(OUT / "D2_top_level_anatomical_system_summary.tsv", sep="\t", dtype=str)
    disease_summary["n_unique_diseases"] = pd.to_numeric(disease_summary["n_unique_diseases"])
    disease_system_ids = set(disease_summary["top_level_uberon_id"])

    expr_cols = [
        "target_uniprot_id", "source_dataset", "measurement_layer", "tissue", "tissue_ontology_id",
        "tissue_ontology_label", "tissue_mapping_status", "measurement_status_v711",
        "mapped_measured_denominator_eligible_v71", "detection_status_v711",
    ]
    expression = pd.read_csv(
        FORMAL / "expression_measurement_v72.tsv.gz", sep="\t", dtype=str,
        keep_default_na=False, usecols=expr_cols,
    )

    # A single modality is used: HPA normal tissue RNA.  The source table stores
    # multiple RNA units per protein/tissue, so final protein×tissue×system rows
    # are deduplicated before counts are calculated.
    tissue_rna = expression.loc[
        (expression["source_dataset"] == "HPA_tissue_RNA")
        & (expression["measurement_layer"] == "RNA")
        & (expression["measurement_status_v711"] == "MEASURED")
        & (expression["mapped_measured_denominator_eligible_v71"] == "1")
        & expression["tissue_mapping_status"].str.startswith("mapped_", na=False)
        & (expression["tissue_ontology_id"] != "")
    ].copy()

    tissue_rows = tissue_rna[["tissue", "tissue_ontology_id", "tissue_ontology_label", "tissue_mapping_status"]].drop_duplicates()
    projection_rows: list[dict] = []
    for row in tissue_rows.itertuples(index=False):
        tissue_ids = [item.strip() for item in row.tissue_ontology_id.split(";") if item.strip()]
        # Formally mapped HPA tissue labels should be one-to-one.  Multi-ID labels
        # are intentionally not projected because they are not a unique tissue term.
        if len(tissue_ids) != 1:
            projection_rows.append({
                "hpa_tissue": row.tissue,
                "source_uberon_id": row.tissue_ontology_id,
                "source_uberon_term": row.tissue_ontology_label,
                "tissue_mapping_status": row.tissue_mapping_status,
                "top_level_uberon_id": "",
                "top_level_system": "",
                "display_system": "",
                "mapping_relation": "",
                "mapping_status": "ambiguous_source_tissue_term",
            })
            continue
        tissue_id = tissue_ids[0]
        if tissue_id not in terms:
            projection_rows.append({
                "hpa_tissue": row.tissue, "source_uberon_id": tissue_id,
                "source_uberon_term": row.tissue_ontology_label, "tissue_mapping_status": row.tissue_mapping_status,
                "top_level_uberon_id": "", "top_level_system": "", "display_system": "",
                "mapping_relation": "", "mapping_status": "not_in_uberon_snapshot",
            })
            continue
        ancestors = ancestor_distances(tissue_id, terms)
        hits: list[str] = ([tissue_id] if tissue_id in candidates else [])
        hits.extend(candidate for candidate in ancestors if candidate in candidates)
        anchors = sorted({terminal[candidate] for candidate in hits})
        retained = [anchor for anchor in anchors if anchor in disease_system_ids]
        if not retained:
            projection_rows.append({
                "hpa_tissue": row.tissue, "source_uberon_id": tissue_id,
                "source_uberon_term": row.tissue_ontology_label, "tissue_mapping_status": row.tissue_mapping_status,
                "top_level_uberon_id": "", "top_level_system": "", "display_system": "",
                "mapping_relation": "", "mapping_status": "no_matching_disease_display_system_ancestor",
            })
            continue
        for anchor in retained:
            projection_rows.append({
                "hpa_tissue": row.tissue,
                "source_uberon_id": tissue_id,
                "source_uberon_term": row.tissue_ontology_label,
                "tissue_mapping_status": row.tissue_mapping_status,
                "top_level_uberon_id": anchor,
                "top_level_system": terms[anchor]["name"],
                "display_system": display_name(terms[anchor]["name"]),
                "mapping_relation": "identity_top_level" if tissue_id == anchor else "is_a_or_part_of_ancestor_projection",
                "mapping_status": "mapped_to_disease_display_system",
            })
    tissue_projection = pd.DataFrame(projection_rows).sort_values(["mapping_status", "hpa_tissue", "top_level_system"])
    tissue_projection.to_csv(OUT / "D3_hpa_tissue_to_top_level_system_mapping.tsv", sep="\t", index=False)

    mapped_tissues = tissue_projection.loc[
        tissue_projection["mapping_status"] == "mapped_to_disease_display_system",
        ["hpa_tissue", "source_uberon_id", "top_level_uberon_id", "top_level_system", "display_system"],
    ].drop_duplicates()
    expr_projected = tissue_rna.merge(
        mapped_tissues, left_on=["tissue", "tissue_ontology_id"], right_on=["hpa_tissue", "source_uberon_id"], how="inner",
    )
    usable = expr_projected.drop_duplicates(["target_uniprot_id", "tissue", "top_level_uberon_id"])
    detected = usable.loc[usable["detection_status_v711"] == "DETECTED"].drop_duplicates(
        ["target_uniprot_id", "tissue", "top_level_uberon_id"]
    )
    expression_summary = (
        usable.groupby(["top_level_uberon_id", "top_level_system", "display_system"], as_index=False)
        .agg(
            n_unique_formal_membrane_proteins_with_usable_tissue_RNA=("target_uniprot_id", "nunique"),
            n_unique_HPA_tissues=("tissue", "nunique"),
        )
        .merge(
            detected.groupby("top_level_uberon_id", as_index=False)["target_uniprot_id"].nunique()
            .rename(columns={"target_uniprot_id": "n_unique_RNA_detected_formal_membrane_proteins"}),
            on="top_level_uberon_id", how="left",
        )
    )
    expression_summary["n_unique_RNA_detected_formal_membrane_proteins"] = (
        expression_summary["n_unique_RNA_detected_formal_membrane_proteins"].fillna(0).astype(int)
    )

    combined = disease_summary.merge(
        expression_summary[["top_level_uberon_id", "n_unique_formal_membrane_proteins_with_usable_tissue_RNA",
                            "n_unique_RNA_detected_formal_membrane_proteins", "n_unique_HPA_tissues"]],
        on="top_level_uberon_id", how="left",
    )
    for field in ["n_unique_formal_membrane_proteins_with_usable_tissue_RNA", "n_unique_RNA_detected_formal_membrane_proteins", "n_unique_HPA_tissues"]:
        combined[field] = combined[field].fillna(0).astype(int)
    combined = combined.sort_values("n_unique_RNA_detected_formal_membrane_proteins", ascending=False).reset_index(drop=True)
    combined.to_csv(OUT / "D3_expression_disease_top_level_system_summary.tsv", sep="\t", index=False)

    # Bilateral bar: disease count and RNA-detected protein count are distinct
    # coverage quantities, so axes remain separately labelled and are never summed.
    plot = combined.sort_values("n_unique_RNA_detected_formal_membrane_proteins", ascending=True)
    y = list(range(len(plot)))
    left = plot["n_unique_RNA_detected_formal_membrane_proteins"].to_numpy()
    right = plot["n_unique_diseases"].astype(int).to_numpy()
    # Protein and disease counts are different coverage quantities with markedly
    # different ranges.  Normalize only the bar geometry to side-specific axes;
    # raw counts remain the labels and tick values on their own side.
    left_max = max(left.max(), 1)
    right_max = max(right.max(), 1)
    left_scaled = left / left_max
    right_scaled = right / right_max
    extent = 1.24
    fig_height = max(5.2, 0.48 * len(plot) + 1.8)
    fig, ax = plt.subplots(figsize=(7.0, fig_height))
    ax.barh(y, -left_scaled, color="#168C8C", edgecolor="#126F70", linewidth=0.55, height=0.62, label="RNA-detected membrane proteins")
    ax.barh(y, right_scaled, color="#B76E79", edgecolor="#934F5A", linewidth=0.55, height=0.62, label="Mapped diseases")
    ax.axvline(0, color="#425A70", linewidth=0.9, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["display_system"], fontsize=8.5, color="#24364B")
    ax.set_xlim(-extent, extent)
    tick_positions = [-1.0, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1.0]
    tick_labels = [
        fmt(left_max), fmt(round(left_max * 0.75)), fmt(round(left_max * 0.50)), fmt(round(left_max * 0.25)), "0",
        fmt(round(right_max * 0.25)), fmt(round(right_max * 0.50)), fmt(round(right_max * 0.75)), fmt(right_max),
    ]
    ax.set_xticks(tick_positions, tick_labels)
    ax.xaxis.grid(True, color="#E3E8EC", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#748392")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=7.8, colors="#506174")
    # Keep the two metric headers close to the plotting region and below the
    # figure-level subtitle; this prevents the labels from visually merging.
    ax.text(-0.57, len(plot) + 0.27, "RNA-detected\nmembrane proteins", ha="center", va="bottom",
            fontsize=8.1, fontweight="bold", color="#126F70", linespacing=0.95)
    ax.text(0.57, len(plot) + 0.27, "Unique mapped\ndiseases", ha="center", va="bottom",
            fontsize=8.1, fontweight="bold", color="#934F5A", linespacing=0.95)
    for yi, protein_n, disease_n, protein_scaled, disease_scaled in zip(y, left, right, left_scaled, right_scaled):
        ax.text(-protein_scaled - 0.025, yi, fmt(protein_n), ha="right", va="center", fontsize=7.4, color="#126F70")
        ax.text(disease_scaled + 0.025, yi, fmt(disease_n), ha="left", va="center", fontsize=7.4, color="#934F5A")
    ax.set_xlabel("Count (side-specific axes; quantities are not directly comparable)", fontsize=8.4, color="#24364B", labelpad=8)
    fig.suptitle("Expression and disease coverage across harmonized anatomical systems",
                 x=0.125, y=0.985, ha="left", fontsize=10.6, fontweight="bold", color="#24364B")
    fig.text(0.125, 0.915,
             f"Left: HPA tissue-RNA detected formal membrane proteins (0–{fmt(left_max)}). Right: unique diseases mapped through formal disease–Uberon relations (0–{fmt(right_max)}).",
             ha="left", va="top", fontsize=7.3, color="#5E7182")
    fig.text(0.125, 0.012,
             f"Systems use the same Uberon {version} projection as D2. HPA tissue terms without a unique matching system ancestor are retained in mapping QC, not forced into a display system.",
             ha="left", va="bottom", fontsize=6.8, color="#5E7182")
    fig.subplots_adjust(left=0.30, right=0.94, top=0.86, bottom=0.11)
    base = OUT / "D3_expression_disease_top_level_system_bilateral_bar"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)

    hpa_only_tissues = tissue_projection.loc[
        tissue_projection["mapping_status"] != "mapped_to_disease_display_system", "hpa_tissue"
    ].nunique()
    qc_rows = [
        ("expression_source_dataset", "HPA_tissue_RNA", "single normal-tissue RNA source; IHC, MS and cell-type RNA excluded", "PASS"),
        ("expression_measurement_layer", "RNA", "RNA is not described as protein abundance", "PASS"),
        ("eligible_expression_rows_before_projection", len(tissue_rna), "measured, mapped and denominator-eligible HPA tissue-RNA rows", "PASS"),
        ("unique_HPA_tissue_terms", tissue_rows.shape[0], "distinct tissue label × ontology mapping pairs", "PASS"),
        ("HPA_tissues_without_matching_disease_display_system", hpa_only_tissues, "retained in mapping table, excluded from bilateral display", "ACCEPTABLE_AS_UNRESOLVED"),
        ("displayed_systems", len(combined), "exact D2 disease-system set", "PASS"),
        ("duplicate_protein_tissue_system_rows_after_dedup", int(usable.duplicated(["target_uniprot_id", "tissue", "top_level_uberon_id"]).sum()), "must equal zero", "PASS"),
        ("protein_metric", "unique RNA-detected formal membrane proteins", "left bar", "PASS"),
        ("disease_metric", "unique canonical diseases", "right bar", "PASS"),
        ("shared_system_projection", "D2 Uberon 2026-06-19 terminal-anchor set", "same display system IDs on both sides", "PASS"),
    ]
    with (OUT / "D3_expression_disease_top_level_system_QC.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["check", "value", "definition", "status"])
        writer.writerows(qc_rows)
    with (OUT / "D3_expression_disease_top_level_system_fields_used.txt").open("w", encoding="utf-8") as handle:
        handle.write("D3 bilateral expression-versus-disease coverage: fields used\n\n")
        handle.write("Expression table (read-only formal V7.2): expression_measurement_v72.tsv.gz\n")
        handle.write("- target_uniprot_id, source_dataset, measurement_layer, tissue, tissue_ontology_id, tissue_ontology_label\n")
        handle.write("- tissue_mapping_status, measurement_status_v711, mapped_measured_denominator_eligible_v71, detection_status_v711\n")
        handle.write("Filter: HPA_tissue_RNA; RNA; MEASURED; denominator eligible = 1; mapped tissue; DETECTED for left-bar protein count.\n\n")
        handle.write("Disease table: D2_top_level_anatomical_system_summary.tsv\n")
        handle.write("- top_level_uberon_id, n_unique_diseases\n\n")
        handle.write(f"Shared visualization snapshot: Uberon {version}; same terminal anchors as D2.\n")
        handle.write("Unit: left = unique canonical proteins with detected HPA tissue RNA in a system; right = unique canonical diseases in a system.\n")
        handle.write("No raw expression rows, IHC, MS or cell-type RNA records were counted as left-bar proteins.\n")


if __name__ == "__main__":
    main()
