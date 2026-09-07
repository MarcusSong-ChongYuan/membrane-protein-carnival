"""P4 only: membrane-protein role × normal-tissue RNA coverage bubble matrix.

No biological threshold is created here.  The frozen HPA `detection_status`
field is used exactly as supplied.  nTPM is used only for the supplementary
median-expression column because it is the explicitly normalised RNA unit.
"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize


ROOT = Path(r"D:\finale\42_MemPro_protein_module_core_figures_20260825")
FORMAL = Path(r"D:\finale\FORMAL")
OUT = ROOT / "results" / "protein_figure_pool"
OUT.mkdir(parents=True, exist_ok=True)

MASTER = FORMAL / "01_core_v72" / "01_release_tables" / "protein_master_v72.tsv.gz"
ROLE = FORMAL / "03_publication_repairs" / "protein_membrane_role_FORMAL.tsv"
EXPRESSION = FORMAL / "01_core_v72" / "01_release_tables" / "expression_measurement_v72.tsv.gz"

PALETTE = ["#7b95c6", "#49c2d9", "#a1d8e8", "#67a583", "#a2c986", "#d0e2c0", "#fded95", "#ffc1a6", "#f59c7c", "#f47254", "#c85e62"]
TEXT = "#24364B"
GREY = "#E8EDF0"
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7.2, "axes.linewidth": 0.8, "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.facecolor": "white", "axes.facecolor": "white", "text.color": TEXT,
})


def clean(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def role_table() -> pd.DataFrame:
    p = pd.read_csv(MASTER, sep="\t", compression="gzip", usecols=["canonical_uniprot_accession", "public_inclusion_v7"])
    p = p.loc[p.public_inclusion_v7.eq(1), ["canonical_uniprot_accession"]].rename(columns={"canonical_uniprot_accession": "target_uniprot_id"})
    if p.target_uniprot_id.nunique() != 7800:
        raise ValueError("P4 expected 7,800 formal canonical proteins.")
    role = pd.read_csv(ROLE, sep="\t")
    idcol = next(c for c in role.columns if c in {"target_uniprot_id", "canonical_uniprot_accession"})
    rolecol = next(c for c in role.columns if "membrane_role" in c and "primary" in c)
    role = role[[idcol, rolecol]].drop_duplicates(idcol).rename(columns={idcol: "target_uniprot_id", rolecol: "primary_role"})
    out = p.merge(role, on="target_uniprot_id", how="left")
    out["primary_role"] = clean(out.primary_role).replace("", "family_defined_membrane_role_unresolved")
    return out


def main() -> None:
    roles = role_table()
    role_count = roles.primary_role.value_counts()
    role_order = role_count.index.tolist()

    usecols = ["target_uniprot_id", "hpa_ensembl_gene_id", "source_dataset", "tissue", "value", "unit", "detection_status", "measurement_status_v71", "mapped_measured_denominator_eligible_v71"]
    # The frozen expression table has >3 million rows; stream it to avoid
    # materialising non-RNA and non-formal-protein records in memory.
    formal_ids = set(roles.target_uniprot_id)
    kept = []
    before_rows = 0
    for chunk in pd.read_csv(EXPRESSION, sep="\t", compression="gzip", usecols=usecols, low_memory=False, chunksize=200_000):
        before_rows += len(chunk)
        eligible = chunk.mapped_measured_denominator_eligible_v71.astype(str).eq("1")
        keep = (
            chunk.source_dataset.eq("HPA_tissue_RNA")
            & chunk.unit.eq("nTPM")
            & chunk.measurement_status_v71.eq("MEASURED")
            & eligible
            & chunk.target_uniprot_id.isin(formal_ids)
            & chunk.tissue.notna()
        )
        if keep.any():
            kept.append(chunk.loc[keep].copy())
    e = pd.concat(kept, ignore_index=True)
    # This is the non-consensus HPA normal-tissue RNA layer.  The frozen table
    # defines both MAPPED_AND_MEASURED_ONLY and DETECTED/NOT_DETECTED states.
    e["tissue"] = clean(e.tissue)
    e["detection_status"] = clean(e.detection_status).str.lower()
    if not set(e.detection_status).issubset({"detected", "not_detected"}):
        raise ValueError(f"Unexpected frozen HPA detection statuses: {sorted(e.detection_status.unique())}")
    # One canonical protein can have several HPA Ensembl mappings.  Collapse
    # only when their frozen detection calls agree; disagreement is retained as
    # an explicit analysis-specific unresolved state, never settled by max/min.
    e["value_numeric"] = pd.to_numeric(e.value, errors="coerce")
    e = e.groupby(["target_uniprot_id", "tissue"], as_index=False, sort=False).agg(
        n_hpa_gene_mappings=("value_numeric", "size"),
        n_distinct_hpa_gene_mappings=("hpa_ensembl_gene_id", "nunique"),
        detection_nunique=("detection_status", "nunique"),
        frozen_detection_status=("detection_status", "first"),
        median_expression_nTPM_protein_tissue=("value_numeric", "median"),
    )
    e["detection_status_resolved"] = np.where(
        e.detection_nunique.eq(1), e.frozen_detection_status, "multi_mapping_detection_unresolved"
    )
    e = e.merge(roles, on="target_uniprot_id", how="inner", validate="many_to_one")

    rows = []
    tissue_order_all = e.groupby("tissue").target_uniprot_id.nunique().sort_values(ascending=False)
    # Stable alphabetical tie-break is automatic: all retained tissues happen
    # to have identical coverage in the HPA normal-tissue RNA layer.
    tissue_order_all = tissue_order_all.sort_index().sort_values(ascending=False, kind="stable")
    main_tissues = tissue_order_all.head(20).index.tolist()
    for role in role_order:
        rtotal = int(role_count[role])
        for tissue in tissue_order_all.index:
            d = e.loc[(e.primary_role == role) & (e.tissue == tissue)]
            usable_d = d.loc[~d.detection_status_resolved.eq("multi_mapping_detection_unresolved")]
            usable = int(usable_d.target_uniprot_id.nunique())
            detected = int(usable_d.loc[usable_d.detection_status_resolved.eq("detected"), "target_uniprot_id"].nunique())
            mapping_conflict = int(d.loc[d.detection_status_resolved.eq("multi_mapping_detection_unresolved"), "target_uniprot_id"].nunique())
            rows.append({
                "primary_role": role,
                "tissue": tissue,
                "n_role_total": rtotal,
                "n_role_with_expression_data": usable,
                "n_detected": detected,
                "detected_fraction": detected / usable if usable else np.nan,
                "median_expression_nTPM": float(usable_d.median_expression_nTPM_protein_tissue.median()) if usable else np.nan,
                "n_multi_mapping_detection_unresolved": mapping_conflict,
                "hpa_source_dataset": "HPA_tissue_RNA",
                "expression_unit": "nTPM",
                "detection_status_field": "detection_status",
                "measurement_policy": "MAPPED_AND_MEASURED_ONLY",
                "included_in_main_bubble_matrix": tissue in main_tissues,
            })
    source = pd.DataFrame(rows)
    source.to_csv(OUT / "P4_role_tissue_bubble.tsv", sep="\t", index=False)

    stats = pd.DataFrame([{
        "protein_id_field": "target_uniprot_id",
        "role_field": "primary_role (from protein_membrane_role_FORMAL)",
        "tissue_field": "tissue",
        "expression_field": "value (nTPM)",
        "detected_field": "detection_status",
        "input_expression_rows_all_layers": before_rows,
        "analysis_protein_tissue_pairs_after_canonical_collapse": len(e),
        "formal_proteins_total": 7800,
        "formal_proteins_with_usable_HPA_tissue_RNA": int(e.target_uniprot_id.nunique()),
        "analysis_specific_excluded_formal_proteins_no_usable_HPA_tissue_RNA": int(7800 - e.target_uniprot_id.nunique()),
        "all_tissues_in_source_data": int(tissue_order_all.size),
        "main_figure_tissues": int(len(main_tissues)),
        "tissue_selection_rule": "top 20 by unique protein coverage; stable alphabetical tie-break",
        "detection_rule": "frozen HPA detection_status == detected; no new threshold created",
        "multi_mapping_rule": "A protein-tissue call is usable only when all retained HPA mappings agree on detected/not_detected; conflicts remain explicit and are excluded from that tissue's detected-fraction denominator.",
    }])
    stats.to_csv(OUT / "P4_role_tissue_statistics.tsv", sep="\t", index=False)

    # Persist role colours now; later tasks read exactly the same mapping.
    colour_map = {role: PALETTE[i % len(PALETTE)] for i, role in enumerate(role_order)}
    (OUT / "category_colors.json").write_text(json.dumps({"primary_role": colour_map}, indent=2), encoding="utf-8")

    plot = source.loc[source.included_in_main_bubble_matrix].copy()
    x_order = main_tissues
    y_order = role_order
    xpos = {t: i for i, t in enumerate(x_order)}
    ypos = {r: i for i, r in enumerate(y_order)}
    fig = plt.figure(figsize=(7.01, 6.55), constrained_layout=False)
    ax = fig.add_axes([0.24, 0.20, 0.62, 0.64])
    cmap = LinearSegmentedColormap.from_list("hpa_detected_fraction", ["#a1d8e8", "#49c2d9", "#67a583"])
    norm = Normalize(vmin=0, vmax=1)
    # Size maps usable protein count (not detection rate), maintaining separate encodings.
    sizes = 18 + 170 * (plot.n_role_with_expression_data / plot.n_role_with_expression_data.max())
    ax.scatter(plot.tissue.map(xpos), plot.primary_role.map(ypos), s=sizes, c=plot.detected_fraction, cmap=cmap, norm=norm,
               edgecolor="white", linewidth=0.45, alpha=0.94)
    ax.set_xlim(-.65, len(x_order) - .35)
    ax.set_ylim(len(y_order) - .45, -.55)
    ax.set_xticks(range(len(x_order)))
    ax.set_xticklabels([textwrap.fill(t.title(), 13) for t in x_order], rotation=60, ha="right", rotation_mode="anchor", fontsize=6.4)
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels([textwrap.fill(r.replace("_", " ").title(), 27) for r in y_order], fontsize=7.0)
    ax.tick_params(axis="x", length=0, pad=4)
    ax.tick_params(axis="y", length=0, pad=3)
    for x in range(len(x_order)):
        ax.axvline(x, color="#EDF1F3", lw=.55, zorder=0)
    for y in range(len(y_order)):
        ax.axhline(y, color="#EDF1F3", lw=.55, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cax = fig.add_axes([0.89, .32, .018, .42])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[0, .25, .5, .75, 1])
    cb.ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
    cb.ax.tick_params(labelsize=6.5, length=2)
    cb.set_label("Detected fraction", fontsize=7, labelpad=5)
    # Bubble-size key is deliberately a separate variable from bubble colour.
    sample_sizes = [int(plot.n_role_with_expression_data.quantile(q)) for q in [.25, .75]]
    handles = [plt.scatter([], [], s=18 + 170*(s / plot.n_role_with_expression_data.max()), facecolor="#D5DDE1", edgecolor="white", label=f"n = {s:,}") for s in sample_sizes]
    ax.legend(handles=handles, title="Usable proteins", title_fontsize=7, fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(1.01, .10), labelspacing=.8)
    fig.text(.05, .955, "Membrane-protein role × tissue RNA detection", fontsize=11, weight="bold", color=TEXT)
    fig.text(.05, .922, "Bubble colour: frozen HPA detected fraction; bubble area: proteins with usable nTPM information. Top 20 tissues by coverage are shown.", fontsize=6.8, color="#5E7182")
    fig.text(.24, .045, "Analysis-specific denominator: 7,531 / 7,800 formal proteins with mapped and measured HPA_tissue_RNA. Missing is not interpreted as zero.", fontsize=6.7, color="#5E7182")
    for ext, kw in [("svg", {}), ("pdf", {}), ("png", {"dpi": 600})]:
        fig.savefig(OUT / f"P4_role_tissue_bubble.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print("P4 complete: only role × tissue bubble matrix was created.")


if __name__ == "__main__":
    main()
