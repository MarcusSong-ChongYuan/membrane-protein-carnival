"""C3: Exact Bemis–Murcko scaffold rank–frequency for MemPro.

One canonical compound is counted once.  The script uses only the existing
formal parent-compound registry and its existing exact Murcko scaffold field;
it never regenerates scaffolds or changes compound identities.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compound_figure_pool"
FORMAL = Path(r"D:\finale\FORMAL")
MASTER = FORMAL / "01_core_v72" / "01_release_tables" / "compound_master_v72.tsv.gz"
PAIR = FORMAL / "01_core_v72" / "01_release_tables" / "protein_compound_pair_v72.tsv.gz"
CHEM = FORMAL / "02_companion_v721" / "02_companion_tables" / "compound_chemical_classification_v721.tsv.gz"
QUARANTINE = FORMAL / "03_publication_repairs" / "COMPOUND_STRUCTURE_REPAIR_AND_QUARANTINE.tsv"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.linewidth": 0.8,
})

INK = "#24364B"
SLATE = "#5E7182"
LINE = "#168C8C"
ACCENT = "#D6A23A"


def gini(values: np.ndarray) -> float:
    values = np.sort(np.asarray(values, dtype=float))
    if values.size == 0 or values.sum() == 0:
        return float("nan")
    index = np.arange(1, values.size + 1)
    return float((2 * np.sum(index * values) / (values.size * values.sum())) - (values.size + 1) / values.size)


def label_count(value: int) -> str:
    if value >= 1000:
        return f"{value / 1000:.1f}k".replace(".0k", "k")
    return str(value)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(
        MASTER,
        sep="\t",
        compression="gzip",
        usecols=["compound_internal_id", "compound_scope_class", "molecule_types", "record_qc_status"],
        low_memory=False,
    )
    chem = pd.read_csv(
        CHEM,
        sep="\t",
        compression="gzip",
        usecols=["compound_internal_id", "exact_murcko_scaffold_smiles", "scaffold_status"],
        low_memory=False,
    )
    pair = pd.read_csv(
        PAIR,
        sep="\t",
        compression="gzip",
        usecols=["pair_id", "target_uniprot_id", "compound_internal_id"],
        low_memory=False,
    ).drop_duplicates("pair_id")

    # Existing formal scope classes, not a guessed PubChem/CID proxy, define
    # the small-molecule population. Extended chemical entities are reported
    # but excluded from this chemical-structure analysis.
    small_scope = {"core_small_molecule", "structure_resolved_small_molecule"}
    selected = master.loc[master["compound_scope_class"].isin(small_scope)].copy()
    excluded_non_small_scope = master.loc[~master["compound_scope_class"].isin(small_scope)].copy()

    quarantine = pd.read_csv(QUARANTINE, sep="\t", low_memory=False)
    eligibility = quarantine[["compound_internal_id", "eligible_for_scaffold_fingerprint_pagtn_similarity"]].drop_duplicates("compound_internal_id")
    selected = selected.merge(eligibility, on="compound_internal_id", how="left")
    selected["eligible_for_scaffold_fingerprint_pagtn_similarity"] = selected[
        "eligible_for_scaffold_fingerprint_pagtn_similarity"
    ].fillna(True).astype(str).str.lower().isin({"true", "1", "yes"})
    selected_before_quarantine = len(selected)
    selected = selected.loc[selected["eligible_for_scaffold_fingerprint_pagtn_similarity"]].copy()

    data = selected.merge(chem, on="compound_internal_id", how="left", validate="one_to_one")
    valid = data.loc[
        data["scaffold_status"].eq("CYCLIC_SCAFFOLD") & data["exact_murcko_scaffold_smiles"].notna()
    ].copy()
    valid["exact_murcko_scaffold_smiles"] = valid["exact_murcko_scaffold_smiles"].astype(str)

    pair_to_scaffold = pair.merge(
        valid[["compound_internal_id", "exact_murcko_scaffold_smiles"]],
        on="compound_internal_id",
        how="inner",
        validate="many_to_one",
    )
    scaffold = (
        valid.groupby("exact_murcko_scaffold_smiles", as_index=False)["compound_internal_id"]
        .nunique()
        .rename(columns={"compound_internal_id": "unique_compound_count"})
    )
    protein_n = (
        pair_to_scaffold.groupby("exact_murcko_scaffold_smiles", as_index=False)["target_uniprot_id"]
        .nunique()
        .rename(columns={"target_uniprot_id": "unique_protein_count"})
    )
    pair_n = (
        pair_to_scaffold.groupby("exact_murcko_scaffold_smiles", as_index=False)["pair_id"]
        .nunique()
        .rename(columns={"pair_id": "unique_protein_compound_pair_count"})
    )
    scaffold = scaffold.merge(protein_n, on="exact_murcko_scaffold_smiles", how="left").merge(
        pair_n, on="exact_murcko_scaffold_smiles", how="left"
    )
    scaffold[["unique_protein_count", "unique_protein_compound_pair_count"]] = scaffold[
        ["unique_protein_count", "unique_protein_compound_pair_count"]
    ].fillna(0).astype(int)
    scaffold = scaffold.sort_values(
        ["unique_compound_count", "exact_murcko_scaffold_smiles"], ascending=[False, True], kind="stable"
    ).reset_index(drop=True)
    scaffold.insert(0, "scaffold_rank", np.arange(1, len(scaffold) + 1))
    scaffold["cumulative_unique_compounds"] = scaffold["unique_compound_count"].cumsum()
    scaffold["cumulative_compound_coverage_fraction"] = scaffold["cumulative_unique_compounds"] / len(valid)
    scaffold["included_in_exact_scaffold_analysis"] = True
    scaffold.to_csv(OUT / "C3_scaffold_rank_frequency.tsv", sep="\t", index=False)

    checkpoints = [1, 5, 10, 100]
    checkpoint_rows = []
    for rank in checkpoints:
        row = scaffold.iloc[min(rank, len(scaffold)) - 1]
        checkpoint_rows.append({
            "checkpoint": f"Top {rank}",
            "scaffold_rank": int(rank),
            "cumulative_unique_compounds": int(row["cumulative_unique_compounds"]),
            "cumulative_compound_coverage_fraction": float(row["cumulative_compound_coverage_fraction"]),
        })
    summary = pd.DataFrame([
        {"metric": "formal_registry_canonical_compounds", "value": int(master["compound_internal_id"].nunique()), "note": "All formal canonical compound identities."},
        {"metric": "included_formal_small_molecules_before_quarantine", "value": int(selected_before_quarantine), "note": "Existing formal scope classes core_small_molecule or structure_resolved_small_molecule."},
        {"metric": "excluded_extended_chemical_entities", "value": int(excluded_non_small_scope["compound_internal_id"].nunique()), "note": "Existing formal extended_chemical_entity scope; not treated as small molecules for structure statistics."},
        {"metric": "excluded_structure_quarantine", "value": int(selected_before_quarantine - len(selected)), "note": "Existing formal structure identity quarantine; no automatic repair applied."},
        {"metric": "acyclic_or_no_exact_murcko_scaffold", "value": int(len(selected) - len(valid)), "note": "Legitimate scaffold status, retained in registry but not rankable by exact Murcko scaffold."},
        {"metric": "canonical_compounds_with_exact_cyclic_scaffold", "value": int(len(valid)), "note": "Analysis N for C3."},
        {"metric": "unique_exact_scaffolds", "value": int(len(scaffold)), "note": "One unique exact Murcko scaffold string."},
        {"metric": "scaffold_frequency_gini", "value": gini(scaffold["unique_compound_count"].to_numpy()), "note": "Inequality of compound counts across exact scaffolds."},
        *[
            {"metric": f"{r['checkpoint']}_compound_coverage_fraction", "value": r["cumulative_compound_coverage_fraction"], "note": f"Cumulative coverage among {len(valid):,} scaffold-bearing canonical small molecules."}
            for r in checkpoint_rows
        ],
    ])
    summary.to_csv(OUT / "C3_scaffold_diversity_summary.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(178 / 25.4, 96 / 25.4), dpi=600)
    ax = fig.add_axes([0.13, 0.20, 0.76, 0.68])
    ax.plot(scaffold["scaffold_rank"], scaffold["unique_compound_count"], color=LINE, linewidth=1.45)
    ax.fill_between(scaffold["scaffold_rank"], scaffold["unique_compound_count"], 1, color=LINE, alpha=0.10, linewidth=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, len(scaffold))
    ax.set_ylim(0.8, max(2, scaffold["unique_compound_count"].max() * 1.32))
    ax.set_xlabel("Exact Bemis–Murcko scaffold rank", fontsize=8, color=INK)
    ax.set_ylabel("Unique canonical compounds per scaffold", fontsize=8, color=INK)
    # Log tick exponents are rendered at a reduced mathtext size; 8 pt keeps
    # their rendered glyphs above the 5 pt publication floor.
    ax.tick_params(axis="both", labelsize=8, colors=INK)
    ax.grid(which="major", color="#E4EAED", linewidth=0.65)
    ax.grid(which="minor", color="#F1F4F5", linewidth=0.35)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    # Fixed, collision-free offsets are used only for the four predeclared
    # checkpoints. Coordinates remain on the actual log-scale data axes.
    annotation_positions = {
        1: (2.4, 6000),
        5: (2.6, 1650),
        10: (13.0, 1250),
        100: (110.0, 360),
    }
    for row in checkpoint_rows:
        rank = row["scaffold_rank"]
        y = scaffold.iloc[rank - 1]["unique_compound_count"]
        ax.scatter(rank, y, s=25, color=ACCENT, edgecolor="white", linewidth=0.7, zorder=4)
        tx, ty = annotation_positions[rank]
        ax.annotate(
            f"Top {rank}: {row['cumulative_compound_coverage_fraction']:.1%}",
            xy=(rank, y),
            xytext=(tx, ty),
            fontsize=7.1,
            color=INK,
            arrowprops={"arrowstyle": "-", "color": SLATE, "linewidth": 0.65},
        )
    ax.text(
        0,
        1.03,
        f"Exact scaffolds; analysis N = {len(valid):,} canonical small molecules; {len(scaffold):,} unique scaffolds",
        transform=ax.transAxes,
        fontsize=7.2,
        color=SLATE,
    )
    ax.text(
        0,
        -0.23,
        "Both axes are logarithmic to show the high-frequency head and the scaffold long tail. Acyclic compounds are not rankable by Murcko scaffold.",
        transform=ax.transAxes,
        fontsize=7.0,
        color=SLATE,
    )
    for suffix, kwargs in {".svg": {}, ".pdf": {}, ".png": {"dpi": 600}}.items():
        fig.savefig(OUT / f"C3_scaffold_rank_frequency{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)

    stats_text = [
        "Figure\tC3 — Bemis–Murcko scaffold rank–frequency",
        f"Analysis N\t{len(valid):,} canonical formal small molecules with an existing exact cyclic Murcko scaffold",
        f"Unique exact scaffolds\t{len(scaffold):,}",
        f"Gini\t{gini(scaffold['unique_compound_count'].to_numpy()):.6f}",
        "Population rule\tExisting formal scope classes core_small_molecule and structure_resolved_small_molecule; extended_chemical_entity excluded.",
        "Scaffold rule\tExisting exact_murcko_scaffold_smiles only; no scaffold regeneration.",
        "Quarantine rule\tExisting ineligible structure records excluded from structure/scaffold analysis and reported in summary TSV.",
        "Pair mapping\tUnique proteins and pairs were left-joined from formal canonical protein–compound pair records.",
    ]
    (OUT / "C3_scaffold_rank_frequency_statistics.txt").write_text("\n".join(stats_text) + "\n", encoding="utf-8")
    (OUT / "C3_reproducibility_manifest.json").write_text(json.dumps({
        "analysis_unit": "one canonical compound",
        "figure": "C3_scaffold_rank_frequency",
        "inputs": [str(MASTER), str(CHEM), str(PAIR), str(QUARANTINE)],
        "exact_scaffold_field": "exact_murcko_scaffold_smiles",
        "outputs": ["C3_scaffold_rank_frequency.svg", "C3_scaffold_rank_frequency.pdf", "C3_scaffold_rank_frequency.png", "C3_scaffold_rank_frequency.tsv", "C3_scaffold_diversity_summary.tsv"],
        "log_axes": ["scaffold_rank", "unique_compound_count"],
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
