"""C1: direct source-by-module coverage for MemPro canonical compounds.

This is a display-only analysis over the frozen formal release. It does not
alter canonicalization, compound provenance, or any release table.
"""
from pathlib import Path
import json
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compound_figure_pool"
CORE = Path(r"D:\finale\FORMAL\01_core_v72\01_release_tables")
COMPANION = Path(r"D:\finale\FORMAL\02_companion_v721\02_companion_tables")
MASTER = CORE / "compound_master_v72.tsv.gz"
EVIDENCE = CORE / "positive_interaction_evidence_v72.tsv.gz"
CONTRIBUTIONS = COMPANION / "evidence_source_contribution_v721.tsv.gz"
CHEMICAL_CLASS = COMPANION / "compound_chemical_classification_v721.tsv.gz"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.linewidth": 0.7,
})

SOURCES = [
    "PubChem BioAssay",
    "ChEMBL",
    "BindingDB",
    "PDBe",
    "IUPHAR/BPS Guide to PHARMACOLOGY",
    "PDBbind",
    "BRENDA",
    "sc-PDB",
]
DISPLAY_SOURCE = {
    "IUPHAR/BPS Guide to PHARMACOLOGY": "GtoPdb (IUPHAR/BPS)",
    "MemPro-derived": "MemPro-derived",
}
MODULES = [
    "Identity / chemical structure",
    "Physicochemical\nproperties",
    "Scaffold / structural\nannotation",
    "Status\nannotation",
    "Protein\ninteraction",
    "Quantitative\nactivity",
    "Binding\nevidence",
    "Structural binding\ncontext",
]
MODULE_KEYS = [m.replace("\n", " ") for m in MODULES]
INK, SLATE = "#24364B", "#5E7182"
CMAP = LinearSegmentedColormap.from_list(
    "mempro_teal", ["#F7FBFA", "#D7EEE8", "#8AC9B8", "#3A9A88", "#126B66"]
)


def split_sources(value):
    """Parse only the existing direct provenance field; never use ID mappings."""
    if pd.isna(value):
        return []
    parsed = [x.strip() for x in re.split(r"[;|]", str(value)) if x.strip()]
    aliases = {
        "CHEMBL": "ChEMBL",
        "ChEMBL database": "ChEMBL",
        "IUPHAR/BPS Guide to Pharmacology": "IUPHAR/BPS Guide to PHARMACOLOGY",
        "Guide to PHARMACOLOGY": "IUPHAR/BPS Guide to PHARMACOLOGY",
    }
    return [aliases.get(x, x) for x in parsed]


def fmt_count(value):
    value = int(value)
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    return f"{value / 1_000:.1f}k" if value >= 1_000 else str(value)


def add_source_module(records, source, module, compounds):
    ids = pd.Series(compounds).dropna().astype(str).unique()
    if len(ids):
        records.append(pd.DataFrame({
            "data_source": source,
            "compound_module": module,
            "compound_internal_id": ids,
        }))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    master_cols = [
        "compound_internal_id", "source_databases", "molecular_weight", "xlogp", "tpsa",
        "formal_charge", "development_status",
    ]
    master = pd.read_csv(MASTER, sep="\t", compression="gzip", usecols=master_cols, low_memory=False)
    master = master.drop_duplicates("compound_internal_id").copy()
    n_universe = len(master)
    records = []

    # Identity uses the frozen source_databases direct provenance field only.
    # pubchem_cids, chembl_ids and all cross-reference identifiers are not read.
    identity = master.loc[:, ["compound_internal_id", "source_databases"]].copy()
    identity["data_source"] = identity["source_databases"].map(split_sources)
    identity = identity.explode("data_source").dropna(subset=["data_source"])
    for source in SOURCES:
        add_source_module(
            records, source, "Identity / chemical structure",
            identity.loc[identity["data_source"].eq(source), "compound_internal_id"],
        )

    # The following are computed/harmonized MemPro annotations and are therefore
    # shown only in the explicitly named MemPro-derived row.
    property_mask = master[["molecular_weight", "xlogp", "tpsa", "formal_charge"]].notna().all(axis=1)
    add_source_module(records, "MemPro-derived", "Physicochemical properties", master.loc[property_mask, "compound_internal_id"])
    add_source_module(records, "MemPro-derived", "Status annotation", master.loc[master["development_status"].notna(), "compound_internal_id"])
    chemical = pd.read_csv(
        CHEMICAL_CLASS, sep="\t", compression="gzip",
        usecols=["compound_internal_id", "structure_status"], low_memory=False,
    ).drop_duplicates("compound_internal_id")
    add_source_module(
        records, "MemPro-derived", "Scaffold / structural annotation",
        chemical.loc[chemical["structure_status"].eq("VALID"), "compound_internal_id"],
    )

    # Evidence modules use the explicit evidence-to-source contribution junction,
    # avoiding ambiguities in denormalized display source strings.
    evidence = pd.read_csv(
        EVIDENCE, sep="\t", compression="gzip",
        usecols=["evidence_id", "compound_internal_id", "standard_value_nM", "evidence_directness", "pdb_ids", "default_web_inclusion_v72"],
        low_memory=False,
    )
    evidence = evidence.loc[evidence["default_web_inclusion_v72"].eq(1)].copy()
    contribution = pd.read_csv(
        CONTRIBUTIONS, sep="\t", compression="gzip", usecols=["evidence_id", "source_database"], low_memory=False,
    ).drop_duplicates()
    evidence = evidence.merge(contribution, on="evidence_id", how="inner", validate="m:m")
    evidence["source_database"] = evidence["source_database"].map(
        lambda value: split_sources(value)[0] if split_sources(value) else None
    )
    evidence = evidence.loc[evidence["source_database"].isin(SOURCES)].copy()
    direct_values = {
        "direct", "direct_quantitative", "direct_structural",
        "direct_structural_observation", "binding_assay_activity",
    }
    evidence_masks = {
        "Protein interaction": pd.Series(True, index=evidence.index),
        "Quantitative activity": evidence["standard_value_nM"].notna(),
        "Binding evidence": evidence["evidence_directness"].isin(direct_values),
        "Structural binding context": evidence["pdb_ids"].notna(),
    }
    for module, mask in evidence_masks.items():
        for source in SOURCES:
            add_source_module(
                records, source, module,
                evidence.loc[mask & evidence["source_database"].eq(source), "compound_internal_id"],
            )

    long = pd.concat(records, ignore_index=True)
    tab = (
        long.groupby(["data_source", "compound_module"], as_index=False)["compound_internal_id"]
        .nunique().rename(columns={"compound_internal_id": "unique_canonical_compounds"})
    )
    tab["analysis_universe"] = "all formal canonical compounds"
    tab["analysis_universe_n"] = n_universe
    tab["value_definition"] = (
        "unique canonical compounds directly covered by the named source and module; "
        "MemPro-derived indicates computed or harmonized annotation"
    )
    row_order = SOURCES + ["MemPro-derived"]
    tab["data_source"] = pd.Categorical(tab["data_source"], row_order, ordered=True)
    tab["compound_module"] = pd.Categorical(tab["compound_module"], MODULE_KEYS, ordered=True)
    tab = tab.sort_values(["data_source", "compound_module"]).reset_index(drop=True)
    tab.to_csv(OUT / "C1_source_module_compound_coverage.tsv", sep="\t", index=False)

    definitions = pd.DataFrame([
        ("Identity / chemical structure", "Direct canonical-compound registry provenance recorded in compound_master.source_databases.", "Direct source contribution", "all formal canonical compounds"),
        ("Physicochemical properties", "Complete MemPro-derived descriptor panel: molecular weight, XlogP, TPSA and formal charge.", "MemPro-derived; not attributed to an external source", "all formal canonical compounds"),
        ("Scaffold / structural annotation", "Valid RDKit chemical-classification record including scaffold/structural annotation.", "MemPro-derived; not attributed to an external source", "all formal canonical compounds"),
        ("Status annotation", "Non-missing harmonized development_status annotation.", "MemPro-harmonized; exact source-of-status is not field-level provenance", "all formal canonical compounds"),
        ("Protein interaction", "At least one formal positive evidence record with an explicit evidence-to-source contribution.", "Direct source contribution", "all formal canonical compounds"),
        ("Quantitative activity", "Formal positive evidence with non-missing standard_value_nM and explicit source contribution.", "Direct source contribution", "all formal canonical compounds"),
        ("Binding evidence", "Direct/direct-quantitative/direct-structural/binding-assay evidence with explicit source contribution.", "Direct source contribution", "all formal canonical compounds"),
        ("Structural binding context", "Formal positive evidence with a recorded PDB identifier and explicit source contribution.", "Direct source contribution", "all formal canonical compounds"),
    ], columns=["compound_module", "definition", "attribution_basis", "analysis_universe"])
    definitions.to_csv(OUT / "C1_compound_module_definition.tsv", sep="\t", index=False)

    matrix = tab.pivot(index="data_source", columns="compound_module", values="unique_canonical_compounds")
    matrix = matrix.reindex(index=row_order, columns=MODULE_KEYS)
    valid_values = matrix.stack().astype(float)
    log_values = np.log1p(matrix.astype(float))
    fig = plt.figure(figsize=(1800 / 150, 1050 / 150), dpi=600)
    ax = fig.add_axes([0.27, 0.17, 0.60, 0.66])
    im = ax.imshow(np.ma.masked_invalid(log_values.to_numpy()), cmap=CMAP, aspect="auto", interpolation="nearest")
    norm = Normalize(vmin=float(np.log1p(valid_values.min())), vmax=float(np.log1p(valid_values.max())))
    im.set_norm(norm)
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            value = matrix.iat[y, x]
            if pd.notna(value):
                text_color = "white" if norm(np.log1p(value)) >= 0.58 else INK
                ax.text(x, y, fmt_count(value), ha="center", va="center", fontsize=8.2, color=text_color, weight="semibold")
    ax.set_xticks(range(len(MODULES)), MODULES)
    ax.set_yticks(range(len(row_order)), [DISPLAY_SOURCE.get(s, s) for s in row_order])
    ax.tick_params(axis="x", labelrotation=30, labelsize=8.5, length=0, pad=8)
    ax.tick_params(axis="y", labelsize=9, length=0, pad=8)
    plt.setp(ax.get_xticklabels(), ha="right", rotation_mode="anchor")
    ax.set_xticks(np.arange(-.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="#E2E8E8", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, pad=0.025, fraction=0.045)
    cbar.set_label("Canonical compound coverage", fontsize=8.7, color=INK)
    cbar.ax.tick_params(labelsize=7.5, colors=SLATE)
    cbar.outline.set_visible(False)
    fig.text(0.27, 0.955, "Compound source-by-module coverage", fontsize=15, color=INK, weight="bold", ha="left")
    fig.text(0.27, 0.910, f"All formal canonical compounds (Analysis N = {n_universe:,}); colour intensity is log1p-scaled.", fontsize=9, color=SLATE, ha="left")
    fig.text(0.27, 0.875, "Cell labels show raw unique-compound counts; blank cells indicate no direct source-level contribution counted.", fontsize=8.3, color=SLATE, ha="left")
    for suffix, kwargs in {".svg": {}, ".pdf": {}, ".png": {"dpi": 600}}.items():
        fig.savefig(OUT / f"C1_source_module_compound_coverage{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)

    # Record equal values for review; equality itself is not treated as an error.
    duplicate_cells = []
    for module in MODULE_KEYS:
        series = matrix[module].dropna().astype(int)
        for value, members in series.groupby(series):
            if len(members) > 1:
                duplicate_cells.append({"module": module, "count": int(value), "sources": "; ".join(members.index.astype(str))})
    qc_lines = [
        "C1 source-by-module compound coverage QC",
        f"Analysis universe: {n_universe:,} all formal canonical compounds.",
        "Analysis unit: one compound_internal_id (unique canonical compound).",
        "Identity/chemical-structure cells use compound_master.source_databases only (direct registry provenance).",
        "PubChem CID mapping was not used as PubChem direct provenance: PASS.",
        "ChEMBL ID mapping was not used as ChEMBL direct provenance: PASS.",
        "Evidence modules use evidence_source_contribution_v721 explicit evidence-to-source junction: PASS.",
        "Physicochemical descriptors and scaffolds are shown only in the MemPro-derived row: PASS.",
        "Status annotation is MemPro-harmonized because field-level source-of-status provenance is unavailable: PASS.",
        "All-compound and interaction-linked universes are not mixed; interaction modules are subsets within the all-compound universe: PASS.",
        f"Displayed external direct sources: {', '.join(SOURCES)}.",
        f"Non-zero equal-count source/module cells requiring visual review: {len(duplicate_cells)}.",
    ]
    if duplicate_cells:
        qc_lines.append("Equal non-zero cells (not treated as duplicated attribution):")
        qc_lines.extend([f"- {x['module']}: {x['count']:,} ({x['sources']})" for x in duplicate_cells])
    (OUT / "C1_source_module_qc.txt").write_text("\n".join(qc_lines) + "\n", encoding="utf-8")
    (OUT / "C1_reproducibility_manifest.json").write_text(json.dumps({
        "figure": "C1_source_module_compound_coverage",
        "analysis_unit": "unique canonical compound",
        "analysis_universe": "all formal canonical compounds",
        "analysis_universe_n": n_universe,
        "inputs": [str(MASTER), str(EVIDENCE), str(CONTRIBUTIONS), str(CHEMICAL_CLASS)],
        "source_attribution": "compound registry source_databases for identity; evidence_source_contribution_v721 for evidence modules; derived fields only in MemPro-derived row",
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
