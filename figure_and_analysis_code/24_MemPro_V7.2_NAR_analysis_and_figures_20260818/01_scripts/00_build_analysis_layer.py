from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE = Path(os.environ["MEMPRO_DATA_ROOT"])
ROOT = Path(os.environ["MEMPRO_FIGURE_OUTPUT_ROOT"])
TABLES = SOURCE / "01_core_v72" / "01_release_tables"
ANALYSIS = ROOT / "02_analysis_data"
SOURCE_DATA = ROOT / "03_source_data"
QA = ROOT / "08_QA"

PLACEHOLDERS = {
    "unclassified",
    "no_specialist_classification",
    "family_defined_membrane_role_unresolved",
    "",
}


def read_tsv(name: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(TABLES / name, sep="\t", compression="infer", dtype=str, **kwargs)


def write_tsv(df: pd.DataFrame, name: str, gzip_output: bool = False) -> None:
    path = ANALYSIS / name
    if gzip_output and not path.name.endswith(".gz"):
        path = path.with_suffix(path.suffix + ".gz")
    df.to_csv(path, sep="\t", index=False, compression="gzip" if gzip_output else None)


def split_values(value: str | float | None) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    return [x.strip() for x in str(value).replace("|", ";").split(";") if x.strip()]


def broad_structure_class(row: pd.Series) -> str:
    """Deterministic broad class from standardised structure-derived fields.

    The purpose is a reproducible publication grouping, not a chemical ontology.
    Exact ChEBI classes remain available in the frozen source table.
    """
    formula = str(row.get("molecular_formula", "") or "")
    smiles = str(row.get("standard_smiles", "") or "")
    try:
        rings = int(float(row.get("ring_count", 0) or 0))
    except ValueError:
        rings = 0
    try:
        max_ring = int(float(row.get("max_ring_size", 0) or 0))
    except ValueError:
        max_ring = 0
    try:
        amides = int(float(row.get("amide_bond_count", 0) or 0))
    except ValueError:
        amides = 0
    try:
        mw = float(row.get("molecular_weight", 0) or 0)
    except ValueError:
        mw = 0
    try:
        tpsa = float(row.get("tpsa", 0) or 0)
    except ValueError:
        tpsa = 0

    has_carbon = bool(re.search(r"(^|[^A-Za-z])C(?:\d|$)", formula)) or bool(re.search(r"(^|[\[\(=#+\-])c", smiles))
    metal_tokens = ("[Fe", "[Zn", "[Mg", "[Mn", "[Cu", "[Co", "[Ni", "[Ca", "[Na", "[K", "[Pt", "[Pd")
    if any(t in smiles for t in metal_tokens):
        return "metal-containing"
    if not has_carbon:
        return "inorganic"
    if amides >= 3:
        return "peptidic/peptidomimetic"
    if max_ring >= 8:
        return "macrocyclic"
    if rings >= 2:
        return "polycyclic organic"
    if rings == 1 and mw > 0 and tpsa / mw >= 0.30 and smiles.count("O") >= 4:
        return "carbohydrate-like"
    if rings == 1:
        return "monocyclic organic"
    return "acyclic organic"


def build_protein_layer() -> tuple[pd.DataFrame, dict[str, set[str]]]:
    proteins = read_tsv("protein_master_v72.tsv.gz")
    assertions = read_tsv("protein_cross_classification_v72.tsv.gz")

    primary = assertions[assertions["primary_flag"].isin(["1", "true", "TRUE"])].copy()
    # Resolve accidental duplicate primary labels deterministically while retaining
    # the complete assertion table in the frozen source release.
    primary = primary.sort_values(
        ["target_uniprot_id", "classification_axis", "is_inferred", "classification_label"],
        ascending=[True, True, True, True],
    ).drop_duplicates(["target_uniprot_id", "classification_axis"], keep="first")
    wide = primary.pivot(index="target_uniprot_id", columns="classification_axis", values="classification_label")
    wide.columns = [f"primary_{c}" for c in wide.columns]
    wide = wide.reset_index()

    protein = proteins.merge(wide, left_on="canonical_uniprot_accession", right_on="target_uniprot_id", how="left")
    protein = protein.drop(columns=["target_uniprot_id"], errors="ignore")
    for axis in ["structural_family", "molecular_function", "biological_process", "membrane_role", "specialist_classification"]:
        col = f"primary_{axis}"
        protein[f"{axis}_informative"] = (~protein[col].fillna("").isin(PLACEHOLDERS)).astype(int)

    module_sets: dict[str, set[str]] = {
        "formal_protein": set(protein["canonical_uniprot_accession"]),
    }
    write_tsv(protein, "protein_analysis.tsv.gz", gzip_output=True)

    completeness = []
    for axis in ["structural_family", "molecular_function", "biological_process", "membrane_role", "specialist_classification"]:
        informative = int(protein[f"{axis}_informative"].sum())
        completeness.append({"classification_axis": axis, "informative": informative, "placeholder_or_missing": len(protein) - informative, "denominator": len(protein)})
    pd.DataFrame(completeness).to_csv(SOURCE_DATA / "F2E_classification_completeness.tsv", sep="\t", index=False)
    return protein, module_sets


def build_pair_and_coverage_layer(protein: pd.DataFrame, module_sets: dict[str, set[str]]) -> pd.DataFrame:
    pair_cols = [
        "pair_id", "target_uniprot_id", "compound_internal_id", "evidence_count", "best_evidence_tier",
        "distinct_database_count", "source_databases", "independent_experiment_count",
        "independent_structure_count", "independent_pubchem_assay_count", "independent_evidence_modality_count",
        "protein_release_role_v72",
    ]
    pairs = read_tsv("protein_compound_pair_v72.tsv.gz", usecols=pair_cols)
    pairs.to_csv(ANALYSIS / "pair_analysis.tsv.gz", sep="\t", index=False, compression="gzip")
    module_sets["compound_interaction"] = set(pairs["target_uniprot_id"])
    return pairs


def build_expression_layer(protein: pd.DataFrame, module_sets: dict[str, set[str]]) -> pd.DataFrame:
    usecols = [
        "target_uniprot_id", "source_dataset", "measurement_layer", "measurement_type", "tissue",
        "tissue_mapping_status", "cell_type", "cell_type_mapping_status", "value", "unit",
        "v7_mapping_state", "measurement_status_v71", "detection_status_v71",
        "mapped_measured_denominator_eligible_v71", "hpa_missingness_status_v72",
    ]
    layer_counter = Counter()
    rows_total = 0
    eligible_total = 0
    mapped_proteins: set[str] = set()
    consensus_chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(TABLES / "expression_measurement_v72.tsv.gz", sep="\t", compression="gzip", dtype=str, usecols=usecols, chunksize=250_000):
        rows_total += len(chunk)
        layer_counter.update(chunk["source_dataset"].fillna("UNSPECIFIED"))
        eligible = chunk["mapped_measured_denominator_eligible_v71"].eq("1")
        eligible_total += int(eligible.sum())
        mapped_proteins.update(chunk.loc[eligible, "target_uniprot_id"].dropna())
        take = chunk[
            chunk["source_dataset"].eq("HPA_consensus_tissue_RNA")
            & eligible
            & chunk["measurement_status_v71"].eq("MEASURED")
            & chunk["tissue"].notna()
        ].copy()
        if not take.empty:
            take["value"] = pd.to_numeric(take["value"], errors="coerce")
            take = take.dropna(subset=["value"])
            consensus_chunks.append(take[["target_uniprot_id", "tissue", "value", "detection_status_v71"]])

    consensus = pd.concat(consensus_chunks, ignore_index=True)
    # Multiple entries per protein/tissue are reduced by median to avoid one tissue
    # receiving extra weight from technical duplication.
    consensus = consensus.groupby(["target_uniprot_id", "tissue"], as_index=False)["value"].median()
    tissue_count = consensus["tissue"].nunique()

    tau_rows = []
    for protein_id, grp in consensus.groupby("target_uniprot_id", sort=False):
        values = grp.set_index("tissue")["value"].reindex(sorted(consensus["tissue"].unique()), fill_value=0).to_numpy(float)
        maximum = float(np.nanmax(values)) if len(values) else 0.0
        if maximum <= 0 or len(values) <= 1:
            tau = np.nan
        else:
            tau = float(np.sum(1.0 - values / maximum) / (len(values) - 1))
        tau_rows.append({
            "target_uniprot_id": protein_id,
            "tau_consensus_tissue_rna": tau,
            "tissues_in_denominator": len(values),
            "maximum_nTPM": maximum,
            "top_tissue": grp.loc[grp["value"].idxmax(), "tissue"],
        })
    tau_df = pd.DataFrame(tau_rows).merge(
        protein[["canonical_uniprot_accession", "primary_membrane_role", "membrane_class_v7", "membrane_evidence_level_v7"]],
        left_on="target_uniprot_id", right_on="canonical_uniprot_accession", how="left",
    ).drop(columns=["canonical_uniprot_accession"])
    tau_df.to_csv(ANALYSIS / "expression_tau.tsv.gz", sep="\t", index=False, compression="gzip")
    tau_df.to_csv(SOURCE_DATA / "F2D_tissue_specificity_tau.tsv", sep="\t", index=False)
    module_sets["expression_mapped_measured"] = mapped_proteins

    summary = pd.DataFrame([
        {"metric": "all_expression_records", "count": rows_total},
        {"metric": "mapped_measured_denominator_eligible", "count": eligible_total},
        {"metric": "source_unresolved_missing_or_not_measured", "count": 103709},
        {"metric": "unmapped_or_ambiguous", "count": 22701},
        {"metric": "consensus_tissue_rna_tissues", "count": tissue_count},
        {"metric": "proteins_with_tau", "count": int(tau_df["tau_consensus_tissue_rna"].notna().sum())},
    ])
    summary.to_csv(SOURCE_DATA / "F1_expression_denominator_audit.tsv", sep="\t", index=False)
    pd.DataFrame(layer_counter.most_common(), columns=["source_dataset", "record_count"]).to_csv(SOURCE_DATA / "S_expression_source_layers.tsv", sep="\t", index=False)
    return tau_df


def build_disease_layer(module_sets: dict[str, set[str]]) -> pd.DataFrame:
    relations = read_tsv("protein_disease_relation_v72.tsv.gz")
    anatomy = read_tsv("disease_anatomical_system_v72.tsv.gz")
    therapeutic = read_tsv("disease_therapeutic_area_v72.tsv.gz")
    master_ids = set(read_tsv("disease_entity_master_v72.tsv.gz", usecols=["canonical_disease_id"])["canonical_disease_id"])
    relation_ids = set(relations["canonical_disease_id"])

    anatomy = anatomy[anatomy["canonical_disease_id"].isin(relation_ids & master_ids)].copy()
    therapeutic = therapeutic[therapeutic["canonical_disease_id"].isin(relation_ids & master_ids)].copy()
    anatomy_map = anatomy.groupby("canonical_disease_id")["anatomical_system_name"].agg(lambda x: ";".join(sorted(set(x)))).to_dict()
    therapeutic_map = therapeutic.groupby("canonical_disease_id")["therapeutic_area_name"].agg(lambda x: ";".join(sorted(set(x)))).to_dict()
    relations["anatomical_systems"] = relations["canonical_disease_id"].map(anatomy_map).fillna("")
    relations["therapeutic_areas"] = relations["canonical_disease_id"].map(therapeutic_map).fillna("")
    relations["anatomy_mapped"] = relations["anatomical_systems"].ne("").astype(int)
    relations["therapeutic_area_mapped"] = relations["therapeutic_areas"].ne("").astype(int)
    relations.to_csv(ANALYSIS / "disease_relation_analysis.tsv.gz", sep="\t", index=False, compression="gzip")
    module_sets["disease"] = set(relations["target_uniprot_id"])

    coverage = pd.DataFrame([
        {"metric": "formal_diseases", "count": len(relation_ids)},
        {"metric": "diseases_with_anatomy", "count": len(relation_ids & set(anatomy["canonical_disease_id"]))},
        {"metric": "diseases_without_anatomy", "count": len(relation_ids - set(anatomy["canonical_disease_id"]))},
        {"metric": "diseases_with_therapeutic_area", "count": len(relation_ids & set(therapeutic["canonical_disease_id"]))},
        {"metric": "diseases_without_therapeutic_area", "count": len(relation_ids - set(therapeutic["canonical_disease_id"]))},
        {"metric": "historical_anatomy_ids_outside_formal_master", "count": len(set(read_tsv("disease_anatomical_system_v72.tsv.gz", usecols=["canonical_disease_id"])["canonical_disease_id"]) - master_ids)},
        {"metric": "historical_therapeutic_ids_outside_formal_master", "count": len(set(read_tsv("disease_therapeutic_area_v72.tsv.gz", usecols=["canonical_disease_id"])["canonical_disease_id"]) - master_ids)},
    ])
    coverage.to_csv(SOURCE_DATA / "F5_disease_mapping_coverage.tsv", sep="\t", index=False)
    anatomy.to_csv(ANALYSIS / "disease_anatomy_formal_intersection.tsv.gz", sep="\t", index=False, compression="gzip")
    therapeutic.to_csv(ANALYSIS / "disease_therapeutic_formal_intersection.tsv.gz", sep="\t", index=False, compression="gzip")
    return relations


def build_site_and_localization_layers(module_sets: dict[str, set[str]]) -> None:
    site_cols = [
        "site_record_id", "evidence_id", "target_uniprot_id", "compound_internal_id", "source_database",
        "mapping_fraction", "site_residue_status", "coordinate_mapping_status", "coordinate_confidence",
        "chain_assignment_status", "membrane_side", "coordinate_docking_eligible", "site_tier_v72",
        "structure_ligand_relationship",
    ]
    sites = read_tsv("interaction_site_v72.tsv.gz", usecols=site_cols)
    sites.to_csv(ANALYSIS / "site_analysis.tsv.gz", sep="\t", index=False, compression="gzip")
    module_sets["binding_site"] = set(sites["target_uniprot_id"])

    localization = read_tsv("subcellular_localization_v72.tsv.gz", usecols=["target_uniprot_id", "location_term", "go_id", "location_role", "direct_observation_flag", "source_dataset", "hpa_reliability"])
    localization.to_csv(ANALYSIS / "localization_analysis.tsv.gz", sep="\t", index=False, compression="gzip")
    module_sets["subcellular_localization"] = set(localization["target_uniprot_id"])


def build_compound_layer(pairs: pd.DataFrame) -> None:
    linked = set(pairs["compound_internal_id"])
    usecols = [
        "compound_internal_id", "preferred_name", "standard_smiles", "standard_inchikey", "molecular_formula",
        "molecular_weight", "xlogp", "tpsa", "hbond_donor_count", "hbond_acceptor_count",
        "rotatable_bond_count", "ring_count", "max_ring_size", "amide_bond_count", "formal_charge",
        "development_status", "is_approved_drug", "is_clinical_candidate", "is_endogenous_ligand",
        "is_natural_product", "is_chemical_probe", "protein_target_count", "identity_confidence",
    ]
    out_path = ANALYSIS / "compound_features.tsv.gz"
    first = True
    registry = 0
    linked_rows = 0
    class_counter = Counter()
    for chunk in pd.read_csv(TABLES / "compound_master_v72.tsv.gz", sep="\t", compression="gzip", dtype=str, usecols=usecols, chunksize=50_000):
        registry += len(chunk)
        chunk["interaction_linked"] = chunk["compound_internal_id"].isin(linked).astype(int)
        linked_rows += int(chunk["interaction_linked"].sum())
        chunk["broad_structure_class_recomputed"] = chunk.apply(broad_structure_class, axis=1)
        class_counter.update(chunk["broad_structure_class_recomputed"])
        chunk.to_csv(out_path, sep="\t", index=False, compression="gzip", mode="wt" if first else "at", header=first)
        first = False

    audit = pd.DataFrame([
        {"metric": "compound_registry", "count": registry},
        {"metric": "interaction_linked_compounds", "count": linked_rows},
        {"metric": "registry_without_formal_pair", "count": registry - linked_rows},
    ])
    audit.to_csv(SOURCE_DATA / "F3_compound_universe_audit.tsv", sep="\t", index=False)
    pd.DataFrame(class_counter.most_common(), columns=["broad_structure_class", "count"]).to_csv(SOURCE_DATA / "F3_structure_class_registry.tsv", sep="\t", index=False)


def write_module_coverage(protein: pd.DataFrame, module_sets: dict[str, set[str]]) -> None:
    all_ids = set(protein["canonical_uniprot_accession"])
    rows = []
    for module, ids in module_sets.items():
        rows.append({"module": module, "protein_count": len(ids & all_ids), "formal_protein_denominator": len(all_ids), "coverage_fraction": len(ids & all_ids) / len(all_ids)})
    pd.DataFrame(rows).sort_values("protein_count", ascending=False).to_csv(SOURCE_DATA / "F1D_module_protein_coverage.tsv", sep="\t", index=False)


def main() -> None:
    for directory in [ANALYSIS, SOURCE_DATA, QA]:
        directory.mkdir(parents=True, exist_ok=True)
    protein, module_sets = build_protein_layer()
    pairs = build_pair_and_coverage_layer(protein, module_sets)
    build_expression_layer(protein, module_sets)
    build_disease_layer(module_sets)
    build_site_and_localization_layers(module_sets)
    build_compound_layer(pairs)
    write_module_coverage(protein, module_sets)

    report = {
        "source_release": "MemPro_V7.2",
        "source_modified": False,
        "protein_rows": len(protein),
        "pair_rows": len(pairs),
        "analysis_status": "BUILT",
    }
    (QA / "analysis_layer_build.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
