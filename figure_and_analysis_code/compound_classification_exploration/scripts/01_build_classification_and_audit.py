from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


FORMAL = Path(r"D:\finale\FORMAL")
OUT = Path(r"D:\finale\compound_classification_exploration")
for sub in [
    "01_data_audit", "02_classification_tables", "03_scaffold", "04_heatmaps",
    "05_chemical_space", "06_target_breadth", "07_optional_3d",
    "08_contact_sheet", "09_review", "scripts", "logs", "qa",
]:
    (OUT / sub).mkdir(parents=True, exist_ok=True)

CM = FORMAL / "01_core_v72/01_release_tables/compound_master_v72.tsv.gz"
PAIR = FORMAL / "01_core_v72/01_release_tables/protein_compound_pair_v72.tsv.gz"
CC = FORMAL / "02_companion_v721/02_companion_tables/compound_chemical_classification_v721.tsv.gz"
ROLE = FORMAL / "03_publication_repairs/protein_membrane_role_FORMAL.tsv"
REPAIR = FORMAL / "03_publication_repairs/COMPOUND_STRUCTURE_REPAIR_AND_QUARANTINE.tsv"

SEED = 42


def sha_id(prefix: str, value: str) -> str:
    if not value:
        return f"{prefix}-ACYCLIC"
    return f"{prefix}-{hashlib.sha256(value.encode('utf-8')).hexdigest()[:16].upper()}"


def first_token(value: object) -> str:
    if pd.isna(value) or not str(value).strip():
        return "NOT_AVAILABLE"
    return str(value).split("|")[0].strip()


METALS = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 56,
    57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
    74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
}
PEPTIDE_BOND = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
NUCLEOBASE = Chem.MolFromSmarts("[nH0,nH1]1~[c,n]~[n,c]~[c,n]~[c,n]1")
PHOSPHATE = Chem.MolFromSmarts("P(=O)(O)O")
LONG_CHAIN = Chem.MolFromSmarts("[C;!R]~[C;!R]~[C;!R]~[C;!R]~[C;!R]~[C;!R]~[C;!R]~[C;!R]")


def chemical_regime(mol: Chem.Mol | None, desc: dict[str, float]) -> tuple[str, str]:
    """Transparent local regime; never presented as ClassyFire/ChEBI ontology."""
    if mol is None:
        return "undefined / unsupported structure", "LOCAL_RULE_UNSUPPORTED"
    atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    carbon = atoms.count(6)
    if carbon == 0:
        return "inorganic / no-carbon compound", "LOCAL_RULE_NO_CARBON"
    metal_atoms = {i for i, z in enumerate(atoms) if z in METALS}
    if metal_atoms:
        metal_carbon_bond = any(
            (b.GetBeginAtomIdx() in metal_atoms and b.GetEndAtom().GetAtomicNum() == 6)
            or (b.GetEndAtomIdx() in metal_atoms and b.GetBeginAtom().GetAtomicNum() == 6)
            for b in mol.GetBonds()
        )
        label = "organometallic compound" if metal_carbon_bond else "metal-containing compound"
        return label, "LOCAL_RULE_METAL"
    if desc["max_ring_size"] >= 12:
        return "macrocyclic compound", "LOCAL_RULE_MACROCYCLE_RING_GE_12"
    peptide_bonds = len(mol.GetSubstructMatches(PEPTIDE_BOND))
    if peptide_bonds >= 3 and desc["hba"] >= 3 and desc["molecular_weight"] >= 250:
        return "peptide / peptide-like", "LOCAL_RULE_PEPTIDE_3_AMIDE"
    if mol.HasSubstructMatch(PHOSPHATE) and mol.HasSubstructMatch(NUCLEOBASE) and desc["oxygen_count"] >= 3:
        return "nucleotide / nucleoside-like", "LOCAL_RULE_PHOSPHATE_NUCLEOBASE"
    if (
        desc["oxygen_count"] >= 4 and desc["fraction_csp3"] >= 0.35
        and desc["ring_count"] >= 1 and desc["molecular_weight"] <= 1200
    ):
        return "carbohydrate / glycoside-like", "LOCAL_RULE_O_RICH_SP3_RING"
    if mol.HasSubstructMatch(LONG_CHAIN) and desc["ring_count"] <= 1 and desc["heteroatom_count"] <= 8:
        return "lipid / lipid-like", "LOCAL_RULE_LONG_ALIPHATIC_CHAIN"
    return "organic small molecule", "LOCAL_RULE_DEFAULT_ORGANIC"


def mol_descriptors(mol: Chem.Mol) -> dict[str, float]:
    rings = mol.GetRingInfo().AtomRings()
    return {
        "molecular_weight": float(Descriptors.MolWt(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
        "aromatic_ring_count": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "fraction_csp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
        "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
        "max_ring_size": max((len(r) for r in rings), default=0),
        "oxygen_count": sum(a.GetAtomicNum() == 8 for a in mol.GetAtoms()),
        "nitrogen_count": sum(a.GetAtomicNum() == 7 for a in mol.GetAtoms()),
        "heteroatom_count": sum(a.GetAtomicNum() not in (1, 6) for a in mol.GetAtoms()),
        "amide_bond_count": len(mol.GetSubstructMatches(PEPTIDE_BOND)),
    }


def canonical_generic_scaffold(mol: Chem.Mol) -> tuple[str, str, str]:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        return "", "", "ACYCLIC"
    exact = Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)
    try:
        generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        generic_smiles = Chem.MolToSmiles(generic, canonical=True, isomericSmiles=False)
        return exact, generic_smiles, "CURRENT_RDKIT"
    except Exception as exc:
        return exact, "", f"CURRENT_RDKIT_GENERIC_FAILED:{type(exc).__name__}"


def write_md(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def main() -> None:
    cm_cols = [
        "compound_internal_id", "preferred_name", "standard_smiles", "standard_inchi",
        "standard_inchikey", "molecular_formula", "molecular_weight", "xlogp", "tpsa",
        "hbond_donor_count", "hbond_acceptor_count", "rotatable_bond_count",
        "formal_charge", "heavy_atom_count", "ring_count", "max_ring_size",
        "computed_structural_class", "development_status", "is_approved_drug",
        "is_clinical_candidate", "is_endogenous_ligand", "is_natural_product",
        "is_chemical_probe", "chebi_ids", "chebi_direct_classes", "chebi_roles",
    ]
    cm = pd.read_csv(CM, sep="\t", usecols=cm_cols, low_memory=False)
    pairs = pd.read_csv(PAIR, sep="\t", usecols=["target_uniprot_id", "compound_internal_id"])
    pairs = pairs.drop_duplicates()
    linked_ids = pd.Index(pairs["compound_internal_id"].unique())
    linked = cm[cm["compound_internal_id"].isin(linked_ids)].copy()
    if len(cm) != 646_670 or len(linked_ids) != 240_055 or len(pairs) != 529_168:
        raise RuntimeError(f"Frozen denominator mismatch: registry={len(cm)}, linked={len(linked_ids)}, pairs={len(pairs)}")

    cc = pd.read_csv(CC, sep="\t", low_memory=False)
    linked = linked.merge(cc, on="compound_internal_id", how="left", validate="one_to_one")
    linked = linked.rename(columns={
        "molecular_weight": "molecular_weight_formal",
        "tpsa": "tpsa_formal",
        "formal_charge": "formal_charge_formal",
        "heavy_atom_count": "heavy_atom_count_formal",
        "ring_count": "ring_count_formal",
        "max_ring_size": "max_ring_size_formal",
    })
    repair = pd.read_csv(REPAIR, sep="\t")
    repair = repair.set_index("compound_internal_id")

    # Apply publication repair/quarantine as an analysis overlay only.
    repaired_id = "HMPD-CMPD-0090098"
    conflict_id = "HMPD-CMPD-0192136"
    linked["analysis_structure_status"] = linked["structure_status"].fillna("MISSING_STRUCTURE_STATUS")
    linked["analysis_smiles"] = linked["canonical_smiles_rdkit"].fillna(linked["standard_smiles"])
    if repaired_id in set(linked["compound_internal_id"]):
        linked.loc[linked.compound_internal_id == repaired_id, "analysis_smiles"] = repair.loc[repaired_id, "reconstructed_canonical_smiles"]
        linked.loc[linked.compound_internal_id == repaired_id, "analysis_structure_status"] = "VALID_REPAIRED_FROM_INCHI"
    if conflict_id in set(linked["compound_internal_id"]):
        linked.loc[linked.compound_internal_id == conflict_id, "analysis_structure_status"] = "STRUCTURE_IDENTITY_CONFLICT_QUARANTINED"
    linked["eligible_structure_analysis"] = linked["analysis_structure_status"].isin(["VALID", "VALID_REPAIRED_FROM_INCHI"])

    results: list[dict[str, object]] = []
    eligible = linked[linked["eligible_structure_analysis"]].copy().reset_index(drop=True)
    print(f"Computing RDKit descriptors/scaffolds for {len(eligible):,} compounds", flush=True)
    for i, row in eligible.iterrows():
        smi = row["analysis_smiles"]
        mol = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
        if mol is None:
            results.append({"compound_internal_id": row.compound_internal_id, "rdkit_parse_status": "FAILED"})
            continue
        can = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        desc = mol_descriptors(mol)
        exact, generic, generic_status = canonical_generic_scaffold(mol)
        if not generic and exact and pd.notna(row.get("generic_murcko_scaffold_smiles")):
            generic = str(row["generic_murcko_scaffold_smiles"])
            generic_status += ";V721_FALLBACK"
        regime, rule = chemical_regime(mol, desc)
        results.append({
            "compound_internal_id": row.compound_internal_id,
            "rdkit_parse_status": "VALID",
            "canonical_smiles_recalculated": can,
            "exact_murcko_smiles_recalculated": exact,
            "exact_scaffold_id": sha_id("EMS", exact),
            "generic_murcko_smiles_recalculated": generic,
            "generic_scaffold_id": sha_id("GMS", generic),
            "generic_scaffold_recalc_status": generic_status,
            "chemical_regime_local": regime,
            "chemical_regime_rule": rule,
            **desc,
        })
        if (i + 1) % 25_000 == 0:
            print(f"  processed {i+1:,}/{len(eligible):,}", flush=True)

    calc = pd.DataFrame(results)
    eligible = eligible.merge(calc, on="compound_internal_id", how="left", validate="one_to_one")
    parse_failed = int((eligible["rdkit_parse_status"] != "VALID").sum())
    if parse_failed:
        raise RuntimeError(f"Unexpected RDKit parse failures after eligibility filter: {parse_failed}")

    eligible["authoritative_chebi_direct_class"] = eligible["chebi_direct_classes"].map(first_token)
    eligible["chemical_superclass"] = "NOT_AVAILABLE_NEEDS_EXTERNAL_MAPPING"
    eligible["chemical_class"] = eligible["authoritative_chebi_direct_class"]
    eligible["chemical_subclass"] = "NOT_AVAILABLE_NEEDS_EXTERNAL_MAPPING"
    eligible["direct_parent"] = "NOT_AVAILABLE_NEEDS_EXTERNAL_MAPPING"
    eligible["taxonomy_assignment_status"] = np.where(
        eligible["authoritative_chebi_direct_class"] != "NOT_AVAILABLE",
        "AUTHORITATIVE_CHEBI_DIRECT_CLASS_ONLY",
        "NEEDS_EXTERNAL_ONTOLOGY_MAPPING",
    )

    # Descriptor comparison: report, never overwrite the formal fields.
    comparisons = {
        "molecular_weight": "molecular_weight_formal",
        "logp": "xlogp",
        "tpsa": "tpsa_formal",
        "hbd": "hbond_donor_count",
        "hba": "hbond_acceptor_count",
        "rotatable_bonds": "rotatable_bond_count",
        "ring_count": "ring_count_formal",
        "formal_charge": "formal_charge_formal",
        "heavy_atom_count": "heavy_atom_count_formal",
    }
    discrepancy = {}
    for new, old in comparisons.items():
        a = pd.to_numeric(eligible[new], errors="coerce")
        b = pd.to_numeric(eligible[old], errors="coerce")
        valid = a.notna() & b.notna()
        tol = 0.05 if new in {"molecular_weight", "logp", "tpsa"} else 0.0
        discrepancy[new] = {
            "compared": int(valid.sum()),
            "missing_formal": int(b.isna().sum()),
            "different_beyond_tolerance": int(((a - b).abs() > tol)[valid].sum()),
            "tolerance": tol,
        }

    # Stable classification/scaffold tables.
    class_cols = [
        "compound_internal_id", "analysis_structure_status", "eligible_structure_analysis",
        "taxonomy_assignment_status", "authoritative_chebi_direct_class", "chemical_superclass",
        "chemical_class", "chemical_subclass", "direct_parent", "chemical_regime_local",
        "chemical_regime_rule", "computed_structural_class", "structural_tags",
    ]
    eligible[class_cols].to_csv(
        OUT / "02_classification_tables/COMPOUND_MULTI_LAYER_CLASSIFICATION.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    scaffold_cols = [
        "compound_internal_id", "canonical_smiles_recalculated",
        "exact_murcko_smiles_recalculated", "exact_scaffold_id",
        "generic_murcko_smiles_recalculated", "generic_scaffold_id",
        "generic_scaffold_recalc_status",
        "ring_count", "aromatic_ring_count", "max_ring_size",
    ]
    eligible[scaffold_cols].to_csv(
        OUT / "03_scaffold/COMPOUND_SCAFFOLD_TAXONOMY.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    desc_cols = [
        "compound_internal_id", "molecular_weight", "logp", "tpsa", "hbd", "hba",
        "rotatable_bonds", "ring_count", "aromatic_ring_count", "fraction_csp3",
        "formal_charge", "heavy_atom_count",
    ]
    eligible[desc_cols].to_csv(
        OUT / "02_classification_tables/COMPOUND_PHYSICOCHEMICAL_DESCRIPTORS.tsv",
        sep="\t", index=False,
    )

    # Target-role linkage from unique formal pairs; unresolved stays explicit.
    roles = pd.read_csv(ROLE, sep="\t", usecols=["canonical_uniprot_accession", "formal_primary_membrane_role"])
    pair_role = pairs.merge(
        roles, left_on="target_uniprot_id", right_on="canonical_uniprot_accession",
        how="left", validate="many_to_one",
    )
    if pair_role["formal_primary_membrane_role"].isna().any():
        raise RuntimeError("Pair targets missing FORMAL membrane roles")
    pair_role = pair_role[pair_role["compound_internal_id"].isin(set(eligible["compound_internal_id"]))].copy()
    role_counts = (
        pair_role.groupby(["compound_internal_id", "formal_primary_membrane_role"]).size()
        .rename("protein_target_count_in_role").reset_index()
    )
    target_summary = pair_role.groupby("compound_internal_id").agg(
        unique_protein_target_count=("target_uniprot_id", "nunique"),
        unique_primary_membrane_role_count=("formal_primary_membrane_role", "nunique"),
    ).reset_index()
    dominant = role_counts.sort_values(
        ["compound_internal_id", "protein_target_count_in_role", "formal_primary_membrane_role"],
        ascending=[True, False, True],
    ).groupby("compound_internal_id", as_index=False).first()[
        ["compound_internal_id", "formal_primary_membrane_role"]
    ].rename(columns={"formal_primary_membrane_role": "dominant_membrane_role_descriptive"})
    target_summary = target_summary.merge(dominant, on="compound_internal_id", how="left")
    target_summary.to_csv(
        OUT / "06_target_breadth/COMPOUND_TARGET_BREADTH_SUMMARY.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    role_counts.to_csv(OUT / "04_heatmaps/COMPOUND_ROLE_COUNTS.tsv.gz", sep="\t", index=False, compression="gzip")

    analysis = eligible[[
        "compound_internal_id", "preferred_name", "canonical_smiles_recalculated",
        "chemical_regime_local", "authoritative_chebi_direct_class",
        "exact_scaffold_id", "exact_murcko_smiles_recalculated",
        "generic_scaffold_id", "generic_murcko_smiles_recalculated",
        *[c for c in desc_cols if c != "compound_internal_id"],
        "is_approved_drug", "is_clinical_candidate", "is_endogenous_ligand",
        "is_natural_product", "is_chemical_probe",
    ]].merge(target_summary, on="compound_internal_id", how="left", validate="one_to_one")
    analysis.to_csv(OUT / "02_classification_tables/INTERACTION_LINKED_COMPOUND_ANALYSIS.tsv.gz", sep="\t", index=False, compression="gzip")

    # Status is explicitly multi-label.
    status_cols = ["is_approved_drug", "is_clinical_candidate", "is_endogenous_ligand", "is_natural_product", "is_chemical_probe"]
    status_summary = []
    for col in status_cols:
        status_summary.append({"status_label": col, "compound_count": int(analysis[col].fillna(0).astype(int).sum()), "denominator": len(analysis)})
    any_status = analysis[status_cols].fillna(0).astype(int).sum(axis=1) > 0
    status_summary.append({"status_label": "no_captured_status", "compound_count": int((~any_status).sum()), "denominator": len(analysis)})
    pd.DataFrame(status_summary).to_csv(OUT / "02_classification_tables/COMPOUND_STATUS_SUMMARY.tsv", sep="\t", index=False)

    # Audit counts.
    duplicate_smiles_groups = int(eligible.groupby("canonical_smiles_recalculated").size().gt(1).sum())
    duplicate_smiles_rows = int(eligible.duplicated("canonical_smiles_recalculated", keep=False).sum())
    duplicate_ik_groups = int(eligible.dropna(subset=["standard_inchikey"]).groupby("standard_inchikey").size().gt(1).sum())
    duplicate_ik_rows = int(eligible.duplicated("standard_inchikey", keep=False).sum())
    scaffold_compare_exact = (eligible["exact_murcko_smiles_recalculated"].fillna("") == eligible["exact_murcko_scaffold_smiles"].fillna(""))
    scaffold_compare_generic = (eligible["generic_murcko_smiles_recalculated"].fillna("") == eligible["generic_murcko_scaffold_smiles"].fillna(""))
    audit = {
        "compound_registry": len(cm),
        "interaction_linked_compounds": len(linked_ids),
        "formal_pairs": len(pairs),
        "structure_valid_interaction_linked_compounds": len(eligible),
        "excluded_or_quarantined_interaction_linked_compounds": int((~linked["eligible_structure_analysis"]).sum()),
        "missing_smiles_interaction_linked": int(linked["analysis_smiles"].isna().sum()),
        "missing_inchikey_interaction_linked": int(linked["standard_inchikey"].isna().sum()),
        "duplicate_canonical_smiles_groups": duplicate_smiles_groups,
        "duplicate_canonical_smiles_rows": duplicate_smiles_rows,
        "duplicate_inchikey_groups": duplicate_ik_groups,
        "duplicate_inchikey_rows": duplicate_ik_rows,
        "chebi_direct_class_available": int((eligible["authoritative_chebi_direct_class"] != "NOT_AVAILABLE").sum()),
        "chebi_hierarchy_complete": False,
        "exact_scaffold_match_to_v721": int(scaffold_compare_exact.sum()),
        "exact_scaffold_mismatch_to_v721": int((~scaffold_compare_exact).sum()),
        "generic_scaffold_match_to_v721": int(scaffold_compare_generic.sum()),
        "generic_scaffold_mismatch_to_v721": int((~scaffold_compare_generic).sum()),
        "rdkit_version_current": rdBase.rdkitVersion,
        "rdkit_version_v721": str(eligible["rdkit_version"].dropna().iloc[0]),
        "descriptor_discrepancies": discrepancy,
        "downstream_denominators": {
            "V1_taxonomy": len(eligible), "V2_scaffold": len(eligible), "V3_gallery": len(eligible),
            "V4_class_role": len(pair_role), "V5_scaffold_role": len(pair_role), "V6_physchem": len(eligible),
            "V7_PCA": "deterministic sample from eligible; exact n written by visualization script",
            "V8_UMAP": "deterministic sample from eligible; exact n written by visualization script",
            "V9_target_breadth": len(eligible), "V10_3D": "deterministic sample from eligible",
        },
    }
    (OUT / "01_data_audit/COMPOUND_CLASSIFICATION_DATA_AUDIT.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    regime_counts = analysis["chemical_regime_local"].value_counts()
    md = f"""# Compound classification data audit

## Frozen universes

- Compound registry: **{len(cm):,}** canonical compounds.
- Interaction-linked compounds: **{len(linked_ids):,}** canonical compounds.
- Formal protein-compound pairs: **{len(pairs):,}** unique pairs.
- Structure-valid interaction-linked compounds: **{len(eligible):,}**.
- Excluded/quarantined from structural analysis: **{int((~linked['eligible_structure_analysis']).sum()):,}**; registry identity retained.

## Identity and structure checks

- Missing SMILES in interaction-linked set: **{int(linked['analysis_smiles'].isna().sum()):,}**.
- Missing InChIKey in interaction-linked set: **{int(linked['standard_inchikey'].isna().sum()):,}**.
- Duplicate recalculated canonical-SMILES groups: **{duplicate_smiles_groups:,}** ({duplicate_smiles_rows:,} rows).
- Duplicate InChIKey groups: **{duplicate_ik_groups:,}** ({duplicate_ik_rows:,} rows).
- ChEBI direct-class annotation available: **{int((eligible['authoritative_chebi_direct_class'] != 'NOT_AVAILABLE').sum()):,} / {len(eligible):,}**.

Duplicate structures are reported rather than silently merged because canonical registry identity and source provenance may differ. The quarantined ferrocenyl name--structure conflict is excluded from every structure-derived analysis.

## Taxonomy boundary

The release does not contain a complete ClassyFire/ChEBI hierarchy. Existing ChEBI direct classes are retained as authoritative direct annotations, but superclass, subclass and direct-parent levels remain `NOT_AVAILABLE_NEEDS_EXTERNAL_MAPPING`. The broad `chemical_regime_local` field is an explicitly local, deterministic rule layer and is not presented as ClassyFire.

## Local chemical-regime counts

{regime_counts.to_markdown()}

## RDKit recalculation

All eligible interaction-linked compounds were reparsed with RDKit **{rdBase.rdkitVersion}**. Murcko scaffolds and descriptors were recalculated; stored values were compared and never silently overwritten. Full discrepancy counts are in `COMPOUND_CLASSIFICATION_DATA_AUDIT.json`.

## Downstream denominators

V1/V2/V3/V6/V9 use {len(eligible):,} eligible interaction-linked compounds. V4/V5 use {len(pair_role):,} eligible compound--protein pairs after exact pair deduplication. PCA, UMAP and 3D plots use deterministic samples whose exact sizes and selection rules are written by the visualization script.
"""
    write_md(OUT / "01_data_audit/COMPOUND_CLASSIFICATION_DATA_AUDIT.md", md)
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
