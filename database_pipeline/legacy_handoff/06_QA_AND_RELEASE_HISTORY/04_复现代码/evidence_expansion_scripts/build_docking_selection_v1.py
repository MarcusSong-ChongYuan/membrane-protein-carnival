#!/usr/bin/env python3
"""Build a reproducible docking-priority preview from the frozen MemPro V6.0 release.

This script does not alter release files. It assigns every protein-compound pair
to one mutually exclusive tier and emits both the full audit table and the
recommended R0-D3 subset.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


RELEASE_DIR = Path(
    r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_0_20260727"
)
RUN_DIR = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\docking_selection_v1_20260728"
)
STAGING_DIR = RUN_DIR / "staging"
QA_DIR = RUN_DIR / "qa"

PAIR_FILE = RELEASE_DIR / "protein_compound_summary_v6_0.tsv.gz"
SITE_FILE = RELEASE_DIR / "binding_site_instances_v6_0.tsv.gz"
PROTEIN_FILE = RELEASE_DIR / "human_membrane_protein_master_v6_0.tsv"
COMPOUND_FILE = RELEASE_DIR / "small_molecule_master_v1_1.tsv"

ALL_OUTPUT = STAGING_DIR / "docking_pair_priority_v1_v60_preview.tsv.gz"
RECOMMENDED_OUTPUT = STAGING_DIR / "docking_recommended_pairs_v1_v60_preview.tsv.gz"
TARGET_OUTPUT = STAGING_DIR / "docking_target_summary_v1_v60_preview.tsv"
VALIDATION_OUTPUT = QA_DIR / "DOCKING_SELECTION_V1_VALIDATION.json"
DICTIONARY_OUTPUT = RUN_DIR / "DOCKING_SELECTION_V1_DATA_DICTIONARY.md"

EXPECTED_TIER_COUNTS = {
    "R0": 5062,
    "D1": 19018,
    "D2": 15347,
    "D3": 9727,
    "X": 526884,
}
RELEASE_BASIS = "MemPro V6.0; preview pending V6.1 refresh"
SELECTION_BASIS_TAG = "MemPro_V6.0_preview_before_V6.1_refresh"
DICTIONARY_TITLE = "Docking selection V1 (V6.0 preview)"


def text(value: str | None) -> str:
    return (value or "").strip()


def integer(value: str | None) -> int:
    try:
        return int(float(text(value)))
    except (TypeError, ValueError):
        return 0


def number(value: str | None) -> float | None:
    try:
        value = text(value)
        return float(value) if value else None
    except (TypeError, ValueError):
        return None


def true_value(value: str | None) -> bool:
    return text(value).lower() in {"1", "true", "yes", "y"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pair_key(target_id: str, compound_id: str) -> tuple[str, str]:
    return text(target_id), text(compound_id)


def collect_pair_compounds() -> set[str]:
    compound_ids: set[str] = set()
    with gzip.open(PAIR_FILE, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            compound_ids.add(text(row["compound_internal_id"]))
    return compound_ids


def load_pair_specific_sites() -> dict[tuple[str, str], set[str]]:
    sites: dict[tuple[str, str], set[str]] = defaultdict(set)
    with gzip.open(SITE_FILE, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pdb_ids = text(row["pdb_ids"])
            if not pdb_ids:
                continue
            key = pair_key(row["target_uniprot_id"], row["compound_internal_id"])
            if not all(key):
                continue
            for pdb_id in pdb_ids.replace(",", ";").split(";"):
                if text(pdb_id):
                    sites[key].add(text(pdb_id).upper())
    return sites


def load_proteins() -> dict[str, dict[str, str | bool]]:
    proteins: dict[str, dict[str, str | bool]] = {}
    with PROTEIN_FILE.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target_id = text(row["target_uniprot_id"])
            has_membrane_structure = any(
                (
                    true_value(row["opm_present_v5"]),
                    true_value(row["pdbtm_present_v5"]),
                    bool(text(row["opm_pdb_ids_v5"])),
                    bool(text(row["pdbtm_pdb_ids_v5"])),
                )
            )
            has_pdb = bool(text(row["pdb_ids"]))
            has_af = bool(text(row["alphafolddb_ids"]))
            if has_membrane_structure:
                base_structure_grade = "S2_membrane_experimental_OPM_PDBTM"
            elif has_pdb:
                base_structure_grade = "S3_other_experimental_PDB"
            elif has_af:
                base_structure_grade = "S4_AlphaFold_only"
            else:
                base_structure_grade = "S0_no_structure"
            evidence = text(row["evidence_level_v52"])
            membrane_class = text(row["membrane_class_v52"])
            website_default = true_value(row["website_default_v52"])
            eligible = (
                evidence in {"E1", "E2"}
                and membrane_class in {"A", "B", "C"}
                and website_default
            )
            exclusion_reasons = []
            if evidence not in {"E1", "E2"}:
                exclusion_reasons.append(f"protein_evidence_{evidence or 'missing'}")
            if membrane_class not in {"A", "B", "C"}:
                exclusion_reasons.append(
                    f"membrane_class_{membrane_class or 'missing'}"
                )
            if not website_default:
                exclusion_reasons.append("not_website_default")
            proteins[target_id] = {
                "approved_symbol": text(row["approved_symbol"]),
                "protein_name": text(row["protein_name"]),
                "protein_evidence_level": evidence,
                "membrane_class": membrane_class,
                "functional_primary_class": text(row["functional_primary_class_v5"]),
                "transmembrane_count": text(row["transmembrane_count_v5"]),
                "sequence_length": text(row["sequence_length_v53"]),
                "has_membrane_experimental_structure": has_membrane_structure,
                "has_any_experimental_pdb": has_pdb,
                "has_alphafold_model": has_af,
                "base_structure_grade": base_structure_grade,
                "opm_pdb_ids": text(row["opm_pdb_ids_v5"]),
                "pdbtm_pdb_ids": text(row["pdbtm_pdb_ids_v5"]),
                "all_pdb_ids": text(row["pdb_ids"]),
                "alphafold_ids": text(row["alphafolddb_ids"]),
                "protein_eligible": eligible,
                "protein_exclusion_reasons": ";".join(exclusion_reasons),
            }
    return proteins


def load_compounds(pair_compound_ids: set[str]) -> dict[str, dict[str, str | bool | float]]:
    compounds: dict[str, dict[str, str | bool | float]] = {}
    with COMPOUND_FILE.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            compound_id = text(row["compound_internal_id"])
            if compound_id not in pair_compound_ids:
                continue
            mw = number(row["molecular_weight"])
            heavy_atoms = number(row["heavy_atom_count"])
            charge = number(row["formal_charge"])
            rotatable = number(row["rotatable_bond_count"])
            structural_class = text(row["computed_structural_class"])
            scope_status = text(row["compound_scope_status"])
            qc_status = text(row["record_qc_status"])
            exclusions = []
            if not text(row["standard_smiles"]):
                exclusions.append("missing_standard_smiles")
            if not text(row["standard_inchikey"]):
                exclusions.append("missing_standard_inchikey")
            if scope_status != "core":
                exclusions.append(f"compound_scope_{scope_status or 'missing'}")
            if qc_status != "ok":
                exclusions.append(f"compound_qc_{qc_status or 'missing'}")
            if mw is None or not (100 <= mw <= 700):
                exclusions.append("molecular_weight_outside_100_700")
            if heavy_atoms is None or not (6 <= heavy_atoms <= 60):
                exclusions.append("heavy_atoms_outside_6_60")
            if charge is None or abs(charge) > 2:
                exclusions.append("absolute_formal_charge_above_2")
            if rotatable is None or rotatable > 20:
                exclusions.append("rotatable_bonds_above_20")
            structural_lower = structural_class.lower()
            if "inorganic" in structural_lower or "metal" in structural_lower:
                exclusions.append(f"excluded_structural_class_{structural_class}")
            eligible = not exclusions
            approved = true_value(row["is_approved_drug"])
            clinical = true_value(row["is_clinical_candidate"])
            endogenous = true_value(row["is_endogenous_ligand"])
            natural_product = true_value(row["is_natural_product"])
            probe = true_value(row["is_chemical_probe"])
            compounds[compound_id] = {
                "compound_preferred_name": text(row["preferred_name"]),
                "compound_scope_status": scope_status,
                "identity_confidence": text(row["identity_confidence"]),
                "standard_smiles": text(row["standard_smiles"]),
                "standard_inchikey": text(row["standard_inchikey"]),
                "molecular_formula": text(row["molecular_formula"]),
                "molecular_weight": "" if mw is None else mw,
                "xlogp": text(row["xlogp"]),
                "tpsa": text(row["tpsa"]),
                "hbond_donor_count": text(row["hbond_donor_count"]),
                "hbond_acceptor_count": text(row["hbond_acceptor_count"]),
                "rotatable_bond_count": "" if rotatable is None else rotatable,
                "formal_charge": "" if charge is None else charge,
                "heavy_atom_count": "" if heavy_atoms is None else heavy_atoms,
                "ring_count": text(row["ring_count"]),
                "computed_structural_class": structural_class,
                "record_qc_status": qc_status,
                "is_approved_drug": approved,
                "is_clinical_candidate": clinical,
                "is_endogenous_ligand": endogenous,
                "is_natural_product": natural_product,
                "is_chemical_probe": probe,
                "special_compound": approved or clinical or probe or endogenous,
                "chemical_eligible": eligible,
                "chemical_exclusion_reasons": ";".join(exclusions),
            }
    return compounds


def assign_tier(
    pair: dict[str, str],
    protein: dict[str, str | bool],
    compound: dict[str, str | bool | float],
    site_pdb_ids: set[str],
) -> tuple[str, str, str, int]:
    protein_ok = bool(protein["protein_eligible"])
    chemical_ok = bool(compound["chemical_eligible"])
    has_pair_site = bool(site_pdb_ids)
    best_evidence = text(pair["best_binding_evidence_level"])
    be1 = integer(pair["BE1_evidence_count"])
    be2 = integer(pair["BE2_evidence_count"])
    be3 = integer(pair["BE3_evidence_count"])
    sources = integer(pair["independent_source_count"])
    potency = number(pair["best_standard_value_nM"])
    special = bool(compound["special_compound"])
    approved_clinical_probe = any(
        (
            bool(compound["is_approved_drug"]),
            bool(compound["is_clinical_candidate"]),
            bool(compound["is_chemical_probe"]),
        )
    )
    endogenous = bool(compound["is_endogenous_ligand"])
    membrane_structure = bool(protein["has_membrane_experimental_structure"])
    any_pdb = bool(protein["has_any_experimental_pdb"])

    if not (protein_ok and chemical_ok):
        reasons = []
        if not protein_ok:
            reasons.append("protein_not_eligible")
        if not chemical_ok:
            reasons.append("compound_not_eligible")
        return "X", "exclude_default", ";".join(reasons), 0

    if has_pair_site:
        return (
            "R0",
            "redocking_protocol_validation",
            "pair_specific_experimental_PDB_binding_site",
            100,
        )

    if (
        membrane_structure
        and be2 > 0
        and potency is not None
        and potency <= 1000
        and (sources >= 2 or special)
    ):
        return (
            "D1",
            "primary_prospective_docking",
            "membrane_experimental_structure;direct_binding_le_1uM;"
            + ("multi_source" if sources >= 2 else "priority_compound"),
            90,
        )

    d2_any_pdb = (
        any_pdb
        and be2 > 0
        and potency is not None
        and potency <= 100
        and (sources >= 2 or special)
    )
    d2_membrane = (
        membrane_structure
        and be2 > 0
        and potency is not None
        and potency <= 10000
        and (sources >= 2 or special)
    )
    if d2_any_pdb or d2_membrane:
        reason = (
            "experimental_PDB;direct_binding_le_100nM"
            if d2_any_pdb
            else "membrane_experimental_structure;direct_binding_le_10uM"
        )
        reason += ";" + ("multi_source" if sources >= 2 else "priority_compound")
        return "D2", "secondary_prospective_docking", reason, 70

    d3_drug_like = (
        any_pdb
        and approved_clinical_probe
        and (
            (be2 > 0 and potency is not None and potency <= 10000)
            or sources >= 2
        )
    )
    d3_endogenous = (
        membrane_structure
        and endogenous
        and (
            best_evidence in {"BE1", "BE2"}
            or (best_evidence == "BE3" and sources >= 2)
            or be1 > 0
            or be2 > 0
            or (be3 > 0 and sources >= 2)
        )
    )
    if d3_drug_like or d3_endogenous:
        reason = (
            "experimental_PDB;approved_clinical_or_probe;potency_or_multisource"
            if d3_drug_like
            else "membrane_experimental_structure;endogenous_ligand;supported_relation"
        )
        return "D3", "exploratory_mechanism_docking", reason, 50

    reasons = []
    if not any_pdb:
        reasons.append("no_experimental_PDB")
    elif not membrane_structure:
        reasons.append("no_membrane_curated_structure")
    if best_evidence == "BE3" or (be1 == 0 and be2 == 0):
        reasons.append("pharmacology_only_or_no_direct_binding")
    if potency is None:
        reasons.append("missing_quantitative_potency")
    elif potency > 10000:
        reasons.append("potency_above_10uM")
    if sources < 2:
        reasons.append("single_source")
    if not special:
        reasons.append("not_priority_compound")
    if not reasons:
        reasons.append("does_not_meet_D1_D3_combination_rules")
    return "X", "exclude_default", ";".join(reasons), 0


def build_outputs() -> dict:
    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)

    pair_compound_ids = collect_pair_compounds()
    sites = load_pair_specific_sites()
    proteins = load_proteins()
    compounds = load_compounds(pair_compound_ids)

    output_fields = [
        "target_uniprot_id",
        "approved_symbol",
        "protein_name",
        "compound_internal_id",
        "compound_preferred_name",
        "docking_tier",
        "priority_score",
        "docking_role",
        "tier_reason",
        "structure_grade",
        "pair_site_pdb_ids",
        "protein_evidence_level",
        "membrane_class",
        "functional_primary_class",
        "transmembrane_count",
        "sequence_length",
        "has_membrane_experimental_structure",
        "has_any_experimental_pdb",
        "has_alphafold_model",
        "opm_pdb_ids",
        "pdbtm_pdb_ids",
        "all_pdb_ids",
        "alphafold_ids",
        "best_binding_evidence_level",
        "BE1_evidence_count",
        "BE2_evidence_count",
        "BE3_evidence_count",
        "binding_evidence_count",
        "independent_source_count",
        "independent_sources",
        "best_standard_value_nM",
        "compound_scope_status",
        "identity_confidence",
        "standard_smiles",
        "standard_inchikey",
        "molecular_formula",
        "molecular_weight",
        "xlogp",
        "tpsa",
        "hbond_donor_count",
        "hbond_acceptor_count",
        "rotatable_bond_count",
        "formal_charge",
        "heavy_atom_count",
        "ring_count",
        "computed_structural_class",
        "record_qc_status",
        "is_approved_drug",
        "is_clinical_candidate",
        "is_endogenous_ligand",
        "is_natural_product",
        "is_chemical_probe",
        "protein_eligible",
        "chemical_eligible",
        "protein_exclusion_reasons",
        "chemical_exclusion_reasons",
        "selection_policy_version",
        "selection_release_basis",
    ]

    tier_counts: Counter[str] = Counter()
    tier_targets: dict[str, set[str]] = defaultdict(set)
    tier_compounds: dict[str, set[str]] = defaultdict(set)
    target_tiers: dict[str, Counter[str]] = defaultdict(Counter)
    unique_keys: set[tuple[str, str]] = set()
    duplicate_keys = 0
    missing_proteins = 0
    missing_compounds = 0

    with (
        gzip.open(ALL_OUTPUT, "wt", encoding="utf-8", newline="") as all_handle,
        gzip.open(RECOMMENDED_OUTPUT, "wt", encoding="utf-8", newline="") as rec_handle,
        gzip.open(PAIR_FILE, "rt", encoding="utf-8-sig", newline="") as pair_handle,
    ):
        all_writer = csv.DictWriter(
            all_handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        rec_writer = csv.DictWriter(
            rec_handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        all_writer.writeheader()
        rec_writer.writeheader()

        for pair in csv.DictReader(pair_handle, delimiter="\t"):
            target_id = text(pair["target_uniprot_id"])
            compound_id = text(pair["compound_internal_id"])
            key = pair_key(target_id, compound_id)
            if key in unique_keys:
                duplicate_keys += 1
            unique_keys.add(key)

            protein = proteins.get(target_id)
            compound = compounds.get(compound_id)
            if protein is None:
                missing_proteins += 1
                protein = {
                    "approved_symbol": text(pair["approved_symbol"]),
                    "protein_name": "",
                    "protein_evidence_level": "",
                    "membrane_class": "",
                    "functional_primary_class": "",
                    "transmembrane_count": "",
                    "sequence_length": "",
                    "has_membrane_experimental_structure": False,
                    "has_any_experimental_pdb": False,
                    "has_alphafold_model": False,
                    "base_structure_grade": "S0_no_structure",
                    "opm_pdb_ids": "",
                    "pdbtm_pdb_ids": "",
                    "all_pdb_ids": "",
                    "alphafold_ids": "",
                    "protein_eligible": False,
                    "protein_exclusion_reasons": "missing_protein_master_record",
                }
            if compound is None:
                missing_compounds += 1
                compound = {
                    "compound_preferred_name": text(pair["preferred_name"]),
                    "compound_scope_status": "",
                    "identity_confidence": "",
                    "standard_smiles": "",
                    "standard_inchikey": "",
                    "molecular_formula": "",
                    "molecular_weight": "",
                    "xlogp": "",
                    "tpsa": "",
                    "hbond_donor_count": "",
                    "hbond_acceptor_count": "",
                    "rotatable_bond_count": "",
                    "formal_charge": "",
                    "heavy_atom_count": "",
                    "ring_count": "",
                    "computed_structural_class": "",
                    "record_qc_status": "",
                    "is_approved_drug": False,
                    "is_clinical_candidate": False,
                    "is_endogenous_ligand": False,
                    "is_natural_product": False,
                    "is_chemical_probe": False,
                    "special_compound": False,
                    "chemical_eligible": False,
                    "chemical_exclusion_reasons": "missing_compound_master_record",
                }

            site_pdb_ids = sites.get(key, set())
            tier, role, reason, priority_score = assign_tier(
                pair, protein, compound, site_pdb_ids
            )
            structure_grade = (
                "S1_pair_specific_experimental_site"
                if site_pdb_ids
                else text(str(protein["base_structure_grade"]))
            )

            row = {
                "target_uniprot_id": target_id,
                "approved_symbol": text(str(protein["approved_symbol"])),
                "protein_name": text(str(protein["protein_name"])),
                "compound_internal_id": compound_id,
                "compound_preferred_name": text(
                    str(compound["compound_preferred_name"])
                ),
                "docking_tier": tier,
                "priority_score": priority_score,
                "docking_role": role,
                "tier_reason": reason,
                "structure_grade": structure_grade,
                "pair_site_pdb_ids": ";".join(sorted(site_pdb_ids)),
                "protein_evidence_level": protein["protein_evidence_level"],
                "membrane_class": protein["membrane_class"],
                "functional_primary_class": protein["functional_primary_class"],
                "transmembrane_count": protein["transmembrane_count"],
                "sequence_length": protein["sequence_length"],
                "has_membrane_experimental_structure": int(
                    bool(protein["has_membrane_experimental_structure"])
                ),
                "has_any_experimental_pdb": int(
                    bool(protein["has_any_experimental_pdb"])
                ),
                "has_alphafold_model": int(bool(protein["has_alphafold_model"])),
                "opm_pdb_ids": protein["opm_pdb_ids"],
                "pdbtm_pdb_ids": protein["pdbtm_pdb_ids"],
                "all_pdb_ids": protein["all_pdb_ids"],
                "alphafold_ids": protein["alphafold_ids"],
                "best_binding_evidence_level": pair["best_binding_evidence_level"],
                "BE1_evidence_count": pair["BE1_evidence_count"],
                "BE2_evidence_count": pair["BE2_evidence_count"],
                "BE3_evidence_count": pair["BE3_evidence_count"],
                "binding_evidence_count": pair["binding_evidence_count"],
                "independent_source_count": pair["independent_source_count"],
                "independent_sources": pair["independent_sources"],
                "best_standard_value_nM": pair["best_standard_value_nM"],
                "compound_scope_status": compound["compound_scope_status"],
                "identity_confidence": compound["identity_confidence"],
                "standard_smiles": compound["standard_smiles"],
                "standard_inchikey": compound["standard_inchikey"],
                "molecular_formula": compound["molecular_formula"],
                "molecular_weight": compound["molecular_weight"],
                "xlogp": compound["xlogp"],
                "tpsa": compound["tpsa"],
                "hbond_donor_count": compound["hbond_donor_count"],
                "hbond_acceptor_count": compound["hbond_acceptor_count"],
                "rotatable_bond_count": compound["rotatable_bond_count"],
                "formal_charge": compound["formal_charge"],
                "heavy_atom_count": compound["heavy_atom_count"],
                "ring_count": compound["ring_count"],
                "computed_structural_class": compound["computed_structural_class"],
                "record_qc_status": compound["record_qc_status"],
                "is_approved_drug": int(bool(compound["is_approved_drug"])),
                "is_clinical_candidate": int(
                    bool(compound["is_clinical_candidate"])
                ),
                "is_endogenous_ligand": int(bool(compound["is_endogenous_ligand"])),
                "is_natural_product": int(bool(compound["is_natural_product"])),
                "is_chemical_probe": int(bool(compound["is_chemical_probe"])),
                "protein_eligible": int(bool(protein["protein_eligible"])),
                "chemical_eligible": int(bool(compound["chemical_eligible"])),
                "protein_exclusion_reasons": protein["protein_exclusion_reasons"],
                "chemical_exclusion_reasons": compound[
                    "chemical_exclusion_reasons"
                ],
                "selection_policy_version": "docking_selection_v1",
                "selection_release_basis": SELECTION_BASIS_TAG,
            }
            all_writer.writerow(row)
            if tier != "X":
                rec_writer.writerow(row)

            tier_counts[tier] += 1
            tier_targets[tier].add(target_id)
            tier_compounds[tier].add(compound_id)
            target_tiers[target_id][tier] += 1

    target_fields = [
        "target_uniprot_id",
        "approved_symbol",
        "protein_name",
        "protein_evidence_level",
        "membrane_class",
        "functional_primary_class",
        "structure_grade",
        "R0_pair_count",
        "D1_pair_count",
        "D2_pair_count",
        "D3_pair_count",
        "recommended_pair_count",
        "X_pair_count",
        "recommended_for_docking",
    ]
    with TARGET_OUTPUT.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=target_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for target_id in sorted(target_tiers):
            protein = proteins[target_id]
            counts = target_tiers[target_id]
            recommended = sum(counts[tier] for tier in ("R0", "D1", "D2", "D3"))
            writer.writerow(
                {
                    "target_uniprot_id": target_id,
                    "approved_symbol": protein["approved_symbol"],
                    "protein_name": protein["protein_name"],
                    "protein_evidence_level": protein["protein_evidence_level"],
                    "membrane_class": protein["membrane_class"],
                    "functional_primary_class": protein["functional_primary_class"],
                    "structure_grade": protein["base_structure_grade"],
                    "R0_pair_count": counts["R0"],
                    "D1_pair_count": counts["D1"],
                    "D2_pair_count": counts["D2"],
                    "D3_pair_count": counts["D3"],
                    "recommended_pair_count": recommended,
                    "X_pair_count": counts["X"],
                    "recommended_for_docking": int(recommended > 0),
                }
            )

    observed_counts = {tier: tier_counts[tier] for tier in ("R0", "D1", "D2", "D3", "X")}
    recommended_count = sum(observed_counts[tier] for tier in ("R0", "D1", "D2", "D3"))
    prospective_count = sum(observed_counts[tier] for tier in ("D1", "D2", "D3"))
    validation = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS"
        if (
            (EXPECTED_TIER_COUNTS is None or observed_counts == EXPECTED_TIER_COUNTS)
            and duplicate_keys == 0
            and missing_proteins == 0
            and missing_compounds == 0
        )
        else "FAIL",
        "release_basis": RELEASE_BASIS,
        "selection_policy_version": "docking_selection_v1",
        "pair_count_total": sum(observed_counts.values()),
        "pair_count_unique": len(unique_keys),
        "duplicate_pair_keys": duplicate_keys,
        "recommended_pair_count_including_R0": recommended_count,
        "prospective_pair_count_D1_D3": prospective_count,
        "tier_counts": observed_counts,
        "expected_tier_counts": (
            EXPECTED_TIER_COUNTS
            if EXPECTED_TIER_COUNTS is not None
            else "not_predeclared_for_release_refresh"
        ),
        "tier_target_counts": {
            tier: len(tier_targets[tier]) for tier in ("R0", "D1", "D2", "D3", "X")
        },
        "tier_compound_counts": {
            tier: len(tier_compounds[tier]) for tier in ("R0", "D1", "D2", "D3", "X")
        },
        "pair_specific_site_pair_count_loaded": len(sites),
        "pair_compound_count_loaded": len(pair_compound_ids),
        "protein_master_count_loaded": len(proteins),
        "compound_master_pair_subset_count_loaded": len(compounds),
        "missing_protein_master_records": missing_proteins,
        "missing_compound_master_records": missing_compounds,
        "input_files": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in (PAIR_FILE, SITE_FILE, PROTEIN_FILE, COMPOUND_FILE)
        },
        "output_files": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in (ALL_OUTPUT, RECOMMENDED_OUTPUT, TARGET_OUTPUT)
        },
    }
    VALIDATION_OUTPUT.write_text(
        json.dumps(validation, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return validation


def write_dictionary(validation: dict) -> None:
    counts = validation["tier_counts"]
    text_body = f"""# {DICTIONARY_TITLE}

This is a reproducible docking-oriented partition based on {RELEASE_BASIS}.

## Mutually exclusive tiers

| Tier | Intended use | Pair count |
|---|---|---:|
| R0 | Redocking and protocol validation using a pair-specific experimental PDB site | {counts['R0']:,} |
| D1 | Primary prospective docking: membrane-curated experimental structure, direct binding <=1 uM, plus multi-source support or a priority compound | {counts['D1']:,} |
| D2 | Secondary prospective docking: very strong direct binding with any experimental PDB, or <=10 uM with a membrane-curated structure, plus support | {counts['D2']:,} |
| D3 | Exploratory/mechanistic docking for approved drugs, clinical candidates, probes, and supported endogenous ligands | {counts['D3']:,} |
| X | Retained in the database but excluded from default docking | {counts['X']:,} |

R0 is a validation set, not a prospective discovery set. D1-D3 contain
{validation['prospective_pair_count_D1_D3']:,} prospective pairs. Including all
R0 controls produces {validation['recommended_pair_count_including_R0']:,}
recommended rows.

## Eligibility gates

Protein: evidence E1/E2, membrane class A/B/C, and website-default inclusion.

Compound: canonical structure and InChIKey present; core scope; QC=ok; molecular
weight 100-700; heavy atoms 6-60; absolute formal charge <=2; rotatable bonds
<=20; not inorganic or metal-containing.

## Structure grades

- S1: pair-specific experimental binding site.
- S2: membrane-curated experimental structure from OPM/PDBTM.
- S3: another experimental PDB structure.
- S4: AlphaFold-only model.
- S0: no structure currently assigned.

## Important interpretation

Docking scores are not binding affinities. Each selected row still requires
receptor-state selection, ligand microstate preparation, binding-box definition,
and protocol validation before HPC execution.
"""
    DICTIONARY_OUTPUT.write_text(text_body, encoding="utf-8")


if __name__ == "__main__":
    report = build_outputs()
    write_dictionary(report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
