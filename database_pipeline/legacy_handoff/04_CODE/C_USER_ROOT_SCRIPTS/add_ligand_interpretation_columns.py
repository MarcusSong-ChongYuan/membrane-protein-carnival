from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")

CLASS_PRIORITY = [
    "structure_affinity_ligand",
    "approved_or_clinical_target_ligand",
    "cofactor_or_functional_ligand",
    "membrane_lipid_or_sterol",
    "curated_target_ligand",
    "structure_bound_ligand",
    "bioassay_active_ligand",
    "ion_or_metal",
    "buffer_salt_solvent",
    "unknown_structure_ligand",
]

PRIORITY_RANK = {value: idx for idx, value in enumerate(CLASS_PRIORITY)}


def split_values(value: str) -> set[str]:
    return {part.strip() for part in value.split(";") if part.strip()}


def classify_binary_row(row: dict[str, str], status_by_drug: dict[str, str]) -> str:
    role = row.get("compound_biological_role", "")
    compound_source = split_values(row.get("compound_source", ""))
    source_database = split_values(row.get("source_database", ""))
    row_status = split_values(row.get("clinical_or_approval_status", ""))
    molecule_status = status_by_drug.get(row.get("drug_id", ""), "")

    if role == "buffer_salt_solvent":
        return "buffer_salt_solvent"
    if role == "ion_or_metal":
        return "ion_or_metal"
    if role == "endogenous_ligand_or_cofactor":
        return "cofactor_or_functional_ligand"
    if role == "membrane_lipid_or_sterol":
        return "membrane_lipid_or_sterol"
    if role == "structure_affinity_ligand" or "PDBbind" in source_database:
        return "structure_affinity_ligand"
    if row_status & {"approved_or_listed", "clinical_trial"} or molecule_status in {
        "approved_drug",
        "clinical_trial_compound",
    }:
        return "approved_or_clinical_target_ligand"
    if "curated_target_relation" in compound_source:
        return "curated_target_ligand"
    if role == "unknown_structure_ligand":
        return "unknown_structure_ligand"
    if "structure_ligand_context" in compound_source:
        return "structure_bound_ligand"
    if "bioassay_active_relation" in compound_source:
        return "bioassay_active_ligand"
    return "bioassay_active_ligand"


def classify_small_molecule_fallback(row: dict[str, str]) -> str:
    role = row.get("compound_biological_role", "")
    compound_source = split_values(row.get("compound_source", ""))
    status = row.get("compound_evidence_status", "")

    if role == "buffer_salt_solvent":
        return "buffer_salt_solvent"
    if role == "ion_or_metal":
        return "ion_or_metal"
    if role == "endogenous_ligand_or_cofactor":
        return "cofactor_or_functional_ligand"
    if role == "membrane_lipid_or_sterol":
        return "membrane_lipid_or_sterol"
    if role == "structure_affinity_ligand":
        return "structure_affinity_ligand"
    if status in {"approved_drug", "clinical_trial_compound"}:
        return "approved_or_clinical_target_ligand"
    if "curated_target_relation" in compound_source:
        return "curated_target_ligand"
    if role == "unknown_structure_ligand":
        return "unknown_structure_ligand"
    if "structure_ligand_context" in compound_source:
        return "structure_bound_ligand"
    if "bioassay_active_relation" in compound_source:
        return "bioassay_active_ligand"
    return "bioassay_active_ligand"


def sorted_classes(values: set[str]) -> list[str]:
    return sorted(values, key=lambda value: PRIORITY_RANK.get(value, len(PRIORITY_RANK)))


def add_binary_column(status_by_drug: dict[str, str]) -> dict[str, set[str]]:
    path = BASE / "drug_protein_binary_relationships.tsv"
    temp = path.with_suffix(".tmp")
    classes_by_drug: dict[str, set[str]] = defaultdict(set)

    with path.open("r", encoding="utf-8", newline="") as src, temp.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        old_fields = reader.fieldnames or []
        insert_after = "compound_biological_role"
        fields: list[str] = []
        for field in old_fields:
            fields.append(field)
            if field == insert_after:
                fields.append("ligand_interpretation_class")

        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            cls = classify_binary_row(row, status_by_drug)
            classes_by_drug[row["drug_id"]].add(cls)
            out: dict[str, str] = {}
            for field in old_fields:
                out[field] = row.get(field, "")
                if field == insert_after:
                    out["ligand_interpretation_class"] = cls
            writer.writerow(out)

    temp.replace(path)
    return classes_by_drug


def add_small_molecule_columns(classes_by_drug: dict[str, set[str]]) -> None:
    path = BASE / "small_molecule_database.tsv"
    temp = path.with_suffix(".tmp")

    with path.open("r", encoding="utf-8", newline="") as src, temp.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        old_fields = reader.fieldnames or []
        insert_after = "compound_biological_role"
        fields: list[str] = []
        for field in old_fields:
            fields.append(field)
            if field == insert_after:
                fields.extend(["ligand_interpretation_classes", "best_ligand_interpretation_class"])

        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            row_classes = classes_by_drug.get(row["drug_id"], set())
            if not row_classes:
                row_classes = {classify_small_molecule_fallback(row)}
            classes = sorted_classes(row_classes)
            class_text = ";".join(classes)
            best = classes[0]
            out: dict[str, str] = {}
            for field in old_fields:
                out[field] = row.get(field, "")
                if field == insert_after:
                    out["ligand_interpretation_classes"] = class_text
                    out["best_ligand_interpretation_class"] = best
            writer.writerow(out)

    temp.replace(path)


def main() -> None:
    status_by_drug: dict[str, str] = {}
    small_path = BASE / "small_molecule_database.tsv"
    with small_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            status_by_drug[row["drug_id"]] = row.get("compound_evidence_status", "")

    classes_by_drug = add_binary_column(status_by_drug)
    add_small_molecule_columns(classes_by_drug)
    print(
        {
            "binary_ligand_interpretation_classes": len(classes_by_drug),
            "priority": CLASS_PRIORITY,
        }
    )


if __name__ == "__main__":
    main()
