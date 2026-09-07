from __future__ import annotations

import csv
import re
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")

SOURCE_MAP = {
    "self-screening": "curated_target_relation",
    "non-screening": "bioassay_active_relation",
    "structure_context": "structure_ligand_context",
}

HIGH_MISSING_FIELDS = [
    "ghs_hazard_classes",
    "ghs_signal_words",
    "toxicity_summary",
    "livertox",
    "drug_classes",
    "pharmacodynamics",
    "drug_indication",
]

NOT_FOUND = "not_found_in_current_sources"

ION_OR_METAL_TERMS = {
    "sodium",
    "potassium",
    "calcium",
    "magnesium",
    "zinc",
    "iron",
    "copper",
    "manganese",
    "chloride",
    "bromide",
    "iodide",
    "fluoride",
    "lithium",
    "cesium",
    "cadmium",
    "nickel",
    "cobalt",
    "mercury",
    "lead",
    "aluminum",
    "aluminium",
    "chromium",
    "selenium",
}

BUFFER_SALT_SOLVENT_TERMS = {
    "water",
    "glycerol",
    "ethylene glycol",
    "propylene glycol",
    "dimethyl sulfoxide",
    "dmso",
    "ethanol",
    "methanol",
    "isopropanol",
    "acetate",
    "formate",
    "sulfate",
    "sulphate",
    "phosphate",
    "nitrate",
    "carbonate",
    "bicarbonate",
    "citrate",
    "tartrate",
    "tris",
    "hepes",
    "mes",
    "mops",
    "pipes",
    "peg",
    "polyethylene glycol",
    "imidazole",
    "azide",
}

SIMPLE_BUFFER_SALT_SOLVENT_NAMES = {
    "water",
    "glycerol",
    "ethylene glycol",
    "propylene glycol",
    "dimethyl sulfoxide",
    "dmso",
    "ethanol",
    "methanol",
    "isopropanol",
    "acetate",
    "formate",
    "sulfate",
    "sulphate",
    "phosphate",
    "nitrate",
    "carbonate",
    "bicarbonate",
    "citrate",
    "tartrate",
    "tris",
    "hepes",
    "mes",
    "mops",
    "pipes",
    "polyethylene glycol",
    "imidazole",
    "azide",
}

COFACTOR_TERMS = {
    "atp",
    "adp",
    "amp",
    "gtp",
    "gdp",
    "nad",
    "nadh",
    "nadp",
    "nadph",
    "fad",
    "fmn",
    "heme",
    "haem",
    "retinal",
    "retinol",
    "coa",
    "coenzyme",
    "biotin",
    "thiamine",
    "pyridoxal",
    "flavin",
    "folate",
    "glutathione",
    "porphyrin",
}

LIPID_STEROL_TERMS = {
    "cholesterol",
    "cholesteryl",
    "ergosterol",
    "sterol",
    "lipid",
    "phospholipid",
    "phosphatidyl",
    "sphingosine",
    "ceramide",
    "fatty acid",
    "palmitate",
    "palmitic",
    "oleate",
    "oleic",
    "stearate",
    "stearic",
    "arachidonic",
    "linoleic",
    "myristate",
    "laurate",
}


def map_compound_source(value: str) -> str:
    if not value:
        return value
    parts = [SOURCE_MAP.get(part.strip(), part.strip()) for part in value.split(";") if part.strip()]
    return ";".join(dict.fromkeys(parts))


def compound_evidence_status(clinical_statuses: str) -> str:
    statuses = {part.strip() for part in clinical_statuses.split(";") if part.strip()}
    if "approved_or_listed" in statuses:
        return "approved_drug"
    if "clinical_trial" in statuses:
        return "clinical_trial_compound"
    return "binding_evidence_only"


def compound_name_category(old_type: str) -> str:
    if old_type in {"approved_drug", "literature_drug", "cofactor_or_functional_ligand"}:
        return "common_or_drug_name"
    if old_type == "dev_code":
        return "research_code"
    if old_type == "chemical":
        return "systematic_or_chemical_name"
    if old_type in {"named_compound", "chembl_mapped_small_molecule", "organic_or_named_small_molecule"}:
        return "pubchem_synonym_name"
    if old_type in {"pdb_ligand_code_unknown", "pdbbind_affinity_ligand"}:
        return "pdb_ligand_code_or_name"
    return "unknown_name_type"


def contains_any_substring(text: str, terms: set[str]) -> bool:
    lower = text.lower()
    return any(term in lower for term in terms)


def contains_any_token(text: str, terms: set[str]) -> bool:
    tokens = set(re.findall(r"[a-z0-9]+", text.lower()))
    return any(term in tokens for term in terms)


def is_simple_buffer_salt_solvent(drug_name: str) -> bool:
    normalized = re.sub(r"[^a-z0-9]+", " ", drug_name.lower()).strip()
    if normalized in SIMPLE_BUFFER_SALT_SOLVENT_NAMES:
        return True
    return normalized.startswith("polyethylene glycol")


def compound_biological_role(row: dict[str, str], old_type: str, evidence_status: str) -> str:
    name = " ".join(
        [
            row.get("drug_name", ""),
            row.get("iupac_name", ""),
            row.get("molecular_formula", ""),
        ]
    )
    compound_id_type = row.get("compound_id_type", "")
    is_structure = "structure_ligand_context" in row.get("compound_source", "")

    if compound_id_type in {"ChEMBL", "DrugCentral"}:
        return "biologic_or_large_molecule"
    if evidence_status in {"approved_drug", "clinical_trial_compound"}:
        return "therapeutic_or_clinical_compound"
    if old_type == "pdbbind_affinity_ligand":
        return "structure_affinity_ligand"

    if is_structure or old_type in {"pdb_ligand_code_unknown", "cofactor_or_functional_ligand"}:
        if contains_any_token(name, ION_OR_METAL_TERMS):
            return "ion_or_metal"
        if old_type == "cofactor_or_functional_ligand" or contains_any_token(name, COFACTOR_TERMS):
            return "endogenous_ligand_or_cofactor"
        if contains_any_substring(name, LIPID_STEROL_TERMS):
            return "membrane_lipid_or_sterol"
        if is_simple_buffer_salt_solvent(row.get("drug_name", "")):
            return "buffer_salt_solvent"
        if compound_id_type == "unmapped_structure_ligand" or old_type == "pdb_ligand_code_unknown":
            return "unknown_structure_ligand"

    return "bioactive_research_ligand"


def rewrite_binary() -> None:
    path = BASE / "drug_protein_binary_relationships.tsv"
    temp = path.with_suffix(".tmp")
    with path.open("r", encoding="utf-8", newline="") as src, temp.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            row["compound_source"] = map_compound_source(row.get("compound_source", ""))
            writer.writerow(row)
    temp.replace(path)


def rewrite_small_molecules() -> None:
    path = BASE / "small_molecule_database.tsv"
    temp = path.with_suffix(".tmp")
    with path.open("r", encoding="utf-8", newline="") as src, temp.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        old_fields = reader.fieldnames or []
        fields: list[str] = []
        for field in old_fields:
            if field == "drug_name_type":
                fields.extend(["compound_name_category", "compound_biological_role"])
            elif field == "clinical_statuses":
                fields.append("compound_evidence_status")
            else:
                fields.append(field)

        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            old_type = row.get("drug_name_type", "")
            evidence_status = compound_evidence_status(row.get("clinical_statuses", ""))
            row["compound_source"] = map_compound_source(row.get("compound_source", ""))

            out: dict[str, str] = {}
            for field in old_fields:
                if field == "drug_name_type":
                    out["compound_name_category"] = compound_name_category(old_type)
                    out["compound_biological_role"] = compound_biological_role(row, old_type, evidence_status)
                elif field == "clinical_statuses":
                    out["compound_evidence_status"] = evidence_status
                else:
                    value = row.get(field, "")
                    if field in HIGH_MISSING_FIELDS and value == "":
                        value = NOT_FOUND
                    out[field] = value
            writer.writerow(out)
    temp.replace(path)


def main() -> None:
    rewrite_binary()
    rewrite_small_molecules()


if __name__ == "__main__":
    main()
