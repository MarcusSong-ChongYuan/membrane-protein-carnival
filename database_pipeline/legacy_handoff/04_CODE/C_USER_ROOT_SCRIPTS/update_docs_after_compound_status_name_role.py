from __future__ import annotations

import hashlib
import json
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")


SMALL_MOLECULE_DICTIONARY_SECTION = """## `small_molecule_database.tsv`

### Identity and Role (1-6)
| # | Column | Meaning |
|---|---|---|
| 1 | `drug_id` | Primary compound identifier. Usually PubChem CID; may be ChEMBL/DrugCentral for unmapped biologics or `UNMAPPED_PDB_LIGAND:*` for structure ligands with no CID. |
| 2 | `compound_id_type` | `PubChem CID`, original type for biologics, or `unmapped_structure_ligand`. |
| 3 | `drug_name` | Preferred display name from the source table or PubChem/PDB mapping. |
| 4 | `compound_name_category` | Type of display name: `common_or_drug_name`, `research_code`, `systematic_or_chemical_name`, `pubchem_synonym_name`, `pdb_ligand_code_or_name`, or `unknown_name_type`. |
| 5 | `compound_biological_role` | Biological/use role for filtering: `therapeutic_or_clinical_compound`, `bioactive_research_ligand`, `structure_affinity_ligand`, `endogenous_ligand_or_cofactor`, `membrane_lipid_or_sterol`, `ion_or_metal`, `buffer_salt_solvent`, `biologic_or_large_molecule`, or `unknown_structure_ligand`. |
| 6 | `compound_source` | Relationship evidence origin: `curated_target_relation`, `bioassay_active_relation`, and/or `structure_ligand_context`. |

### Source and Evidence Summary (7-10)
| # | Column | Meaning |
|---|---|---|
| 7 | `source_databases` | Semicolon-separated source databases. |
| 8 | `source_row_count` | Number of Sheet 1 relationship rows for this compound. |
| 9 | `unique_target_count` | Number of unique proteins this compound targets. |
| 10 | `compound_evidence_status` | Compound-level status: `approved_drug`, `clinical_trial_compound`, or `binding_evidence_only`. |

### Activity Summary (11-13)
| # | Column | Meaning |
|---|---|---|
| 11 | `activity_types` | Activity measurement types present. |
| 12 | `activity_value_uM_min` | Minimum activity value in uM. |
| 13 | `activity_value_uM_median` | Median activity value in uM. |

### PubChem Chemical Properties (14-26)
From PubChem PUG REST. Populated for PubChem CID rows when available.
| # | Column | Meaning |
|---|---|---|
| 14 | `molecular_formula` | Molecular formula (e.g., `C9H8O4`). |
| 15 | `molecular_weight` | Molecular weight (g/mol). |
| 16 | `canonical_smiles` | Canonical SMILES string. |
| 17 | `isomeric_smiles` | Isomeric SMILES string. |
| 18 | `inchikey` | InChIKey (standard). |
| 19 | `iupac_name` | IUPAC systematic name. |
| 20 | `xlogp` | Computed logP (partition coefficient). |
| 21 | `tpsa` | Topological polar surface area. |
| 22 | `hbond_donor_count` | Hydrogen bond donor count. |
| 23 | `hbond_acceptor_count` | Hydrogen bond acceptor count. |
| 24 | `rotatable_bond_count` | Rotatable bond count. |
| 25 | `complexity` | Molecular complexity score. |
| 26 | `formal_charge` | Formal charge. |

### Safety/Toxicity (27-33)
From PubChem PUG-View. Limited coverage; `not_found_in_current_sources` means the current sources did not provide that field and does not prove absence of toxicity, indication, or clinical use.
| # | Column | Meaning |
|---|---|---|
| 27 | `ghs_hazard_classes` | GHS hazard class codes. |
| 28 | `ghs_signal_words` | GHS signal words (`Warning` / `Danger`) when available. |
| 29 | `toxicity_summary` | Human-readable toxicity summary. |
| 30 | `livertox` | LiverTox summary text when available. |
| 31 | `drug_classes` | Drug class/category text when available. |
| 32 | `pharmacodynamics` | Pharmacodynamic description when available. |
| 33 | `drug_indication` | Indication/use description when available. |
"""


SHEET3_GUIDE_SECTION = """## Sheet 3: `3_Small_Molecules`

Current shape: 87,552 rows x 33 columns. One row is one compound/drug identifier.

### Compound Identity and Role

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 1 | `drug_id` | Primary compound identifier used to join Sheet 1. | Usually a PubChem CID. Non-CID examples include ChEMBL/DrugCentral IDs for biologics and `UNMAPPED_PDB_LIGAND:*` for structure ligands with no reliable CID. Unique in this table. |
| 2 | `compound_id_type` | Type of `drug_id`. | Current values: `PubChem CID`, `ChEMBL`, `DrugCentral`, `unmapped_structure_ligand`. `PubChem CID` is the normal small-molecule case. |
| 3 | `drug_name` | Preferred display name for the compound. | Always filled. May be a readable drug name, research code, PDB ligand name, or systematic chemical name. |
| 4 | `compound_name_category` | What kind of name is shown in `drug_name`. | Values: `common_or_drug_name`, `research_code`, `systematic_or_chemical_name`, `pubchem_synonym_name`, `pdb_ligand_code_or_name`, `unknown_name_type`. |
| 5 | `compound_biological_role` | How to interpret the compound biologically or structurally. | Values include `therapeutic_or_clinical_compound`, `bioactive_research_ligand`, `structure_affinity_ligand`, `endogenous_ligand_or_cofactor`, `membrane_lipid_or_sterol`, `biologic_or_large_molecule`, and `unknown_structure_ligand`. This is for broad filtering, not a final binding-site relevance score. |
| 6 | `compound_source` | Broad origin class of this compound in the dataset. | Values can be combined with `;`: `curated_target_relation` for ChEMBL/DrugCentral-style curated target relations, `bioassay_active_relation` for PubChem BioAssay active records, `structure_ligand_context` for BioLiP/sc-PDB/PDBbind ligand context. |

### Relationship Summary From Sheet 1

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 7 | `source_databases` | Databases in which this compound has protein relationships. | Examples: `PubChem BioAssay`, `ChEMBL`, `DrugCentral`, `BioLiP`, `sc-PDB`, `PDBbind`, or combinations separated by `;`. |
| 8 | `source_row_count` | Number of Sheet 1 relationship rows involving this compound. | Integer. Larger values mean the compound appears with more relationship/evidence rows. |
| 9 | `unique_target_count` | Number of unique proteins linked to this compound. | Integer. `1` means a single target in the current dataset; larger values indicate multi-target compounds or broad screening hits. |
| 10 | `compound_evidence_status` | Compound-level evidence/status class. | `approved_drug` means approved/listed status was captured; `clinical_trial_compound` means clinical-trial status was captured but approval was not; `binding_evidence_only` means the current table only captured binding/activity/structure evidence. |
| 11 | `activity_types` | Activity measurement types observed for the compound. | Examples: `IC50`, `Ki`, `EC50`, `KD`, `Kd`, or combinations separated by `;`. Empty means no quantitative activity type captured. |
| 12 | `activity_value_uM_min` | Minimum activity value in micromolar across relationships. | Numeric. Lower values generally mean stronger measured activity, but compare only within compatible assay contexts. Empty means no numeric value captured. |
| 13 | `activity_value_uM_median` | Median activity value in micromolar across relationships. | Numeric summary of all available values. More robust than minimum when multiple assays exist. Empty means no numeric value captured. |

### PubChem Chemical Properties

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 14 | `molecular_formula` | Molecular formula from PubChem. | Example: `C9H8O4`. Empty mainly occurs for unmapped biologics, unmapped structure ligands, or CIDs without property retrieval. |
| 15 | `molecular_weight` | Molecular weight in g/mol. | Numeric. Useful for drug-likeness filters and comparing ligand size. |
| 16 | `canonical_smiles` | PubChem canonical SMILES. | Non-isomeric structural string. Useful for cheminformatics but does not fully preserve stereochemistry. |
| 17 | `isomeric_smiles` | PubChem isomeric SMILES. | Structural string with stereochemistry where available. Prefer this for modeling when stereochemistry matters. |
| 18 | `inchikey` | Standard InChIKey. | Useful for cross-database identity matching. Empty means no chemical structure mapping was available. |
| 19 | `iupac_name` | IUPAC systematic name. | Often long. Useful as precise chemical name, less useful for human browsing. |
| 20 | `xlogp` | PubChem XLogP estimate. | Numeric lipophilicity estimate. Higher values are more hydrophobic; empty when not computed. |
| 21 | `tpsa` | Topological polar surface area. | Numeric area in square angstrom-like units. Higher values often imply lower membrane permeability. |
| 22 | `hbond_donor_count` | Hydrogen bond donor count. | Integer. Used in drug-likeness and permeability heuristics. |
| 23 | `hbond_acceptor_count` | Hydrogen bond acceptor count. | Integer. Used in drug-likeness and permeability heuristics. |
| 24 | `rotatable_bond_count` | Number of rotatable bonds. | Integer. Higher values indicate more conformational flexibility. |
| 25 | `complexity` | PubChem molecular complexity score. | Numeric. Higher values generally indicate more structurally complex molecules. |
| 26 | `formal_charge` | Formal charge. | Integer, usually `0`. Nonzero values indicate charged compounds or salts/ions. |

### Safety, Toxicity, and Clinical Text

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 27 | `ghs_hazard_classes` | GHS hazard classifications from PubChem/PUG-View where available. | Examples include `Acute Tox.`, `Repr.`, `Carc.`. `not_found_in_current_sources` means not captured in the current sources. |
| 28 | `ghs_signal_words` | GHS signal word. | Intended values are `Warning` or `Danger`; `not_found_in_current_sources` means not captured. |
| 29 | `toxicity_summary` | Human-readable toxicity summary text from PubChem/PUG-View where available. | Long free text or `not_found_in_current_sources`. Absence of a captured summary does not mean no toxicity. |
| 30 | `livertox` | LiverTox summary text where available. | Long free text or `not_found_in_current_sources`. |
| 31 | `drug_classes` | Drug class or category text. | Coverage is sparse and mainly for known drugs; otherwise `not_found_in_current_sources`. |
| 32 | `pharmacodynamics` | Pharmacodynamic description. | Long free text or `not_found_in_current_sources`. |
| 33 | `drug_indication` | Indication or use description. | Long text/source labels or `not_found_in_current_sources`; absence does not mean no medical use. |
"""


def replace_between(text: str, start: str, end: str, replacement: str) -> str:
    a = text.index(start)
    b = text.index(end, a)
    return text[:a] + replacement.rstrip() + "\n\n" + text[b:]


def update_data_dictionary() -> None:
    path = BASE / "DATA_DICTIONARY.md"
    text = path.read_text(encoding="utf-8")
    text = text.replace(
        "| `small_molecule_database.tsv` | 87,552 | 32 | 42.7 MB | One row per unique compound/drug ID. Identity, activity, PubChem properties, safety/toxicity, clinical fields. |",
        "| `small_molecule_database.tsv` | 87,552 | 33 | 66.1 MB | One row per unique compound/drug ID. Identity, role, activity, PubChem properties, safety/toxicity, clinical fields. |",
    )
    text = text.replace(
        "`self-screening` (ChEMBL curated + DrugCentral), `non-screening` (PubChem BioAssay HTS), or `structure_context`",
        "`curated_target_relation` (ChEMBL curated + DrugCentral), `bioassay_active_relation` (PubChem BioAssay active records), or `structure_ligand_context`",
    )
    text = text.replace("compound_source = structure_context", "compound_source = structure_ligand_context")
    text = replace_between(text, "## `small_molecule_database.tsv`", "---\n\n## `pocket_instances.tsv`", SMALL_MOLECULE_DICTIONARY_SECTION)
    text = text.replace(
        "| `compound_source`: self-screening vs non-screening | Separates curated pharmacological targets (9.7%) from HTS hits (90.3%). |",
        "| `compound_source` evidence-origin labels | Uses `curated_target_relation`, `bioassay_active_relation`, and `structure_ligand_context` to separate curated target relations, active bioassays, and structure ligand context. |",
    )
    path.write_text(text, encoding="utf-8")


def update_sheet_guide() -> None:
    path = BASE / "SHEET2_SHEET3_COLUMN_GUIDE.md"
    text = path.read_text(encoding="utf-8")
    text = replace_between(text, "## Sheet 3: `3_Small_Molecules`", "## Practical Reading Tips", SHEET3_GUIDE_SECTION)
    path.write_text(text, encoding="utf-8")


def update_logs() -> None:
    process = BASE / "PROCESS_LOG.md"
    text = process.read_text(encoding="utf-8")
    text = text.replace(
        "| `upload/normalized_tables/small_molecule_database.tsv` | 87,552 | 32 | 42.7 MB |",
        "| `upload/normalized_tables/small_molecule_database.tsv` | 87,552 | 33 | 66.1 MB |",
    )
    text = text.replace(
        "- Removed `evidence_levels` from `protein_database.tsv` and `small_molecule_database.tsv`; it duplicated broad source/category information and did not provide useful evidence grading for these entity summary tables.",
        "- Removed `evidence_levels` from `protein_database.tsv` and `small_molecule_database.tsv`; it duplicated broad source/category information and did not provide useful evidence grading for these entity summary tables.\n- Renamed current `compound_source` values to `curated_target_relation`, `bioassay_active_relation`, and `structure_ligand_context`.\n- Replaced small-molecule `clinical_statuses` with `compound_evidence_status`: `approved_drug`, `clinical_trial_compound`, or `binding_evidence_only`.\n- Split `drug_name_type` into `compound_name_category` and `compound_biological_role`.\n- Filled high-missing safety/clinical text fields with `not_found_in_current_sources` when no value was captured.",
    )
    text = text.replace("| `small_molecule_database.tsv` | 87,552 | 32 |", "| `small_molecule_database.tsv` | 87,552 | 33 |")
    process.write_text(text, encoding="utf-8")

    full = BASE / "FULL_WORK_LOG.md"
    text = full.read_text(encoding="utf-8")
    text = text.replace("small_molecule_database.tsv             87,552 rows 脳 32 cols", "small_molecule_database.tsv             87,552 rows 脳 33 cols")
    text = text.replace(
        "- Removed `evidence_levels` from `protein_database.tsv` and `small_molecule_database.tsv`; it was too coarse for the current entity-summary tables and duplicated source/context information.",
        "- Removed `evidence_levels` from `protein_database.tsv` and `small_molecule_database.tsv`; it was too coarse for the current entity-summary tables and duplicated source/context information.\n- Renamed current `compound_source` values to `curated_target_relation`, `bioassay_active_relation`, and `structure_ligand_context`.\n- Replaced small-molecule `clinical_statuses` with `compound_evidence_status` and split `drug_name_type` into `compound_name_category` plus `compound_biological_role`.\n- Filled high-missing safety/clinical text fields with `not_found_in_current_sources` to distinguish missing source coverage from biological absence.",
    )
    text = text.replace("| `small_molecule_database.tsv` | 87,552 | 32 |", "| `small_molecule_database.tsv` | 87,552 | 33 |")
    full.write_text(text, encoding="utf-8")


def update_manifest() -> None:
    path = BASE / "MANIFEST.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    tracked = [
        "drug_protein_binary_relationships.tsv",
        "protein_database.tsv",
        "protein_source_summaries.tsv",
        "normalized_tables_v3.xlsx",
        "small_molecule_database.tsv",
        "pocket_instances.tsv",
        "DATA_DICTIONARY.md",
        "PROCESS_LOG.md",
        "FULL_WORK_LOG.md",
        "SHEET2_SHEET3_COLUMN_GUIDE.md",
        "drug_name_normalization_report.tsv",
        "protein_structure_ligand_context_review.tsv",
        "protein_structure_ligand_context_review.xlsx",
        "protein_structure_ligand_context_review_summary.json",
        "protein_structure_ligand_context_binary_format_review.tsv",
        "protein_structure_ligand_context_binary_format_review.xlsx",
        "structure_context_merge_report.json",
        "pubchem_structure_context_properties_cache.json",
        "stitch_cid_only_normalization_report.json",
    ]
    items = []
    for name in tracked:
        p = BASE / name
        if p.exists():
            items.append(
                {
                    "path": f"upload/normalized_tables/{name}",
                    "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                    "bytes": p.stat().st_size,
                }
            )
    data["files"] = items
    data["outputs"] = [
        item
        for item in items
        if Path(item["path"]).name
        in {
            "drug_protein_binary_relationships.tsv",
            "protein_database.tsv",
            "small_molecule_database.tsv",
            "pocket_instances.tsv",
            "normalized_tables_v3.xlsx",
        }
    ]
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    update_data_dictionary()
    update_sheet_guide()
    update_logs()
    update_manifest()


if __name__ == "__main__":
    main()
