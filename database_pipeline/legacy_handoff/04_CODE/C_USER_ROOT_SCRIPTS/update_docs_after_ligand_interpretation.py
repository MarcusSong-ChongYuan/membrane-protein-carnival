from __future__ import annotations

import hashlib
import json
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")

BINARY_SECTION = """## `drug_protein_binary_relationships.tsv`

Column order: protein identity -> compound identity/role -> ligand interpretation -> relationship/evidence -> binding-site details.

### Protein Identity (1-4)
| # | Column | Meaning |
|---|---|---|
| 1 | `target_uniprot_id` | UniProt accession; primary protein identifier. |
| 2 | `approved_symbol` | HGNC-approved gene symbol. |
| 3 | `ncbi_gene_id` | NCBI Gene identifier mapped for the protein. |
| 4 | `protein_name` | Human-readable protein name. |

### Compound Identity, Role, and Interpretation (5-11)
| # | Column | Meaning |
|---|---|---|
| 5 | `drug_id` | Compound identifier. Usually PubChem CID; unmapped biologics retain ChEMBL/DrugCentral ID; unmapped structure ligands use stable IDs such as `UNMAPPED_PDB_LIGAND:ITC`. |
| 6 | `compound_id_type` | `PubChem CID`, original type for unmapped biologics, or `unmapped_structure_ligand`. |
| 7 | `drug_name` | Preferred compound display name. |
| 8 | `compound_name_category` | Type of display name: `common_or_drug_name`, `research_code`, `systematic_or_chemical_name`, `pubchem_synonym_name`, `pdb_ligand_code_or_name`, or `unknown_name_type`. |
| 9 | `compound_biological_role` | Broad role of the compound itself, such as `therapeutic_or_clinical_compound`, `bioactive_research_ligand`, `structure_affinity_ligand`, `endogenous_ligand_or_cofactor`, `membrane_lipid_or_sterol`, `biologic_or_large_molecule`, or `unknown_structure_ligand`. |
| 10 | `ligand_interpretation_class` | Row-level interpretation of why this ligand is useful in this compound-protein relationship. Values include `structure_affinity_ligand`, `approved_or_clinical_target_ligand`, `cofactor_or_functional_ligand`, `membrane_lipid_or_sterol`, `curated_target_ligand`, `structure_bound_ligand`, `bioassay_active_ligand`, `ion_or_metal`, `buffer_salt_solvent`, and `unknown_structure_ligand`. |
| 11 | `compound_source` | Relationship evidence origin: `curated_target_relation`, `bioassay_active_relation`, or `structure_ligand_context`. |

### Relationship Evidence (12-16)
| # | Column | Meaning |
|---|---|---|
| 12 | `source_database` | Originating database(s) for this compound-protein row. |
| 13 | `assay_or_mechanism` | Assay type, mechanism-of-action text, or structure-context descriptor. |
| 14 | `activity_type` | `IC50`, `Ki`, `Kd`, `EC50`, etc. |
| 15 | `activity_value_uM` | Activity value in uM when available. |
| 16 | `clinical_or_approval_status` | Row-level captured status from source relationship data: `approved_or_listed`, `clinical_trial`, combined value, or empty. |

### Matched Binding Evidence and Details (17-31)
Boolean flags use `1` = present and `0` = absent.
| # | Column | Meaning |
|---|---|---|
| 17 | `has_pdb_biolip_matched` | BioLiP co-crystal/context evidence for this row. |
| 18 | `pdb_biolip_matched_sites` | BioLiP site detail, usually `LIGAND(CODE)@PDB:CHAIN=RESIDUES`. |
| 19 | `has_pdb_scpdb_matched` | sc-PDB co-crystal/context evidence for this row. |
| 20 | `pdb_scpdb_matched_sites` | sc-PDB site detail. |
| 21 | `has_pdbbind_matched` | PDBbind affinity/structure context for this row. |
| 22 | `pdbbind_matched_sites` | PDBbind site/affinity detail. |
| 23 | `has_stitch_matched` | STITCH evidence for this compound-protein pair. |
| 24 | `stitch_matched_compounds` | `CIDxxx:confidence`; ChEMBL prefixes from the original STITCH mapping are removed. |
| 25 | `has_exp_binding_site` | UniProt literature binding-site annotation exists for this protein. |
| 26 | `exp_binding_sites` | Protein-level residue binding-site annotations. |
| 27 | `exp_ligands` | Known ligands from UniProt experimental annotation. |
| 28 | `exp_evidence_quality` | UniProt/ECO evidence quality tiers. |
| 29 | `exp_pubmed` | PubMed IDs for binding-site evidence. |
| 30 | `exp_pdb` | Cross-referenced PDB entries when captured. |
| 31 | `target_group_id` | Group ID preserving source rows that originally listed multiple UniProt targets. |
"""

SMALL_SECTION = """## `small_molecule_database.tsv`

### Identity, Role, and Interpretation (1-8)
| # | Column | Meaning |
|---|---|---|
| 1 | `drug_id` | Primary compound identifier. Usually PubChem CID; may be ChEMBL/DrugCentral for unmapped biologics or `UNMAPPED_PDB_LIGAND:*` for structure ligands with no CID. |
| 2 | `compound_id_type` | `PubChem CID`, original type for biologics, or `unmapped_structure_ligand`. |
| 3 | `drug_name` | Preferred display name from the source table or PubChem/PDB mapping. |
| 4 | `compound_name_category` | Type of display name: `common_or_drug_name`, `research_code`, `systematic_or_chemical_name`, `pubchem_synonym_name`, `pdb_ligand_code_or_name`, or `unknown_name_type`. |
| 5 | `compound_biological_role` | Biological/use role for filtering, such as `therapeutic_or_clinical_compound`, `bioactive_research_ligand`, `structure_affinity_ligand`, `endogenous_ligand_or_cofactor`, `membrane_lipid_or_sterol`, `biologic_or_large_molecule`, or `unknown_structure_ligand`. |
| 6 | `ligand_interpretation_classes` | Semicolon-separated set of all row-level ligand interpretation classes observed for this compound in the binary table. |
| 7 | `best_ligand_interpretation_class` | Single highest-priority interpretation class for quick compound-level filtering. Priority: structure affinity > approved/clinical > cofactor > membrane lipid/sterol > curated target > structure-bound > bioassay active > ion/metal > buffer/salt/solvent > unknown structure ligand. |
| 8 | `compound_source` | Relationship evidence origin: `curated_target_relation`, `bioassay_active_relation`, and/or `structure_ligand_context`. |

### Source and Evidence Summary (9-12)
| # | Column | Meaning |
|---|---|---|
| 9 | `source_databases` | Semicolon-separated source databases. |
| 10 | `source_row_count` | Number of Sheet 1 relationship rows for this compound. |
| 11 | `unique_target_count` | Number of unique proteins this compound targets. |
| 12 | `compound_evidence_status` | Compound-level status: `approved_drug`, `clinical_trial_compound`, or `binding_evidence_only`. |

### Activity Summary (13-15)
| # | Column | Meaning |
|---|---|---|
| 13 | `activity_types` | Activity measurement types present. |
| 14 | `activity_value_uM_min` | Minimum activity value in uM. |
| 15 | `activity_value_uM_median` | Median activity value in uM. |

### PubChem Chemical Properties (16-28)
From PubChem PUG REST. Populated for PubChem CID rows when available.
| # | Column | Meaning |
|---|---|---|
| 16 | `molecular_formula` | Molecular formula (e.g., `C9H8O4`). |
| 17 | `molecular_weight` | Molecular weight (g/mol). |
| 18 | `canonical_smiles` | Canonical SMILES string. |
| 19 | `isomeric_smiles` | Isomeric SMILES string. |
| 20 | `inchikey` | InChIKey (standard). |
| 21 | `iupac_name` | IUPAC systematic name. |
| 22 | `xlogp` | Computed logP (partition coefficient). |
| 23 | `tpsa` | Topological polar surface area. |
| 24 | `hbond_donor_count` | Hydrogen bond donor count. |
| 25 | `hbond_acceptor_count` | Hydrogen bond acceptor count. |
| 26 | `rotatable_bond_count` | Rotatable bond count. |
| 27 | `complexity` | Molecular complexity score. |
| 28 | `formal_charge` | Formal charge. |

### Safety/Toxicity (29-35)
From PubChem PUG-View. Limited coverage; `not_found_in_current_sources` means the current sources did not provide that field and does not prove absence of toxicity, indication, or clinical use.
| # | Column | Meaning |
|---|---|---|
| 29 | `ghs_hazard_classes` | GHS hazard class codes. |
| 30 | `ghs_signal_words` | GHS signal words (`Warning` / `Danger`) when available. |
| 31 | `toxicity_summary` | Human-readable toxicity summary. |
| 32 | `livertox` | LiverTox summary text when available. |
| 33 | `drug_classes` | Drug class/category text when available. |
| 34 | `pharmacodynamics` | Pharmacodynamic description when available. |
| 35 | `drug_indication` | Indication/use description when available. |
"""

SHEET3_SECTION = """## Sheet 3: `3_Small_Molecules`

Current shape: 87,552 rows x 35 columns. One row is one compound/drug identifier.

### Compound Identity, Role, and Interpretation

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 1 | `drug_id` | Primary compound identifier used to join Sheet 1. | Usually a PubChem CID. Non-CID examples include ChEMBL/DrugCentral IDs for biologics and `UNMAPPED_PDB_LIGAND:*` for structure ligands with no reliable CID. Unique in this table. |
| 2 | `compound_id_type` | Type of `drug_id`. | Current values: `PubChem CID`, `ChEMBL`, `DrugCentral`, `unmapped_structure_ligand`. |
| 3 | `drug_name` | Preferred display name for the compound. | Always filled. May be a readable drug name, research code, PDB ligand name, or systematic chemical name. |
| 4 | `compound_name_category` | What kind of name is shown in `drug_name`. | `common_or_drug_name`, `research_code`, `systematic_or_chemical_name`, `pubchem_synonym_name`, `pdb_ligand_code_or_name`, `unknown_name_type`. |
| 5 | `compound_biological_role` | What the compound itself broadly is. | Examples: `therapeutic_or_clinical_compound`, `bioactive_research_ligand`, `structure_affinity_ligand`, `endogenous_ligand_or_cofactor`, `membrane_lipid_or_sterol`, `biologic_or_large_molecule`, `unknown_structure_ligand`. |
| 6 | `ligand_interpretation_classes` | All ligand interpretation classes observed for this compound in Sheet 1. | Semicolon-separated; preserves multiple reasons a compound may matter. |
| 7 | `best_ligand_interpretation_class` | Single highest-priority interpretation class. | Use for fast filtering. Priority is structure affinity > approved/clinical > cofactor > membrane lipid/sterol > curated target > structure-bound > bioassay active > ion/metal > buffer/salt/solvent > unknown structure ligand. |
| 8 | `compound_source` | Broad origin class of this compound in the dataset. | Values can be combined with `;`: `curated_target_relation`, `bioassay_active_relation`, `structure_ligand_context`. |

### Relationship Summary From Sheet 1

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 9 | `source_databases` | Databases in which this compound has protein relationships. | Examples: `PubChem BioAssay`, `ChEMBL`, `DrugCentral`, `BioLiP`, `sc-PDB`, `PDBbind`, or combinations separated by `;`. |
| 10 | `source_row_count` | Number of Sheet 1 relationship rows involving this compound. | Integer. Larger values mean the compound appears with more relationship/evidence rows. |
| 11 | `unique_target_count` | Number of unique proteins linked to this compound. | Integer. `1` means a single target in the current dataset; larger values indicate multi-target compounds or broad screening hits. |
| 12 | `compound_evidence_status` | Compound-level evidence/status class. | `approved_drug`, `clinical_trial_compound`, or `binding_evidence_only`. |
| 13 | `activity_types` | Activity measurement types observed for the compound. | Examples: `IC50`, `Ki`, `EC50`, `KD`, `Kd`, or combinations separated by `;`. |
| 14 | `activity_value_uM_min` | Minimum activity value in micromolar across relationships. | Numeric. Lower values generally mean stronger measured activity, but compare only within compatible assay contexts. |
| 15 | `activity_value_uM_median` | Median activity value in micromolar across relationships. | Numeric summary of all available values. |

### PubChem Chemical Properties

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 16 | `molecular_formula` | Molecular formula from PubChem. | Example: `C9H8O4`. Empty mainly occurs for unmapped biologics, unmapped structure ligands, or CIDs without property retrieval. |
| 17 | `molecular_weight` | Molecular weight in g/mol. | Numeric. Useful for drug-likeness filters and comparing ligand size. |
| 18 | `canonical_smiles` | PubChem canonical SMILES. | Non-isomeric structural string. |
| 19 | `isomeric_smiles` | PubChem isomeric SMILES. | Structural string with stereochemistry where available. |
| 20 | `inchikey` | Standard InChIKey. | Useful for cross-database identity matching. |
| 21 | `iupac_name` | IUPAC systematic name. | Often long. |
| 22 | `xlogp` | PubChem XLogP estimate. | Numeric lipophilicity estimate. |
| 23 | `tpsa` | Topological polar surface area. | Numeric area in square angstrom-like units. |
| 24 | `hbond_donor_count` | Hydrogen bond donor count. | Integer. |
| 25 | `hbond_acceptor_count` | Hydrogen bond acceptor count. | Integer. |
| 26 | `rotatable_bond_count` | Number of rotatable bonds. | Integer. |
| 27 | `complexity` | PubChem molecular complexity score. | Numeric. |
| 28 | `formal_charge` | Formal charge. | Integer, usually `0`. |

### Safety, Toxicity, and Clinical Text

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 29 | `ghs_hazard_classes` | GHS hazard classifications from PubChem/PUG-View where available. | `not_found_in_current_sources` means not captured in the current sources. |
| 30 | `ghs_signal_words` | GHS signal word. | Intended values are `Warning` or `Danger`; `not_found_in_current_sources` means not captured. |
| 31 | `toxicity_summary` | Human-readable toxicity summary text. | Long free text or `not_found_in_current_sources`. |
| 32 | `livertox` | LiverTox summary text where available. | Long free text or `not_found_in_current_sources`. |
| 33 | `drug_classes` | Drug class or category text. | Sparse coverage; otherwise `not_found_in_current_sources`. |
| 34 | `pharmacodynamics` | Pharmacodynamic description. | Long free text or `not_found_in_current_sources`. |
| 35 | `drug_indication` | Indication or use description. | Long text/source labels or `not_found_in_current_sources`. |
"""


def replace_between(text: str, start: str, end: str, replacement: str) -> str:
    a = text.index(start)
    b = text.index(end, a)
    return text[:a] + replacement.rstrip() + "\n\n" + text[b:]


def refresh_manifest() -> None:
    path = BASE / "MANIFEST.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    tracked = [Path(item["path"]).name for item in data.get("files", [])]
    if "MANIFEST.json" not in tracked:
        tracked.append("MANIFEST.json")
    items = []
    for name in tracked:
        p = BASE / name
        if p.exists() and name != "MANIFEST.json":
            items.append({"path": f"upload/normalized_tables/{name}", "sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "bytes": p.stat().st_size})
    data["files"] = items
    data["outputs"] = [item for item in items if Path(item["path"]).name in {
        "drug_protein_binary_relationships.tsv",
        "protein_database.tsv",
        "small_molecule_database.tsv",
        "pocket_instances.tsv",
        "normalized_tables_v3.xlsx",
    }]
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    data_dict = BASE / "DATA_DICTIONARY.md"
    text = data_dict.read_text(encoding="utf-8")
    text = text.replace(
        "| `drug_protein_binary_relationships.tsv` | 101,171 | 30 | 53.0 MB | One row per compound-protein relationship/context row, including curated structure ligand context rows. |",
        "| `drug_protein_binary_relationships.tsv` | 101,171 | 31 | 60.5 MB | One row per compound-protein relationship/context row, including curated structure ligand context rows. |",
    )
    text = text.replace(
        "| `small_molecule_database.tsv` | 87,552 | 33 | 66.1 MB | One row per unique compound/drug ID. Identity, role, activity, PubChem properties, safety/toxicity, clinical fields. |",
        "| `small_molecule_database.tsv` | 87,552 | 35 | 70.2 MB | One row per unique compound/drug ID. Identity, role, interpretation, activity, PubChem properties, safety/toxicity, clinical fields. |",
    )
    text = replace_between(text, "## `drug_protein_binary_relationships.tsv`", "---\n\n## `protein_database.tsv`", BINARY_SECTION)
    text = replace_between(text, "## `small_molecule_database.tsv`", "---\n\n## `pocket_instances.tsv`", SMALL_SECTION)
    data_dict.write_text(text, encoding="utf-8")

    guide = BASE / "SHEET2_SHEET3_COLUMN_GUIDE.md"
    text = guide.read_text(encoding="utf-8")
    text = replace_between(text, "## Sheet 3: `3_Small_Molecules`", "## Practical Reading Tips", SHEET3_SECTION)
    guide.write_text(text, encoding="utf-8")

    for name in ["PROCESS_LOG.md", "FULL_WORK_LOG.md"]:
        p = BASE / name
        text = p.read_text(encoding="utf-8")
        text = text.replace("| `upload/normalized_tables/drug_protein_binary_relationships.tsv` | 101,171 | 30 | 53.0 MB |", "| `upload/normalized_tables/drug_protein_binary_relationships.tsv` | 101,171 | 31 | 60.5 MB |")
        text = text.replace("| `upload/normalized_tables/small_molecule_database.tsv` | 87,552 | 33 | 66.1 MB |", "| `upload/normalized_tables/small_molecule_database.tsv` | 87,552 | 35 | 70.2 MB |")
        text = text.replace("| `drug_protein_binary_relationships.tsv` | 101,171 | 30 |", "| `drug_protein_binary_relationships.tsv` | 101,171 | 31 |")
        text = text.replace("| `small_molecule_database.tsv` | 87,552 | 33 |", "| `small_molecule_database.tsv` | 87,552 | 35 |")
        text = text.replace("drug_protein_binary_relationships.tsv   101,171 rows 脳 30 cols", "drug_protein_binary_relationships.tsv   101,171 rows 脳 31 cols")
        text = text.replace("small_molecule_database.tsv             87,552 rows 脳 33 cols", "small_molecule_database.tsv             87,552 rows 脳 35 cols")
        marker = "- Filled high-missing safety/clinical text fields with `not_found_in_current_sources`"
        if marker in text and "ligand_interpretation_class" not in text:
            text = text.replace(
                marker,
                "- Added row-level `ligand_interpretation_class` to the binary table and compound-level `ligand_interpretation_classes` plus `best_ligand_interpretation_class` to the small-molecule table.\n" + marker,
            )
        p.write_text(text, encoding="utf-8")

    refresh_manifest()


if __name__ == "__main__":
    main()
