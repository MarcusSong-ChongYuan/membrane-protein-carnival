# Sheet 2 and Sheet 3 Column Guide

Date: 2026-07-09
Workbook: `normalized_tables_v3.xlsx`

This guide explains the current Sheet 2 (`2_Protein_Database`) and Sheet 3 (`3_Small_Molecules`) columns after cleanup.

General conventions:
- Empty cell means no value was available from the current sources, not necessarily biological absence.
- Semicolon `;` separates multiple values inside one cell.
- Counts such as `0` mean the feature was explicitly counted as zero in the current processed table.
- For complete one-to-many records, use the normalized long tables: drug-protein relationships in Sheet 1 and AlphaFold pocket instances in Sheet 4.

## Sheet 2: `2_Protein_Database`

Current shape: 1,773 rows x 25 columns. One row is one UniProt protein.

### Protein Identity

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 1 | `target_uniprot_id` | Primary protein ID. This is the UniProt accession and the main key for joining to Sheet 1 and Sheet 4. | Always filled and unique in this table. Example: `P08172`. |
| 2 | `approved_symbol` | HGNC-style gene symbol used as a readable protein/gene label. | Always filled. Example: `CHRM2`, `KCNH2`, `PTGS1`. Treat as a label, not the primary key. |
| 3 | `ncbi_gene_id` | NCBI Gene identifier mapped for the protein. | Usually a numeric ID. Some rows have multiple IDs separated by `;`, such as `801;805;808`. Empty means no NCBI Gene ID was mapped. |
| 4 | `protein_name` | Full protein name, generally from UniProt or source annotation. | Always filled. This is descriptive text, not a controlled identifier. |

### AlphaFold Pocket Summary

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 5 | `n_pockets_alphafold` | Number of predicted AlphaFold pockets for this protein. | Integer from the pocket pipeline. `0` means the protein is retained in the protein table but has no pocket instance in Sheet 4. Full pocket details are in Sheet 4. |

### STITCH Context

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 6 | `stitch_compounds_original_count` | Number of original STITCH compounds found for the protein before high-confidence filtering. | Integer. `0` means no original STITCH compound list was available for this protein. This is a pre-cleaning count. |
| 7 | `stitch_compounds_cleaned` | STITCH compounds retained after filtering to high-confidence CID-protein pairs. | Format is `CID:score`, separated by `;`, for example `CID5494449:996`. Score is STITCH combined score. Empty means no retained STITCH context after filtering. |
| 8 | `stitch_compounds_cleaned_count` | Number of retained STITCH compounds in `stitch_compounds_cleaned`. | Integer. `0` means no retained cleaned STITCH compounds. This should match the number of entries in column 7. |

### UniProt Experimental Binding Site Annotations

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 9 | `exp_binding_sites` | UniProt literature/curated binding site annotations at residue level. | Format resembles `ligand@start-end(quality)`, separated by `;`. Example: `ATP@219-226(medium)`. Empty means no UniProt binding-site annotation was captured. |
| 10 | `exp_ligands` | Ligand names extracted from UniProt experimental binding-site annotations. | Examples: `ATP`, `Ca(2+)`, `L-glutamate`, `heme`. Multiple ligands are separated by `;`. Empty means no ligand annotation was captured. |
| 11 | `exp_evidence_quality` | Residue ranges and evidence quality for experimental annotations. | Format resembles `219-226:medium` or `697-697:high`. `high` usually reflects stronger literature evidence; `medium` is still useful but less direct. |
| 12 | `exp_pubmed` | PubMed IDs supporting experimental binding site annotations. | Semicolon-separated PubMed IDs. Empty means no PubMed ID was captured for the annotation. |
| 13 | `exp_pdb` | PDB IDs linked from experimental annotations. | Currently empty in this dataset. Keep as a reserved field unless later evidence fills it. |

### Relationship Summary From Sheet 1

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 14 | `unique_drug_count` | Count of distinct `drug_id` values linked to this protein in Sheet 1. | Integer. Includes PubChem CID, retained biologic IDs, and `UNMAPPED_PDB_LIGAND:*` structure-context IDs. `0` means no drug/small-molecule relationship in Sheet 1. |
| 15 | `unique_drug_cid_sample` | Preview of up to 50 distinct `drug_id` values linked to the protein. | Semicolon-separated sample, not a complete list when `unique_drug_count > 50`. Use Sheet 1 for complete relationships. |
| 16 | `source_databases` | Aggregated databases contributing relationships for this protein. | Examples: `ChEMBL`, `DrugCentral`, `PubChem BioAssay`, `BioLiP`, `sc-PDB`, `PDBbind`. Multiple values are separated by `;`. Empty means no Sheet 1 relationship. |
| 17 | `has_matched_evidence` | Whether any matched evidence exists for this protein in the relationship/evidence system. | `1` means at least one matched evidence flag exists, such as STITCH, UniProt experimental site, or structure context. `0` means none in the current table. |

### UniProt Annotation

| # | Column | Meaning | Values and interpretation |
|---:|---|---|---|
| 18 | `reviewed` | UniProt official reviewed status. | Value `reviewed` means the UniProt entry is Swiss-Prot reviewed. Empty usually means unreviewed or not filled from the enrichment step. This replaces the removed internal `annotation_status` field. |
| 19 | `go_cellular_component` | GO cellular component annotations. | Semicolon-separated GO terms with IDs, for example `plasma membrane [GO:0005886]`. Empty means no GO cellular component was captured. |
| 20 | `go_molecular_function` | GO molecular function annotations. | Semicolon-separated GO terms with IDs, for example `ATP binding [GO:0005524]`. Empty means no GO molecular function was captured. |
| 21 | `pdb_structures` | Count of linked PDB structures from UniProt enrichment. | Integer. `0` means UniProt did not list linked PDB structures in the enrichment data. |
| 22 | `disease_association` | UniProt disease annotation text. | Long free text. Empty means no disease association was captured, not necessarily no disease relevance. |
| 23 | `transmembrane_count` | Count of UniProt TRANSMEM features. | Integer. `0` means no UniProt transmembrane feature was counted. Higher values usually indicate multi-pass membrane proteins. |
| 24 | `subcellular_location` | UniProt subcellular location comment. | Long text such as membrane, secreted, nucleus, extracellular matrix. Empty means no location comment was captured. |
| 25 | `protein_function` | UniProt function comment. | Long functional description. Empty means no function comment was captured. |

## Sheet 3: `3_Small_Molecules`

Current shape: 87,611 rows x 35 columns. One row is one compound/drug identifier.

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

## Practical Reading Tips

- For protein-level browsing, use Sheet 2.
- For compound-level browsing, use Sheet 3.
- For exact molecule-protein relationships, always use Sheet 1.
- For AlphaFold pocket geometry, always use Sheet 4.
- If a Sheet 2 or Sheet 3 field contains semicolon-separated values, it is a summary field. Use Sheet 1 or Sheet 4 for row-level detail.
