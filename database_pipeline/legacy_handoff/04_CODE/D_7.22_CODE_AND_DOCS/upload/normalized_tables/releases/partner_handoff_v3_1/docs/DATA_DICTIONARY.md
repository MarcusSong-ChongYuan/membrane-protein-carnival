# Normalized Tables Data Dictionary — v3 (CID Edition)

Date: 2026-07-07
Scope: 4-table normalized split from mempro1 v1.0 master table, with full CID conversion,
       STITCH cleaning, compound_source classification, and AlphaFold column suffix.
       Protein-first column ordering for the binary relationships table.

## Files

| File | Rows | Cols | Size | Meaning |
|---|---|---|---|---|
| `drug_protein_binary_relationships.tsv` | 101,802 | 31 | 60.5 MB | One row per compound-protein relationship/context row, including curated structure ligand context rows. |
| `protein_database.tsv` | 1,773 | 25 | 2.8 MB | One row per UniProt protein. Drug counts and source summaries are recomputed from the binary table. |
| `small_molecule_database.tsv` | 87,611 | 35 | 70.2 MB | One row per unique compound/drug ID. Identity, role, interpretation, activity, PubChem properties, safety/toxicity, clinical fields. |
| `pocket_instances.tsv` | 6,815 | 10 | 1.4 MB | One row per predicted AlphaFold pocket. Reused from Mintong v2 with debug columns removed. |
| `protein_source_summaries.tsv` | 811 | 13 | 0.2 MB | Audit-only archived source-level protein summary records from the pre-collapse protein table; not included as a result sheet in `normalized_tables_v3.xlsx`. |

---

## `drug_protein_binary_relationships.tsv`

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
| 6 | `compound_id_type` | `PubChem CID`, original type for unmapped biologics, `unmapped_structure_ligand`, or `unmapped_uniprot_ligand`. |
| 7 | `drug_name` | Preferred compound display name. |
| 8 | `compound_name_category` | Type of display name: `common_or_drug_name`, `research_code`, `systematic_or_chemical_name`, `pubchem_synonym_name`, `pdb_ligand_code_or_name`, or `unknown_name_type`. |
| 9 | `compound_biological_role` | Broad role of the compound itself, such as `therapeutic_or_clinical_compound`, `bioactive_research_ligand`, `structure_affinity_ligand`, `endogenous_ligand_or_cofactor`, `membrane_lipid_or_sterol`, `biologic_or_large_molecule`, or `unknown_structure_ligand`. |
| 10 | `ligand_interpretation_class` | Row-level interpretation of why this ligand is useful in this compound-protein relationship. Values include `structure_affinity_ligand`, `approved_or_clinical_target_ligand`, `cofactor_or_functional_ligand`, `membrane_lipid_or_sterol`, `uniprot_experimental_site_ligand`, `curated_target_ligand`, `structure_bound_ligand`, `bioassay_active_ligand`, `ion_or_metal`, `buffer_salt_solvent`, and `unknown_structure_ligand`. |
| 11 | `compound_source` | Relationship evidence origin: `curated_target_relation`, `bioassay_active_relation`, `structure_ligand_context`, or `uniprot_binding_site_annotation`. |

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

---

## `protein_database.tsv`

### Protein Identity (1-4)
| # | Column | Meaning |
|---|---|---|
| 1 | `target_uniprot_id` | UniProt accession (primary key). |
| 2 | `approved_symbol` | Gene symbol. |
| 3 | `ncbi_gene_id` | NCBI Gene ID. |
| 4 | `protein_name` | Full protein name. |

Dataset-level constants removed from this table: all rows are human (`Homo sapiens`) membrane proteins.

### AlphaFold Predicted Pocket Summary (5)
Detailed AlphaFold pocket instances live in `pocket_instances.tsv`; the protein entity table keeps only the protein-level count.
| # | Column | Meaning |
|---|---|---|
| 5 | `n_pockets_alphafold` | Number of predicted pockets. |

Pocket IDs, centers, scores, residues, and pocket batch provenance are stored in `pocket_instances.tsv`, not in the protein entity table.

### STITCH - Cleaned (6-8)
Only compounds with high-confidence (CID, UniProt) pairs in the master table.
| # | Column | Meaning |
|---|---|---|
| 6 | `stitch_compounds_original_count` | Total STITCH entries before cleaning. |
| 7 | `stitch_compounds_cleaned` | Semicolon-separated `CIDxxx:confidence` entries retained after cleaning. |
| 8 | `stitch_compounds_cleaned_count` | Number of entries retained. |

Cleaning rule: 83,372 -> 26,288 (31% retained).

PDB/sc-PDB/BioLiP/PDBbind same-protein other-ligand strings were removed from this protein entity table after they were curated into `protein_structure_ligand_context_review.tsv` and merged into the binary table as `compound_source = structure_ligand_context`.

### UniProt Experimental Sites (9-13)
| # | Column | Meaning |
|---|---|---|
| 9 | `exp_binding_sites` | Deduplicated binding site descriptions. |
| 10 | `exp_ligands` | Deduplicated known ligands. |
| 11 | `exp_evidence_quality` | Evidence quality tiers. |
| 12 | `exp_pubmed` | Deduplicated PubMed IDs. |
| 13 | `exp_pdb` | Deduplicated PDB cross-references. |

### Relationship Summary (14-17)
| # | Column | Meaning | Typical Range |
|---|---|---|---|
| 14 | `unique_drug_count` | Count of distinct `drug_id` values linked to this protein in `drug_protein_binary_relationships.tsv`. | |
| 15 | `unique_drug_cid_sample` | Preview list of up to 50 distinct `drug_id` values from the binary table. Complete relationships live in the binary table. | |
| 16 | `source_databases` | Semicolon-separated source databases aggregated from the binary table. | |
| 17 | `has_matched_evidence` | `1` if any matched PDB/STITCH/structure-context evidence exists. | `0` or `1` |

### UniProt Enrichment (18-25)
| # | Column | Meaning |
|---|---|---|
| 18 | `reviewed` | UniProt reviewed status when available. |
| 19 | `go_cellular_component` | GO cellular component annotations. |
| 20 | `go_molecular_function` | GO molecular function annotations. |
| 21 | `pdb_structures` | Count of linked PDB structures from UniProt enrichment. |
| 22 | `disease_association` | UniProt disease annotation text. |
| 23 | `transmembrane_count` | Count of UniProt TRANSMEM features. |
| 24 | `subcellular_location` | UniProt subcellular location comment. |
| 25 | `protein_function` | UniProt function comment. |

---

## `small_molecule_database.tsv`

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

---

## `pocket_instances.tsv` (Mintong v2, Reused)

| # | Column | Meaning |
|---|---|---|
| 1 | `pocket_batch` | Batch label: `1` or `2`. |
| 2 | `target_uniprot_id` | UniProt accession. |
| 3 | `pocket_id` | Unique pocket identifier (`{UniProt}_p{number}`). |
| 4 | `center_x` | Pocket center X coordinate. |
| 5 | `center_y` | Pocket center Y coordinate. |
| 6 | `center_z` | Pocket center Z coordinate. |
| 7 | `n_residues` | Number of residues in the pocket. |
| 8 | `score` | Heuristic pocket concavity score. |
| 9 | `pocket_residues` | Semicolon-separated residue list (`CHAIN:RESNUM`). |
| 10 | `source_file` | Original pocket CSV filename. |

---

## `protein_source_summaries.tsv`

This audit-only auxiliary table preserves the source-level protein summary rows that previously caused duplicate
`target_uniprot_id` values in `protein_database.tsv`. It is not a protein entity table and is intentionally kept outside `normalized_tables_v3.xlsx`.

| # | Column | Meaning |
|---|---|---|
| 1 | `target_uniprot_id` | UniProt accession. |
| 2 | `source_summary_id` | Stable local row ID for the archived source summary record. |
| 3 | `approved_symbol` | Gene symbol copied from the source record. |
| 4 | `ncbi_gene_id` | NCBI Gene ID copied from the source record. |
| 5 | `protein_name` | Protein name copied from the source record. |
| 6 | `source_databases` | Source database string from the pre-collapse protein record. |
| 7 | `unique_drug_count_in_source_record` | Original source-record drug count, retained for audit only. |
| 8 | `unique_drug_cid_sample_in_source_record` | Original source-record CID sample, retained for audit only. |
| 9 | `n_pockets_alphafold` | Pocket count from the source record. |
| 10 | `pocket_batch` | Pocket batch from the source record. |
| 11 | `has_matched_evidence` | Original source-record matched evidence flag. |
| 12 | `reviewed` | Original source-record UniProt reviewed flag. |
| 13 | `annotation_status` | Original source-record annotation status. |

---

## Key Design Decisions

| Decision | Rationale |
|---|---|
| Protein-first column order in binary table | Follows meeting decision: protein-centric design. |
| `*_other` columns removed from binary table | Same-protein/different-ligand context, not direct pair evidence. |
| STITCH cleaned to high-confidence only | 83,372 → 26,288 entries (31%). Eliminates unvalidated compound-protein noise. |
| AlphaFold columns suffixed `_alphafold` | Distinguishes computational predictions from experimental PDB data. |
| `compound_source` evidence-origin labels | Uses `curated_target_relation`, `bioassay_active_relation`, and `structure_ligand_context` to separate curated target relations, active bioassays, and structure ligand context. |
| CID-first molecule identity | 98.8% CID coverage. 985 biologics retain ChEMBL/DrugCentral ID. |
| Binary table stores 101,171 relationship/context rows | Unique drug-protein pairs plus curated structure ligand context rows are retained; source databases are merged with `;` where exact duplicates existed. |
| Protein table collapsed to one UniProt per row | `target_uniprot_id` is the primary key; `unique_drug_count`, `unique_drug_cid_sample`, and `source_databases` are derived from the binary table. |
| Source-level protein summaries archived separately | Pre-collapse source summary records are retained in `protein_source_summaries.tsv` for auditability. |
| Structure ligand context merged into binary table | BioLiP/sc-PDB/PDBbind `*_other` ligands with `keep_for_binding_site_context = 1` are merged as `compound_source = structure_ligand_context`; empty ligand labels and buffer/salt-only rows are excluded. |

## Rebuilding

```powershell
$env:PYTHONPATH = ""

# CID mapping (5 layers)
python bindingsite\map_all_to_cid.py

# 4-table split
python bindingsite\build_four_tables.py
```

`pocket_instances.tsv` is reused from Mintong v2:
```powershell
python memprot-drug-to-pro\scripts\build_pocket_subtables.py
```
