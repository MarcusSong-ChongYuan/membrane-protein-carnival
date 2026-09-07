# Normalized Tables Data Dictionary — v3 (CID Edition)

Date: 2026-07-06
Scope: 4-table normalized split from mempro1 v1.0 master table, with full CID conversion,
       STITCH cleaning, compound_source classification, and AlphaFold column suffix.

## Files

| File | Rows | Cols | Meaning |
|---|---|---|---|
| `drug_protein_binary_relationships.tsv` | 117,149 | 28 | One row per drug-protein pair. Only matched (same-drug) binding evidence retained. `*_other` columns removed. |
| `protein_database.tsv` | 1,189 | 27 | One row per protein. Pocket predictions, PDB co-crystal (other), cleaned STITCH, UniProt experimental sites. |
| `small_molecule_database.tsv` | 84,299 | 29 | One row per unique CID. Compound identity, source aggregation, activity summary, PubChem chemical properties, safety/toxicity. |
| `pocket_instances.tsv` | 6,815 | 13 | One row per predicted pocket. Geometry-based concavity detection output (Mintong v2, reused). |

---

## `drug_protein_binary_relationships.tsv`

### Drug Identity (1–4)
| Column | Meaning |
|---|---|
| `drug_id` | PubChem CID (98.8% coverage). Unmapped biologics keep ChEMBL/DrugCentral ID. |
| `drug_name` | Preferred drug/compound name from source table. |
| `compound_id_type` | `PubChem CID` for mapped compounds; original type for unmapped biologics. |
| `compound_source` | `self-screening` (ChEMBL curated + DrugCentral pharmacological targets) or `non-screening` (PubChem BioAssay HTS). |

### Protein Identity (5–7)
| Column | Meaning |
|---|---|
| `target_uniprot_id` | UniProt accession. |
| `approved_symbol` | HGNC-approved gene symbol. |
| `protein_name` | Human-readable protein name. |

### Evidence (8–14)
| Column | Meaning |
|---|---|
| `source_database` | Originating database(s) for this drug-protein pair. |
| `evidence_level` | Original evidence level: `high_confidence`, `medium_confidence`, `low_confidence`, `active_no_numeric_value`. |
| `relationship_confidence` | Normalized confidence tier derived from evidence_level: `high`, `medium`, `low`. |
| `assay_or_mechanism` | Assay type or mechanism-of-action description. |
| `activity_type` | Measured activity type: `IC50`, `Ki`, `Kd`, `EC50`, etc. |
| `activity_value_uM` | Activity value in µM. ≤10 µM considered active. |
| `clinical_or_approval_status` | Clinical development or FDA approval status where available. |

### Matched Binding Evidence Flags (15–19)
Boolean flags (`1`/`0`) for quick filtering. Full detail in the site columns below.
| Column | Meaning |
|---|---|
| `has_pdb_biolip_matched` | BioLiP co-crystal evidence for the exact same drug-protein pair. |
| `has_pdb_scpdb_matched` | sc-PDB co-crystal evidence for the exact same drug-protein pair. |
| `has_pdbbind_matched` | PDBbind affinity data for the exact same drug-protein pair. |
| `has_stitch_matched` | STITCH interaction evidence (combined_score ≥ 700) for the exact same drug-protein pair. |
| `has_exp_binding_site` | UniProt literature-curated experimental binding site for this protein. |

### Matched Binding Site Detail (20–24)
| Column | Meaning |
|---|---|
| `pdb_biolip_matched_sites` | Format: `LIGAND(CODE)@PDB:CHAIN=RESIDUES` |
| `pdb_scpdb_matched_sites` | Format: `LIGAND(CODE)@PDB:CHAIN=RESIDUES` |
| `pdbbind_matched_sites` | Format: `LIGAND(CODE)@PDB:Ki=3.2nM=CHAIN:RESIDUES` |
| `stitch_matched_compounds` | Format: `CHEMBLxxx(CIDxxx):confidence` or `CIDxxx:confidence` |

### UniProt Experimental Sites (25–28)
Protein-level literature-validated binding site annotations.
| Column | Meaning |
|---|---|
| `exp_binding_sites` | Residue-level binding site descriptions. |
| `exp_ligands` | Known ligands from experimental annotation. |
| `exp_evidence_quality` | ECO evidence codes (`high`: ECO:0000269). |
| `exp_pubmed` | PubMed IDs for binding site evidence. |
| `exp_pdb` | Cross-referenced PDB entries. |

---

## `protein_database.tsv`

### Protein Identity (1–5)
| Column | Meaning |
|---|---|
| `target_uniprot_id` | UniProt accession (primary key). |
| `approved_symbol` | Gene symbol. |
| `protein_name` | Full protein name. |
| `organism` | `Homo sapiens`. |
| `membrane_protein` | `true` for all entries. |

### AlphaFold Predicted Pockets (6–11)
All columns suffixed `_alphafold` to distinguish from experimental data.
| Column | Meaning |
|---|---|
| `n_pockets_alphafold` | Number of predicted pockets. |
| `pocket_ids_alphafold` | Semicolon-separated pocket IDs (e.g., `P12345_p1;P12345_p2`). |
| `pocket_centers_alphafold` | Pipe-separated 3D coordinates (e.g., `x,y,z\|x,y,z`). |
| `pocket_scores_alphafold` | Semicolon-separated heuristic pocket scores. |
| `pocket_residues_all_alphafold` | Semicolon-separated unique residues across all pockets. |
| `pocket_batch` | Batch label: `1`, `2`, or `1;2`. |

### PDB Co-crystal (Other Ligands) (12–14)
Other-ligand evidence: same protein, different compound. Not direct evidence for any specific drug-protein pair.
| Column | Meaning |
|---|---|
| `pdb_biolip_other` | Deduplicated BioLiP entries for all ligands on this protein. |
| `pdb_scpdb_other` | Deduplicated sc-PDB entries. |
| `pdbbind_other` | Deduplicated PDBbind entries with Ki/Kd/IC50. |

### STITCH (Cleaned) (15–17)
| Column | Meaning |
|---|---|
| `stitch_compounds_original_count` | Total STITCH `*_other` entries before cleaning. |
| `stitch_compounds_cleaned` | Semicolon-separated STITCH entries retained after cleaning. Only compounds whose (CID, UniProt) pair exists as `high_confidence` in the master table. |
| `stitch_compounds_cleaned_count` | Number of entries after cleaning. |

**Cleaning rule**: 83,372 entries → 26,288 retained (31%). Removed entries describe compounds without high-confidence evidence for this protein.

### UniProt Experimental Sites (18–22)
| Column | Meaning |
|---|---|
| `exp_binding_sites` | Deduplicated binding site descriptions. |
| `exp_ligands` | Deduplicated known ligands. |
| `exp_evidence_quality` | Evidence quality tiers. |
| `exp_pubmed` | Deduplicated PubMed IDs. |
| `exp_pdb` | Deduplicated PDB cross-references. |

### Summary Statistics (23–27)
| Column | Meaning |
|---|---|
| `unique_drug_count` | Number of unique CIDs linked to this protein. |
| `unique_drug_cid_sample` | Up to 50 sample CIDs (semicolon-separated). |
| `source_databases` | All source databases represented. |
| `evidence_levels` | All evidence levels present. |
| `has_matched_evidence` | `1` if any matched PDB/STITCH evidence exists for this protein. |

---

## `small_molecule_database.tsv`

### Identity (1–5)
| Column | Meaning |
|---|---|
| `drug_id` | PubChem CID (primary key for mapped compounds). |
| `drug_name` | Preferred name from source table. |
| `compound_id_type` | `PubChem CID` or original type for biologics. |
| `compound_source` | `self-screening` or `non-screening`. |
| `cid_mapping_method` | `mapped` (CID converted) or `unmapped_biologic` (no PubChem CID). |

### Source Evidence Summary (6–10)
| Column | Meaning |
|---|---|
| `source_databases` | Semicolon-separated source databases. |
| `source_row_count` | Number of master-table rows for this compound. |
| `unique_target_count` | Number of unique proteins this compound targets. |
| `evidence_levels` | Evidence levels present. |
| `clinical_statuses` | Clinical/approval statuses if any. |

### Activity Summary (11–13)
| Column | Meaning |
|---|---|
| `activity_types` | Activity measurement types present. |
| `activity_value_uM_min` | Minimum activity value in µM. |
| `activity_value_uM_median` | Median activity value in µM. |

### PubChem Chemical Properties (14–26)
From PubChem PUG REST. Populated for 80,364 CIDs (Mintong v2 cache, reused).
| Column | Meaning |
|---|---|
| `molecular_formula` | Molecular formula. |
| `molecular_weight` | Molecular weight (g/mol). |
| `canonical_smiles` | Canonical SMILES string. |
| `isomeric_smiles` | Isomeric SMILES string. |
| `inchikey` | InChIKey (standard). |
| `iupac_name` | IUPAC systematic name. |
| `xlogp` | Computed logP (octanol-water partition coefficient). |
| `tpsa` | Topological polar surface area (Å²). |
| `hbond_donor_count` | Hydrogen bond donor count. |
| `hbond_acceptor_count` | Hydrogen bond acceptor count. |
| `rotatable_bond_count` | Rotatable bond count. |
| `complexity` | Molecular complexity score. |
| `formal_charge` | Formal charge. |

### Safety/Toxicity (27–29)
From PubChem PUG-View. Populated for a limited subset (not a bulk endpoint).
| Column | Meaning |
|---|---|
| `ghs_hazard_classes` | GHS hazard class codes. |
| `ghs_signal_words` | GHS signal words (Warning / Danger). |
| `toxicity_summary` | Human-readable toxicity summary. |

---

## `pocket_instances.tsv` (Mintong v2, Reused)

| Column | Meaning |
|---|---|
| `pocket_batch` | Batch label: `1` or `2`. |
| `target_uniprot_id` | UniProt accession. |
| `pocket_id` | Unique pocket identifier (`{UniProt}_p{number}`). |
| `center_x`, `center_y`, `center_z` | 3D pocket center coordinates. |
| `n_atoms`, `n_residues` | Atom and residue counts. |
| `score` | Heuristic pocket score. |
| `pocket_residues` | Semicolon-separated residue list. |
| `source_file` | Original pocket CSV. |
| `source_row_number` | 1-based row number in the source file. |
| `residue_count_from_list` | Recomputed residue count (QC check). |

---

## Key Design Decisions (v3 vs v2)

| Decision | Rationale |
|---|---|
| `*_other` columns removed from binary relationships | These describe same-protein/different-ligand context, not direct evidence for the row's drug-protein pair. Archived separately by Mintong v2. |
| STITCH cleaned to high-confidence only | 83,372 → 26,288 entries. Eliminates noise from unvalidated compound-protein interactions. |
| AlphaFold columns suffixed `_alphafold` | Distinguishes computational predictions from experimental PDB data. |
| `compound_source` added | Separates curated pharmacological targets (self-screening) from HTS hits (non-screening). |
| CID-first molecule identity | 98.8% of drugs mapped to PubChem CID for cross-database interoperability. |
| Binary table preserves 117K rows | Multi-source cross-validation retained; duplicates are different evidence sources for the same pair. |

## Rebuilding

```powershell
$env:PYTHONPATH = ""

# CID mapping
python bindingsite\map_all_to_cid.py

# 4-table split
python bindingsite\build_four_tables.py
```

The `pocket_instances.tsv` file is reused from Mintong v2 and was built by:
```powershell
python memprot-drug-to-pro\scripts\build_pocket_subtables.py
```
