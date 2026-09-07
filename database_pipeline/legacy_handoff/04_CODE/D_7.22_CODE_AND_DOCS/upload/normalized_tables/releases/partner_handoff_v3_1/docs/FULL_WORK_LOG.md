# Full Work Log — mempro1 v3 Normalized Tables

Date: 2026-06-15 to 2026-07-08
Team: Marcus Song, Mintong Yu

---

## Phase 1: Master Table Assembly (June 15-24)

### 1.1 Drug → Membrane Protein (Mintong, Step 1)
- FDA/clinical drugs → curated MoA targets → UniProt membrane verification
- Output: 511 unique membrane proteins, 2,520 drug-target pairs
- Evidence levels: A1 (FDA+single protein) 2234, A2 (FDA+complex) 218, B1 (clinical+single) 60, B2 (clinical+complex) 8
- Filter: excluded docking/bioassay-only/pathway/family-level/fuzzy targets

### 1.2 Protein → Drug Reverse Search (Mintong, Step 2)
- ChEMBL drug_mechanism, DrugCentral interactions, PubChem BioAssay bounded by UniProt accessions
- 1,703 total target accessions
- Output: 27,039 rows at checkpoint

### 1.3 AlphaFold Pocket Prediction (Marcus, Step 3a)
- 1,757 membrane proteins → AlphaFold DB v6 structure download
- Geometry-based concavity detection: CA atoms → NeighborSearch → surface residues (<40th percentile) → concavity (density > 1.15× mean) → spatial clustering (10 Å radius)
- Output: 1,690 proteins with pockets, 6,815 total pockets (Batch1: 503/2,177 + Batch2: 1,187/4,638)

### 1.4 PDB Co-crystal Evidence (Marcus, Step 3b)
- BioLiP: 389 PDB ligands, 299 proteins, 56,996 binding sites
- sc-PDB: 17,594 binding sites, 187 proteins, 44,256 rows
- PDBbind v2020: 19,037 complexes, 269 proteins, 3,805 entries with Ki/Kd/IC50
- Ligand 3-letter codes → PDBe API → InChIKey → UniChem → ChEMBL ID mapping
- 3,400 codes → 1,913 mapped to ChEMBL

### 1.5 STITCH + UniProt (Marcus, Step 3c)
- STITCH v5: 9606.protein_chemical.links (73MB), combined_score ≥ 700
- 746 proteins, 83,374 rows with CID→ChEMBL mapping (1,365/3,274)
- UniProt experimental sites: 858 proteins, 5,800 annotations, 69,767 rows covered

### 1.6 PubChem Completion (Marcus)
- 17 accessions at checkpoint → 1,154 production complete
- PUG REST API → AID lists → curl.exe concise CSV → filter active + IC50/Ki/Kd/EC50 ≤ 10µM
- 27,972 AIDs, 651,820 raw records → 105,817 final rows
- ~500 accessions genuinely have no PubChem data (HTTP 404)

### 1.7 Final Merge
- Dedup, ChEMBL ID normalization, matched/other split
- **Master table: 117,149 rows × 33 columns, 1,189 proteins, 85,584 drugs**

---

## Phase 2: Mintong's v2 Normalized Split (June 29 - July 1)

Mintong's design decisions:
- `*_other` columns removed from binary table (same protein, different ligand — not direct evidence)
- Pocket granularity chosen as the fundamental unit (not drug-protein pair)
- External API enrichments: IntAct PSICQUIC (1,587 proteins), BioGRID REST (1,599), UniProt REST (1,690), PubChem PUG REST (80,364 CIDs)

Outputs (at `mempro1-codex-pocket-subtables/upload/`):
- `pocket_instances.tsv`: 6,815 rows × 13 cols
- `pocket_protein_database.tsv`: 1,690 rows × 64 cols (with PPI networks)
- `small_molecule_database.tsv`: 85,447 rows × 47 cols (with PubChem properties + safety)
- `drug_target_binary_relationships.tsv.gz`: 96,211 rows (deduplicated)
- `removed_other_small_molecule_interactions.tsv`: 18,306 other-ligand records archived
- Documentation: DATA_DICTIONARY.md, MANIFEST.json (SHA-256), PROCESS_LOG × 3

---

## Phase 3: Marcus v3 — CID Edition (July 5-8)

### 3.1 Full CID Conversion

**Goal**: Convert ALL drug IDs to PubChem CID format.

**5-layer mapping strategy**:

| Layer | Method | Mapped |
|---|---|---|
| 1 | Direct PubChem CID already in source | 80,364 |
| 2 | Embedded CID extracted from name (e.g., CHEMBL66654(CID3025961)) | +292 |
| 3 | Mintong's small_molecule_database mapping (source_ids → CID) | +137 |
| 4 | ChEMBL SQLite → InChIKey → PubChem (parallel 8 workers) | +3,002 |
| 5 | PubChem name lookup (parallel 2 workers) | +1 |
| **Total** | | **84,584 / 85,584 (98.8%)** |

**Unmapped**: ~1,000 biologics (monoclonal antibodies, fusion proteins, peptide hormones) — PubChem does not index biologics. Kept ChEMBL IDs with `unmapped_biologic` flag.

### 3.2 drug_name Standardization

**Problem**: 69.4% of rows had "PubChem CID xxxxx" as drug_name placeholder. 63% were IUPAC chemical formulas, not usable drug names.

**Fix pipeline**:
| Step | Method | Result |
|---|---|---|
| 1 | Load 85,583 source names from master table | Baseline |
| 2 | Batch PubChem Title query (63,696 CIDs, 200/batch, 8 workers) | 63,296 titles found |
| 3 | Replace "PubChem CID xxx" placeholders | 82,063 → 509 remaining |
| 4 | Batch PubChem synonyms query (53,127 CIDs) | 52,753 found better names |
| 5 | Unify: one name per drug_id (pick most frequent) | 0 multi-name remnants |
| 6 | Classify into drug_name_type (see below) | 6 types |

### 3.3 drug_name_type Classification

| Type | Count | Definition |
|---|---|---|
| `approved_drug` | 6,470 | FDA `approved_or_listed` confirmed |
| `literature_drug` | 3,471 | Self-screening source, no FDA record |
| `dev_code` | 1,814 | Research code (e.g., BMS-794833) |
| `named_compound` | 24,381 | PubChem HTS hit, has readable name |
| `chemical` | 59,098 | IUPAC/systematic name only |
| `unknown` | 0 | (was 396, all resolved by re-querying PubChem) |

`biologic` type abolished — biologics merged into `approved_drug`.

**Bug fix**: 58,000 compounds originally misclassified as `named_compound` or `dev_code` when they were actually chemical/IUPAC names. Reclassified properly using IUPAC pattern detection.

### 3.4 Multi-UniProt Explosion

**Problem**: 355 DrugCentral rows had `|`-separated multi-UniProt values (e.g., `P03372|Q92731` for fosfestrol targeting both ESR1 and ESR2). Violated 1NF.

**Fix**: Exploded into 1,106 additional rows. Added `target_group_id` column (e.g., `DC_GROUP_1`) for traceability.

### 3.5 ENSG → UniProt Mapping

**Problem**: 44 proteins used Ensembl Gene IDs (ENSG...) instead of UniProt accessions.

**Fix**: 
- 42 mapped to UniProt via REST API (`xref:ensembl-{ENSG}`)
- 2 were miRNAs (MIR155, MIR21) — deleted from all tables
- 9 new UniProt IDs emerged from multi-UniProt explosion — added to protein table

### 3.6 Activity Data Gap-filling

**Problem**: 12,434 rows had activity_type but no activity_value_uM.

**Sources used**:
| Source | Method | Filled |
|---|---|---|
| DrugCentral raw file (`drug.target.interaction.tsv`) | Match by (drug_name, UniProt) | 2,167 |
| DrugCentral gene-name fallback | Match by (drug_name, gene_symbol) | 12 |
| ChEMBL SQLite (`chembl_37.db`) | molregno → assay_id → tid join | 15 |
| PubChem assaysummary API | Per-CID, match by NCBI Gene ID | 132 (IC50/Ki/Kd/EC50 only) |

**Rule**: Only IC50, Ki, Kd, EC50 accepted. "Potency" and other vague labels rejected.

### 3.7 Column Restructuring

**Binary Relationships**:
- Removed `evidence_level`, `relationship_confidence` (always `high_confidence` — no discriminatory power)
- Removed `drug_name_standardized` (CID already unifies identity)
- Moved `compound_id_type` next to `drug_id`
- Reorganized binding site columns into paired structure: `has_X_matched` + `X_matched_sites` together
- Merged 22,989 duplicate rows → 95,234 unique (drug, protein) pairs
- `source_database` merged with `;` separator for duplicates
- `assay_or_mechanism` keeps most informative version (ChEMBL > DrugCentral > PubChem)

**Protein Database**:
- Removed `evidence_levels` (always `high_confidence`)
- Added 7 UniProt enrichment columns (see 3.9)

**Small Molecule Database**:
- Removed `drug_name_standardized`
- Moved `compound_id_type` next to `drug_id`
- Added 4 clinical columns (see 3.10)

**Pocket Instances**:
- Removed `n_atoms` (always 0 — hardcoded, not recorded)
- Removed `source_row_number` (debugging artifact)
- Removed `residue_count_from_list` (QC check, already verified)

### 3.8 Orphan Cleanup

- 28 rows with empty drug_id → deleted
- 4 rows targeting deleted miRNAs → deleted
- 3 molecules only linked to deleted miRNA rows → deleted
- Tables verified: binary ↔ protein (1,169 == 1,169), binary ↔ molecule (84,294 == 84,294)

### 3.9 Protein Enrichment (UniProt REST)

**Source**: `rest.uniprot.org/uniprotkb/search` — batch query (50/batch, 0.2s delay)

**New columns**:
| Column | Fill Rate | Source Field |
|---|---|---|
| `protein_function` | 51% | `cc_function` |
| `subcellular_location` | 42% | `cc_subcellular_location` |
| `transmembrane_count` | 27% | `ft_transmem` (count of TRANSMEM features) |
| `disease_association` | 27% | `cc_disease` |
| `pdb_structures` | 42% | `xref_pdb` (count of PDB entries) |
| `go_molecular_function` | ~40% | `go_f` |
| `go_cellular_component` | ~40% | `go_c` |

### 3.10 Small Molecule Enrichment (PubChem PUG-View)

**Source**: `pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{CID}/JSON` — 4,598 priority compounds (8 parallel workers)

**New columns**:
| Column | Filled | Source Section |
|---|---|---|
| `drug_indication` | 244 | Drug and Medication Information → Drug Indication |
| `pharmacodynamics` | 187 | Pharmacology → Pharmacodynamics |
| `drug_classes` | 165 | Drug and Medication Information → Drug Classes |
| `livertox` | 96 | Drug and Medication Information → LiverTox Summary |

### 3.11 Documentation

**New files created**:
- `DATA_DICTIONARY.md`: column-by-column definitions for all 4 tables
- `MANIFEST.json`: SHA-256 hashes of all output files
- `PROCESS_LOG.md`: step-by-step processing log
- `drug_name_normalization_report.tsv`: 84,295-row audit trail (original → standardized)
- `FULL_WORK_LOG.md`: this file

**Cache files** (for reproducibility):
- `bindingsite/cache/complete_cid_map.json`: 84,584 CID mappings
- `bindingsite/cache/cid_name_cache.json`: 63,696 PubChem Titles
- `bindingsite/cache/pubchem_synonyms.json`: 53,127 PubChem synonyms
- `bindingsite/cache/pubchem_assay_cache.json`: 2,581 PubChem assay summaries
- `bindingsite/cache/pubchem_clinical.json`: 4,598 PubChem clinical data
- `bindingsite/cache/ensg_to_uniprot.json`: 44 ENSG → UniProt mappings
- `bindingsite/cache/ncbi_gene_cache.json`: 1,773 NCBI Gene IDs
- `bindingsite/cache/uniprot_enrich.json`: 1,773 UniProt annotations
- `bindingsite/cache/uniprot_annotations.json`: 429 UniProt gene names

**Build scripts**:
- `bindingsite/map_all_to_cid.py`: CID conversion pipeline
- `bindingsite/build_four_tables.py`: 4-table split
- `bindingsite/normalize_drug_names.py`: drug name standardization
- `bindingsite/fill_chembl_activity.py`: ChEMBL SQLite activity lookup
- `bindingsite/fill_chembl_from_pubchem.py`: PubChem assay summary lookup
- `bindingsite/fill_protein_annotations.py`: UniProt annotation queries
- `bindingsite/tsv_to_xlsx.py`: Excel generation (4 sheets, formatted)

---

## Final Deliverables

```
C:\github-repos\upload\normalized_tables\
├── drug_protein_binary_relationships.tsv   101,802 rows × 29 cols   55 MB
├── protein_database.tsv                     1,773 rows × 26 cols   5.4 MB
├── protein_source_summaries.tsv            811 rows × 13 cols     0.2 MB
├── small_molecule_database.tsv             87,611 rows × 33 cols   37 MB
├── pocket_instances.tsv                     6,815 rows × 10 cols  1.4 MB
├── normalized_tables_v3.xlsx               4-sheet workbook        31 MB
├── drug_name_normalization_report.tsv      84,295 row audit trail
├── DATA_DICTIONARY.md                      Column definitions
├── MANIFEST.json                            SHA-256 hashes
├── PROCESS_LOG.md                           Processing log
└── FULL_WORK_LOG.md                         This file
```

## Key Stats

| Metric | Value |
|---|---|
| Binary relationships | 101,171 relationship/context rows |
| Unique proteins | 1,773 (1,134 with binary drug relationships + 639 pocket-only/no-binary retained) |
| Unique molecules | 84,295 |
| CID coverage | 98.8% |
| Binding site coverage | BioLiP 0.04%, sc-PDB 0.01%, PDBbind 0%, STITCH 0.8%, UniProt 59% |
| Matched evidence | 948 rows (0.8% — data nature) |
| Self-screening vs non-screening | 9.7% vs 90.3% |
| STITCH cleaned | 83,372 → 26,288 (31% retained) |
| Protein function annotated | 51% |
| Clinical data | 244 indications, 187 MOAs |

## Post-log Cleanup on 2026-07-09

- Collapsed `protein_database.tsv` from 2,242 source-summary rows to 1,773 unique UniProt rows.
- Recomputed `unique_drug_count`, `unique_drug_cid_sample`, and `source_databases` from `drug_protein_binary_relationships.tsv`.
- Preserved the 811 pre-collapse duplicate/source-summary records in `protein_source_summaries.tsv`.
- Rebuilt `normalized_tables_v3.xlsx` as a four-sheet result workbook; kept `protein_source_summaries.tsv` outside the workbook as audit-only provenance.
- Removed low-information/provenance columns from `protein_database.tsv`: `organism`, `membrane_protein`, and protein-level `pocket_batch`. Pocket batch remains in `pocket_instances.tsv`.
- Curated BioLiP/sc-PDB/PDBbind `*_other` ligands into structure ligand context review files, mapped reliable ligands to PubChem CID, and left unresolved non-empty ligands as stable `UNMAPPED_PDB_LIGAND:*` IDs.
- Merged 5,937 structure ligand context rows into `drug_protein_binary_relationships.tsv` as `compound_source = structure_context`: 5,809 mapped PubChem CID rows and 128 unmapped structure-ligand rows.
- Synchronized `small_molecule_database.tsv` with 3,257 added entities: 3,158 PubChem CID rows and 99 unmapped structure-ligand rows.
- Recomputed protein and small-molecule summary fields from the merged binary table and rebuilt `normalized_tables_v3.xlsx`.
- Removed `pdb_biolip_other`, `pdb_scpdb_other`, and `pdbbind_other` from `protein_database.tsv`; structure ligand context is now represented in the binary table and retained in the review TSV for audit.
- Standardized STITCH evidence strings to CID-only format in binary/protein tables: `CHEMBLxxx(CIDyyy):score` -> `CIDyyy:score`.
- Removed AlphaFold pocket detail string columns from `protein_database.tsv`: `pocket_ids_alphafold`, `pocket_centers_alphafold`, `pocket_scores_alphafold`, and `pocket_residues_all_alphafold`.
- Kept `n_pockets_alphafold` as the protein-level summary; complete AlphaFold pocket details remain in `pocket_instances.tsv`.
- Removed `annotation_status` from `protein_database.tsv`; `reviewed` remains as the standard UniProt review-status column.
- Removed `cid_mapping_method` from `small_molecule_database.tsv`; it is internal CID-mapping provenance rather than a final result field.
- Removed `evidence_levels` from `protein_database.tsv` and `small_molecule_database.tsv`; it was too coarse for the current entity-summary tables and duplicated source/context information.
- Renamed current `compound_source` values to `curated_target_relation`, `bioassay_active_relation`, and `structure_ligand_context`.
- Replaced small-molecule `clinical_statuses` with `compound_evidence_status` and split `drug_name_type` into `compound_name_category` plus `compound_biological_role` in both the binary and small-molecule tables.
- Added row-level `ligand_interpretation_class` to the binary table and compound-level `ligand_interpretation_classes` plus `best_ligand_interpretation_class` to the small-molecule table.
- Filled high-missing safety/clinical text fields with `not_found_in_current_sources` to distinguish missing source coverage from biological absence.

Updated final row counts after structure ligand context merge:

| Table | Rows | Columns |
|---|---:|---:|
| `drug_protein_binary_relationships.tsv` | 101,802 | 31 |
| `protein_database.tsv` | 1,773 | 25 |
| `small_molecule_database.tsv` | 87,611 | 35 |
| `pocket_instances.tsv` | 6,815 | 10 |

- Split explicit UniProt experimental binding-site ligands into binary rows with `source_database = UniProt` and `compound_source = uniprot_binding_site_annotation`; generic/non-specific terms were excluded to `uniprot_exp_ligand_excluded_terms.tsv`.

## Known Limitations

1. 985 biologics (mAbs, fusion proteins) have no PubChem CID — kept ChEMBL IDs
2. DrugCentral PostgreSQL dump not downloaded — 643 activity values unfilled
3. ChEMBL MoA records are qualitative — 8,373 rows have no IC50 by design
4. PDB co-crystal matched rate <1% — research tool ligands ≠ marketed drugs
5. AlphaFold pockets are computational predictions, not experimentally validated
6. IntAct/BioGRID PPI data not included (meeting decision: future work)
7. Molecular docking not performed
