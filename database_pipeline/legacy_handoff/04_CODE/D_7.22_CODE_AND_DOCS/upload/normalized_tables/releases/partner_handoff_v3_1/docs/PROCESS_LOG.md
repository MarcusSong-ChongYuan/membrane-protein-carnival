# Normalized Tables Process Log — v3 (CID Edition)

Generated at UTC: 2026-07-06

## Purpose

Build 4 normalized tables from the v1.0 master table with full CID conversion, STITCH cleaning,
compound_source classification, and AlphaFold column suffix renaming. Follows Mintong v2 design
philosophy (remove *_other from binary table, enrich protein/molecule tables) with meeting-driven
enhancements.

## Inputs

| File | Source | Rows |
|---|---|---|
| `upload/drug_target_pocket.tsv.gz` | Marcus v1.0 master table | 117,149 |
| `upload/binding_pockets_batch1.csv` | AlphaFold pocket prediction batch 1 | 2,177 |
| `upload/binding_pockets_batch2.csv` | AlphaFold pocket prediction batch 2 | 4,638 |
| `upload/pocket_subtables/small_molecule_database.tsv` | Mintong v2 (PubChem properties) | 85,447 |
| `upload/pocket_subtables/pocket_instances.tsv` | Mintong v2 (pocket instances, reused as-is) | 6,815 |
| `bindingsite/cache/complete_cid_map.json` | CID mapping (5-layer, see below) | 84,584 mappings |

## Outputs

| File | Rows | Cols | Size |
|---|---|---|---|
| `upload/normalized_tables/drug_protein_binary_relationships.tsv` | 101,802 | 31 | 60.5 MB |
| `upload/normalized_tables/protein_database.tsv` | 1,773 | 25 | 2.8 MB |
| `upload/normalized_tables/small_molecule_database.tsv` | 87,611 | 35 | 70.2 MB |
| `upload/normalized_tables/pocket_instances.tsv` | 6,815 | 10 | 1.4 MB |
| `upload/normalized_tables/protein_source_summaries.tsv` | 811 | 13 | 0.2 MB |
| `upload/normalized_tables/DATA_DICTIONARY.md` | — | — | — |
| `upload/normalized_tables/MANIFEST.json` | — | — | — |
| `upload/normalized_tables/PROCESS_LOG.md` | — | — | — |

## Row Counts

| Metric | Count |
|---|---|
| Source master table rows | 117,149 |
| Binary relationship rows | 101,171 |
| Unique proteins (binary table) | 1,134 |
| Unique proteins (protein database) | 1,773 |
| Pocket-only/no-binary proteins retained | 639 |
| Unique molecules (small molecule database) | 87,552 |
| Unique molecule CIDs | 86,472 |
| Unique molecule CIDs with PubChem properties | 80,364 |
| CID-mapped molecules | 86,472 (98.8%) |
| Unmapped biologics | 981 (1.1%) |
| Unmapped structure ligands | 99 (0.1%) |
| Pocket instances (reused) | 6,815 |

## Processing Steps

### Phase 1: CID Mapping

1. Load Mintong v2 `small_molecule_database.tsv` to extract 80,501 source-ID → CID mappings.
2. Collect 85,584 unique drug IDs from master table.
3. Apply 5-layer mapping:
   - **Layer 1** (direct CID): 80,364 IDs already PubChem CID → use directly.
   - **Layer 2** (embedded CID): 292 IDs with embedded CID in name (e.g., `CHEMBL66654(CID3025961)`) → extract CID.
   - **Layer 3** (Mintong mapping): 137 IDs mapped via Mintong's table.
   - **Layer 4** (ChEMBL SQLite → InChIKey → PubChem): 2,165 unique InChIKeys queried via PubChem PUG REST with 8 parallel workers. 3,002 ChEMBL IDs mapped.
   - **Layer 5** (PubChem name lookup): 1,211 drug names queried with 2 parallel workers. 1 additional ID mapped.
4. Final mapping: 84,584 / 85,584 (98.8%). 1,000 unmapped are biologics (monoclonal antibodies, fusion proteins, peptide hormones) without PubChem CIDs.

Script: `bindingsite/map_all_to_cid.py`

### Phase 2: Master Table Processing

5. Read `drug_target_pocket.tsv.gz` with gzip streaming.
6. Substitute all drug IDs with CIDs using the complete mapping.
7. Classify `compound_source`: `self-screening` for ChEMBL + DrugCentral rows, `non-screening` for PubChem BioAssay rows.
8. Derive `relationship_confidence` from `evidence_level`: high / medium / low.
9. Create boolean `has_*_matched` flags for quick filtering.
10. Aggregate protein-level data: pocket info, PDB other, STITCH other, UniProt experimental sites.
11. Aggregate molecule-level data: source databases, target counts, activity ranges, evidence levels.

### Phase 3: STITCH Cleaning

12. Identify high-confidence pairs: 94,936 (CID, UniProt) pairs with evidence_level = high_confidence.
13. For each protein, parse `stitch_other` entries and extract CID.
14. Retain only entries whose (CID, UniProt) pair exists in the high-confidence set.
15. Result: 83,372 → 26,288 entries retained (31%). 746 → 147 proteins with cleaned STITCH data.

### Phase 4: Output

16. Merge Mintong v2 PubChem properties (chemical + safety/toxicity) into molecule table by CID.
17. Rename AlphaFold pocket columns with `_alphafold` suffix.
18. Write 3 TSV files with UTF-8 encoding, tab delimiter.
19. Reuse Mintong v2 `pocket_instances.tsv` unchanged.

Script: `bindingsite/build_four_tables.py`

### Phase 5: Documentation

20. Compute SHA-256 hashes for all output files.
21. Write DATA_DICTIONARY.md, PROCESS_LOG.md, MANIFEST.json.

## CID Mapping Detail

| Layer | Method | Source | Mapped | Cumulative |
|---|---|---|---|---|
| 1 | Direct PubChem CID | `compound_id_type = PubChem CID` | 80,364 | 80,364 |
| 2 | Embedded CID | Regex `(CID\d+)` in drug_id | +292 | 80,656 |
| 3 | Mintong v2 mapping | `small_molecule_database.tsv` source_ids | +137 | 80,793 |
| 4 | ChEMBL SQLite → InChIKey → PubChem | `chembl_37.db` molecule_dictionary JOIN compound_structures | +3,002 | 83,795 |
| 5 | PubChem name lookup | Drug names for remaining unmapped | +1 | 83,796 |

PubChem queries: 3,418 total (2,165 InChIKey + 1,253 name), 3,208 mapped (93.9% query success rate).
Layer 4 used 8 parallel workers (ThreadPoolExecutor). Layer 5 used 2 workers to avoid rate limiting.

## STITCH Cleaning Detail

| Metric | Before | After | Reduction |
|---|---|---|---|
| Proteins with STITCH other | 746 | 147 | — |
| Total STITCH entries | 83,372 | 26,288 | 68.5% |
| High-confidence pair filter | — | (CID, UniProt) ∈ high_conf_pairs | — |

Rule: Compound must have at least one high_confidence relationship with the target protein.

## compound_source Classification

| Source | Rows | Criteria |
|---|---|---|
| self-screening | 11,332 (9.7%) | source_database contains ChEMBL or DrugCentral |
| non-screening | 105,817 (90.3%) | source_database contains PubChem |

## External Source Status

### PubChem PUG REST
- Status: `ok`
- InChIKey queries: 2,165 (8 parallel workers, ~4.4 queries/sec, 761s)
- Name queries: 1,253 (2 parallel workers, ~5.0 queries/sec, 244s)
- Rate limiting observed: switched from 8 to 2 workers for name queries
- Cache files: `bindingsite/cache/inchikey_to_cid_v2.json`, `bindingsite/cache/name_to_cid_v2.json`

### ChEMBL SQLite
- Status: `ok`
- Database: `data/raw/chembl/chembl_37.db`
- Tables used: `molecule_dictionary`, `compound_structures`
- 2,165 InChIKeys retrieved for 3,166 ChEMBL IDs (1,001 ChEMBL IDs had no InChIKey)

### Mintong v2 Small Molecule Properties
- Status: `reused`
- 80,364 CIDs with PubChem chemical properties
- 1,495 CIDs with limited safety/toxicity data
- Cache reused: `provenance/pubchem_property_cache.tsv`, `provenance/pubchem_safety_toxicity_cache.tsv`

## Changes From Mintong v2

| Dimension | Mintong v2 | v3 |
|---|---|---|
| Drug ID format | Mixed CID + ChEMBL + DrugCentral | 98.8% PubChem CID |
| STITCH in protein table | 83,372 entries (all) | 26,288 entries (high-confidence only) |
| compound_source | None | self-screening / non-screening |
| AlphaFold column names | `pocket_ids`, `pocket_scores`... | `pocket_ids_alphafold`, `pocket_scores_alphafold`... |
| Binary table rows | 96,211 (deduplicated) | 117,149 (preserves multi-source) |
| Protein table rows | 1,690 (incl. pocket-only) | 1,189 (only master-table proteins) |
| has_matched_* flags | None | 5 boolean flags |
| relationship_confidence | None | high / medium / low |
| PPI networks | IntAct + BioGRID | Not included (meeting decision: future work) |

## Interpretation Notes

- `pocket_instances.tsv` is computational output from AlphaFold predictions, not experimental binding proof.
- `protein_database.tsv` is protein-level and has exactly one row per `target_uniprot_id`.
- `unique_drug_count`, `unique_drug_cid_sample`, and `source_databases` in `protein_database.tsv` are recomputed from `drug_protein_binary_relationships.tsv`.
- Dataset-level constants (`organism = Homo sapiens`, `membrane_protein = true`) were removed from `protein_database.tsv`; `pocket_batch` remains in `pocket_instances.tsv`.
- Historical source-level protein summary rows are retained in `protein_source_summaries.tsv` for auditability.
- The source master table remains authoritative for drug-target evidence; these files are derived normalized views.
- STITCH cleaning reduces noise but may remove genuine low-confidence interactions. Original STITCH data is available in Mintong v2.
- 985 biologics (1.2%) lack PubChem CIDs because PubChem does not index biologics.

## Post-log Structure Ligand Context Merge on 2026-07-09

- Extracted BioLiP/sc-PDB/PDBbind `*_other` same-protein ligand context into a curated review table before merging.
- Mapped structure ligands to PubChem CID where reliable identifiers were available: local small-molecule table, ChEMBL/PubChem name lookup, RCSB CCD PubChem xref, and UniChem InChIKey cross-reference.
- Merged only rows with `keep_for_binding_site_context = 1` into `drug_protein_binary_relationships.tsv`.
- Excluded 309 empty ligand labels `()` and 1 buffer/salt-only row from the main binary table.
- Added 5,937 `structure_context` binary rows: 5,809 with PubChem CID and 128 with stable `UNMAPPED_PDB_LIGAND:*` IDs.
- Updated `small_molecule_database.tsv` with 3,257 new entities: 3,158 PubChem CID rows and 99 unmapped structure-ligand rows.
- Recomputed protein-level `unique_drug_count`, `unique_drug_cid_sample`, `source_databases`, and `has_matched_evidence` from the merged binary table.
- Rebuilt `normalized_tables_v3.xlsx` as a four-sheet workbook.
- Removed `pdb_biolip_other`, `pdb_scpdb_other`, and `pdbbind_other` from `protein_database.tsv`; the structured result now lives in the binary table and the review TSV.
- Standardized STITCH fields to CID-only evidence strings: `CHEMBLxxx(CIDyyy):score` -> `CIDyyy:score` in both binary and protein tables.
- Removed AlphaFold pocket detail string columns from `protein_database.tsv`: `pocket_ids_alphafold`, `pocket_centers_alphafold`, `pocket_scores_alphafold`, and `pocket_residues_all_alphafold`.
- Kept only `n_pockets_alphafold` in the protein entity table; complete pocket instances remain authoritative in `pocket_instances.tsv`.
- Removed `annotation_status` from `protein_database.tsv`; it only duplicated a low-information internal `ok` status, while UniProt's standard `reviewed` field is retained.
- Removed `cid_mapping_method` from `small_molecule_database.tsv`; mapping method is internal provenance and remains documented in process/audit artifacts rather than the final result table.
- Removed `evidence_levels` from `protein_database.tsv` and `small_molecule_database.tsv`; it duplicated broad source/category information and did not provide useful evidence grading for these entity summary tables.
- Renamed current `compound_source` values to `curated_target_relation`, `bioassay_active_relation`, and `structure_ligand_context`.
- Replaced small-molecule `clinical_statuses` with `compound_evidence_status`: `approved_drug`, `clinical_trial_compound`, or `binding_evidence_only`.
- Split `drug_name_type` into `compound_name_category` and `compound_biological_role` in both the binary and small-molecule tables.
- Added row-level `ligand_interpretation_class` to the binary table and compound-level `ligand_interpretation_classes` plus `best_ligand_interpretation_class` to the small-molecule table.
- Filled high-missing safety/clinical text fields with `not_found_in_current_sources` when no value was captured.

Current row counts after this merge:

| Table | Rows | Columns |
|---|---:|---:|
| `drug_protein_binary_relationships.tsv` | 101,802 | 31 |
| `protein_database.tsv` | 1,773 | 25 |
| `small_molecule_database.tsv` | 87,611 | 35 |
| `pocket_instances.tsv` | 6,815 | 10 |

- Split explicit UniProt experimental binding-site ligands into binary rows with `source_database = UniProt` and `compound_source = uniprot_binding_site_annotation`; generic/non-specific terms were excluded to `uniprot_exp_ligand_excluded_terms.tsv`.

## Rebuild Command

```powershell
$env:PYTHONPATH = ""
python bindingsite\map_all_to_cid.py
python bindingsite\build_four_tables.py
```
