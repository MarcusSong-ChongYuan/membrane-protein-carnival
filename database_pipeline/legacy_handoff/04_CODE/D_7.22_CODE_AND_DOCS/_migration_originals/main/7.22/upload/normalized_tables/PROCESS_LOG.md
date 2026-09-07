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
| `upload/normalized_tables/drug_protein_binary_relationships.tsv` | 117,149 | 28 | 54.2 MB |
| `upload/normalized_tables/protein_database.tsv` | 1,189 | 27 | 17.6 MB |
| `upload/normalized_tables/small_molecule_database.tsv` | 84,299 | 29 | 37.0 MB |
| `upload/normalized_tables/pocket_instances.tsv` | 6,815 | 13 | 1.4 MB |
| `upload/normalized_tables/DATA_DICTIONARY.md` | — | — | — |
| `upload/normalized_tables/MANIFEST.json` | — | — | — |
| `upload/normalized_tables/PROCESS_LOG.md` | — | — | — |

## Row Counts

| Metric | Count |
|---|---|
| Source master table rows | 117,149 |
| Binary relationship rows | 117,149 |
| Unique proteins (binary table) | 1,189 |
| Unique proteins (protein database) | 1,189 |
| Unique molecule CIDs | 84,299 |
| Unique molecule CIDs with PubChem properties | 80,364 |
| CID-mapped molecules | 83,314 (98.8%) |
| Unmapped biologics | 985 (1.2%) |
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
- `protein_database.tsv` is protein-level; interaction counts are source-specific summaries.
- The source master table remains authoritative for drug-target evidence; these files are derived normalized views.
- STITCH cleaning reduces noise but may remove genuine low-confidence interactions. Original STITCH data is available in Mintong v2.
- 985 biologics (1.2%) lack PubChem CIDs because PubChem does not index biologics.

## Rebuild Command

```powershell
$env:PYTHONPATH = ""
python bindingsite\map_all_to_cid.py
python bindingsite\build_four_tables.py
```
