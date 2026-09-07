# Partner Handoff — normalized_tables v3.1

## Read this first

`v3_1/normalized_tables_v3_1.xlsx` and the four TSV files beside it are the current handoff deliverables. They are **v3.1**, created from the original v3 release in `baseline_v3/`; do not treat the baseline folder as the latest result.

The four-table layout is intentionally retained:

1. `1_Binary_Relationships` — compound–protein relationship/context records.
2. `2_Protein_Database` — one row per UniProt accession.
3. `3_Small_Molecules` — one row per compound identifier.
4. `4_Pocket_Instances` — one row per AlphaFold-derived predicted pocket.

The original source snapshots and several historical scripts referenced by older logs are not available in this workspace. This package therefore supports continuation from v3/v3.1, not full reconstruction from raw ChEMBL, DrugCentral, PubChem, BioLiP, sc-PDB, PDBbind, STITCH, and AlphaFold inputs.

## Package layout

```text
partner_handoff_v3_1/
├── PARTNER_HANDOFF.md              # this guide
├── HANDOFF_MANIFEST.json           # package inventory and SHA-256 hashes
├── baseline_v3/                    # frozen pre-v3.1 input tables/workbook
├── v3_1/                           # current deliverables and QC reports
├── docs/                           # original workflow and data dictionaries
├── audit/                          # provenance/audit tables for prior merges
└── scripts/                        # scripts used to create and validate v3.1
```

## What changed in v3.1

### Protein table

All 1,773 `target_uniprot_id` values were queried from UniProt REST on 2026-07-10. The following fields were added or refreshed:

- `sequence_length`
- `uniprot_family_annotation`, `protein_family`, `protein_subfamily`
- `protein_class`, `protein_class_rule`
- `membrane_protein_type`, `membrane_type_rule`
- `transmembrane_topology_summary`, `transmembrane_annotation_status`
- `signal_peptide_summary`
- `alphafold_model_id`, `experimental_pdb_ids`
- `uniprot_keywords`, `annotation_source`, `annotation_retrieved_at_utc`, `record_qc_status`

`transmembrane_count`, `pdb_structures`, and `reviewed` were refreshed from the returned UniProt record. Protein classes are deliberately conservative keyword/family rules; use `protein_class_rule` when interpreting them as they are not manually curated ontology assignments.

Results: 1,182 proteins have a filled family annotation; 1,764 have an AlphaFold cross-reference. 620 records have `record_qc_status = review_membrane_scope`, meaning their current UniProt record contains no explicit transmembrane feature and needs scope review rather than being silently interpreted as a verified non-membrane protein.

### Relationship table

Added:

- `relationship_context_id` — stable v3.1 row identifier.
- `evidence_type` — `curated_target_assertion`, `assay_activity`, `structure_ligand_observation`, or `uniprot_binding_site_annotation`.
- `evidence_scope` and `evidence_directness` — context interpretation fields.
- `activity_relation` — separates `<`, `>`, `=`, etc. from the endpoint.
- `activity_value_unit` — `uM` when a standardized value exists.
- `record_qc_status`.

Normalizations:

- `KI` → `Ki` (310 records)
- `KD` → `Kd` (1,225 records)
- `IC50<` → `activity_type = IC50`, `activity_relation = <` (1 record)

Nine rows with `activity_type = A2` or `D2` were retained and marked `review_suspected_activity_type`; do not use them as standard affinity endpoints without source review.

### Compound table

Compound-level relationship counts and target counts were recomputed from the v3.1 relationship table. CID `126480232` (COBOMARSEN) has no current relationship row; it was retained rather than deleted, set to zero relationship/target counts, and marked `orphan_no_current_relationship` with its prior summary values retained in `record_qc_note`.

### Pocket table

Added reproducibility metadata:

- `structure_model_id`, `structure_model_source`, `structure_model_version`
- `pocket_method`, `pocket_method_version`, `pocket_qc_status`

The source log identifies the models as AlphaFold DB v6 and the method as geometry-based concavity detection. Pocket pLDDT/PAE and pocket physicochemical descriptors are not yet included.

## Recommended workflow for continuation

1. Start from `v3_1/` and preserve it as an immutable release snapshot.
2. Read `docs/DATA_DICTIONARY.md`, `docs/SHEET2_SHEET3_COLUMN_GUIDE.md`, `docs/PROCESS_LOG.md`, then `docs/FULL_WORK_LOG.md`.
3. Read `v3_1/V3_1_QC_REPORT.json`, `V3_1_VALIDATION_REPORT.json`, and `V3_1_XLSX_VALIDATION_REPORT.json` before changing any rows.
4. Resolve the 620 `review_membrane_scope` proteins using an explicit membrane-inclusion rule and source evidence.
5. Review the 9 flagged `A2/D2` activity rows against their original source records.
6. Decide whether the structure-ligand and UniProt site records remain in the same export table or move to dedicated evidence tables in a later v4 design.
7. If adding more data, preserve `relationship_context_id` for existing rows and create new IDs only for new relationship contexts.
8. Rebuild the TSV files first, run validation, then regenerate the XLSX.

## Scripts

`scripts/` contains the exact v3.1 build and validation scripts:

- `build_v31_data.py` — UniProt enrichment and four TSV transformations.
- `validate_v31.py` — TSV shape, ID, foreign-key, and pocket-count validation.
- `build_v31_xlsx_openpyxl.py` — memory-efficient four-Sheet XLSX generation.
- `validate_v31_xlsx_zip.py` — XLSX ZIP/XML validation.

They retain the original absolute working paths. Before reuse outside this machine, update `ROOT`, `BASE`, and `OUT` to the recipient's copy of the package or project directory. `build_v31_data.py` requires network access to UniProt REST and will cache responses in `uniprot_v3_1_enrichment_cache.json`.

## Important semantic cautions

- A relationship row is a **relationship/context record**, not always a direct pharmacological drug–target assertion.
- `structure_ligand_observation` means structural context; it is not automatically direct pharmacological evidence.
- `assay_activity` does not by itself prove direct binding.
- `uniprot_binding_site_annotation` is curated site evidence and should not be propagated to every compound linked to that protein.
- `target_group_id` tracks prior multi-UniProt expansion and must not be repurposed as a protein-family identifier.
- Empty, `0`, and `not_found_in_current_sources` have different meanings; maintain that distinction.

## Known limitations / next priorities

1. Define and audit membrane-protein inclusion criteria for the 620 review-flagged records.
2. Inspect the 9 `A2/D2` activity records.
3. Preserve original assay values, units, comparator, assay ID, and citation in future enrichment; many v3 rows only retain standardized `uM` values.
4. Add AlphaFold pLDDT/PAE, model file checksum, and pocket descriptors if the next use case involves docking or structure-aware ML.
5. For a future v4, split source records, evidence, assay/activity, structure-ligand observations, and literature references into dedicated tables; retain the current four tables as export views.

## Validation status at handoff

- TSV rows: Binary 101,802; Protein 1,773; Small Molecules 87,611; Pocket 6,815.
- All four TSVs: no malformed rows.
- Entity keys and `relationship_context_id`: unique.
- Relationship-to-protein, relationship-to-compound, and pocket-to-protein links: closed.
- Protein pocket counts: consistent with pocket rows.
- XLSX: ZIP integrity passed; all four Sheets have expected rows/columns, frozen header row, auto-filter, dark-blue bold header styling, and no formula cells.
