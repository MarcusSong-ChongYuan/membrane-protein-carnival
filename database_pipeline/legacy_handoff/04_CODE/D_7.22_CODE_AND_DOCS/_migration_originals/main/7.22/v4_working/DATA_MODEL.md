# v4 Data Model

## Reference entities

- `human_membrane_protein_master`: one row per canonical human protein.
- `protein_isoforms`: canonical-to-isoform relationships.
- `protein_identifier_history`: current and historical identifiers.
- `membrane_annotation_evidence`: source-level membrane assertions.
- `protein_classification`: functional, topology, location, and family axes.
- `small_molecule_master`: one row per normalized chemical entity.
- `compound_identifier_crosswalk`: source and normalized compound identifiers.

## Evidence entities

- `sources`: database releases and retrieval metadata.
- `publications`: PMID, DOI, patent, and source citations.
- `assays`: assay definition, target type, organism, and method.
- `activity_measurements`: original and standardized quantitative measurements.
- `compound_protein_assertions`: curated or source-level target assertions.
- `structure_ligand_observations`: structure-context records.
- `binding_site_annotations`: curated residue/site annotations.

## Derived views

- `human_membrane_protein_small_molecule_status`: one derived status per protein.
- `protein_small_molecule_summary`: protein-level counts and best evidence.
- `compound_protein_evidence_export`: denormalized exchange view.

