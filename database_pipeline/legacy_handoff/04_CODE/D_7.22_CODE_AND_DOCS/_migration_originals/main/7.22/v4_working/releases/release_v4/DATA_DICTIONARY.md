# v4 Data Dictionary

## `human_membrane_protein_master_v4.tsv`

One row per reviewed canonical human membrane or membrane-associated protein.

- `membrane_protein_id`: stable v4 internal identifier.
- `target_uniprot_id`: canonical UniProt accession.
- `membrane_scope`: integral, monotopic, lipid-anchored, or peripheral-associated.
- `membrane_topology`: single-pass type, multi-pass, beta-barrel, intramembrane,
  lipid-anchor, or no explicit transmembrane feature.
- `membrane_inclusion_status`: accepted core, accepted extended, or review.
- `functional_primary_class`, `functional_subclass`: conservative functional class.
- `present_in_v3_1`: whether the accession was in the v3.1 protein release.

## `human_membrane_protein_small_molecule_status_v4.tsv`

One derived status row per membrane reference protein.

- `has_small_molecule_evidence`: 1 only when detailed qualifying evidence exists.
- `has_small_molecule_source_index_hit`: BindingDB or ChEMBL target-index hit.
- `small_molecule_evidence_status`: confirmed, index-only candidate, or none found.
- `small_molecule_evidence_level`: final detailed-evidence level SM0–SM5.
- `has_direct_quantitative_binding`: at least one detailed SM4/SM5 record.
- `small_molecule_count`: normalized compounds in detailed evidence.
- `independent_sources`: detailed evidence source labels.

## `small_molecule_master_v4.tsv`

Current v3.1 compounds plus accepted IUPHAR/BPS additions.

- `small_molecule_scope_status`: accepted core, excluded, review, or separate.
- `small_molecule_scope_class`: normalized chemical scope class.
- `is_core_small_molecule`: whether the record contributes to core coverage.
- `record_qc_status`: record-level quality or enrichment status.

## `compound_protein_evidence_v4.tsv.gz`

Denormalized v3.1 evidence exchange view.

- `relationship_context_id`: retained unique v3.1 evidence-context identifier.
- `protein_in_reviewed_membrane_reference`: reference-membership flag.
- `qualifies_for_core_small_molecule_flag`: protein and compound scope flag.
- `candidate_sm_level`: row-level evidence category.
- `assay_id`: derived assay identifier where applicable.

## Evidence child tables

- `activity_measurements_v4.tsv`: activity values and conservative direct-binding flag.
- `compound_protein_assertions_v4.tsv`: curated target assertions.
- `structure_ligand_observations_v4.tsv`: structure-context evidence.
- `binding_site_annotations_v4.tsv`: UniProt binding-site evidence.
- `assays_v4.tsv`: derived assay contexts.
- `publications_v4.tsv`: retained publication identifiers.
- `sources_v4.tsv`: source registry extracted from v3.1.

