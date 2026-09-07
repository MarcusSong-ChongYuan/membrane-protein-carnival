# v5.4 Data Dictionary

## `human_membrane_audit_master_v5_4.tsv`

The complete 10,997-row v5.3 protein master plus disease and small-molecule
summary fields.

- `primary_gene_id`: Ensembl gene ID used by the selected primary relation.
- `has_supported_disease`: at least one included disease relation.
- `supported_disease_count`: unique source disease identifiers.
- `primary_disease_id`, `primary_disease_name`: selected primary relation.
- `best_disease_evidence_level`: `medium`, `high`, or `very_high`.
- `disease_relation_ids`: stable IDs linking to the detailed relation table.
- `has_BE1_structural_binding`: at least one default-included BE1 record.
- `has_BE2_direct_binding`: at least one default-included BE2 record.
- `has_BE3_pharmacology`: at least one default-included BE3 record.
- `experimental_small_molecule_count`: unique included compounds.
- `best_binding_evidence_level`: strongest included BE tier.
- `binding_evidence_count`: included source-observation count.
- `has_predicted_pocket`, `predicted_pocket_count`: BP prediction status; not
  experimental interaction evidence.

## `protein_gene_disease_relations_v5_4.tsv`

One protein–gene–disease source relation per row. The table preserves v3 Rule-B
fields needed for evidence interpretation and adds current UniProt/OMIM
relations as an independent source layer.

## `binding_evidence_master_v4_1.tsv`

Canonical full evidence table. One row is one source observation in a
protein–compound context.

- `evidence_id`: stable observation identifier.
- `target_uniprot_id`: canonical protein accession.
- `compound_id`: normalized source-aware compound key.
- `evidence_tier`: BE1, BE2 or BE3.
- `evidence_type`, `evidence_directness`: interpretation of the source row.
- `target_assignment_status`: single-protein or complex/ambiguous assignment.
- `activity_type`, `activity_relation`, `activity_value`, `activity_unit`:
  original/retained measurement representation.
- `standard_value_nM`: standardized value only when valid.
- `pdb_ids`, `binding_site_residues`, `ligand_het_id`: structure/site context.
- `default_release_inclusion`: whether the row contributes to derived website
  summaries.
- `originating_dataset`: migration or supplementation layer.

## Derived and child tables

- `binding_site_instances_v4_1.tsv`: experimental structure/site child records.
- `protein_compound_summary_v4_1.tsv`: one included protein–compound pair.
- `protein_binding_summary_v4_1.tsv`: one binding-status row per protein.
- `small_molecule_index_v4_1.tsv`: compound identifiers, structures, scope and
  source coverage.
- `binding_evidence_workbook_view_v4_1.tsv`: Excel-sized presentation view; all
  non-ChEMBL evidence plus default-included ChEMBL Ki/Kd at 10 µM or better.
- `v3_pgd_relations_not_in_v5_3_audit.tsv`: legacy disease relations whose
  targets are outside the v5.3 protein universe.
- `v3_binding_evidence_not_in_v5_3_audit.tsv`: analogous legacy binding audit.

