# Data dictionary

## Canonical tables

- `human_membrane_protein_master_v6_0.tsv`: protein identity, sequence,
  membrane evidence, classification, disease summary and refreshed canonical
  small-molecule counts.
- `protein_gene_disease_relations_v6_0.tsv`: normalized protein–gene–disease
  relations with evidence tiers and source provenance.
- `binding_evidence_master_v6_0.tsv.gz`: baseline and accepted incremental
  positive evidence. V6 fields record outcome, incremental source, duplicate,
  identity and conflict status.
- `small_molecule_master_v1_1.tsv`: canonical parent compounds, standardized
  structures, properties and refreshed membrane-protein evidence summaries.

## Supporting tables

- `compound_form_hierarchy_v1_1.tsv`: exact forms linked to canonical parents.
- `protein_compound_summary_v6_0.tsv.gz`: one canonical protein–compound pair
  per row.
- `binding_site_instances_v6_0.tsv.gz`: compound-specific sites with PDB and
  residue-coordinate provenance.
- `negative_binding_evidence_v1_0.tsv.gz`: reported inactive outcomes.
- `binding_evidence_review_queue_v6_0.tsv.gz`: records not eligible for
  canonical merging.
- `incremental_source_coverage_v6_0.tsv`: source-level reconciliation.

## Key fields

- `target_uniprot_id`: canonical protein foreign key.
- `compound_internal_id`: canonical parent-compound foreign key.
- `compound_form_id`: exact-form foreign key when available.
- `evidence_tier`: BE1 structural, BE2 direct quantitative, or BE3 functional
  and pharmacological evidence.
- `default_release_inclusion`: one only for the default high-confidence layer.
- `duplicate_status_v60`: exact mirror or unique/complementary status.
- `activity_outcome_v60`: source outcome retained independently of evidence
  tier.
