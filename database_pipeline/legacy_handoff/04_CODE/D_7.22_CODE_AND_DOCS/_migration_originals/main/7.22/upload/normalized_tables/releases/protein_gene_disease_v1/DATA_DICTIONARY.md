# Data dictionary

## Entity tables

- `protein_gene_mapping_v1.*`: one protein-gene mapping; protein names and UniProt disease annotations are stored here once rather than repeated across every association. `uniprot_disease_objects` preserves the structured annotation payload as JSON text.
- `disease_dictionary_v1.*`: one disease/phenotype/measurement entity; descriptions, categories and therapeutic areas are stored here.

## Core relation

`protein_gene_disease_direct_core_v1.parquet` contains identifiers, direct association score, project tier, primary flag, evidence-category booleans, nullable evidence count, mapping status, source release and retrieval date. It deliberately excludes repeated protein/disease descriptions and score JSON.

## Evidence summary

`protein_disease_evidence_summary_v1.parquet` contains one row per Ensembl gene, disease and datasource, with the datasource score, inferred contributing datatype, project evidence category and source release. `evidence_count` is null when the association API does not expose an exact count.

`project_evidence_tier` and `primary_disease_flag` remain mempro1 internal ranking outputs, not official Open Targets or UniProt ratings and not causality assertions.

<!-- ABCE_SCORING_V1_1:START -->
## v1.1 scoring fields and outputs

`rule_a_pass`, `rule_b_pass`, `rule_c_pass` and `rule_e_pass` are boolean rule membership fields. `highest_rule_passed` is A/B/C/E. `evidence_level` is `low`, `medium`, `high` or `very_high`; `evidence_level_rank` is 1 through 4. `default_release_inclusion` equals `rule_b_pass`. `scoring_policy_version` is `pgd_evidence_level_abce_v1`.

`protein_gene_disease_abce_scored_v1.*` is archive-only and contains all Rule A relationships. `protein_gene_disease_rule_b_scored_v1.*` is the default main release and contains only medium/high/very-high rows. `protein_primary_disease_rule_b_v1.tsv`, `protein_gene_disease_rule_b_top10_v1.tsv` and `protein_gene_disease_rule_b_summary_v1.xlsx` are derived exclusively from Rule B relationships. Independent Rule A/C/E and level-split Parquets are archive-only.
<!-- ABCE_SCORING_V1_1:END -->
