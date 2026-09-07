# MemPro FORMAL — website development handoff

## Purpose and status

This package is the **read-only website-development dataset** for MemPro. It contains the frozen public-positive core release and its V7.2.1 companion annotations. It is intended to let a web team design the data model, import data locally, and build search/detail/download views without access to historical project folders.

**Do not change the source TSV files in this package.** Create derived database tables, indexes and search caches in a separate application workspace.

## Package scope

Included:

- `data/core/` — V7.2 core release tables, release rules, manifest and QA.
- `data/companion/` — V7.2.1 additive annotations: B/C membrane compartment/orientation, chemical classifications, structure-model metadata and source contribution.
- `metadata/` — source registry, formal manifest, table schema catalogue and data dictionary.
- `docs/` — earlier web data-model proposal, filter taxonomy and five-axis labels.
- `qa/` — formal release checks.

Deliberately excluded:

- historical V3–V6/V7 candidate directories;
- raw downloaded upstream databases;
- external-source audit candidates (DrugCentral, TTD and current ChEBI review), which are **not part of FORMAL**;
- docking files/results;
- `04_archive_negative_evidence` (a large frozen audit archive, not a default public website dataset);
- legacy pipeline scripts.

## Core public scale

| Entity / relation | Release scale |
|---|---:|
| Formal membrane proteins | 7,800 |
| Default website proteins (`default_web_inclusion_v72=1`) | 7,374 |
| Retained positive evidence records | 942,455 |
| All formal protein–canonical compound pairs | 529,168 |
| Default website pairs (`default_web_inclusion_v72=1`) | 501,242 |
| Registered canonical compounds | 646,670 |

Always show a table-specific denominator and inclusion field. Do not mix all-formal counts with default-website counts.

## Suggested import order

1. **Reference/source tables**
   - `source_registry_v72.tsv.gz`
   - `source_registry_v721.tsv.gz`
   - disease master/hierarchy/source mapping tables.
2. **Protein and gene layer**
   - `protein_master_v72.tsv.gz`
   - `canonical_protein_membrane_class_v72.tsv.gz`
   - `gene_membrane_summary_v72.tsv.gz`
   - `protein_membrane_evidence_and_risk_v72.tsv.gz`
   - V7.2.1 `protein_membrane_orientation_v721.tsv.gz`, `protein_membrane_role_inference_v721.tsv.gz`, and `protein_structure_model_v721.tsv.gz`.
3. **Compound layer**
   - `compound_master_v72.tsv.gz`
   - `compound_form_v72.tsv.gz`
   - V7.2.1 `compound_chemical_classification_v721.tsv.gz` and `compound_gtopdb_exact_xref_v721.tsv.gz`.
4. **Interaction layer**
   - `protein_compound_pair_v72.tsv.gz`
   - `positive_interaction_evidence_v72.tsv.gz`
   - `evidence_lineage_v72.tsv.gz`
   - `evidence_target_entity_assignment_v72.tsv.gz`
   - `evidence_conflict_v72.tsv.gz`
   - V7.2.1 `evidence_source_contribution_v721.tsv.gz`.
5. **Disease, expression, structure and complex views**
   - disease, expression, subcellular localization, binding-site and complex tables in `data/core/01_release_tables/`.

## Key website relations

```text
gene
  └─ canonical protein ──< protein–compound pair >── canonical compound
         │                       │                         └─ compound form
         │                       └─ positive evidence ──< evidence lineage/source
         ├─ protein–disease relation ── canonical disease ── anatomy / therapeutic area
         ├─ expression / localization
         ├─ binding-site assertion
         ├─ protein complex ── components / state / stoichiometry
         └─ isoform / processed-form / membrane orientation annotations
```

## Required UI semantics

- Display `source_databases` as contributing databases, **not** as a count of independently confirmed experiments.
- Keep evidence tier (`BE1/BE2/BE3`) and membrane-evidence tier (`E1/E2/E3`) distinct.
- Show missing/unmapped expression separately from detected/not detected. RNA, IHC and localization must not share a denominator.
- Binding-site absence does not mean absence of a protein–compound interaction.
- B/C orientation fields distinguish `DIRECT`, `MECHANISM_INFERRED`, and `UNRESOLVED`; do not collapse them.
- A compound registry entry is not necessarily interaction-linked. Use pair/evidence tables to assert an interaction.
- Complex rows may be `complex-specific`, `source-asserted context`, or `component-context-only`; label the context visibly.

## Documentation and validation

- Field meanings: `metadata/release_metadata/V72_TABLE_SCHEMA_CATALOG.tsv` and `metadata/release_metadata/V721_DATA_DICTIONARY.tsv`.
- Formal core routing/inclusion rules: `data/core/02_metadata/V72_TABLE_ROUTING.tsv` and `V72_RELEASE_RULES.tsv`.
- SHA-256 manifests: `metadata/release_metadata/FORMAL_MANIFEST_SHA256.tsv`, `data/core/02_metadata/SHA256SUMS.tsv`, and V7.2.1 QA manifests.
- Existing Chinese web planning proposal: `docs/web_planning_files/MemPro_V7.2_WEB_DATABASE_PLAN_CN.md`.

## Basic PostgreSQL import pattern

Use a staging schema with all columns imported as text first; validate keys and enumerations before typing/indexing. TSV files are UTF-8 and most large tables are gzip-compressed.

```sql
CREATE TABLE staging.protein_compound_pair (...text columns...);
-- Decompress outside the database process, then COPY from STDIN or a controlled file loader.
```

Recommended initial indexes:

- protein: `canonical_uniprot_accession`, `approved_symbol`, `gene_entity_id`;
- compound: `compound_internal_id`, `standard_inchikey`, `pubchem_cids`, `chembl_ids`;
- pair: `(target_uniprot_id, compound_internal_id)`;
- evidence: pair identifiers, `source_database`, evidence tier;
- disease: canonical disease ID and target accession.

Do not expose every column in the first UI release. Begin with protein, compound, pair, evidence, disease, expression summary, binding-site and complex-detail endpoints.
