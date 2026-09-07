# Entity and evidence schema

## Entity hierarchy

```text
gene_entity
  └── canonical_protein_entity (one or more; never merged by symbol)
        └── protein_isoform (zero or more; exact UniProt isoform accession)

complex_target
  └── complex_target_component
        ├── canonical_protein_entity
        ├── protein_isoform
        ├── external_complex_protein_entity
        └── external protein_isoform
```

## Evidence target semantics

| `target_resolution_level` | Meaning |
|---|---|
| `isoform` | Source explicitly provided an isoform accession and it resolved exactly |
| `canonical_protein` | Source record points to a canonical UniProt accession only |
| `gene_product_unspecified` | Source operates at gene level; no isoform assertion is allowed |
| `complex_context_ambiguous` | Historical single-protein projection came from an unresolved complex context |
| `external_isoform` | Exact isoform belongs to a non-core protein used only as a complex component |
| `external_canonical_protein` | Non-core canonical protein used only as a complex component |

## Primary and foreign keys

| Table | Primary key | Important foreign key |
|---|---|---|
| gene entity | `gene_entity_id` | — |
| canonical protein | `canonical_protein_entity_id` | `gene_entity_id` |
| protein isoform | `protein_isoform_entity_id` | `canonical_protein_entity_id` |
| binding target bridge | `evidence_id` | gene/canonical/isoform entity fields |
| disease target bridge | `disease_relation_id` | gene/canonical entity fields |
| complex component | `component_assertion_id` | canonical and/or isoform entity fields |

## Important non-equivalences

- Same gene symbol does not mean same canonical protein record.
- Canonical accession evidence does not mean canonical-isoform-specific evidence.
- A complex external component is not automatically a membrane protein.
- A canonical protein sequence and an isoform sequence may be identical; their identifiers and evidence scopes remain distinct.
