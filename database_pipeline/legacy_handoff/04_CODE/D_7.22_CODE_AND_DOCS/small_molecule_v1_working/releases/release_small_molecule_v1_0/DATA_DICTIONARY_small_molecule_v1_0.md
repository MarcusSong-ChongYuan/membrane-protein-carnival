# Data Dictionary: Small-Molecule v1.0

## Identity fields

- `compound_internal_id`: stable parent-level identifier.
- `compound_form_id`: exact salt, charge, stereochemical or multicomponent
  form identifier.
- `standard_inchikey`: full Standard InChIKey used for parent identity.
- `connectivity_key`: first InChIKey block, used only for related-structure
  searching.
- `identity_confidence`: `high` for cross-source or ChEBI-supported identity;
  otherwise `medium` for a valid single-source structure.

## Scope fields

- `compound_scope_status`: `core` or `extended` in the master.
- `final_scope_status`: `core`, `extended`, `unresolved`, or `excluded` in the
  source cross-reference table.
- `mapping_status`: mapping outcome for each source record.

## Evidence summary fields

- `protein_target_count`: distinct membrane proteins supported by default
  v4.2 evidence.
- `BE1_evidence_count`: compound-specific structure or site observations.
- `BE2_evidence_count`: direct quantitative binding observations.
- `BE3_evidence_count`: pharmacology, functional, or curated target
  observations.
- `best_standard_value_nM`: minimum positive normalized Ki/Kd value among
  default observations; it is not a potency average.

## Quality fields

- `record_qc_status`: `ok` or `review`.
- `qc_notes`: explicit representation, parent-mapping, or merge cautions.
- `structure_qc_note`: source-level parsing or identifier discrepancies.

