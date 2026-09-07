# Complex-target module data dictionary

## Entity model

### `complex_target_master_v0_2.tsv`

- `complex_target_id`: stable HuMemLigDB complex identifier.
- `complex_name`: preferred source name.
- `protein_component_ids`: canonical UniProt accessions asserted by sources.
- `membrane_component_ids`: components found in the 10,997-protein V6.2 audit.
- `default_membrane_component_ids`: components in the E1/E2 V6.2 default layer.
- `complex_class`: heteromer, homomer/single-component, or unresolved.
- `curated_source_count`: distinct curated contributing databases.
- `predicted_source_count`: distinct predicted contributing databases.
- `component_set_conflict_flag`: sources assigned to one entity disagree on
  component membership.
- `required_optional_subunit_status`: remains `not_yet_curated`; absence of a
  flag must never be interpreted as optional or required.
- `default_release_inclusion`: curated, has at least one E1/E2 membrane
  component, and has no component-set conflict.

### `complex_target_components_v0_2.tsv`

One row is a source component assertion, not necessarily one unique biological
subunit.

- `component_assertion_id`: assertion primary key.
- `complex_target_id`: foreign key to the complex master.
- `component_entity_type`: protein, subcomplex, small molecule, RNA or other.
- `component_uniprot_id`: canonical accession for protein components.
- `component_uniprot_isoform_id`: exact isoform when present in the source.
- `copy_count`: explicit count; zero means unknown, not absent.
- `component_relation`: direct source participant or expansion of a nested
  complex.
- `required_subunit_status` / `optional_subunit_status`: unknown until a source
  makes this distinction or it is manually curated.
- `protein_membrane_evidence_level`: V6.2 E1/E2/E3 status.
- `protein_website_default_flag`: membership in the V6.2 default protein layer.

### `complex_binding_evidence_candidates_v0_2.tsv`

- `evidence_scope`: complex-level annotation or component-only context.
- `candidate_entity_type`: small-molecule candidate, protein ligand, multiple
  identifiers, or ambiguous free text.
- `source_compound_id_type` / `source_compound_id`: parsed explicit identifier.
- `compound_internal_id`: populated only for a unique V6.2 identifier match.
- `compound_mapping_status`: exact, absent, ambiguous or excluded.
- `complex_relation_directness_status`: whether direct binding remains unknown.
- `default_release_inclusion`: always zero in V0.2.

### `complex_structure_evidence_candidates_v0_2.tsv`

PDB cross-references from Complex Portal are intersected with cached PDBe
biological assemblies. A PDB cross-reference alone does not identify the correct
biological assembly, chain, bound compound or binding subunit, so all rows remain
candidates.

## Identity rule

Records are merged only through an explicit source identity cross-reference.
Identical component sets alone generate a review record because the same proteins
can form different stoichiometries, activation states or functional assemblies.

