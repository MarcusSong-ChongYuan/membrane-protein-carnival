# MemPro V6.1 finalization execution plan

## Operating rule

V6.0 is immutable. Current completed inputs are used to build a refreshable
premerge and delta preview. No preview file is a formal release. Final HMPD
compound and form identifiers are assigned only after PubChem P2, PubChem P3,
and BRENDA are complete.

## Current upstream sequence

1. PubChem P1 fetch and identity mapping: complete.
2. PubChem P2 fetch, then P2 identity mapping: running.
3. PubChem P3 fetch, then P3 identity mapping: automatic after P2.
4. BRENDA structure-backed name resolution: automatic after P3 fetch.

## Finalization stages

### F1. Cross-source global identity merge

- Inputs: V6.0 compound master, PubChem P1/P2/P3, PDBe CCD, PDBbind SDF,
  BRENDA structure-verified mappings.
- Exact full InChIKey is the primary identity key.
- Connectivity-only agreement never triggers an automatic full-identity merge.
- Output: one global exact-structure cluster per full InChIKey and a separate
  connectivity/form audit table.

Estimated runtime after all inputs are ready: 1–2 hours.

### F2. Parent, salt, charge, and stereochemistry hierarchy

- RDKit cleanup is recorded with the original structure retained.
- Fragment parent and charge parent are computed separately.
- Stereochemically specified and unspecified records remain distinct exact
  forms even when they share a connectivity family.
- Ambiguous families remain in the audit queue.

Estimated runtime: 1–2 hours.

### F3. Evidence and master-table rebuild

- Update the small-molecule master and compound-form hierarchy.
- Remap binding evidence to final compound and form identifiers.
- Preserve active, inactive, inconclusive, and unspecified assay outcomes.
- Keep covalent, non-small-molecule, solvent/buffer, and unresolved identities
  in their policy-specific layers.

Estimated runtime: 2–4 hours.

### F4. Conflict and negative-evidence processing

- Reconcile CID/structure conflicts, cross-assay activity conflicts, and
  quantitative heterogeneity.
- Rebuild the negative evidence and review tables without promoting uncertain
  records into the default positive layer.

Estimated runtime: 1–2 hours.

### F5. Referential and source reconciliation

- Validate all primary keys and foreign keys.
- Reconcile every input row to a release, negative, review, excluded, or
  duplicate disposition.
- Compare source totals, distinct objects, and overlap counts with staging
  reports.

Estimated runtime: 1–2 hours.

### F6. Validation and freeze

- Generate the V6.1 validation report, release manifest, SHA-256 checksums,
  methods, data dictionary, source registry, and release notes.
- Freeze only if every blocking validation is zero.

Estimated runtime: about 1 hour.

## Formal release products

- Human membrane-protein master table.
- Protein–gene–disease relation table.
- Small-molecule master table.
- Compound parent/form hierarchy.
- Small-molecule–membrane-protein binding evidence table.
- Binding-site instances.
- Protein–compound summary.
- Negative evidence, review queue, conflict audit, source coverage, and
  reproducibility manifest.
