# Inclusion policy

## Default release evidence

A new evidence record is default-included only when:

- the target is a canonical A/B/C human membrane protein;
- the compound has a unique canonical parent identity;
- the source-specific evidence and target assignment pass QA;
- the outcome supports binding, activity or an experimental structure;
- the ligand is not a solvent, buffer, crystal artifact, incomplete ligand or
  non-small-molecule oligomer;
- the record is not an exact quantitative source mirror.

## Nondefault but retained

Mapped supporting evidence that does not meet the default threshold remains
traceable with `default_release_inclusion=0`.

## Separate audit layers

- Inactive PubChem results: `negative_binding_evidence_v1_0.tsv.gz`
- Unresolved identities, inconclusive outcomes and unsuitable structures:
  `binding_evidence_review_queue_v6_0.tsv.gz`
- Same-assay positive/negative contexts:
  `pubchem_context_conflicts_v1_0.tsv.gz`
- Large cross-source quantitative ranges:
  `quantitative_heterogeneity_audit_v1_0.tsv.gz`

BRENDA names are not merged without an identifier or structure-level identity
check. DrugBank is outside the release scope by project decision.
