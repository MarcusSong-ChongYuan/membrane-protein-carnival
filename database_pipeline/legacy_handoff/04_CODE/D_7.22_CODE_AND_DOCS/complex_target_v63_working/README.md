# HuMemLigDB V6.3 candidate: complex-target module

This working module adds macromolecular complexes as entities separate from the
frozen V6.2 protein table. V6.2 is read-only and has not been modified.

## Current outputs

The hardened preview is in `output_v0_2/`.

- `complex_target_master_v0_2.tsv`: source-anchored complex entities involving
  at least one V6.2 membrane-protein audit accession.
- `complex_target_components_v0_2.tsv`: component assertions, UniProt and
  isoform identifiers, stoichiometry status, and V6.2 membrane evidence status.
- `complex_source_crosswalk_v0_2.tsv`: source record to complex-target mapping.
- `complex_identity_review_v0_2.tsv`: same-component-set candidates that were
  deliberately not auto-merged.
- `complex_structure_evidence_candidates_v0_2.tsv`: PDB/assembly candidates;
  chain-to-component mapping is not yet resolved.
- `complex_binding_evidence_candidates_v0_2.tsv`: raw complex-context ligand or
  drug annotations with conservative compound-identity classification.
- `complex_binding_evidence_v0_2.tsv`: formal release table schema; it is empty
  until direct complex-level binding and component identity are verified.
- `COMPLEX_TARGET_MODULE_V0_2_VALIDATION.json`: counts, policies, checksums and
  blocking-error status.

## Current status

- 19,791 audit-level complex entities.
- 5,298 curated complex candidates contain at least one V6.2 E1/E2 default
  membrane protein and pass the current automatic release gate.
- 13,051 predicted-only entities remain outside the default layer.
- 97,993 source component assertions.
- 425 cross-source same-component-set groups remain in identity review.
- 24,954 complex ligand/drug candidates after exact deduplication.
- 190 candidates map uniquely through an explicit ChEBI/ChEMBL/PubChem ID, but
  none is promoted to formal complex-binding evidence without directness review.
- 4,053 structure candidates await PDB chain, biological assembly and component
  mapping.

## Rebuild

The scripts contain no machine-specific project root. Supply V6.2 input paths
explicitly:

```powershell
python scripts/build_complex_target_module_v0_1.py `
  --protein-master <human_membrane_protein_master_v6_2.tsv> `
  --complexportal-curated raw/complexportal_9606_20260114.tsv `
  --complexportal-predicted raw/complexportal_9606_predicted_20260114.tsv `
  --corum-human raw/corum_5_3_humanComplexes.txt `
  --pdbe-assembly <protein_biological_assembly_v2.tsv.gz> `
  --output output_v0_1

python scripts/finalize_complex_target_module_v0_2.py `
  --protein-master <human_membrane_protein_master_v6_2.tsv> `
  --compound-master <small_molecule_master_v1_3.tsv> `
  --v01 output_v0_1 `
  --output output_v0_2
```

Python 3.12 and pandas are required for V0.2.

## Scientific boundary

Membership of a membrane protein in a curated complex does not prove that every
ligand associated with the protein binds the assembled complex. Consequently,
V6.2 single-protein evidence is never promoted automatically. Complex Portal or
CORUM ligand/drug annotations are candidates until compound identity, relation
directness, binding subunit, assembly state and experimental context are checked.

