# Human Membrane Protein Master v5.3 — sequence enhanced

This release preserves the v5.2 ABC × E0–E3 biological/evidence model and adds
the current UniProtKB canonical amino-acid sequence to every audited record.

## Coverage

- All audited records: 10,997
- E1 experimental core: 3,451
- E2 strongly supported: 4,453
- Default website release (E1+E2): 7,904
- E3 prediction candidates: 2,652
- E0 excluded audit records: 441
- Canonical sequences retrieved: 10,997 (100%)
- Missing sequences: 0
- Sequence-length mismatches against the previous table: 0

## Sequence files

- `human_membrane_all_audit_v5_3.fasta` contains all 10,997 canonical sequences.
- `human_membrane_default_E1_E2_v5_3.fasta` contains the 7,904 default website records.
- Every TSV includes one full `canonical_sequence` field per record.
- Excel stores sequences in `Canonical sequence part 1` and `part 2`. Only Q8WZ42
  requires part 2 because its 34,350-aa sequence exceeds Excel's per-cell limit.

The v5.2 release is not modified.
