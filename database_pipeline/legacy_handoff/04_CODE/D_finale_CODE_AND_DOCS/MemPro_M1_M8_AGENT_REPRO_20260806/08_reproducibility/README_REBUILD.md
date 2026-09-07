# Rebuild runbook

## Inputs

- Frozen V6.2: `D:/7.22/evidence_expansion_v2_working/releases/release_mempro_v6_2_20260730`
- V6.3 candidate: `D:/finale/02_V6.3_candidate_20260803`
- Revision output: `D:/finale/07_M1-M8图表修订_20260806`

## Order

1. Set `MEMPRO_REVISION_ROOT` to the revision output directory.
2. Run `build_revision_figures.py` for the initial M1-M8 set.
3. Run `build_m2_fixed.py` for scoped Spearman/bootstrap M2.
4. Run `classify_binding_site_side.py`, then rebuild M6.
5. Run `fetch_alphafold_confidence.py`.
6. In the isolated PAGTN environment, run `run_pagtn_sample_experiment_v2.py`.
7. Run `build_m5_pagtn_final.py`.
8. Build the deck with the artifact-tool script. The supplied `reference_template.pptx` is the approved template.
9. Run the presentation `slides_test.py` validator and `finalize_revision_release.py`.

## Environment

- Standard plots: the existing `D:/7.22/evidence_expansion_v2_working/tools/plotting_venv`.
- PAGTN: `04_pagtn/env`; versions are frozen in `pagtn_environment_freeze.txt`.
- PPT: Node.js plus `@oai/artifact-tool`; do not replace this step with python-pptx.

The scripts retain explicit source paths to prevent accidental use of a different release. For another device, mount/copy the frozen releases and update only the three documented root constants/environment variables.
