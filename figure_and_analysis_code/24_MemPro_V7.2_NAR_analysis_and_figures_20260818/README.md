# MemPro V7.2 NAR analysis and figure package

This directory contains a read-only analytical layer derived from the frozen
MemPro V7.2 release. The source release is never modified.

## Reproduction order

Use the Python environment recorded in `08_QA/environment.json` and run:

1. `01_scripts/00_build_analysis_layer.py`
2. `01_scripts/01_compute_statistics.py`
3. `01_scripts/01b_compute_evidence_flow.py`
4. `01_scripts/02_compute_chemical_space.py`
5. `01_scripts/03_draw_main_figures.py`
6. `01_scripts/04_draw_graphics.py`
7. `01_scripts/05_run_qa.py`

The presentation is built after the figures pass QA. All figure-level source
data are under `03_source_data`; rendered output is under `04_figures`.

## Interpretation boundaries

- The compound registry and interaction-linked compound universe are separate.
- RNA, IHC and mass-spectrometry layers are never pooled as one expression
  measurement.
- Only mapped and measured records enter expression denominators.
- Contributing database count is not described as independent experiment count.
- Membrane-side unknowns and unmapped diseases remain visible.
- Cross-module prioritisation is hypothesis-generating and is not causal proof.
