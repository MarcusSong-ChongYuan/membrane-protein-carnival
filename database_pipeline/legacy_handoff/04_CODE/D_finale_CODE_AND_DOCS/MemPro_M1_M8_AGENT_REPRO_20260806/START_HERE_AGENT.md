# START HERE — MemPro M1–M8 Agent reproduction bundle

This package contains the exact input subset, scripts, cached computational augmentations, figures, PowerPoint, speaker script and QA needed to reproduce the 2026-08-06 MemPro M1–M8 scientific figure revision.

## Fast rebuild

Open PowerShell in the extracted folder:

```powershell
.\RUN_REPRODUCE.ps1 -Python python -Install -SkipPpt
```

This redraws M1–M8 from the included frozen inputs and the included PAGTN, AlphaFold and membrane-side cached outputs. The final check is written to `09_qa/AGENT_REPRODUCTION_VALIDATION.json`.

## Full computational rerun

```powershell
.\RUN_REPRODUCE.ps1 -Python python -Install -Full -SkipPpt
```

Full mode additionally retrains the property-supervised PAGTN experiment with seed 42, rebuilds the membrane-side site classification, and refreshes AlphaFold API metadata. The AlphaFold step needs internet. DGL wheels are platform-specific; on non-Windows systems install a compatible PyTorch/DGL pair before running full mode.

## Rebuild the PowerPoint

The verified 18-slide PowerPoint is already included in `07_ppt`. To rebuild it, use a Codex environment, set `CODEX_NODE_MODULES` to the bundled runtime `node_modules` that contains `@oai/artifact-tool`, and omit `-SkipPpt`.

## Package structure

- `inputs/`: exact read-only input subset used by the figures.
- `01_source_audit/` and `02_data/`: source/field audit and intermediate statistics.
- `03_figures/`: PNG, SVG and PDF outputs.
- `04_pagtn/`, `05_alphafold/`, `06_binding_site_side/`: method outputs and QA.
- `07_ppt/`: verified PPT and detailed script.
- `08_reproducibility/scripts/`: reproduction code.
- `09_qa/`: validation reports and all slide renders.
- `BUNDLE_MANIFEST.json` and `BUNDLE_SHA256SUMS.tsv`: provenance and integrity.

## Scientific boundaries

- PAGTN is an encoder and HDBSCAN performs clustering; this experiment is physicochemical-property supervised, not a universal pretrained chemical embedding.
- AlphaFold pLDDT is model confidence, not experimental proof.
- Membrane-side assignments are protein-topology-relative and are not atom-coordinate OPM z-axis calculations.
- Frozen files under `inputs/` must not be edited.
