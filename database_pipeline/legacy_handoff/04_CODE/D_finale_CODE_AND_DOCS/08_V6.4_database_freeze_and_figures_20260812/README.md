# MemPro V6.4 internal release and refreshed scientific figures

Freeze date: 2026-08-12

Status: **FROZEN_INTERNAL**. All blocking scientific/relational QA gates passed. A public redistributable export is **conditional**, because BRENDA- and PDBbind-derived raw fields require source-specific permission review.

## What is included

- `01_release_overlay`: changed or newly normalized V6.4 tables. This is a thin release overlay and does not duplicate unchanged V6.2/V6.3 baselines.
- `02_statistics`: frozen tabular inputs used by the refreshed figures.
- `03_figures`: M1-M6 in 600-dpi PNG, editable SVG and PDF, with captions.
- `04_QA`: validation reports, source/version/license registry, contribution baseline and the review workbook.
- `05_review_queues`: records intentionally retained for manual review; none are silently discarded.
- `06_reproducibility/scripts`: executable scripts used for audits, candidate construction, statistics and figures.

## Required immutable baselines

- `../01_濮濓絽绱￠弫鐗堝祦_V6.2`
- `../02_V6.3_candidate_20260803`
- `../04_V6.3.1_isoform_accession_audit_20260803`

Their referenced files and SHA-256 values are recorded in `RELEASE_MANIFEST_V6_4.tsv`. V6.2 and V6.3 were kept read-only.

## Main release definitions

- All evidence: 3,003,306 evidence records and 1,639,146 unique protein-compound pairs.
- Default release: 991,847 evidence records and 576,038 unique pairs (`default_release_inclusion=1`).
- These scopes are intentionally separate in tables and figures.

## Controlled corrections

- PDB identity: current RCSB holdings snapshot (2026-08-12) used to remove concentration/unit contamination from PDB fields; 1,272 evidence rows retained with partial cleanup and 2 all-invalid structural rows excluded from the default structural layer.
- HPA IHC: 35,351 exact duplicate rows removed; 9,506 distinct records sharing legacy IDs were re-keyed deterministically.
- Complex accessions: V6.3.1 patch applied; 20 high-confidence replacements, 9 ambiguous assertions retained with uncertainty, and 1 conflicting assertion excluded.

## Reproduction order

Run scripts in `06_reproducibility/scripts` in the order documented by their numeric/functional names. The final figure entrypoint is `build_v64_figures.py`; it reads only frozen TSV/JSON statistics from `02_statistics` after paths are configured.

The manuscript stage and docking production were deliberately not started in this phase.