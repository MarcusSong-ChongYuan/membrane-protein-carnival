# Completeness and residual review report

## What is complete in this release

- All 20,416 reviewed human reference-proteome entries were eligible for the external-source union.
- All 9,823 TrEMBL transmembrane candidates have canonical-resolution rows.
- The exact 82 no-reviewed-gene-symbol special rows are isolated and adjudicated by additional identifiers/external evidence.
- All 377 v3.1 legacy rows have identifier and membrane-scope outcomes.
- All 680 unresolved single-pass rows have a v5 topology outcome.
- All v5 master rows have a release tier, membrane decision, evidence matrix row and classification status.

## Residual review, intentionally not hidden

- Tier R contains **1,109** reviewed proteins.
- TrEMBL rows without one unique reviewed canonical accession: **232**.
- Single-pass exact type still unresolved: **165**; another **73** have orientation but remain type I/III ambiguous.
- Strict unclassified: **1,426** in the full union, **574** in Tier A+B, and **1,132** in Tier A+B+C.
- Cross-source conflict/review rows: **1,531**.
- OPM, Membranome, UniTmp/PDBTM, GPCRdb and TCDB provide public downloads/APIs, but explicit redistribution terms were not uniformly stated. Confirm those terms before mirroring third-party annotations in a public production website.

## Coverage is versioned

The union is complete only relative to the sources and dates in `SOURCE_REGISTRY_v5.tsv`. A release refresh must rerun downloads, mapping, tiering and QA; counts should not be expected to remain constant.
