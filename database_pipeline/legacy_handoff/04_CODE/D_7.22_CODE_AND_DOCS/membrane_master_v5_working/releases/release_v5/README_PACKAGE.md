# Human Membrane Protein Master v5 release package

Release date: 2026-07-24

## Entry points

- `human_membrane_protein_master_v5.xlsx`: human-readable review workbook.
- `human_membrane_protein_master_v5.tsv`: full machine-readable master table.
- `human_integral_membrane_view_v5.tsv`: Tier A.
- `human_core_membrane_view_v5.tsv`: Tiers A+B.
- `human_extended_membrane_view_v5.tsv`: Tiers A+B+C.
- `review/human_membrane_review_queue_v5.tsv`: Tier R records.

## Supporting folders

- `crosswalks/`: TrEMBL canonicalization, legacy-record, topology,
  source-evidence, and classification audit tables.
- `normalized/`: normalized source-specific evidence tables used by the build.
- `review/`: manual-review and source-conflict queues.
- `reports/`: machine-readable build statistics and validation results.
- `reproducibility/`: scripts, tier policy, source-download manifest, and
  baseline hashes.

The raw third-party downloads are not redistributed in this package. Their
official URLs, release identifiers, roles, and license notes are recorded in
`SOURCE_REGISTRY_v5.tsv` and the reproducibility manifest.
