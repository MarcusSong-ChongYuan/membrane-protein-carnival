# Human membrane-protein master database v5

This release completes the first-table rebuild as a **reviewed human membrane-protein union with explicit evidence tiers**, while keeping unreviewed TrEMBL candidates and uncertain external-only entries in separate review layers.

## Main result

- Reviewed human reference proteome examined: **20,416**
- v4 membrane/membrane-associated baseline: **9,813**
- v5 reviewed union: **10,997**
- New reviewed proteins surfaced by independent sources: **1,184**
- Tier A integral core: **5,496**
- Tier B monotopic/lipid-anchored extension: **438**
- Tier C peripheral membrane-associated extension: **3,954**
- Tier R evidence/review queue: **1,109**

The three publication views are:

1. `human_integral_membrane_view_v5.tsv` — Tier A.
2. `human_core_membrane_view_v5.tsv` — Tiers A+B.
3. `human_extended_membrane_view_v5.tsv` — Tiers A+B+C.

Tier R is not mixed into the default website result set. It is preserved in `human_membrane_review_queue_v5.tsv` so that coverage is broad without presenting predictions or location-only evidence as settled integral-membrane biology.

## Previously open issues

- **9,823 TrEMBL candidates:** all received a deterministic canonical-resolution status. 9,591 map to one reviewed canonical accession; the remainder are retained in auditable ambiguous or review states.
- **82 special TrEMBL rows:** the exact original 82 `no_reviewed_gene_symbol_match` records are isolated; 28 were rescued by HGNC/GeneID/Ensembl evidence, 36 have external membrane support without a reviewed canonical match, and 18 remain unconfirmed.
- **377 legacy v3.1 proteins:** all identifiers were resolved: 376 remain current reviewed accessions and one was mapped by approved symbol. Membrane-scope decisions are {'not_supported_as_membrane_in_v5': 328, 'external_support_review': 47, 'included_core_or_extended': 2}.
- **680 single-pass unresolved:** 442 received exact I/II/III assignments, 73 received orientation but remain I-versus-III ambiguous, and 165 remain unresolved.
- **4,423 v4 unclassified:** strict residual unclassified among retained v4 rows is now **1,111**. Family-level classifications are explicitly labelled and not misrepresented as mature functional classes.
- **GPCR/ion-channel/transporter hierarchy:** GPCRdb, GtoPdb 2026.2 and TCDB identifiers are stored in dedicated columns.
- **Peripheral proteins:** Tier C is a separate extension and is excluded from integral/core defaults.

## Important interpretation

“Complete” here means the union of the named, versioned sources under the documented policy. It does not mean that biology has a permanently closed list. New isoforms, new reviewed accessions, revised topology predictions and database updates will change future releases.
