# Compound classification data audit

## Frozen statistical universes

- Compound registry: **646,670** canonical compounds.
- Interaction-linked compounds: **240,055** canonical compounds.
- Formal protein–compound pairs: **529,168** unique pairs.
- Structure-valid interaction-linked compounds: **240,054**.
- Excluded/quarantined from structural analysis: **1**; registry identity retained.

These denominators remain separate throughout the exploration pack.

## Identity and structure checks

- Missing SMILES among interaction-linked compounds: **0**.
- Missing InChIKey: **0**.
- Duplicate recalculated canonical-SMILES groups: **0** involving **0** rows.
- Duplicate InChIKey groups: **0** involving **0** rows.
- ChEBI direct-class annotations: **6,161 / 240,054**.

Duplicate structures are reported but not silently merged because registry identity and source provenance can differ. The quarantined name–structure conflict remains registered but is excluded from scaffold, descriptor, fingerprint, PCA, UMAP, clustering and similarity analyses.

## Taxonomy boundary

No complete ClassyFire/ChEBI hierarchy exists locally. Existing ChEBI direct classes are preserved as authoritative direct annotations. `chemical_superclass`, `chemical_subclass` and `direct_parent` remain `NOT_AVAILABLE_NEEDS_EXTERNAL_MAPPING`. The broad `chemical_regime_local` field is an explicitly local deterministic rule layer, not ClassyFire.

## Local regime composition

| chemical_regime_local          |   compound_count |   percent |
|:-------------------------------|-----------------:|----------:|
| organic small molecule         |           192110 |     80.03 |
| carbohydrate / glycoside-like  |            29277 |     12.20 |
| peptide / peptide-like         |            11984 |      4.99 |
| macrocyclic compound           |             5416 |      2.26 |
| lipid / lipid-like             |              829 |      0.35 |
| nucleotide / nucleoside-like   |              401 |      0.17 |
| inorganic / no-carbon compound |               26 |      0.01 |
| organometallic compound        |               11 |      0.00 |

## Recalculation policy

Eligible interaction-linked compounds were reparsed with RDKit 2026.03.4. Murcko scaffolds and physicochemical descriptors were recalculated and compared with formal stored fields. Formal fields were not overwritten. Exact discrepancy counts are retained in the JSON audit.

## Downstream denominators

- V1, V2, V3, V6 and V9: **240,054 compounds**.
- V4 and V5: **529,167 eligible unique compound–protein pairs**.
- V7 PCA: deterministic stratified sample; exact n in its coordinate/source table.
- V8 UMAP: deterministic stratified sample; exact n and parameters in its JSON.
- V10 3D: deterministic stratified sample; exploratory only.
