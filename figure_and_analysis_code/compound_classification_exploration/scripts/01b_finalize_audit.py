import json
from pathlib import Path
import pandas as pd

OUT = Path(r"D:\finale\compound_classification_exploration")
audit = json.loads((OUT / "01_data_audit/COMPOUND_CLASSIFICATION_DATA_AUDIT.json").read_text(encoding="utf-8"))
analysis = pd.read_csv(OUT / "02_classification_tables/INTERACTION_LINKED_COMPOUND_ANALYSIS.tsv.gz", sep="\t", low_memory=False)
regime = analysis["chemical_regime_local"].value_counts().rename_axis("chemical_regime_local").reset_index(name="compound_count")
regime["percent"] = regime.compound_count / len(analysis) * 100

text = f"""# Compound classification data audit

## Frozen statistical universes

- Compound registry: **{audit['compound_registry']:,}** canonical compounds.
- Interaction-linked compounds: **{audit['interaction_linked_compounds']:,}** canonical compounds.
- Formal protein–compound pairs: **{audit['formal_pairs']:,}** unique pairs.
- Structure-valid interaction-linked compounds: **{audit['structure_valid_interaction_linked_compounds']:,}**.
- Excluded/quarantined from structural analysis: **{audit['excluded_or_quarantined_interaction_linked_compounds']:,}**; registry identity retained.

These denominators remain separate throughout the exploration pack.

## Identity and structure checks

- Missing SMILES among interaction-linked compounds: **{audit['missing_smiles_interaction_linked']:,}**.
- Missing InChIKey: **{audit['missing_inchikey_interaction_linked']:,}**.
- Duplicate recalculated canonical-SMILES groups: **{audit['duplicate_canonical_smiles_groups']:,}** involving **{audit['duplicate_canonical_smiles_rows']:,}** rows.
- Duplicate InChIKey groups: **{audit['duplicate_inchikey_groups']:,}** involving **{audit['duplicate_inchikey_rows']:,}** rows.
- ChEBI direct-class annotations: **{audit['chebi_direct_class_available']:,} / {audit['structure_valid_interaction_linked_compounds']:,}**.

Duplicate structures are reported but not silently merged because registry identity and source provenance can differ. The quarantined name–structure conflict remains registered but is excluded from scaffold, descriptor, fingerprint, PCA, UMAP, clustering and similarity analyses.

## Taxonomy boundary

No complete ClassyFire/ChEBI hierarchy exists locally. Existing ChEBI direct classes are preserved as authoritative direct annotations. `chemical_superclass`, `chemical_subclass` and `direct_parent` remain `NOT_AVAILABLE_NEEDS_EXTERNAL_MAPPING`. The broad `chemical_regime_local` field is an explicitly local deterministic rule layer, not ClassyFire.

## Local regime composition

{regime.to_markdown(index=False, floatfmt='.2f')}

## Recalculation policy

Eligible interaction-linked compounds were reparsed with RDKit {audit['rdkit_version_current']}. Murcko scaffolds and physicochemical descriptors were recalculated and compared with formal stored fields. Formal fields were not overwritten. Exact discrepancy counts are retained in the JSON audit.

## Downstream denominators

- V1, V2, V3, V6 and V9: **{audit['structure_valid_interaction_linked_compounds']:,} compounds**.
- V4 and V5: **{audit['downstream_denominators']['V4_class_role']:,} eligible unique compound–protein pairs**.
- V7 PCA: deterministic stratified sample; exact n in its coordinate/source table.
- V8 UMAP: deterministic stratified sample; exact n and parameters in its JSON.
- V10 3D: deterministic stratified sample; exploratory only.
"""
(OUT / "01_data_audit/COMPOUND_CLASSIFICATION_DATA_AUDIT.md").write_text(text, encoding="utf-8")
print(text)
