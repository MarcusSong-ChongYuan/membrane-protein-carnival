# Protein-gene-disease direct mapping v1

This release uses a normalized storage model. Protein text is stored once in `protein_gene_mapping_v1.*`, disease text and ontology metadata once in `disease_dictionary_v1.*`, and the complete direct association set in `protein_gene_disease_direct_core_v1.parquet`. The complete data are release assets and are not embedded in Excel or intended for ordinary Git history.

## Python / pandas

```python
import pandas as pd

mapping = pd.read_parquet("protein_gene_mapping_v1.parquet")
diseases = pd.read_parquet("disease_dictionary_v1.parquet")
relations = pd.read_parquet("protein_gene_disease_direct_core_v1.parquet")

result = (
    relations
    .merge(mapping, on=["target_uniprot_id", "ensembl_gene_id"], how="left")
    .merge(diseases, on="disease_id", how="left")
)
```

## DuckDB

```sql
SELECT
    r.target_uniprot_id,
    p.protein_name,
    p.approved_symbol,
    r.ensembl_gene_id,
    r.disease_id,
    d.disease_name,
    d.pathology_category,
    r.project_evidence_tier,
    r.ot_overall_score
FROM read_parquet('protein_gene_disease_direct_core_v1.parquet') r
LEFT JOIN read_parquet('protein_gene_mapping_v1.parquet') p
USING (target_uniprot_id, ensembl_gene_id)
LEFT JOIN read_parquet('disease_dictionary_v1.parquet') d
USING (disease_id);
```

`protein_gene_disease_summary_v1.xlsx` contains summaries only. The complete core is also distributed as stable, protein-grouped TSV.GZ shard(s) under `full_data/`.

<!-- ABCE_SCORING_V1_1:START -->
## v1.1 A/B/C/E evidence scoring

The default `main` release is Rule B and above. Rule A-only relationships are `low` and archive-only; Rule B is the `medium` inclusion threshold; Rule C is `high`; corrected Rule E is `very_high`. The unified Rule B table contains 6,978 relationships and retains `rule_c_pass`, `rule_e_pass`, `highest_rule_passed`, `evidence_level` and `evidence_level_rank` so higher-confidence subsets remain visible without duplicate main-branch datasets.

Historical Rule D had an entity-type defect and is not a formal evidence level. Rule E is recomputed with direct disease scope, specificity exclusions, known sources, non-weak evidence and the Rule C hierarchy. EVA remains `EVA` and is not inferred to be ClinVar. The four levels are evidence-policy categories, not disease-causality probabilities; `medium` is the database inclusion threshold and does not mean poor evidence quality.

Full A/B/C/E scored outputs and independent A/C/E files are stored on `archive/protein-gene-disease-abce-v1.1`.
<!-- ABCE_SCORING_V1_1:END -->
