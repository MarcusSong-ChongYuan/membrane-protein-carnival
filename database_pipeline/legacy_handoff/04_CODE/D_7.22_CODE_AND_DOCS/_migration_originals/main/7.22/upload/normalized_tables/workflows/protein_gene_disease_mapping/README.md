# Protein-gene-disease mapping workflow

This module builds protein-centred UniProt → Ensembl Gene → Open Targets direct disease/phenotype mappings from the immutable partner-handoff v3.1 protein release. It never reads drug, compound or pocket tables into its data flow.

Run from the repository root:

```powershell
python upload\normalized_tables\workflows\protein_gene_disease_mapping\run_pipeline.py --config upload\normalized_tables\workflows\protein_gene_disease_mapping\config\config.yaml --resume
```

Supported options are `--resume`, `--force-refresh`, `--skip-excel`, and `--validate-only`. Successful UniProt and Open Targets responses are cached per accession/target and reused. The cache, Parquet intermediates and log files are intentionally excluded from Git.

The workflow explicitly uses `enableIndirect: false`. `project_evidence_tier` is an internal mempro1 filter/ranking and is not an official Open Targets or UniProt rating. The primary disease is an internal summary; it does not discard other direct associations and does not imply causality.

