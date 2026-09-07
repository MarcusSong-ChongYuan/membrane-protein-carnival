# Authoritative one-command reproduction

Use `scripts/run_pipeline.ps1`. It is the authoritative workflow and supersedes the longer development command examples in `REPRODUCE.md`.

```powershell
powershell -ExecutionPolicy Bypass -File <ROOT>/scripts/run_pipeline.ps1 `
  -Root <ROOT> `
  -V62 <PATH_TO_FROZEN_V6.2> `
  -ComplexRoot <PATH_TO_COMPLEX_MODULE> `
  -Python python
```

The runner:

1. downloads the official UniProt human reference proteome with isoforms only when the frozen FASTA is absent;
2. verifies all explicit complex isoform accessions;
3. builds gene, canonical protein, isoform and evidence-target bridges;
4. adds external complex component entities without adding them to the membrane-protein count;
5. verifies reference-proteome-missing external accessions individually;
6. runs final relational QA and stops with an error unless V0.3 reports `PASS`.

Required upstream inputs are read-only:

- frozen V6.2 directory;
- complex module containing `output_v0_2/complex_target_components_v0_2.tsv` and `complex_target_master_v0_2.tsv`.

The scripts contain no machine-specific C: or D: paths. All locations are supplied as parameters.
