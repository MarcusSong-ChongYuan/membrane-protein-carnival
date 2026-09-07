# Reproduce the V6.3 candidate identity layer

All commands are examples. Replace `<ROOT>`, `<V62>` and `<COMPLEX>` with local paths.

```powershell
python <ROOT>/scripts/fetch_explicit_isoforms.py `
  --complex-components <COMPLEX>/output_v0_2/complex_target_components_v0_2.tsv `
  --output-fasta <ROOT>/raw/uniprot/complex_explicit_isoforms.fasta `
  --output-manifest <ROOT>/raw/uniprot/complex_explicit_isoforms_manifest.json
```

```powershell
python <ROOT>/scripts/build_identity_layer.py `
  --protein-master <V62>/human_membrane_protein_master_v6_2.tsv `
  --binding-evidence <V62>/binding_evidence_master_v6_2.tsv.gz `
  --binding-sites <V62>/binding_site_instances_v6_2.tsv.gz `
  --disease-relations <V62>/protein_gene_disease_relations_v6_2.tsv `
  --complex-components <COMPLEX>/output_v0_2/complex_target_components_v0_2.tsv `
  --uniprot-fasta <ROOT>/raw/uniprot/UP000005640_human_with_isoforms.fasta.gz `
  --uniprot-headers <ROOT>/raw/uniprot/response_headers.txt `
  --output-dir <ROOT>/output_v0_1
```

The checked-in compatibility runner preserves the V0.2 builder logic while adding the canonical-column alias required by the V0.1 table schema.

```powershell
python <ROOT>/scripts/run_extend_complex_identity_v2_fixed.py [the arguments documented by --help]
```

```powershell
python <ROOT>/scripts/fetch_missing_external_canonical.py `
  --external-entities <ROOT>/output_v0_1/external_complex_protein_entity_v0_1.tsv.gz `
  --output-fasta <ROOT>/raw/uniprot/external_canonical_direct.fasta `
  --output-manifest <ROOT>/raw/uniprot/external_canonical_direct_manifest.json
```

```powershell
python <ROOT>/scripts/finalize_identity_layer_v3.py `
  --root <ROOT> `
  --direct-fasta <ROOT>/raw/uniprot/external_canonical_direct.fasta `
  --direct-manifest <ROOT>/raw/uniprot/external_canonical_direct_manifest.json
```

Successful reproduction ends with `qa/IDENTITY_LAYER_V0_3_VALIDATION.json` reporting `"status": "PASS"`.
