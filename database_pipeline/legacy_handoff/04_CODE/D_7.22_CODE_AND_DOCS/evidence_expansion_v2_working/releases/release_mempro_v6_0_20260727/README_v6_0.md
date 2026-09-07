# MemProDB V6.0

MemProDB V6.0 is an incremental supplement to the frozen V5.5 protein,
V5.4 protein–gene–disease, V4.2 binding-evidence and V1.0 small-molecule
baselines. The baselines were never overwritten.

The release contains four linked canonical tables:

1. `human_membrane_protein_master_v6_0.tsv`
2. `protein_gene_disease_relations_v6_0.tsv`
3. `binding_evidence_master_v6_0.tsv.gz`
4. `small_molecule_master_v1_1.tsv`

Supporting tables preserve protein–compound pairs, compound forms, binding
sites, negative evidence, identity-review records, source overlap and
quantitative heterogeneity.

Only records with a canonical A/B/C membrane protein and a uniquely mapped
canonical small molecule can enter the high-confidence binding layer.
Unresolved identities are retained in the review queue. Inactive BioAssay
results are retained in the negative-evidence table.

This candidate becomes a formal release only after primary-key, foreign-key,
row-reconciliation, baseline-integrity and checksum validation all pass.
