# MemPro V6.4.1 final publication triage

This release retains the V6.4 formal layer, performs one conservative rescue pass, and freezes all remaining records as nonpublic R/X archive.

- Public evidence: 1,046,677
- Public unique protein-compound pairs: 588,047
- Newly promoted evidence: 54,830
- Frozen nonpublic evidence: 1,956,629

`public_binding_evidence_master_v6_4_1.tsv.gz` is the only evidence table intended for website/figure/manuscript counts.
`public_protein_compound_pairs_v6_4_1.tsv.gz` is its deduplicated pair table.
`frozen_R_X_archive_v6_4_1.tsv.gz` is internal-only and must not be exposed through the website/API or counted in publication results. It preserves reasons so decisions remain auditable.

Promotion was deliberately conservative: unique protein and canonical compound identities, current default membrane-protein status, nonduplicate evidence, acceptable QC, and a PubMed/DOI/PDB anchor were required.
