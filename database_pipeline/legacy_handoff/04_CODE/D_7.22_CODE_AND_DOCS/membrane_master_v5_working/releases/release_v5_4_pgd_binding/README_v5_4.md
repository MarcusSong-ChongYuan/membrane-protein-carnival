# Human Membrane Protein Database v5.4

This release adds two annotation layers to the frozen v5.3 membrane-protein
universe:

1. protein–gene–disease relations, following the v3 Rule-B structure and
   supplemented with Open Targets Platform and current UniProtKB disease
   annotations (including linked disease identifiers where supplied);
2. a source-level small-molecule binding-evidence database using BE1/BE2/BE3
   evidence classes.

The v5.3 canonical sequence and membrane classification fields are unchanged.

## Main files

- `human_membrane_audit_master_v5_4.tsv`: all 10,997 audited proteins with
  disease and binding summary columns.
- `protein_gene_disease_relations_v5_4.tsv`: detailed ternary relations.
- `protein_primary_disease_v5_4.tsv`: selected primary disease per supported
  protein.
- `binding_evidence_master_v4_1.tsv`: canonical full binding-evidence table.
- `binding_evidence_excel_review_v4_1.tsv`: compact review view capped per
  protein at 5 BE1, 7 BE2 and 3 BE3 observations; it does not replace the
  canonical full evidence table.
- `binding_site_instances_v4_1.tsv`: experimental structure/site child records.
- `protein_compound_summary_v4_1.tsv`: website-ready protein–compound pairs.
- `protein_binding_summary_v4_1.tsv`: one binding summary per protein.
- `small_molecule_index_v4_1.tsv`: compound cross-reference and scope index.

## Audit files

Legacy v3 records outside the v5.3 protein universe are retained in explicit
audit tables and are not silently deleted. ChEMBL and BindingDB ambiguous or
failed-QC rows remain visible in the full evidence table but do not contribute to
default website summaries.

## Evidence axes

Membrane evidence (E1–E3), disease evidence
(`medium`/`high`/`very_high`) and binding evidence (BE1–BE3) are independent.
Predicted pockets are BP annotations and do not establish a protein–compound
interaction.
