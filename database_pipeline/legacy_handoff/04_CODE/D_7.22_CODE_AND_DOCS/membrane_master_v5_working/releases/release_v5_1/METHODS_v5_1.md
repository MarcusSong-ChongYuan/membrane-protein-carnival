# Methods v5.1

## Scope

The v5 Tier R queue contained 1,109 reviewed human accessions assembled from HPA location/prediction, UniTmp high-throughput topology, OPM, and PDBTM. v5.1 audited all 1,109 without changing the immutable v5 release.

## Current evidence refresh

Current UniProt records were retrieved on 2026-07-24. TMbed 2022 human-proteome predictions were accepted only when their stored amino-acid sequence exactly matched the current UniProt sequence. GOA Human 2026-07-08 and GO basic 2026-06-26 supplied evidence-coded cellular-component annotations. InterPro 2026-06-10 supplied domain context.

## Chain-level structural correction

All PDB identifiers attached to former Tier R records were mapped to UniProt accessions through PDBe/SIFTS and compared with PDBTM chain topology. No reviewed target accession in Tier R mapped to a PDBTM chain that itself contained a TM segment. Therefore, PDBTM membership alone was not used to promote any R record; in many cases the membrane chain was a receptor or partner in the same complex.

## Decision process

An initial deterministic evidence matrix integrated current UniProt features, TMbed, hydropathy, GOA, HPA, UniTmp, OPM, InterPro, and target-chain mapping. Refined rules rejected short TM calls, non-reproduced external predictions, and signal-peptide confusion. Nineteen targeted overrides document high-impact edge cases such as retroviral products, lipid anchors, and false signal-peptide assignments. Every record retains the rule, evidence basis, cautions, source links, and review date.
