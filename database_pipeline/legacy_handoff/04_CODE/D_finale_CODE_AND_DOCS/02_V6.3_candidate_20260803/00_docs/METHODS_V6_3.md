# MemPro V6.3 incremental methods

## Protein identity and complex targets

The V6.2 canonical UniProt accession remains the protein primary key. Stable HGNC, Ensembl Gene and NCBI Gene identifiers define the gene layer. Exact UniProt isoform identifiers and sequences define the isoform layer; isoforms are never inferred from a gene symbol or a canonical-only record. Protein complexes are separate entities with component assertions that can reference canonical proteins, explicit isoforms, external context proteins or non-protein components.

## Disease identity

Source disease identifiers are retained verbatim. A canonical disease identity is assigned only for a direct MONDO identifier or a unique official exact/equivalent cross-reference. Name equality, broad/narrow/related mappings, obsolete targets and one-to-many mappings never trigger automatic identity merging. Open Targets Platform 26.06 supplies multi-label therapeutic areas. Disease-anatomy mappings are derived from asserted or inherited DO/MONDO axioms pointing to Uberon, with provenance and inference state retained.

## Five-axis protein classification

All 10,997 reviewed proteins were mapped against UniProtKB release 2026_02. Direct GO annotations were propagated only through the frozen GO 2026-06-15 `is_a` graph to controlled molecular-function and biological-process categories. Structural-family labels use UniProt families and then InterPro/Pfam fallbacks. Reactome release 97 provides pathway hierarchy. Membrane-role assignment uses conservative specialist, GO, existing curated class and selected family rules. Family annotations do not force a specific function.

## Evidence provenance

The legacy source count is audited as a literal distinct contributing database count. Separate hash keys represent structure evidence, chain/site-granular structure assertions, PubChem assays, and literature experiment proxies. Evidence modality is separately controlled. These counts are not interchangeable and are not summed into a claim of independent experiments.

## Docking prioritization

V6.2 shortlist membership and tier are frozen. V6.3 adds a capped score term for structure, PubChem, literature and modality lineages and subtracts a fixed penalty for an explicit positive-negative conflict. Protein classifications are joined for interpretation only. The score prioritizes preparation and computation; it is not a binding-affinity prediction.

## Chemical diversity

All core, QC-passed canonical compounds were parsed with RDKit 2026.03.5. Ring-containing molecules were reduced to Bemis-Murcko frameworks. Acyclic, missing and invalid SMILES were counted separately. No random subset, PCA or UMAP embedding was used for M5E.

## Quality control

Each module preserves input row counts and validates primary/foreign keys, mapping policies and release membership. A machine-selected lineage audit universe is retained, while a fixed source/tier-stratified sample is explicitly distinguished from full manual review. Candidate freezing requires every module gate and figure/workbook validation to pass.
