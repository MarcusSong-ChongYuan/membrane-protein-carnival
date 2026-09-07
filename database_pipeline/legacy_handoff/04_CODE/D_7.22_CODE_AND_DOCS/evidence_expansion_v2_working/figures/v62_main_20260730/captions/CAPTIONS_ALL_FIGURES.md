# MemPro V6.2 figure captions

# Main figure captions — M1–M4

## Figure M1. Database architecture, integration and quality control

Sources are grouped by the data module to which they contribute. The pipeline
distinguishes record-level identity resolution from pair-level aggregation.
The incremental positive review input comprised 3,290,514 records; the negative
identity input comprised 10,246,948 records. Counts in the release snapshot are
from the frozen MemPro V6.2 release dated 2026-07-30. Heatmap values are
log10(count + 1), with raw counts printed in cells.

## Figure M2. Landscape of the human membrane proteome

All panels use one row per UniProt accession (n = 10,997). Membrane
classes A, B and C are biological scope categories; E1–E3 are evidence grades
and should not be interpreted as the same variable. Structure grade separates
membrane-curated experimental structures, other experimental PDB entries,
AlphaFold-only models and proteins without an assigned structure. Oligomeric
states marked unresolved or state-dependent remain explicit.

## Figure M3. Tissue, cell-type and subcellular localization atlas

HPA tissue and cell-type values were log1p transformed and summarized by
functional class. Heatmap colors are within-class z-scores and therefore show
relative patterns rather than absolute expression. The anatomical silhouette is
a navigation schematic; quantitative counts are shown in the adjacent bars.
RNA–IHC concordance is evaluated as detection breadth per protein, not as a
claim of assay equivalence.

## Figure M4. Disease-association landscape

Disease analyses use unique protein–disease pairs. Organ systems are a
transparent keyword-derived grouping of disease names for visualization; they
are not prevalence estimates or a replacement for the original MONDO/OMIM
identifiers. The network is filtered to high-confidence, high-degree nodes for
readability. Expression–disease concordance is descriptive Jaccard overlap and
does not imply tissue causality.

# Main figure captions — M5–M8

## Figure M5. Chemical identity and small-molecule space

Canonical parents and exact forms are separate entities. The chemical-space
embedding uses Morgan fingerprints (radius 2, 1024 bits) and a fixed random seed
on a deterministic sample of QC-passed core compounds. UMAP is descriptive:
distances should not be interpreted as an absolute chemical-similarity scale.
Physicochemical ECDFs retain outliers within the displayed ranges.

## Figure M6. Binding evidence and binding-site landscape

Evidence records are grouped by source and mutually exclusive highest evidence
tier. Independent-source counts and quantitative potency use unique
protein–canonical compound pairs. Standardized nM values are displayed on a log
scale but activity types are not assumed to be thermodynamically equivalent.
Residue composition is calculated from parsable residue mentions in
coordinate-supported site descriptions.

## Figure M7. Negative evidence, conflicts and data reliability

Positive and negative coverage is reported at the unique-pair level. Conflict
rate is the number of positive pairs with at least one negative record divided
by positive pairs within each functional class. Conflicts may reflect assay or
context differences and are flagged rather than automatically adjudicated.
Missing, not-applicable and review states remain explicit.

## Figure M8. Structure-aware docking prioritization

The docking funnel begins with 1,502,456 positive pairs. R0 is a protocol
validation/redocking tier and is not a prospective discovery tier. D1–D3 are
ranked prospective or mechanistic sets. Candidate PDB identifiers indicate
triage readiness only; chain, state, biological assembly, binding box and ligand
microstates still require preparation before HPC submission. Docking scores
must not be interpreted as binding affinities without protocol validation.

# Supplementary figure captions

## Figure S1 | Cross-database integration and source support
Independent-source coverage, common source signatures, source-count distributions and disease/compound provenance. Counts describe V6.2 records and do not imply equal source sensitivity or independence.

## Figure S2 | Membrane topology, formal classification and assembly
Topology and transmembrane-segment distributions are shown together with single-pass resolution, GPCRdb/IUPHAR/TCDB alignment and V6.2 oligomeric/biological-assembly evidence. Unresolved and conflict states remain explicit.

## Figure S3 | Complete expression and subcellular-localization detail
Clustered high-information tissue and cell-type summaries are paired with localization roles, expression breadth by membrane class and mapping/conflict quality controls. HPA expression values are contextual annotations and are not used as membrane-protein inclusion evidence by themselves.

## Figure S4 | Disease ontology, evidence provenance and network sensitivity
Disease namespaces, evidence-level composition, high-degree diseases, evidence-channel prevalence and degree distributions summarize the protein–gene–disease module. Organ-system or ontology aggregation reflects database annotation density, not disease prevalence.

## Figure S5 | Compound identity, structure standardization and QC
Identity-confidence and scope states, parent/form multiplicity, source-record support, physicochemical space and biological status summarize canonical compound standardization. Distribution panels use a deterministic stratified sample for rendering; all categorical totals use the complete 2,016,064-row master table.

## Figure S6 | Docking shortlist robustness and HPC workload controls
The V6.2 standard docking shortlist is decomposed by scientific tier, receptor/ligand/box readiness, per-target workload caps and ranking-score distributions. These panels describe prioritization and preparation readiness, not docking outcomes or predicted affinity.

