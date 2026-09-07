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
