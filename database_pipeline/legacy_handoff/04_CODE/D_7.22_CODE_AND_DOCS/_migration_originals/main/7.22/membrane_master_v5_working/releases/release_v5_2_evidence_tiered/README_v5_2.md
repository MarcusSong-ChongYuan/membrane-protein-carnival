# Human Membrane Protein Master v5.2 — evidence-tiered

This release separates biological membrane class from evidence strength.

## Biological class

- **A** — integral membrane protein.
- **B** — directly membrane-inserted without a TM span, chiefly lipid-anchored.
- **C** — peripheral or stable membrane-associated protein.

## Evidence level

- **E1 / 实验证实** — direct experimental topology, lipid-anchor, or membrane-location evidence.
- **E2 / 强证据支持** — curated UniProt topology/location, reliable HPA localization, structure/family support, or multiple independent sources, without explicit direct evidence in the available record.
- **E3 / 预测候选** — prediction, single-source, by-similarity, or limited localization evidence.
- **E0 / 排除** — not retained as a membrane protein.

The strict experimental core contains 3,451 records. The recommended website default contains E1+E2 (7,904 records). E3 contains 2,652 separately searchable candidates, and E0 contains 441 internal exclusion records.
