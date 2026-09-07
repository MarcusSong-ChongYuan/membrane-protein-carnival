# MemPro expression and localization module v1

## Scope

This module separates four concepts:

1. organ system;
2. tissue expression;
3. cell-type expression;
4. subcellular localization.

Expression does not prove localization. RNA and protein measurements retain
their original units and are never averaged together.

## Evidence codes

- `PEX1`: protein-level tissue or cell-type evidence;
- `PEX2`: reserved for concordant independent expression sources;
- `PEX3`: RNA-only HPA expression evidence;
- `PEX0`: absent or unresolved expression evidence;
- `LOC1`: HPA IF with Enhanced or Supported reliability;
- `LOC2`: HPA Approved or a source GO cellular-component annotation;
- `LOC3`: HPA Uncertain or otherwise low-certainty localization.

## Output design

- `protein_tissue_expression_v1.tsv.gz`: one protein–tissue–measurement row.
- `protein_cell_type_expression_v1.tsv.gz`: one protein–cell-type–measurement row.
- `protein_subcellular_localization_v1.tsv.gz`: one protein–location–source row.
- `anatomy_cell_location_ontology_mapping_v1.tsv`: HPA terms mapped to
  UBERON, Cell Ontology, or GO.
- `expression_location_summary_v1.tsv.gz`: one row per MemPro protein.
- `expression_location_conflict_review_v1.tsv.gz`: identifier, ontology,
  RNA/protein, and localization review items.
- `human_membrane_protein_master_v6_2_expression_preview.tsv.gz`: a preview
  only. It must not become a formal release until the V6.1 frozen protein
  master is used as its base.

## Important limitation

The HPA `proteinatlas.tsv` file contains specificity-focused nTPM/nCPM and
protein-intensity subsets, not the full all-tissue or all-cell expression
matrix. Every detailed row is therefore labelled with its coverage scope.
