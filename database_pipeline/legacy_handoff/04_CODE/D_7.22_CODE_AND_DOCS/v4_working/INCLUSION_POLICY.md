# Human Membrane Protein and Small-Molecule Inclusion Policy

Version: membrane-small-molecule-scope-v1

## Purpose

The database universe is the human membrane proteome. Small-molecule evidence is
an annotation layer over that universe and is not a prerequisite for retaining a
membrane-protein record.

## Protein record unit

The master-table unit is one canonical UniProt accession from human reference
proteome `UP000005640`. Isoforms are retained in a child table and are not counted
as additional proteins. Reviewed and unreviewed entries are both permitted, with
their review status retained.

## Membrane scope

Core records are integral alpha-helical or beta-barrel membrane proteins.
Extended records are lipid-anchored or explicitly supported peripheral
membrane-associated proteins. Records with missing or conflicting membrane
evidence remain in a review queue and are not silently treated as verified
integral membrane proteins.

Membrane evidence must retain source, source version or retrieval date, evidence
type, and the source-specific assertion. Signal peptides alone are not sufficient
to classify a protein as transmembrane.

## Small-molecule scope

Core small molecules are defined, non-polymeric chemical entities including
drugs, clinical compounds, probes, research ligands, endogenous metabolites,
cofactors, lipids, sterols, natural products, and defined structure ligands.

Antibodies, proteins, nucleic acids, undefined polymers, biologics, buffers,
solvents, water, and undefined mixtures are excluded from the core small-molecule
set. Inorganic ions and metals are retained separately and do not contribute to
the core small-molecule coverage rate.

## Evidence interpretation

- `SM0`: no qualifying small-molecule evidence found in the currently integrated sources.
- `SM1`: annotation or structure-context evidence only.
- `SM2`: functional activity evidence.
- `SM3`: curated pharmacological target assertion.
- `SM4`: direct quantitative binding evidence.
- `SM5`: direct evidence supported by multiple independent sources.

`SM0` does not mean that no interaction exists. Functional activity does not
automatically establish direct binding. A structure-bound ligand does not
automatically establish a pharmacological target relationship.

