# MemProDB evidence expansion v2 working area

## Purpose

This working area extends the frozen v5.5 protein, v4.2 binding-evidence and
v1.0 compound releases with additional experimental small-molecule evidence.
Existing releases are read-only inputs and are never overwritten.

## Scientific scope

- Species: Homo sapiens.
- Protein universe: A/B/C records in the frozen membrane-protein audit master.
- Default evidence requires an unambiguous human target and a canonical core
  small-molecule mapping.
- Protein, antibody, nucleic-acid and cell-therapy records remain out of the
  small-molecule default layer.
- Source rows, exclusions, negative observations and ambiguous mappings are
  retained in audit tables.

## Priority sources

1. PubChem BioAssay target-linked concise bioactivity.
2. BRENDA human membrane-enzyme ligands and kinetics.
3. GPCRdb and PDSP Ki ligand bioactivity not already represented by ChEMBL or
   Guide to PHARMACOLOGY.
4. Current PDB/PDBe protein-ligand complexes and binding-site contacts, with
   refreshed structural specialist-source cross-checks.

## Directory policy

- `raw/`: immutable source downloads and API responses.
- `intermediate/`: normalized source-specific tables.
- `reports/`: mapping, overlap, exclusion and validation reports.
- `config/`: frozen paths, schemas and evidence policy.
- `releases/`: only validated integrated releases.
- `scripts/`: reproducible acquisition and normalization code.

