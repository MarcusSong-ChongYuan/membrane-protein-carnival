# MemPro V6.2 — Binding Site Pair Export

Generated: 2026-08-08

## Overview

All experimentally validated small-molecule x membrane-protein binding site pairs from MemPro V6.2.

| Metric | Count |
|--------|-------|
| Total binding site pairs | 92,560 |
| Unique compounds | 14,533 |
| Unique membrane proteins | 2,770 |
| With bio-status flags | 1,339 compounds |
| Pairs with residue + PDB | 91,956 |
| Pairs with residue only (UniProt, no PDB) | 604 |
| Pairs with PDB only (BindingDB, no residue) | 0 |

## Files

| File | Content |
|------|---------|
| `01_compound_table.tsv` | 14,533 unique compounds with chemical properties, bio-status, source DBs |
| `02_protein_table.tsv` | 2,770 unique membrane proteins with class, topology, function |
| `03_binding_site_pairs.tsv` | 92,560 rows: full compound-protein binding site pair data (60 columns) |
| `04_source_summary.tsv` | Per-source database statistics |
| `05_site_type_reference.tsv` | Site type formats with examples |
| `README.md` | This file |

## Binding Site Completeness

Not all pairs have complete binding site residue data:

| Status | Count | site_type |
|--------|-------|-----------|
| Complete (residue + PDB) | 91,956 | residue_contact, complex_pocket, structure_match |
| No PDB (sequence-level only) | 604 | uniprot_curated_binding_site |
| No residue data | 0 | bindingdb_ligand_target_complex |

Even among "complete" pairs, residue coverage varies:
- **experimental_structure_residue_contact**: avg 7 residues, but 25% have only 1-2 residues (sparse)
- **experimental_complex_pocket**: avg 84 residues (full pocket)
- **experimental_structure_match**: only key affinity residues (intentionally incomplete)

Use the `residue_count` column to filter by residue coverage.

## Site Types

1. **experimental_structure_residue_contact** — Atomic contacts within 4A of ligand (PDB structures)
2. **experimental_complex_pocket** — Full binding pocket (50-200+ residues)
3. **experimental_structure_match** — Key affinity residues, inferred from similar compound
4. **uniprot_curated_binding_site** — Literature-curated, sequence positions (no PDB)
5. **bindingdb_ligand_target_complex** — PDB available but no residue annotation

See `05_site_type_reference.tsv` for format examples.

## Membrane Protein Classification

- **Class A**: Multi-pass transmembrane (GPCRs, channels, transporters)
- **Class B**: Lipid-anchored or peripheral membrane
- **Class C**: Single-pass transmembrane

## Column Reference

### Compound identity (01_compound_table.tsv)
- `compound_internal_id`: MemPro internal compound ID
- `preferred_name`: Primary compound name
- `standard_smiles`: Canonical SMILES string
- `standard_inchikey`: InChI Key (standardized)
- `molecular_formula`, `molecular_weight`, `xlogp`, `tpsa`: Computed properties
- `pubchem_cids`, `chembl_ids`, `drugcentral_ids`, `chebi_ids`: Cross-reference IDs
- `is_approved_drug` .. `is_chemical_probe`: Biological status flags (0/1)
- `biological_status_flags`: Active flags as semicolon-separated string
- `compound_scope_class`: M5 classification level
- `identity_confidence`: high / medium / low

### Protein identity (02_protein_table.tsv)
- `target_uniprot_id`: UniProt accession
- `approved_symbol`: HGNC gene symbol
- `membrane_class_v52`: A (multi-pass) / B (lipid-anchored) / C (single-pass)
- `membrane_topology`: e.g. multi_pass, single_pass_type_i, lipid_anchor
- `functional_primary_class`: receptor / enzyme / transporter / ion_channel / ...

### Binding site pairs (03_binding_site_pairs.tsv)
- `site_type`: One of the 5 site types (see reference)
- `residue_or_site_description`: Binding site residues in source-specific format
- `residue_count`: Number of binding site residues (0 = none; filter this column to exclude incomplete sites)
- `pdb_ids_site` / `pdb_ids_evidence`: PDB structure IDs
- `pdb_chain_ids_v60`: Chain identifier(s) for the binding site
- `standard_value_nM`: Binding affinity in nanomolar
- `activity_outcome_v60`: positive_or_observed / negative / ...
- `pubmed_ids` / `doi`: Literature references

## Source Databases

- **PDBe**: 81,103 pairs, 13,630 compounds, 2,666 proteins, 27,330 PDBs
- **PDBbind**: 8,465 pairs, 4,604 compounds, 479 proteins, 8,379 PDBs
- **sc-PDB**: 1,948 pairs, 1,488 compounds, 185 proteins, 1,927 PDBs
- **UniProt**: 604 pairs, 101 compounds, 432 proteins, 0 PDBs
- **BioLiP**: 440 pairs, 376 compounds, 292 proteins, 377 PDBs
