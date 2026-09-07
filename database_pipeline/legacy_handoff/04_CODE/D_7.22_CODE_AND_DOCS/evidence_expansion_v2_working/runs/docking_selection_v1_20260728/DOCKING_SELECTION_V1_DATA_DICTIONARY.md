# Docking selection V1 (V6.0 preview)

This is a reproducible docking-oriented partition of the frozen MemPro V6.0
protein-compound summary. It is a preview and must be refreshed after the V6.1
database freeze.

## Mutually exclusive tiers

| Tier | Intended use | Pair count |
|---|---|---:|
| R0 | Redocking and protocol validation using a pair-specific experimental PDB site | 5,062 |
| D1 | Primary prospective docking: membrane-curated experimental structure, direct binding <=1 uM, plus multi-source support or a priority compound | 19,018 |
| D2 | Secondary prospective docking: very strong direct binding with any experimental PDB, or <=10 uM with a membrane-curated structure, plus support | 15,347 |
| D3 | Exploratory/mechanistic docking for approved drugs, clinical candidates, probes, and supported endogenous ligands | 9,727 |
| X | Retained in the database but excluded from default docking | 526,884 |

R0 is a validation set, not a prospective discovery set. D1-D3 contain
44,092 prospective pairs. Including all
R0 controls produces 49,154
recommended rows.

## Eligibility gates

Protein: evidence E1/E2, membrane class A/B/C, and website-default inclusion.

Compound: canonical structure and InChIKey present; core scope; QC=ok; molecular
weight 100-700; heavy atoms 6-60; absolute formal charge <=2; rotatable bonds
<=20; not inorganic or metal-containing.

## Structure grades

- S1: pair-specific experimental binding site.
- S2: membrane-curated experimental structure from OPM/PDBTM.
- S3: another experimental PDB structure.
- S4: AlphaFold-only model.
- S0: no structure currently assigned.

## Important interpretation

Docking scores are not binding affinities. Each selected row still requires
receptor-state selection, ligand microstate preparation, binding-box definition,
and protocol validation before HPC execution.
