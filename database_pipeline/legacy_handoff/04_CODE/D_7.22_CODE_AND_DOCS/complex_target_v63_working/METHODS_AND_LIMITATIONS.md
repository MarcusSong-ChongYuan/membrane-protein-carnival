# Methods and remaining work

## Sources frozen for this preview

- Complex Portal human ComplexTab, curated snapshot dated 2026-01-14.
- Complex Portal human predicted ComplexTab, snapshot dated 2026-01-14.
- CORUM 5.3 human complexes, release 2026-04-14.
- HuMemLigDB V6.2 protein, compound and PDBe assembly tables.

Only complexes containing at least one accession from the V6.2 membrane-protein
audit are retained. Predicted-only complexes are never included in the default
layer. Curated complexes additionally require at least one E1/E2 V6.2 membrane
component for the V0.2 default candidate layer.

## What is intentionally not inferred

- Required versus optional subunits.
- A physiological stoichiometry from a PDB asymmetric unit.
- A complex identity from equal component sets alone.
- A complex-level ligand relation from a component-level drug annotation.
- Direct complex binding from a ligand, agonist or antagonist annotation.
- A small-molecule identity from name similarity alone.
- A binding subunit, binding chain or binding-site residue without explicit
  structural mapping.

## Next gates before formal V6.3 inclusion

1. Resolve the 425 cross-source identity review groups using source cross-
   references, stoichiometry, function and literature.
2. Add Complex Portal detail/PSI-MI component features and CORUM partial-complex
   flags where they explicitly distinguish required, optional or incomplete
   membership.
3. Map PDB entity and chain identifiers to UniProt through PDBe/SIFTS, then select
   the supported biological assembly.
4. Review the 190 explicit compound-ID matches for name/identifier consistency
   and relation directness.
5. Parse CORUM formal drug-complex relations separately from component drug lists.
6. Join verified complex-level relations to V6.2 compounds and create non-empty
   formal `complex_binding_evidence` and `complex_binding_site` tables.
7. Manually audit high-value GPCR heteromers, ion-channel assemblies,
   transporters and the complexes entering the docking shortlist.

