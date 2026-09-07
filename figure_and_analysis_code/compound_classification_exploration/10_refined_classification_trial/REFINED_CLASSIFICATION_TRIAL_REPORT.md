# Refined compound-classification trial results

## Direct answer

Yes. The refined ring topology is demonstrably more informative than the old acyclic/monocyclic/other-polycyclic/macrocyclic scheme, but ring topology alone remains insufficient. The strongest practical representation is a multi-axis annotation in which ring topology, functional groups, physicochemical bins, charge, Murcko scaffold and target-specific interaction type remain separate fields.

## Full-set comparison (n=240,054)

- Old ring classification: 4 nominal categories, effective category number 1.29; largest category 94.7%; Cramer's V versus dominant target role 0.0746.
- Refined ring topology: 7 categories, effective category number 3.27; largest category 52.5%; Cramer's V 0.0855.
- Therefore the refined topology reduces collapse into one broad polycyclic category and modestly improves role association.

## Held-out defined-role comparison

The diagnostic classifier used 239,182 compounds in six defined roles (role n>=500); 872 unresolved or rare-role compounds were excluded only from this inferential comparison and remain in the database.

- Old ring: balanced accuracy 0.191, macro-F1 0.008.
- Refined ring: balanced accuracy 0.220, macro-F1 0.060.
- Best balanced accuracy: Functional groups, 0.378.
- Functional groups: balanced accuracy 0.378.

This demonstrates that the new scheme is better, but the largest gain comes from orthogonal multi-label chemistry rather than creating more mutually exclusive ring labels.

## Recommended database fields

1. `ring_topology_primary` plus multi-label ring flags.
2. exact and generic Bemis-Murcko scaffold IDs.
3. multi-label functional groups.
4. molecular-weight, LogP, TPSA, flexibility and Fsp3 bins, while retaining continuous values.
5. structure-record charge class, explicitly not physiological-pH charge.
6. recognizable scaffold motifs as local annotations, not ClassyFire/ChEBI ontology.
7. pharmacological action kept at the protein-compound interaction level.

## Important limitations

- The named motif list is curated but not an exhaustive chemical ontology.
- Formal charge reflects the stored structure; it does not predict protonation at pH 7.4.
- Dominant target role depends on current database coverage.
- Functional-group enrichment is descriptive; it does not establish causal activity determinants.
