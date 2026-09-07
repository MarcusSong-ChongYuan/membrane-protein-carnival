# MemPro docking HPC handoff

The manifest is a ranked computational worklist, not evidence that docking has been run.

1. Select one biologically appropriate receptor conformation and assembly for every target.
2. Remove irrelevant crystallization additives while retaining required cofactors and structural ions.
3. Assign protonation, missing side chains/residues and membrane-aware orientation; record every choice.
4. Generate ligand stereoisomers/protomers only under an explicit policy and keep the canonical/form IDs.
5. Define the box from a co-crystallized ligand or curated binding-site residues. Do not blind-dock by default.
6. Redock the reference ligand. Production is allowed only when heavy-atom RMSD is at most 2.0 Å and the pose recovers key contacts.
7. Run the pilot set, inspect score distributions and failure modes, then unlock the standard set.

`docking_array_map_v631.tsv` uses stable zero-based `array_index` values and batches of 250 pairs. Cluster-specific resource requests and module names must be filled in by the HPC operator.
