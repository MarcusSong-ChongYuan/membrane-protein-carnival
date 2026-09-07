# MemPro V6.1 review-queue resolution plan

V6.0 is an immutable input. All new files are written under this V6.1 run directory.

1. Build unique PubChem CID, PDBe CCD, PDBbind complex, and BRENDA name inventories.
2. Resolve PDBbind ligand structures and PDBe CCD identities; classify solvents,
   buffers, ions, polymers, covalent ligands, and genuine small molecules.
3. Fetch official PubChem structures for active BE2 CIDs first with rate limiting,
   retry, checkpoints, and an explicit terminal status for every requested CID.
4. Process other active PubChem CIDs, then lower-priority records and high-value
   BRENDA names.
5. Match exact full structures first, audit connectivity-level parent/form matches,
   and never merge by name alone.
6. Rebuild compound/form identities and evidence foreign keys in a V6.1 candidate
   release. Preserve unresolved records with reason codes.
7. Reconcile all source rows, validate keys and hashes, then freeze V6.1 only when
   every required quality-control test passes.
