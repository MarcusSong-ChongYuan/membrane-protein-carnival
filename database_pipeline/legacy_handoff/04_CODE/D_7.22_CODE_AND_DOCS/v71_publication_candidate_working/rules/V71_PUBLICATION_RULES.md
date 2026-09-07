# MemPro V7.1 publication rules

V7 is immutable input. V7.1 is a clean-room incremental publication candidate.

## Public admission requirements

Every public evidence row must satisfy all of the following:

1. Canonical target identity is a single human UniProt protein, or an explicitly modelled protein complex.
2. Protein belongs to the confirmed membrane-protein layer (E1 or E2) and has an A/B/C membrane class.
3. Compound identity maps to a unique canonical core compound. Parent/form identity and stereochemistry are retained.
4. Source database, source version and source record identifier are specific and non-empty.
5. Evidence has a verified PubMed ID, DOI, or an active PDB entry with primary citation.
6. Evidence semantics are assigned to exactly one primary class: P0 direct binding, P1 functional pharmacology, or P2 structural/site observation.
7. The record is not an exact mirrored experiment or mirrored structure. Contributing databases remain in a crosswalk.
8. Public fields are allowed by source license/terms and required attribution is retained.
9. Record-level QC passes and all public foreign keys resolve.

## P0 direct binding

Accepted measurements include Kd, Ki and explicitly defined direct binding constants (Ka/Kb/Ke) from a direct biochemical or biophysical binding assay. IC50 and EC50 are never automatically treated as direct binding constants.

## P1 functional pharmacology

Includes IC50, EC50, agonism, antagonism, inhibition, activation and confirmatory cell-based or functional assays. Primary screening-only activity is excluded.

## P2 structural/site observation

Requires an active PDB or expert-curated site record, correct chain-to-UniProt mapping and a non-artifactual ligand/site. Water, common buffers, solvents, crystallization additives and ligands contacting only companion chains are excluded.

## Publication identity

`verified_peer_reviewed` requires a valid PubMed ID or DOI. `verified_structure_publication` requires an active PDB and its primary citation. Database-record-only evidence is retained outside the manuscript core layer.

## Nonpublic disposition

Failed or unresolved records are never silently deleted. They are frozen with a machine-readable reason and are excluded from website, API, downloads, figures and manuscript statistics.

