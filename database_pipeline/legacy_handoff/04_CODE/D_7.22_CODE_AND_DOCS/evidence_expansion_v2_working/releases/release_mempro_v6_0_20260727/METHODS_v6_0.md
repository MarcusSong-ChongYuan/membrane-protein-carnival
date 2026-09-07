# Methods

## Incremental design

Published V5.5/V5.4/V4.2/V1.0 files were used as read-only baselines. New
records were normalized in isolated source-specific staging directories and
merged only after source QA passed.

## Identity normalization

Proteins were required to map to one canonical UniProt accession in the A/B/C
membrane-protein universe. Small molecules were mapped to the V1.0 parent and
form hierarchy by a unique PubChem CID, standard InChIKey, curated PDB chemical
component mapping or validated structural identity. Names alone were not used
to merge compounds.

## Evidence processing

PubChem BioAssay retained active, inactive, inconclusive and unspecified
outcomes. Inactive rows were separated from positive evidence. GPCRdb-derived
records retained the original source database and GPCRdb as the originating
curation layer. PDBe contacts and PDBbind affinity complexes retained PDB,
chain and residue-coordinate provenance.

## Deduplication and conflicts

Stable source evidence identifiers remove exact repeated records. An exact
target–compound–activity–value–assay fingerprint flags likely cross-source
measurement mirrors; the provenance row is preserved but excluded from default
pair counting. Same-assay active/inactive BioAssay contexts and cross-source
quantitative ranges of at least 100-fold are reported for review. Values are
never averaged across assay contexts.

## Release summaries

Protein–compound summaries are calculated only from default, uniquely mapped,
non-mirror evidence. Protein and compound masters are then refreshed from the
pair table. Disease relation content is carried forward unchanged from V5.4
and receives a V6 provenance field.
