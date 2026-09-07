# Methods for v5.4 Protein–Disease and Binding Integration

## Inputs

- Human membrane-protein audit master v5.3 with canonical sequences.
- v3.1 normalized protein–compound evidence and binding-site tables.
- v3 protein–gene–disease Rule-B release (Open Targets 26.06).
- UniProtKB 2026_02 disease comments and disease-database cross-references.
- IUPHAR/BPS Guide to PHARMACOLOGY 2026.2 interactions.
- BindingDB 2026-07 measurements curated directly from articles.
- ChEMBL 37 standardized single-protein binding-assay Ki/Kd measurements.

## Disease workflow

The v3 relation schema is preserved. Relations are mapped by canonical UniProt
accession to the v5.3 protein universe. Current UniProt disease comments are
retrieved in restartable API batches. Each `DISEASE:` statement becomes a
separate relation; disease name, OMIM identifier, causal wording and PubMed
identifiers are retained. Primary disease selection first preserves a v3 primary
designation; otherwise the highest-ranked current UniProt relation is selected.
Only compact disease summary fields are merged into the protein master.

## Binding workflow

The v3.1 evidence partition is migrated from activity measurements, curated
target assertions, structure-ligand observations and UniProt binding-site
annotations. Guide to PHARMACOLOGY, BindingDB and ChEMBL are then appended as
source-specific records. Each row receives a BE tier, target-assignment status,
compound-scope status and default-inclusion flag. Protein–compound and
protein-level summaries are derived only from default-included evidence.

Structure evidence, activity values and curated mechanisms remain distinct;
measurement types are never pooled into a single potency value. Original
relations and units are retained, with nM values added only when a valid
standardization is available.

## Reproducibility and validation

All source and normalized flat files are retained on the D drive. Stable
identifiers are SHA-256-derived from source keys. Validation checks cover key
uniqueness, protein foreign-key closure, full 10,997-row summary closure,
evidence-tier counts, source counts and preservation of the v5.3 protein
universe.

