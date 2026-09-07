# MemPro V7.1 approved release rules

V7.0.2 is immutable. V7.1 is an incremental data-model and publication-semantics release.

## Interaction and site independence

A reliable human membrane-protein–small-molecule interaction does not require residue or coordinate data. Site status is represented independently as RESIDUE_COMPLETE, RESIDUE_PARTIAL, RESIDUE_SOURCE_ONLY, NO_RESIDUE_REPORTED, PREDICTED_POCKET, or NO_SITE_INFORMATION. Missing site coordinates never remove an otherwise valid interaction.

## Protein entity granularity

Evidence targets exactly one explicit level when supported: gene product unspecified, canonical protein, protein isoform, processed protein form, protein complex, or unresolved target. No isoform or subunit is selected from a name alone.

## Membrane classes

A = transmembrane or strongly integral; B = directly embedded or lipid/GPI anchored without a conventional transmembrane span; C = peripheral/state- or complex-mediated membrane association. Class is stored per canonical protein, isoform, processed form, and state. The most direct mechanism is primary; other mechanisms remain secondary.

## Structure and chain

Source-reported residues are immutable. Mapped and unmapped residues are separate. SIFTS residue mapping is preferred, followed by sequence-confirmed and segment-safe DBREF mapping. Multiple chains remain a candidate set unless the source or structure context selects one uniquely.

## Negative evidence

Positive evidence is the public core. Negative evidence is public-contextual only when it conflicts with or informs selectivity for a positive pair. Isolated inactive records are archived outside website, API, figure, and manuscript default statistics.

## Expression

Detection denominators use mapped and measured records only. NOT_DETECTED, UNMAPPED, NOT_MEASURED, and MISSING are distinct. RNA abundance and IHC protein staining are separate; normal and disease contexts are separate.

## Provenance

Distinct contributing databases, independent experiments, independent structures, PubChem assays, and evidence modalities are separate quantities. Database mirrors do not increase independent-experiment counts.

## Classification and disease

Five independent classification axes are retained with source/version/method: structural family, molecular function, biological process, membrane role, and specialist classification. Disease identity uses MONDO with official exact/equivalent mappings only; broad/narrow/related links do not merge identities. Anatomy is ontology-derived and multi-label.

## Manual validation

Manual sampling is WAIVED_BY_USER for this release. The release may report automated QA, but may not claim that all records were manually reviewed or quote a manual accuracy estimate.
