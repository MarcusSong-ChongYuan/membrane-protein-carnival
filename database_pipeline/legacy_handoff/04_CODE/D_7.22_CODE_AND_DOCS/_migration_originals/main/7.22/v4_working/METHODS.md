# v4 Methods

## Reference universe

The reference universe was built from reviewed human entries in UniProt reference
proteome `UP000005640`, retrieved on 2026-07-24. Canonical UniProt accessions are
the counting unit. Unreviewed transmembrane candidates were retrieved separately
and placed in a canonical-resolution review queue.

Membrane scope is assigned from UniProt transmembrane and intramembrane features,
beta-barrel annotations, lipid-anchor annotations, subcellular-location text,
keywords, and cellular-component GO annotations. Integral, lipid-anchored, and
peripheral membrane-associated proteins remain distinct.

## Existing evidence migration

The immutable v3.1 relationship release was separated into activity measurements,
curated target assertions, structure-ligand observations, and UniProt binding-site
annotations. Original v3.1 relationship identifiers were retained. Compound scope
was re-evaluated so that biologics and large molecules do not contribute to core
small-molecule coverage; ions and metals are retained separately.

## External source integration

IUPHAR/BPS Guide to PHARMACOLOGY version 2026.2 interaction, ligand, target, and
mapping files were downloaded and filtered to human targets in the reviewed
membrane reference. BindingDB 2026-07 and ChEMBL 37 target indexes were intersected
with the reference set. Index-only hits are discovery candidates and are not
promoted to final quantitative evidence levels until measurement records are
imported.

## Evidence levels

- SM0: no qualifying detailed evidence in the integrated release.
- SM1: structure-context or binding-site annotation.
- SM2: functional or bioactivity assay.
- SM3: curated pharmacological target assertion.
- SM4: direct quantitative Ki or Kd evidence.
- SM5: the same protein-compound direct evidence supported by at least two
  independent source labels.

The protein level is the strongest qualifying detailed evidence. An index-only hit
is stored separately from the final SM level.

## Reproducibility

Raw downloads, source metadata, scripts, build reports, release manifests, and
SHA-256 checksums are retained. The v3.1 baseline is not modified.

