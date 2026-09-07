# Inclusion policy v5.1

## Separation of biology and evidence

`candidate_membrane_class_v51` may be A, B, or C even when `final_membrane_class_v51` is `unknown`. Only records with `release_disposition_v51=included` and a final class A/B/C enter the publication-ready master.

## Evidence thresholds

- **Confirmed**: curated/experimental cellular membrane evidence, or a target chain itself demonstrated as transmembrane by structure-level mapping.
- **Probable**: concordant independent evidence, such as external topology plus a full current-sequence TM prediction, experimentally supported membrane-organelle localization for a non-integral protein, an experimental lipid anchor plus membrane localization, or OPM placement at a defined cellular membrane.
- **Uncertain**: a plausible candidate class with incomplete or ambiguous evidence, including HPA Approved localization alone, isolated hydropathy, or an N-terminal signal-anchor/signal-peptide ambiguity.
- **Excluded**: prediction not reproducible on the current sequence, HPA localization rated Uncertain with no corroboration, secreted/lumen proteins without stable cellular membrane association, or a PDB/PDBTM association belonging to another chain.

InterPro membrane-binding domains are contextual evidence and never sufficient alone. HPA `Approved` does not mean higher confidence than `Supported`; it denotes a different validation situation and is retained conservatively as uncertain when uncorroborated.
