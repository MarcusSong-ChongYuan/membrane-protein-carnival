# Evidence policy v5.2

The fields `membrane_class_v52` and `evidence_level_v52` answer different questions. A/B/C describes how a protein associates with membrane; E1/E2/E3/E0 describes confidence.

## E1

E1 requires an explicit ECO:0000269 experimental qualifier attached to a UniProt TM/intramembrane feature, lipid-anchor statement, or membrane/membrane-organelle localization. Former-R records previously marked `confirmed` also map to E1.

## E2

E2 includes curated UniProt TM/intramembrane/lipid-anchor assignments without explicit direct evidence in the exported text; HPA Enhanced/Supported plasma-membrane localization; OPM/Membranome support; at least two independent external sources; and curated UniProt membrane-location statements without computational/by-similarity evidence codes. Former-R `probable` records map to E2.

## E3

E3 contains single-source predictions, by-similarity statements, HPA Approved-only localization, and unresolved signal-anchor/signal-peptide cases. It is not included in default statistics.

## E0

E0 contains secreted/signal-peptide false positives, non-reproduced predictions, target-chain structural mismatches, and other excluded records. It is retained only for audit and re-import prevention.
