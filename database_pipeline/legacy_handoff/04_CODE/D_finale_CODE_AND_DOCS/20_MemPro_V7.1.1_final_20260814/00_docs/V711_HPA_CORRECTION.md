# V7.1.1 HPA semantic correction

V7.1.1 corrects the inherited display-layer field for HPA tissue mass-spectrometry records. `HPA_tissue_MS` is now `tissue_protein_MS`, never `tissue_RNA`. The original value is preserved in `legacy_v7_expression_layer_original`. Empty MS intensity remains `MISSING_OR_NOT_MEASURED_SOURCE_UNRESOLVED` and is excluded from mapped+measured denominators; it is not interpreted as zero or not detected. RNA, IHC and MS cross-field consistency is now blocking QA.
