# MemPro M1-M8 scientific figure revision — completion report

Generated: 2026-08-06T09:03:03.450148+00:00

## Scope and data protection

- All work was written under `D:\finale\07_M1-M8图表修订_20260806`.
- Frozen V6.2 and V6.3.1 source releases were read only; no release table was overwritten.
- Eight figures were exported as PNG (600 dpi), editable SVG, and PDF.
- The presentation preserves the approved 18-slide two-page-per-module template and includes speaker notes on every slide.

## Completed scientific corrections

1. **Panel provenance:** 48 M1-M8 panels now have a source, field/definition, version policy, statistical unit, and limitation record.
2. **M2 Spearman scope:** A-class, E1/E2, canonical multipass proteins only; 2-40 TM segments and 50-5000 aa; n=2,844; rho=0.445; bootstrap 95% CI 0.410-0.481 (seed 42).
3. **M3 provenance:** HPA 25.1 is explicitly separated into tissue RNA nTPM, cell-type RNA nCPM, IHC protein staining, and mapping/missingness layers.
4. **M4 anatomy:** keyword bins were replaced by ontology-derived, multi-label disease-anatomy summaries using available MONDO/DO/Uberon mappings.
5. **M5C classes:** all 899,591 core-QC compounds were reclassified by one deterministic, mutually exclusive descriptor rule set.
6. **PAGTN experiment:** 2,400 compounds, seed 42, 64-dimensional property-supervised PAGTN embeddings, HDBSCAN clustering, and a Morgan radius-2/1024-bit baseline. PAGTN improved non-noise silhouette (0.360 vs 0.148), while Morgan retained higher 10-neighbor scaffold purity (11.9% vs 9.7%); no claim of universal superiority is made.
7. **AlphaFold:** 10,908/10,908 API metadata requests succeeded; local site pLDDT mapped for 857/859 attempted site targets. pLDDT is presented as model confidence, not experimental validation.
8. **Membrane-side sites:** 9,036 residue-contact sites were classified relative to protein topology: 773 intramembrane, 1,331 interface/mixed, 407 cytoplasmic, 384 non-cytoplasmic, 3 both sides/mixed, and 6,138 extramembrane-side unresolved.
9. **PPT/QA:** 18 slides rendered; no slide overflow was detected by the bundled presentation validator.

## Interpretation boundaries

- PAGTN is the encoder; HDBSCAN is the clustering algorithm. The current encoder is trained on four physicochemical descriptors and is a method experiment, not a universal pretrained chemical embedding.
- The membrane-side classifier is protein-topology-relative. It does not yet calculate ligand/residue atom z-coordinates in an OPM membrane frame. Non-cytoplasmic includes extracellular and organelle-lumen sides.
- AlphaFold coverage and pLDDT cannot replace experimental structures or prove a binding pocket is biologically correct.
- M4 anatomical associations describe ontology mappings, not tissue causality.
- Database source counts are distinct contributing database counts unless an experiment/structure lineage key explicitly supports independence.
