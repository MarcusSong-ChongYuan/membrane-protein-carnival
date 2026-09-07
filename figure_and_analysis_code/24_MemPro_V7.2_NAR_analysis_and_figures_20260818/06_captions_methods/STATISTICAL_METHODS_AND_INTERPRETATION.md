# MemPro V7.2 statistical methods and interpretation contract

## Frozen analysis population

All analyses read the frozen MemPro V7.2 release and write only to this derivative package. The primary protein denominator is 7,800 formal A/B/C proteins. E1/E2 proteins (n=7,374) form the confirmed default layer; E3 proteins (n=426) remain visible as a separate evidence stratum. The compound registry contains 646,670 canonical compounds, whereas 240,055 compounds occur in at least one formal protein–compound pair. These universes are never interchanged.

The formal interaction layer contains 529,168 unique protein–compound pairs and 942,455 public positive evidence records. Counts of source records, contributing databases and putative experiment-lineage keys are reported separately. The disease layer contains 3,739 canonical diseases and 6,753 formal protein–disease relations. Anatomy is available for 2,926 diseases; 813 diseases remain anatomy-unmapped. The expression module contains 3,046,789 records, of which 2,921,389 are mapped and measured. The 103,709 unresolved missing/not-measured records and 22,701 unmapped/ambiguous records are excluded from expression-rate denominators but reported separately.

## Missingness and denominator rules

- RNA, IHC and mass-spectrometry layers are never pooled.
- Tissue specificity uses normal-tissue HPA consensus RNA only.
- Expression proportions use mapped + measured observations only.
- Missing, not measured, unmapped and ambiguous states remain explicit.
- Disease anatomy and therapeutic-area analyses use the intersection with the formal disease master. Historical out-of-release ontology rows are excluded and counted in the audit.
- Binding-site membrane-side analyses show the 48,331 unknown-side assertions explicitly.
- Absence of residue coordinates does not invalidate an interaction relation; it limits structural interpretation only.

## Statistical tests

Categorical associations use Pearson's chi-square test. Effect size is reported as bias-corrected Cramer's V. Cell-level enrichment is reported as standardized Pearson residuals; where individual Fisher tests are used, P values are adjusted by Benjamini–Hochberg false-discovery rate (BH-FDR). The role-by-molecular-function association has bias-corrected Cramer's V=0.619 (n=7,800); role-by-biological-process V=0.336. Structure class by membrane role is weak overall (V=0.046) and is not described as a strong preference.

Tissue specificity is summarized using Tau calculated across 51 HPA consensus tissue-RNA categories. Tau ranges from 0 (broad) to 1 (tissue-specific). Eleven membrane-role groups with at least 30 proteins are compared using Kruskal–Wallis H=457.89, P=4.34×10^-92. Pairwise follow-up uses two-sided Mann–Whitney U tests with BH-FDR. The unit is a protein, not an individual tissue measurement.

Evidence concentration is summarized by a Lorenz curve and Gini coefficient (Gini=0.361). Source-set overlap uses pair-level Jaccard similarity. The source-intersection view is an UpSet representation of unique formal protein–compound pairs.

Disease–expression anatomical concordance is the fraction of ontology-mapped disease anatomy labels for which the linked protein is detected in the corresponding mapped, measured normal-tissue RNA category. The null distribution uses 200 degree-preserving permutations that retain protein degree, disease degree, disease anatomy-label count and the measured-tissue frame. Observed mean concordance is 0.4866 versus null mean 0.3273 (SD=0.00361), z=44.18, empirical P=0.004975. This supports anatomical consistency, not disease causality.

Chemical structures are parsed from V7.2 standard SMILES. RDKit is used to derive physicochemical properties, Bemis–Murcko scaffolds and Morgan fingerprints (radius 2, 1,024 bits). Two of 240,055 interaction-linked structures fail parsing and remain in the audit. Butina clustering is run on a deterministic bounded sample of 4,000 compounds at Tanimoto similarity 0.65; it yields 3,546 clusters and is exploratory. UMAP uses a deterministic status-enriched sample of 12,000 compounds with seed 42 and is supplementary only. UMAP coordinates have no direct chemical interpretation.

The hypothesis-generation subset requires high disease support, low-but-nonzero chemical coverage and at least one coordinate-ready site. It contains 58 proteins under the recorded thresholds (disease-support upper quartile=4; nonzero chemical-pair lower quartile=2). It is a transparent filter, not a validated ranking, causal claim or predicted binding result.

## Reproducibility

All random procedures use seed 42. Every quantitative panel has a TSV source-data file in `03_source_data`. Derived analysis tables are in `02_analysis_data`; frozen source tables are never modified. Vector exports retain editable text where supported. The package manifest records SHA-256 hashes.
