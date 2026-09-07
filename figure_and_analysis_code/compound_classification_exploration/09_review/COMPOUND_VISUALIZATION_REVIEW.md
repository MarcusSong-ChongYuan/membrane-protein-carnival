# MemPro compound-classification visualization review

This is an exploratory selection report, not the final manuscript Figure 3. All structural analyses use 240,054 structure-valid interaction-linked compounds; the quarantined identity-conflict record remains in the registry but is excluded from structure-derived analyses.

## Denominator guardrails

- Compound registry: **646,670** identities.
- Interaction-linked compounds: **240,055** identities.
- Structure-valid interaction-linked compounds: **240,054**.
- Formal protein-compound pairs: **529,168**.
- These denominators are not interchangeable.

## V1 | SUPPLEMENTARY

- **Question:** How are interaction-linked compounds distributed across broad local chemical regimes?
- **Denominator:** 240,054 structure-valid interaction-linked compounds
- **Analysis type:** No; descriptive/exploratory.
- **Strength:** Useful overview, but the regime is a deterministic local taxonomy rather than ClassyFire/ChEBI ontology..
- **Limitation:** Useful overview, but the regime is a deterministic local taxonomy rather than ClassyFire/ChEBI ontology.
- **Recommended placement:** SUPPLEMENTARY.

## V2 | MAIN CANDIDATE

- **Question:** Is scaffold use concentrated or long-tailed?
- **Denominator:** 240,054 structure-valid interaction-linked compounds
- **Analysis type:** No; descriptive/exploratory.
- **Strength:** Directly quantifies scaffold diversity and concentration without treating clusters as classes..
- **Limitation:** Directly quantifies scaffold diversity and concentration without treating clusters as classes.
- **Recommended placement:** MAIN CANDIDATE.

## V3 | MAIN CANDIDATE

- **Question:** Which recurrent Bemis-Murcko scaffolds dominate the linked chemical collection?
- **Denominator:** Top recurrent scaffolds among 240,054 compounds
- **Analysis type:** No; descriptive/exploratory.
- **Strength:** Chemically intuitive.
- **Limitation:** Chemically intuitive; use with rank-abundance rather than as a stand-alone frequency gallery.
- **Recommended placement:** MAIN CANDIDATE.

## V4 | SUPPLEMENTARY

- **Question:** Do broad chemical regimes show different membrane-role profiles?
- **Denominator:** 529,168 formal pairs; row-normalized descriptive profiles
- **Analysis type:** Descriptive association view; no causal inference.
- **Strength:** Shows preference patterns, but broad local classes may hide scaffold-level specificity..
- **Limitation:** Shows preference patterns, but broad local classes may hide scaffold-level specificity.
- **Recommended placement:** SUPPLEMENTARY.

## V5 | MAIN CANDIDATE

- **Question:** Which abundant scaffolds share membrane-target role profiles?
- **Denominator:** Top 50 generic scaffolds across 529,168 formal pairs
- **Analysis type:** Descriptive association view; no causal inference.
- **Strength:** Strong cross-module view.
- **Limitation:** Strong cross-module view; retain support counts and avoid interpreting clustering as causal taxonomy.
- **Recommended placement:** MAIN CANDIDATE.

## V6 | SUPPLEMENTARY

- **Question:** How do physicochemical properties vary across local chemical regimes?
- **Denominator:** 240,054 structure-valid interaction-linked compounds
- **Analysis type:** No; descriptive/exploratory.
- **Strength:** Compact descriptor comparison.
- **Limitation:** Compact descriptor comparison; medians are descriptive and do not establish drug-likeness.
- **Recommended placement:** SUPPLEMENTARY.

## V7 | SUPPLEMENTARY

- **Question:** How much variance in standardized descriptors is captured by a linear projection?
- **Denominator:** Deterministic stratified sample n=50,000
- **Analysis type:** No; descriptive/exploratory.
- **Strength:** Axes are interpretable through loadings, but overlap indicates coarse regimes are not cleanly separable..
- **Limitation:** Axes are interpretable through loadings, but overlap indicates coarse regimes are not cleanly separable.
- **Recommended placement:** SUPPLEMENTARY.

## V8 | EXPLORATORY ONLY

- **Question:** Do Morgan-fingerprint neighborhoods form local scaffold islands?
- **Denominator:** Deterministic stratified sample n=20,000
- **Analysis type:** No; descriptive/exploratory.
- **Strength:** Useful for local-neighborhood inspection.
- **Limitation:** Useful for local-neighborhood inspection; UMAP axes and global distances have no direct chemical meaning.
- **Recommended placement:** EXPLORATORY ONLY.

## V9 | MAIN CANDIDATE

- **Question:** How selective or polypharmacological are linked compounds by observed target breadth?
- **Denominator:** 240,054 interaction-linked compounds; formal-pair target counts
- **Analysis type:** No; descriptive/exploratory.
- **Strength:** Directly supports target-breadth narrative.
- **Limitation:** Directly supports target-breadth narrative; observed breadth remains coverage-dependent.
- **Recommended placement:** MAIN CANDIDATE.

## V10 | EXPLORATORY ONLY

- **Question:** Can users interactively inspect MW-LogP-TPSA space and target breadth?
- **Denominator:** Deterministic stratified sample n=20,000
- **Analysis type:** No; descriptive/exploratory.
- **Strength:** Interactive exploration aid.
- **Limitation:** Interactive exploration aid; perspective and overplotting make it unsuitable as primary evidence.
- **Recommended placement:** EXPLORATORY ONLY.

## Recommended Figure 3 core

The strongest nonredundant manuscript candidates are V2 (scaffold concentration), V3 (representative recurrent scaffolds), V5 (scaffold-by-membrane-role profiles), and V9 (observed target breadth). V1, V4, V6, and V7 are useful supporting context. V8 and V10 should remain exploratory because their geometry is projection- or perspective-dependent.

## Classification boundary

`chemical_regime_local` is a transparent deterministic classification derived locally from structure/descriptor rules. It must not be described as ClassyFire, ChEBI, or an experimentally validated chemical ontology. Bemis-Murcko scaffold identity and Morgan/Butina neighborhoods answer different questions and must not be merged into one 'class' field.
