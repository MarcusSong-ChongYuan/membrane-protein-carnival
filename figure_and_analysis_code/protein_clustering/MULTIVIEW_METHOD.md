# Multi-view protein clustering trial

This trial complements, rather than replaces, the ESM-only analysis.

Input blocks are given equal total variance after within-block standardisation:

1. ESM-2 sequence representation (PCA 50);
2. membrane topology and sequence length;
3. role-agnostic functional annotation: structural family, molecular function, biological process and specialist classification;
4. ligand scaffold/tag profile plus interaction-evidence summaries.

The primary membrane-role overlay is intentionally excluded from the feature matrix. It is used only after clustering to assess interpretability. The output must therefore be described as an **integrative functional/pharmacological map**, not a purely sequence-derived embedding or an independent validation of existing functional annotations.
