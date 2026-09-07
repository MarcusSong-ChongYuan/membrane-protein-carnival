# Refined compound-classification trial contract

Core conclusion: A multi-axis structural classification should distinguish membrane-target role profiles more effectively than the former acyclic/single-ring/multi-ring/macrocycle scheme without treating projection clusters as chemical classes.

Figure archetype: quantitative comparison grid of independent candidate views.

Target/output: exploratory NAR analysis package; Python-only; editable SVG/PDF plus 600-dpi TIFF and PNG preview.

Evidence hierarchy:

1. Full-set coverage and category resolution on 240,054 structure-valid interaction-linked compounds.
2. Association strength and held-out dominant-role classification compared on the same compounds.
3. Role-profile heatmaps for biological interpretability.
4. Morgan/Butina and UMAP remain separate exploratory neighbourhood analyses, not taxonomy ground truth.

Reviewer risks:

- local SMARTS labels could be mistaken for ClassyFire/ChEBI ontology;
- dominant target role is database-coverage dependent;
- rare categories can inflate apparent specificity;
- functional groups are multi-label and cannot be counted as mutually exclusive classes;
- formal charge is structure-record charge, not a physiological-pH microspecies prediction.
