# Methods

The preserved legacy Direct TSV (1,439,756 rows) was profiled before migration. The normalized export was rebuilt from cached Open Targets target responses and the retained protein-gene mappings, with identical evidence classification and primary-disease rules. Parquet uses Zstandard compression, dictionary encoding and statistics. ID columns remain strings, boolean evidence flags are native booleans, and unavailable evidence counts remain nullable integers rather than fabricated values.

The core relation is one row per `target_uniprot_id + ensembl_gene_id + disease_id`. Protein and disease descriptions are not repeated in that table. Datasource-level scores are stored in `protein_disease_evidence_summary_v1.parquet`, keyed by `ensembl_gene_id + disease_id + datasource_id`.

No row-level evidence-detail dataset was generated because the v1 association workflow retrieved aggregated association channels, not individual evidence records. The preserved legacy TSV and workbook were not modified or deleted. Old/new equivalence status: `passed`.

<!-- ABCE_SCORING_V1_1:START -->
## v1.1 scoring method

Source channels are collapsed to independent families before counting: UniProt channels to `UNIPROT`, CRISPR channels to `CRISPR`, and EVA channels to `EVA`. Europe PMC and Clinical Precedence do not qualify for Rule B; Expression Atlas and IMPC cannot pass without a core family.

- Rule A (`low` when A-only): at least two qualifying raw datasource IDs.
- Rule B (`medium`, default inclusion): direct disease, at least two independent qualifying families, at least one core family, known source and not weak-only.
- Rule C (`high`): Rule B plus a clinical/curated genetic family and a different eligible core/human-genetic family.
- Rule E (`very_high`): Rule C plus the original strong multi-source predicate and corrected direct/disease/specificity/known/non-weak constraints.

Counts: A=63,269, B=6,978, C=3,860, E=2,384. Formal levels: low=56,291, medium=3,118, high=1,476, very_high=2,384. Single high-authority sources remain outside Rule B and require raw evidence validation.
<!-- ABCE_SCORING_V1_1:END -->
