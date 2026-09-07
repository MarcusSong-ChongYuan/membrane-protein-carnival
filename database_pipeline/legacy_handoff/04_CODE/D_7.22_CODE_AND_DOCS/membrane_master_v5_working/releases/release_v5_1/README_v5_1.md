# Human Membrane Protein Master v5.1

Release date: 2026-07-24

v5.1 resolves every one of the 1,109 records that v5 placed in workflow Tier R. The biological classes and evidence state are now separate:

- **A**: integral membrane protein with one or more membrane-spanning/intramembrane segments.
- **B**: directly membrane-inserted without a transmembrane span, chiefly lipid-anchored proteins.
- **C**: peripheral or stable membrane-associated protein without a membrane-spanning segment.
- **Evidence status**: confirmed, probable, uncertain, or excluded.

The publication-ready master contains 10,076 included proteins: A 5,505, B 440, and C 4,131. The audit master retains all 10,997 v5 union records. Of the former R records, 188 are included, 480 remain explicit evidence-limited candidates, and 441 are excluded from the publication-ready membrane set.

The candidate set is a completed review outcome, not an unfinished queue. It records proteins for which a candidate biological class can be proposed but current evidence does not justify a final A/B/C assignment.

## Principal files

- `human_membrane_protein_master_v5_1.tsv`: publication-ready included A/B/C set.
- `human_membrane_protein_audit_master_v5_1.tsv`: all v5 union records with v5.1 disposition.
- `human_integral_membrane_view_v5_1.tsv`, `human_lipid_anchored_membrane_view_v5_1.tsv`, `human_peripheral_membrane_view_v5_1.tsv`: class-specific views.
- `human_membrane_candidate_view_v5_1.tsv`: evidence-limited candidates.
- `human_membrane_excluded_audit_v5_1.tsv`: excluded records and reasons.
- `review/r1109_final_review_v51.tsv`: full evidence matrix for the former R queue.
