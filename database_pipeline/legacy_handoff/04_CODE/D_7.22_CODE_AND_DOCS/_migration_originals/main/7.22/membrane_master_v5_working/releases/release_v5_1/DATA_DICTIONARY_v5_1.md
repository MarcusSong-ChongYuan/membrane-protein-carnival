# Data dictionary additions in v5.1

| Field | Meaning |
|---|---|
| `candidate_membrane_class_v51` | Best biological hypothesis (A/B/C), including evidence-limited candidates. |
| `final_membrane_class_v51` | Release-level biological class A/B/C, or `unknown` if evidence is insufficient/excluded. |
| `evidence_status_v51` | `confirmed`, `probable`, `uncertain`, or `excluded`. |
| `release_disposition_v51` | `included`, `candidate`, or `excluded`; controls publication-ready membership. |
| `decision_code_v51` | Machine-readable reason for the final disposition. |
| `decision_basis_v51` | Human-readable evidence summary. |
| `cautions_v51` | Conflicts and limitations that should accompany the record. |
| `v51_decision_layer` | Inherited v5, refined deterministic rule, or targeted manual override. |
| `manual_review_flag_v51` | `1` for targeted edge-case override, otherwise `0`. |

Legacy v5 fields are preserved for provenance. `release_tier_v5=R` is historical workflow state only and must not be interpreted as a biological class in v5.1.
