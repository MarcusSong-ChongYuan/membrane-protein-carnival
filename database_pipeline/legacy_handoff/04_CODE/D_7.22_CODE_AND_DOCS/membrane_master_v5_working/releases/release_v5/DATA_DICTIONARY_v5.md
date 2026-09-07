# Data dictionary (v5 additions)

| Column | Meaning |
|---|---|
| `release_tier_v5` | A integral; B monotopic/lipid-anchored; C peripheral; R review |
| `record_status_v5` | Retained v4 record or new reviewed external-union record |
| `membrane_decision_v5` | Machine-readable reason for tier assignment |
| `membrane_confidence_v5` | High, moderate or review |
| `membrane_scope_v5` | Final v5 scope; preserves the original v4 scope in `membrane_scope` |
| `membrane_topology_v5` | Final v5 topology after conservative single-pass resolution |
| `transmembrane_count_v5` | Final v5 TM count used for browsing/QC |
| `independent_membrane_sources_v5` | Semicolon-separated source evidence present for the record |
| `hpa_*_v5` | HPA v25.1 prediction/localization fields; localization is not insertion proof |
| `htp_*_v5` | UniTmp HTP evidence class, TM count, reliability and terminal sides |
| `membranome_*_v5` | Membranome single-pass family and helix segment |
| `opm_*_v5` | OPM structure and integral-segment status |
| `pdbtm_*_v5` | PDB IDs that overlap PDBTM membrane structures |
| `single_pass_resolution_v5` | Exact, partial or unresolved single-pass topology result |
| `functional_primary_class_v5` | Main website browsing class |
| `functional_subclass_v5` | Specialist hierarchy or family-level fallback |
| `classification_source_v5` | Source/rule that produced the v5 class |
| `classification_status_v5` | Specialist aligned, rule classified, family-level only or unclassified |
| `gpcrdb_*_v5` | GPCRdb class/family/ligand type/subfamily |
| `gtopdb_*_v5` | GtoPdb target type and family |
| `tcdb_tcids_v5` | Semicolon-separated TC classification identifiers |
| `cross_source_conflict_v5` | 1 when the record needs conflict/review attention |
| `review_reason_v5` | Explicit conflict or review reason |
