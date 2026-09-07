import sys
from pathlib import Path

import pandas as pd

WORKFLOW = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKFLOW))

from src.export_abce_release_v1_1 import apply_rules
from src.normalize_disease_entities import normalize_disease
from src.select_primary_disease import primary_sort_key


def test_mondo_and_measurement_normalization():
    disease = normalize_disease({
        "disease_id": "MONDO_0005233", "disease_name": "non-small cell lung carcinoma",
        "therapeutic_areas": '[{"id":"MONDO_0045024","name":"cancer or benign tumor"}]',
        "ancestor_ids": '["MONDO_0004992"]', "is_therapeutic_area": "false",
    })
    assert disease["mondo_id"] == "MONDO:0005233"
    assert disease["pathology_category"] == "cancer"
    measurement = normalize_disease({"disease_id": "OBA_0000001", "disease_name": "measurement", "ancestor_ids": "[]"})
    assert measurement["entity_type"] == "biological_measurement"
    legacy_measurement = normalize_disease({
        "disease_id": "EFO_0004842", "disease_name": "eosinophil count",
        "therapeutic_areas": '[{"id":"EFO_0001444","name":"measurement"}]',
        "ancestor_ids": '["EFO_0001444"]',
    })
    assert legacy_measurement["entity_type"] == "biological_measurement"
    behavioral_trait = normalize_disease({
        "disease_id": "EFO_0004329", "disease_name": "alcohol drinking",
        "therapeutic_areas": '[{"id":"GO_0008150","name":"biological_process"}]',
        "ancestor_ids": '["GO_0008150"]',
    })
    assert behavioral_trait["entity_type"] == "phenotype"
    orphanet = normalize_disease({"disease_id": "Orphanet_101016", "disease_name": "rare disease", "ancestor_ids": '["MONDO_0000001"]'})
    assert orphanet["entity_type"] == "disease"


def test_primary_sort_is_stable():
    base = {"project_evidence_tier": "Tier A", "entity_type": "disease", "clinical_genetic_flag": "true",
            "human_genetic_flag": "true", "specificity_status": "specific", "datasource_count": 2,
            "ot_overall_score": 0.5, "evidence_count": 2}
    assert primary_sort_key({**base, "disease_id": "A"}) < primary_sort_key({**base, "disease_id": "B"})


def test_abce_levels_are_nested_and_rule_e_is_corrected():
    base = {
        "association_scope": "direct",
        "entity_type": "disease",
        "unknown_source_flag": False,
        "weak_sources_only_flag": False,
        "rule_a_qualifying_datasource_count": 2,
        "independent_source_family_count": 2,
        "core_source_count": 1,
        "clinical_or_curated_genetic_source_count": 0,
        "second_core_or_human_genetic_source_count": 0,
        "high_authority_source_count": 0,
        "entity_exclusion_reason": "",
    }
    frame = pd.DataFrame([
        {**base, "independent_source_family_count": 1, "core_source_count": 0},
        base,
        {**base, "clinical_or_curated_genetic_source_count": 1, "second_core_or_human_genetic_source_count": 1},
        {**base, "independent_source_family_count": 3, "core_source_count": 2,
         "clinical_or_curated_genetic_source_count": 1, "second_core_or_human_genetic_source_count": 1},
        {**base, "entity_type": "phenotype", "independent_source_family_count": 3, "core_source_count": 2,
         "clinical_or_curated_genetic_source_count": 1, "second_core_or_human_genetic_source_count": 1},
    ])
    scored = apply_rules(frame)
    assert scored["evidence_level"].tolist() == ["low", "medium", "high", "very_high", "low"]
    assert (scored["rule_e_pass"] <= scored["rule_c_pass"]).all()
    assert (scored["rule_c_pass"] <= scored["rule_b_pass"]).all()
    assert (scored["rule_b_pass"] <= scored["rule_a_pass"]).all()
