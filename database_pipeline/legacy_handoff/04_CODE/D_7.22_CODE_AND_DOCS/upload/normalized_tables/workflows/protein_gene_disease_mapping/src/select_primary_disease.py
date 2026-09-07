from __future__ import annotations

from typing import Any


TIER_RANK = {"Tier A": 0, "Tier B": 1, "Tier C": 2, "Tier D": 3}
ENTITY_RANK = {"disease": 0, "phenotype": 1, "biological_measurement": 2, "other": 3}


def primary_candidate(row: dict[str, Any], config: dict[str, Any]) -> bool:
    primary = config["primary_disease"]
    if primary.get("exclude_measurements", True) and row.get("entity_type") == "biological_measurement":
        return False
    if primary.get("exclude_ontology_roots", True) and row.get("specificity_status") != "specific":
        return False
    if primary.get("exclude_literature_only_from_primary", True):
        contexts = set(str(row.get("relation_context", "")).split(";"))
        if contexts and contexts <= {"literature_association"}:
            return False
    return True


def primary_sort_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        TIER_RANK.get(str(row.get("project_evidence_tier")), 9),
        ENTITY_RANK.get(str(row.get("entity_type")), 9),
        0 if str(row.get("clinical_genetic_flag", "")).lower() == "true" else 1,
        0 if str(row.get("human_genetic_flag", "")).lower() == "true" else 1,
        0 if row.get("specificity_status") == "specific" else 1,
        -int(row.get("datasource_count") or 0),
        -float(row.get("ot_overall_score") or 0),
        -int(row.get("evidence_count") or 0),
        str(row.get("disease_id", "")),
    )


def selection_reason(row: dict[str, Any]) -> str:
    return (
        f"selected by tier={row.get('project_evidence_tier')}; entity_type={row.get('entity_type')}; "
        f"clinical_genetic={row.get('clinical_genetic_flag')}; human_genetic={row.get('human_genetic_flag')}; "
        f"specificity={row.get('specificity_status')}; datasource_count={row.get('datasource_count')}; "
        f"ot_overall_score={row.get('ot_overall_score')}; stable disease_id tie-break"
    )

