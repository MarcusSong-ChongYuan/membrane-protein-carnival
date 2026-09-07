from __future__ import annotations

import json
import re
from typing import Any

from .common import strict_name


CATEGORY_RULES = [
    ("cancer", ("cancer", "neoplasm", "tumor", "carcinoma")),
    ("cardiovascular", ("cardiovascular", "heart", "vascular", "circulatory")),
    ("metabolic", ("metabolic", "diabetes", "endocrine")),
    ("neurological", ("nervous system", "neurolog", "neurodegenerative")),
    ("psychiatric", ("psychiatric", "mental", "behavioral")),
    ("immune", ("immune system", "autoimmune")),
    ("inflammatory", ("inflammatory", "inflammation")),
    ("infectious", ("infectious", "infection")),
    ("respiratory", ("respiratory", "thoracic", "lung")),
    ("gastrointestinal", ("gastrointestinal", "digestive", "intestinal")),
    ("renal", ("renal", "kidney", "urinary")),
    ("musculoskeletal", ("musculoskeletal", "skeletal", "bone", "muscle")),
    ("dermatological", ("skin", "dermatological")),
    ("ophthalmological", ("eye", "visual", "ophthalm")),
    ("reproductive", ("reproductive", "pregnancy")),
    ("hematological", ("hematological", "blood")),
    ("rare_genetic", ("rare", "genetic syndrome", "hereditary")),
]


def normalize_disease(row: dict[str, Any]) -> dict[str, Any]:
    disease_id = str(row.get("disease_id", ""))
    name = str(row.get("disease_name", ""))
    therapeutic_areas = parse_json(row.get("therapeutic_areas"), [])
    ancestor_ids = parse_json(row.get("ancestor_ids"), [])
    is_area = str(row.get("is_therapeutic_area", "")).lower() == "true"
    prefix = disease_id.split("_", 1)[0].upper()
    area_names = [str(x.get("name", "")) for x in therapeutic_areas if isinstance(x, dict)]
    normalized_areas = {strict_name(area) for area in area_names if area}
    measurement_area = any(strict_name(area) in {"measurement", "biological measurement"} for area in area_names)
    measurement_name = bool(re.search(
        r"\b(count|measurement|level|levels|concentration|ratio|percentage|percent|index|pressure|rate|volume|height|weight)$",
        strict_name(name),
    ))
    if prefix == "OBA" or measurement_area or (prefix == "EFO" and measurement_name):
        entity_type = "biological_measurement"
    elif prefix in {"HP", "MP"}:
        entity_type = "phenotype"
    elif prefix == "EFO" and normalized_areas and normalized_areas <= {"phenotype", "biological process"}:
        entity_type = "phenotype"
    elif prefix in {"MONDO", "ORPHA", "ORPHANET", "EFO", "OTAR"}:
        entity_type = "disease"
    else:
        entity_type = "other"

    category = "other"
    source = "unclassified"
    for candidate, needles in CATEGORY_RULES:
        if any(any(needle in strict_name(area) for needle in needles) for area in area_names):
            category = candidate
            source = "Open Targets therapeutic area"
            break
    if category == "other":
        normalized_name = strict_name(name)
        for candidate, needles in CATEGORY_RULES:
            if any(needle in normalized_name for needle in needles):
                category = candidate
                source = "strict disease-name rule"
                break

    root_names = {"disease", "phenotype", "biological process", "measurement", "cancer or benign tumor"}
    if entity_type == "biological_measurement":
        specificity = "biological_measurement"
    elif is_area:
        specificity = "broad_therapeutic_area"
    elif strict_name(name) in root_names or not ancestor_ids:
        specificity = "ontology_root_or_parent"
    else:
        specificity = "specific"

    if disease_id.startswith("MONDO_"):
        mondo_id = disease_id.replace("MONDO_", "MONDO:", 1)
        mondo_name = name
        mondo_status = "direct_open_targets_identifier"
    else:
        mondo_id = ""
        mondo_name = ""
        mondo_status = "not_mapped_v1"
    return {
        "entity_type": entity_type,
        "pathology_category": category,
        "pathology_category_source": source,
        "specificity_status": specificity,
        "mondo_id": mondo_id,
        "mondo_name": mondo_name,
        "mondo_mapping_status": mondo_status,
    }


def parse_json(value: Any, default: Any) -> Any:
    if isinstance(value, (list, dict)):
        return value
    if not value:
        return default
    try:
        return json.loads(str(value))
    except json.JSONDecodeError:
        return default
