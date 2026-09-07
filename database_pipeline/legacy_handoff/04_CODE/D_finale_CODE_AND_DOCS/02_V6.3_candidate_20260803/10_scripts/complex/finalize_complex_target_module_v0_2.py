#!/usr/bin/env python3
"""Harden the V0.1 complex module with release-tier and compound-ID controls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


TRUE_VALUES = {"1", "true", "yes", "y"}
ID_PATTERNS = {
    "ChEBI": re.compile(r"CHEBI:\s*([0-9]+)", re.I),
    "ChEMBL": re.compile(r"\b(CHEMBL[0-9]+)\b", re.I),
    "PubChem CID": re.compile(r"(?:PUBCHEM(?:\s+CID)?|CID)\s*[:=]?\s*([0-9]+)", re.I),
}


def clean(value: object | None) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in {"-", "None", "NULL", "nan"} else text


def truth(value: object | None) -> bool:
    return clean(value).lower() in TRUE_VALUES


def stable_id(prefix: str, value: str) -> str:
    return f"{prefix}-{hashlib.sha256(value.encode('utf-8')).hexdigest()[:20].upper()}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def load_protein_status(path: Path) -> dict[str, dict[str, object]]:
    output: dict[str, dict[str, object]] = {}
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = clean(row.get("target_uniprot_id")).upper()
            if not accession:
                continue
            output[accession] = {
                "website_default": truth(row.get("website_default_v52")),
                "evidence_level": clean(row.get("evidence_level_v52")),
                "release_tier": clean(row.get("release_tier_v5")),
            }
    return output


def load_compound_maps(path: Path) -> tuple[dict[str, dict[str, set[str]]], dict[str, str]]:
    maps: dict[str, dict[str, set[str]]] = {
        "ChEBI": defaultdict(set),
        "ChEMBL": defaultdict(set),
        "PubChem CID": defaultdict(set),
    }
    names: dict[str, str] = {}
    usecols = [
        "compound_internal_id",
        "preferred_name",
        "chebi_ids",
        "chembl_ids",
        "pubchem_cids",
    ]
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        chunksize=150_000,
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    ):
        for row in chunk.itertuples(index=False):
            compound_id = clean(row.compound_internal_id)
            names[compound_id] = clean(row.preferred_name)
            for value in re.split(r"[;|,]", clean(row.chebi_ids)):
                match = re.search(r"([0-9]+)", value)
                if match:
                    maps["ChEBI"][match.group(1)].add(compound_id)
            for value in re.split(r"[;|,]", clean(row.chembl_ids)):
                match = re.search(r"CHEMBL[0-9]+", value, flags=re.I)
                if match:
                    maps["ChEMBL"][match.group().upper()].add(compound_id)
            for value in re.split(r"[;|,]", clean(row.pubchem_cids)):
                match = re.search(r"([0-9]+)", value)
                if match:
                    maps["PubChem CID"][match.group(1)].add(compound_id)
    return maps, names


def extract_identifier(raw: str) -> tuple[str, str]:
    matches: list[tuple[str, str]] = []
    for id_type, pattern in ID_PATTERNS.items():
        for match in pattern.finditer(raw):
            identifier = match.group(1).upper() if id_type == "ChEMBL" else match.group(1)
            matches.append((id_type, identifier))
    unique = sorted(set(matches))
    if len(unique) == 1:
        return unique[0]
    if len(unique) > 1:
        return "multiple", ";".join(f"{kind}:{identifier}" for kind, identifier in unique)
    if re.search(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?\b", raw):
        return "UniProt", ""
    return "", ""


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.casefold())


def source_name(raw: str) -> str:
    return clean(re.sub(r"\([^)]*(?:CHEBI|CHEMBL|CID|PUBCHEM)[^)]*\)", "", raw, flags=re.I))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-master", type=Path, required=True)
    parser.add_argument("--compound-master", type=Path, required=True)
    parser.add_argument("--v01", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    protein_status = load_protein_status(args.protein_master)
    components = read_rows(args.v01 / "complex_target_components_v0_1.tsv")
    master = read_rows(args.v01 / "complex_target_master_v0_1.tsv")
    bindings = read_rows(args.v01 / "complex_binding_evidence_candidates_v0_1.tsv")
    crosswalk = read_rows(args.v01 / "complex_source_crosswalk_v0_1.tsv")
    source_entities = read_rows(args.v01 / "complex_source_entities_v0_1.tsv")
    identity_review = read_rows(args.v01 / "complex_identity_review_v0_1.tsv")
    structure_candidates = read_rows(
        args.v01 / "complex_structure_evidence_candidates_v0_1.tsv"
    )

    default_members: dict[str, set[str]] = defaultdict(set)
    all_members: dict[str, set[str]] = defaultdict(set)
    for row in components:
        accession = clean(row.get("component_uniprot_id")).upper()
        status = protein_status.get(accession, {})
        row["protein_membrane_evidence_level"] = status.get("evidence_level", "")
        row["protein_website_default_flag"] = int(bool(status.get("website_default")))
        if accession:
            all_members[row["complex_target_id"]].add(accession)
        if status.get("website_default"):
            default_members[row["complex_target_id"]].add(accession)

    for row in master:
        target_id = row["complex_target_id"]
        default_ids = sorted(default_members.get(target_id, set()))
        row["default_membrane_component_ids"] = ";".join(default_ids)
        row["default_membrane_component_count"] = len(default_ids)
        curated = int(clean(row.get("curated_source_count")) or 0) > 0
        conflict = truth(row.get("component_set_conflict_flag"))
        default_inclusion = curated and bool(default_ids) and not conflict
        row["default_release_inclusion"] = int(default_inclusion)
        if conflict:
            row["review_status"] = "component_conflict_review"
        elif not curated:
            row["review_status"] = "predicted_only_review"
        elif not default_ids:
            row["review_status"] = "no_E1_E2_membrane_component_review"
        else:
            row["review_status"] = "curated_candidate"
        row["module_version"] = "complex_target_v0.2_v63_candidate"

    compound_maps, compound_names = load_compound_maps(args.compound_master)
    refined_bindings: list[dict[str, object]] = []
    duplicate_keys: set[tuple[str, str, str, str, str]] = set()
    for row in bindings:
        raw = clean(row.get("compound_name_or_id_raw"))
        id_type, identifier = extract_identifier(raw)
        mapped_ids: set[str] = set()
        candidate_entity_type = "ambiguous_name_or_free_text"
        mapping_status = "unmapped_name_or_free_text"
        if id_type in compound_maps:
            candidate_entity_type = "small_molecule_or_chemical_entity"
            mapped_ids = compound_maps[id_type].get(identifier, set())
            if len(mapped_ids) == 1:
                mapping_status = "exact_source_identifier_unique"
            elif len(mapped_ids) > 1:
                mapping_status = "source_identifier_maps_multiple_compounds"
            else:
                mapping_status = "source_identifier_not_in_v62_compound_master"
        elif id_type == "multiple":
            candidate_entity_type = "multiple_identifiers_review"
            mapping_status = "multiple_source_identifiers_not_auto_mapped"
        elif id_type == "UniProt":
            candidate_entity_type = "protein_ligand_not_small_molecule"
            mapping_status = "excluded_from_small_molecule_mapping"

        compound_id = next(iter(mapped_ids)) if len(mapped_ids) == 1 else ""
        preferred_name = compound_names.get(compound_id, "")
        raw_name = source_name(raw)
        raw_norm = normalize_name(raw_name)
        preferred_norm = normalize_name(preferred_name)
        if not compound_id or not raw_norm or not preferred_norm:
            name_consistency = "not_assessed"
        elif raw_norm == preferred_norm or raw_norm in preferred_norm or preferred_norm in raw_norm:
            name_consistency = "compatible_string"
        else:
            name_consistency = "requires_synonym_or_source_error_review"

        scope = clean(row.get("evidence_scope"))
        directness = (
            "complex_context_relation_directness_unverified"
            if scope == "complex_level_annotation"
            else "component_context_not_complex_specific"
        )
        key = (
            row["complex_target_id"],
            row["source_database"],
            row["evidence_type"],
            compound_id,
            raw,
        )
        if key in duplicate_keys:
            continue
        duplicate_keys.add(key)
        row["source_compound_id_type"] = id_type
        row["source_compound_id"] = identifier
        row["candidate_entity_type"] = candidate_entity_type
        row["compound_internal_id"] = compound_id
        row["mapped_compound_preferred_name"] = preferred_name
        row["compound_mapping_status"] = mapping_status
        row["source_name_identifier_consistency"] = name_consistency
        row["complex_relation_directness_status"] = directness
        row["default_release_inclusion"] = 0
        row["module_version"] = "complex_target_v0.2_v63_candidate"
        refined_bindings.append(row)

    master_fields = list(master[0])
    component_fields = list(components[0])
    binding_fields = list(refined_bindings[0]) if refined_bindings else list(bindings[0])
    write_rows(args.output / "complex_target_master_v0_2.tsv", master, master_fields)
    write_rows(
        args.output / "complex_target_components_v0_2.tsv", components, component_fields
    )
    write_rows(
        args.output / "complex_binding_evidence_candidates_v0_2.tsv",
        refined_bindings,
        binding_fields,
    )
    for filename, rows in [
        ("complex_source_crosswalk_v0_2.tsv", crosswalk),
        ("complex_source_entities_v0_2.tsv", source_entities),
        ("complex_identity_review_v0_2.tsv", identity_review),
        ("complex_structure_evidence_candidates_v0_2.tsv", structure_candidates),
    ]:
        write_rows(args.output / filename, rows, list(rows[0]) if rows else [])

    formal_binding_fields = [
        "complex_binding_evidence_id",
        "complex_target_id",
        "compound_internal_id",
        "compound_form_id",
        "evidence_type",
        "evidence_directness",
        "activity_type",
        "standard_value_nM",
        "binding_subunit_uniprot_id",
        "pdb_id",
        "biological_assembly_id",
        "binding_site_residues",
        "source_database",
        "source_record_id",
        "pubmed_ids",
        "record_qc_status",
        "default_release_inclusion",
    ]
    write_rows(
        args.output / "complex_binding_evidence_v0_2.tsv", [], formal_binding_fields
    )

    master_ids = {row["complex_target_id"] for row in master}
    blocking = {
        "duplicate_complex_ids": len(master_ids) != len(master),
        "component_fk_failures": sum(
            row["complex_target_id"] not in master_ids for row in components
        ),
        "binding_candidate_fk_failures": sum(
            row["complex_target_id"] not in master_ids for row in refined_bindings
        ),
        "predicted_only_default_release": sum(
            truth(row["default_release_inclusion"])
            for row in master
            if int(clean(row.get("curated_source_count")) or 0) == 0
        ),
        "complex_binding_candidates_in_default_release": sum(
            truth(row["default_release_inclusion"]) for row in refined_bindings
        ),
    }
    status = "PASS" if not any(bool(value) for value in blocking.values()) else "FAIL"
    mapping_counts = Counter(row["compound_mapping_status"] for row in refined_bindings)
    entity_counts = Counter(row["candidate_entity_type"] for row in refined_bindings)
    report = {
        "status": status,
        "module_version": "complex_target_v0.2_v63_candidate",
        "v62_release_modified": False,
        "counts": {
            "complex_targets": len(master),
            "default_curated_complex_targets_with_E1_E2_component": sum(
                truth(row["default_release_inclusion"]) for row in master
            ),
            "complexes_without_default_membrane_component": sum(
                int(clean(row["default_membrane_component_count"]) or 0) == 0
                for row in master
            ),
            "component_assertions": len(components),
            "binding_candidates_after_exact_dedup": len(refined_bindings),
            "formal_complex_binding_evidence": 0,
        },
        "binding_candidate_mapping_status": dict(mapping_counts),
        "binding_candidate_entity_types": dict(entity_counts),
        "policies": {
            "E1_E2_component_required_for_default_complex_layer": True,
            "protein_ligands_excluded_from_small_molecule_mapping": True,
            "name_only_compounds_auto_mapped": False,
            "identifier_mapped_annotations_promoted_to_direct_binding": False,
            "formal_complex_binding_table_empty_until_directness_verified": True,
        },
        "blocking_errors": blocking,
    }
    output_files = sorted(args.output.glob("*.tsv"))
    report["outputs"] = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in output_files
    }
    (args.output / "COMPLEX_TARGET_MODULE_V0_2_VALIDATION.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
