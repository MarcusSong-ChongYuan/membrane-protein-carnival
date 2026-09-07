#!/usr/bin/env python3
"""Build a conservative complex-target entity module for HuMemLigDB V6.3.

The frozen V6.2 protein table is used only as the membrane-protein reference.
Complexes remain separate entities and are never inserted into that protein table.
Cross-source records are merged only through an explicit identity cross-reference;
component-set similarity is emitted for review instead of being treated as identity.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


UNIPROT_PATTERN = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$"
)
PDB_XREF_PATTERN = re.compile(r"wwpdb:([0-9A-Za-z]{4})\(([^)]+)\)", re.I)
CORUM_XREF_PATTERN = re.compile(r"corum:([0-9]+)\(([^)]+)\)", re.I)


def now_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_id(prefix: str, value: str, length: int = 20) -> str:
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:length].upper()
    return f"{prefix}-{digest}"


def clean(value: object | None) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in {"-", "None", "NULL", "nan"} else text


def split_values(value: str, separators: str = r"[|;]") -> list[str]:
    return [item.strip() for item in re.split(separators, clean(value)) if item.strip()]


def normalize_uniprot(value: str) -> tuple[str, str]:
    identifier = clean(value).replace("uniprotkb:", "").upper()
    if not UNIPROT_PATTERN.fullmatch(identifier):
        return "", ""
    if "-" in identifier:
        return identifier.split("-", 1)[0], identifier
    return identifier, ""


def parse_participant_token(token: str) -> dict[str, object]:
    token = clean(token)
    match = re.match(r"^(.+?)\(([0-9]+)\)$", token)
    identifier = clean(match.group(1) if match else token)
    count = int(match.group(2)) if match else 0
    base_uniprot, isoform = normalize_uniprot(identifier)
    lowered = identifier.lower()
    if base_uniprot:
        entity_type = "protein"
    elif lowered.startswith(("cpx-", "complex portal:")):
        entity_type = "subcomplex"
    elif lowered.startswith(("chebi:", "pubchem:", "chembl:")):
        entity_type = "small_molecule"
    elif lowered.startswith(("rnacentral:", "urs")):
        entity_type = "rna"
    else:
        entity_type = "other"
    return {
        "component_identifier": identifier,
        "component_entity_type": entity_type,
        "component_uniprot_id": base_uniprot,
        "component_uniprot_isoform_id": isoform,
        "copy_count": count,
        "stoichiometry_status": "known" if count > 0 else "unknown",
    }


def read_tsv(path: Path) -> Iterable[dict[str, str]]:
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and reader.fieldnames[0].startswith("#"):
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        yield from reader


def read_gzip_tsv(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def load_membrane_reference(path: Path) -> tuple[set[str], dict[str, str]]:
    identifiers: set[str] = set()
    names: dict[str, str] = {}
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = clean(row.get("target_uniprot_id")).upper()
            if accession:
                identifiers.add(accession)
                names[accession] = clean(row.get("approved_symbol")) or clean(
                    row.get("protein_name")
                )
    return identifiers, names


def parse_complex_portal(
    path: Path,
    dataset_class: str,
    membrane_ids: set[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    entities: list[dict[str, object]] = []
    components: list[dict[str, object]] = []
    binding_candidates: list[dict[str, object]] = []
    for row in read_tsv(path):
        source_id = clean(row.get("Complex ac"))
        direct_tokens = split_values(
            row.get("Identifiers (and stoichiometry) of molecules in complex", "")
        )
        direct = [parse_participant_token(token) for token in direct_tokens]
        has_nested = any(item["component_entity_type"] == "subcomplex" for item in direct)
        assertions = [(item, "direct") for item in direct]
        if has_nested:
            direct_proteins = {
                item["component_uniprot_id"]
                for item in direct
                if item["component_uniprot_id"]
            }
            for token in split_values(row.get("Expanded participant list", "")):
                item = parse_participant_token(token)
                if item["component_uniprot_id"] and item["component_uniprot_id"] not in direct_proteins:
                    assertions.append((item, "expanded_from_nested_complex"))
        protein_ids = sorted(
            {
                str(item["component_uniprot_id"])
                for item, relation in assertions
                if item["component_uniprot_id"]
                and relation in {"direct", "expanded_from_nested_complex"}
            }
        )
        membrane_members = sorted(set(protein_ids) & membrane_ids)
        if not membrane_members:
            continue
        xrefs = clean(row.get("Cross references"))
        corum_identity = sorted(
            {
                match.group(1)
                for match in CORUM_XREF_PATTERN.finditer(xrefs)
                if match.group(2).lower() in {"identity", "complex-primary"}
            }
        )
        source_entity_id = f"CP:{source_id}"
        entity = {
            "source_entity_id": source_entity_id,
            "source_database": "Complex Portal",
            "source_complex_id": source_id,
            "source_dataset_class": dataset_class,
            "recommended_name": clean(row.get("Recommended name")) or source_id,
            "aliases": clean(row.get("Aliases for complex")),
            "organism_taxon_id": clean(row.get("Taxonomy identifier")),
            "complex_assembly": clean(row.get("Complex assembly")),
            "description": clean(row.get("Description")),
            "go_annotations": clean(row.get("Go Annotations")),
            "evidence_code": clean(row.get("Evidence Code")),
            "experimental_evidence": clean(row.get("Experimental evidence")),
            "cross_references": xrefs,
            "pubmed_ids": ";".join(
                sorted(set(re.findall(r"pubmed:([0-9]+)", xrefs, flags=re.I)))
            ),
            "corum_identity_ids": ";".join(corum_identity),
            "protein_component_count": len(protein_ids),
            "membrane_component_count": len(membrane_members),
            "membrane_component_ids": ";".join(membrane_members),
            "component_set_signature": ";".join(protein_ids),
            "curation_status": "predicted" if dataset_class == "predicted" else "curated",
            "default_release_inclusion": 0 if dataset_class == "predicted" else 1,
            "source_file": path.name,
        }
        entities.append(entity)
        for index, (item, relation) in enumerate(assertions, start=1):
            component = dict(item)
            component.update(
                {
                    "component_assertion_id": stable_id(
                        "HMCX-COMP",
                        f"{source_entity_id}|{index}|{item['component_identifier']}|{relation}",
                    ),
                    "source_entity_id": source_entity_id,
                    "source_database": "Complex Portal",
                    "source_complex_id": source_id,
                    "component_index": index,
                    "component_relation": relation,
                    "required_subunit_status": "unknown_not_asserted_by_source",
                    "optional_subunit_status": "unknown_not_asserted_by_source",
                    "membrane_component_flag": int(
                        str(item["component_uniprot_id"]) in membrane_ids
                    ),
                }
            )
            components.append(component)
        for evidence_type, field in [
            ("ligand_annotation", "Ligand"),
            ("agonist_annotation", "Agonist"),
            ("antagonist_annotation", "Antagonist"),
        ]:
            raw_value = clean(row.get(field))
            if not raw_value:
                continue
            for token in split_values(raw_value):
                binding_candidates.append(
                    {
                        "complex_binding_candidate_id": stable_id(
                            "HMCX-BIND", f"{source_entity_id}|{evidence_type}|{token}"
                        ),
                        "source_entity_id": source_entity_id,
                        "source_database": "Complex Portal",
                        "source_complex_id": source_id,
                        "evidence_scope": "complex_level_annotation",
                        "evidence_type": evidence_type,
                        "compound_name_or_id_raw": token,
                        "compound_internal_id": "",
                        "compound_mapping_status": "unmapped_raw_annotation",
                        "complex_specificity_status": "source_asserted_complex_context",
                        "pubmed_ids": entity["pubmed_ids"],
                        "default_release_inclusion": 0,
                        "review_reason": "compound_identity_and_direct_binding_not_yet_verified",
                    }
                )
    return entities, components, binding_candidates


def parse_corum(
    path: Path,
    membrane_ids: set[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    entities: list[dict[str, object]] = []
    components: list[dict[str, object]] = []
    binding_candidates: list[dict[str, object]] = []
    for row in read_tsv(path):
        source_id = clean(row.get("complex_id"))
        accessions = split_values(row.get("subunits_uniprot_id", ""), r";")
        gene_names = split_values(row.get("subunits_gene_name", ""), r";")
        stoichiometry = split_values(row.get("subunits_stoechiometrie", ""), r";")
        parsed = [parse_participant_token(token) for token in accessions]
        protein_ids = sorted(
            {str(item["component_uniprot_id"]) for item in parsed if item["component_uniprot_id"]}
        )
        membrane_members = sorted(set(protein_ids) & membrane_ids)
        if not membrane_members:
            continue
        source_entity_id = f"CORUM:{source_id}"
        entity = {
            "source_entity_id": source_entity_id,
            "source_database": "CORUM",
            "source_complex_id": source_id,
            "source_dataset_class": "curated",
            "recommended_name": clean(row.get("complex_name")) or f"CORUM complex {source_id}",
            "aliases": clean(row.get("synonyms")),
            "organism_taxon_id": "9606",
            "complex_assembly": clean(row.get("subunits_stoechiometrie")),
            "description": clean(row.get("comment_complex")),
            "go_annotations": ";".join(
                value
                for value in [clean(row.get("functions_go_id")), clean(row.get("functions_go_name"))]
                if value
            ),
            "evidence_code": clean(row.get("functions_evi")),
            "experimental_evidence": clean(row.get("purification_methods")),
            "cross_references": "",
            "pubmed_ids": ";".join(
                sorted(
                    set(
                        split_values(row.get("pmid", ""), r";")
                        + split_values(row.get("functions_pmid", ""), r";")
                    )
                )
            ),
            "corum_identity_ids": source_id,
            "protein_component_count": len(protein_ids),
            "membrane_component_count": len(membrane_members),
            "membrane_component_ids": ";".join(membrane_members),
            "component_set_signature": ";".join(protein_ids),
            "curation_status": "curated",
            "default_release_inclusion": 1,
            "source_file": path.name,
        }
        entities.append(entity)
        for index, item in enumerate(parsed, start=1):
            count = 0
            if index <= len(stoichiometry):
                match = re.search(r"[0-9]+", stoichiometry[index - 1])
                count = int(match.group()) if match else 0
            item["copy_count"] = count
            item["stoichiometry_status"] = "known" if count > 0 else "unknown"
            component = dict(item)
            component.update(
                {
                    "component_assertion_id": stable_id(
                        "HMCX-COMP", f"{source_entity_id}|{index}|{item['component_identifier']}"
                    ),
                    "source_entity_id": source_entity_id,
                    "source_database": "CORUM",
                    "source_complex_id": source_id,
                    "component_index": index,
                    "component_relation": "direct",
                    "component_gene_symbol": gene_names[index - 1] if index <= len(gene_names) else "",
                    "required_subunit_status": "unknown_not_asserted_by_source",
                    "optional_subunit_status": "unknown_not_asserted_by_source",
                    "membrane_component_flag": int(
                        str(item["component_uniprot_id"]) in membrane_ids
                    ),
                }
            )
            components.append(component)
        drug_fields = [
            ("complex_drug_comment", "comment_drug", "complex_level_annotation"),
            ("formal_complex_drug_relation", "comment_drug_formal", "complex_level_annotation"),
            ("component_drug_annotation", "subunits_drugs", "component_context_only"),
        ]
        for evidence_type, field, scope in drug_fields:
            raw_value = clean(row.get(field))
            if not raw_value:
                continue
            tokens = split_values(raw_value, r";") if field == "subunits_drugs" else [raw_value]
            for token in tokens:
                binding_candidates.append(
                    {
                        "complex_binding_candidate_id": stable_id(
                            "HMCX-BIND", f"{source_entity_id}|{evidence_type}|{token}"
                        ),
                        "source_entity_id": source_entity_id,
                        "source_database": "CORUM",
                        "source_complex_id": source_id,
                        "evidence_scope": scope,
                        "evidence_type": evidence_type,
                        "compound_name_or_id_raw": token,
                        "compound_internal_id": "",
                        "compound_mapping_status": "unmapped_name_or_free_text",
                        "complex_specificity_status": (
                            "source_asserted_complex_context"
                            if scope == "complex_level_annotation"
                            else "not_complex_specific"
                        ),
                        "pubmed_ids": entity["pubmed_ids"],
                        "default_release_inclusion": 0,
                        "review_reason": "compound_identity_and_relation_directness_not_yet_verified",
                    }
                )
    return entities, components, binding_candidates


def assign_canonical_ids(
    entities: list[dict[str, object]],
    components: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    corum_to_cpx: dict[str, set[str]] = defaultdict(set)
    for entity in entities:
        if entity["source_database"] != "Complex Portal":
            continue
        for corum_id in split_values(str(entity["corum_identity_ids"]), r";"):
            corum_to_cpx[corum_id].add(str(entity["source_complex_id"]))

    crosswalk: list[dict[str, object]] = []
    source_to_target: dict[str, str] = {}
    for entity in entities:
        source = str(entity["source_database"])
        source_id = str(entity["source_complex_id"])
        if source == "Complex Portal":
            target_id = f"HMCX-{source_id}"
            method = "complex_portal_primary_accession"
            status = "exact_source_identity"
        elif source == "CORUM" and len(corum_to_cpx.get(source_id, set())) == 1:
            cpx = next(iter(corum_to_cpx[source_id]))
            target_id = f"HMCX-{cpx}"
            method = "complex_portal_official_corum_identity_xref"
            status = "exact_official_cross_reference"
        elif source == "CORUM" and len(corum_to_cpx.get(source_id, set())) > 1:
            target_id = f"HMCX-CORUM-{source_id}"
            method = "ambiguous_complex_portal_corum_xref_not_merged"
            status = "review"
        else:
            target_id = f"HMCX-CORUM-{source_id}"
            method = "corum_primary_accession_unmerged"
            status = "source_specific"
        source_entity_id = str(entity["source_entity_id"])
        source_to_target[source_entity_id] = target_id
        crosswalk.append(
            {
                "source_entity_id": source_entity_id,
                "complex_target_id": target_id,
                "source_database": source,
                "source_complex_id": source_id,
                "mapping_method": method,
                "mapping_status": status,
            }
        )
        entity["complex_target_id"] = target_id

    for component in components:
        component["complex_target_id"] = source_to_target[str(component["source_entity_id"])]

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for entity in entities:
        grouped[str(entity["complex_target_id"])].append(entity)
    master: list[dict[str, object]] = []
    for target_id, records in sorted(grouped.items()):
        records.sort(
            key=lambda row: (
                row["source_dataset_class"] == "predicted",
                row["source_database"] != "Complex Portal",
                row["source_complex_id"],
            )
        )
        preferred = records[0]
        protein_sets = [set(split_values(str(row["component_set_signature"]), r";")) for row in records]
        union = sorted(set().union(*protein_sets)) if protein_sets else []
        membrane_union = sorted(
            set().union(
                *[set(split_values(str(row["membrane_component_ids"]), r";")) for row in records]
            )
        )
        component_conflict = len({tuple(sorted(values)) for values in protein_sets}) > 1
        curated_sources = {
            str(row["source_database"])
            for row in records
            if row["curation_status"] == "curated"
        }
        predicted_sources = {
            str(row["source_database"])
            for row in records
            if row["curation_status"] == "predicted"
        }
        if len(union) == 1:
            complex_class = "homomer_or_single_component_complex"
        elif len(union) > 1:
            complex_class = "heteromeric_complex"
        else:
            complex_class = "component_unresolved"
        master.append(
            {
                "complex_target_id": target_id,
                "complex_name": preferred["recommended_name"],
                "aliases": preferred["aliases"],
                "organism_taxon_id": "9606",
                "complex_class": complex_class,
                "protein_component_ids": ";".join(union),
                "protein_component_count": len(union),
                "membrane_component_ids": ";".join(membrane_union),
                "membrane_component_count": len(membrane_union),
                "membrane_relevance": (
                    "all_protein_components_membrane"
                    if union and len(membrane_union) == len(union)
                    else "contains_membrane_protein"
                ),
                "curated_source_count": len(curated_sources),
                "predicted_source_count": len(predicted_sources),
                "contributing_databases": ";".join(
                    sorted({str(row["source_database"]) for row in records})
                ),
                "source_complex_ids": ";".join(
                    sorted(
                        f"{row['source_database']}:{row['source_complex_id']}"
                        for row in records
                    )
                ),
                "component_set_conflict_flag": int(component_conflict),
                "stoichiometry_status": (
                    "source_annotations_available"
                    if any(clean(row["complex_assembly"]) for row in records)
                    else "unknown"
                ),
                "required_optional_subunit_status": "not_yet_curated",
                "default_release_inclusion": int(bool(curated_sources) and not component_conflict),
                "review_status": (
                    "component_conflict_review"
                    if component_conflict
                    else "predicted_only_review"
                    if not curated_sources
                    else "curated_candidate"
                ),
                "module_version": "complex_target_v0.1_v63_candidate",
            }
        )

    by_signature: dict[str, list[dict[str, object]]] = defaultdict(list)
    for entity in entities:
        signature = str(entity["component_set_signature"])
        if signature:
            by_signature[signature].append(entity)
    similarity_review: list[dict[str, object]] = []
    for signature, records in by_signature.items():
        target_ids = sorted({str(row["complex_target_id"]) for row in records})
        if len(target_ids) < 2:
            continue
        databases = sorted({str(row["source_database"]) for row in records})
        if len(databases) < 2:
            continue
        similarity_review.append(
            {
                "review_id": stable_id("HMCX-XREF-REV", signature),
                "component_set_signature": signature,
                "candidate_complex_target_ids": ";".join(target_ids),
                "source_entities": ";".join(
                    sorted(str(row["source_entity_id"]) for row in records)
                ),
                "contributing_databases": ";".join(databases),
                "mapping_status": "same_component_set_not_auto_merged",
                "review_reason": "same_components_do_not_prove_same_functional_complex_or_stoichiometry",
            }
        )
    return master, components, crosswalk, similarity_review


def build_structure_candidates(
    entities: list[dict[str, object]],
    assembly_path: Path,
) -> list[dict[str, object]]:
    assembly_by_pdb: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_gzip_tsv(assembly_path):
        assembly_by_pdb[clean(row.get("pdb_id")).upper()].append(row)
    output: list[dict[str, object]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for entity in entities:
        if entity["source_database"] != "Complex Portal":
            continue
        for match in PDB_XREF_PATTERN.finditer(str(entity["cross_references"])):
            pdb_id = match.group(1).upper()
            qualifier = match.group(2)
            assemblies = assembly_by_pdb.get(pdb_id, [])
            if not assemblies:
                assemblies = [{"biological_assembly_id": ""}]
            assembly_ids = sorted(
                {clean(row.get("biological_assembly_id")) for row in assemblies}
            )
            for assembly_id in assembly_ids:
                key = (str(entity["complex_target_id"]), pdb_id, assembly_id, qualifier)
                if key in seen:
                    continue
                seen.add(key)
                output.append(
                    {
                        "complex_structure_evidence_id": stable_id(
                            "HMCX-STRUCT", "|".join(key)
                        ),
                        "complex_target_id": entity["complex_target_id"],
                        "source_entity_id": entity["source_entity_id"],
                        "pdb_id": pdb_id,
                        "biological_assembly_id": assembly_id,
                        "source_xref_qualifier": qualifier,
                        "candidate_assembly_count_for_pdb": len(assembly_ids),
                        "structure_mapping_status": (
                            "pdb_xref_single_candidate_assembly"
                            if len(assembly_ids) == 1
                            else "pdb_xref_multiple_candidate_assemblies"
                            if assembly_ids
                            else "pdb_xref_no_cached_assembly"
                        ),
                        "chain_to_component_mapping_status": "not_yet_resolved",
                        "default_release_inclusion": 0,
                        "review_reason": "PDB_xref_does_not_identify_biological_assembly_or_binding_chain",
                    }
                )
    return output


def fields_for(rows: list[dict[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-master", type=Path, required=True)
    parser.add_argument("--complexportal-curated", type=Path, required=True)
    parser.add_argument("--complexportal-predicted", type=Path, required=True)
    parser.add_argument("--corum-human", type=Path, required=True)
    parser.add_argument("--pdbe-assembly", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    membrane_ids, _ = load_membrane_reference(args.protein_master)
    cp_curated = parse_complex_portal(
        args.complexportal_curated, "curated", membrane_ids
    )
    cp_predicted = parse_complex_portal(
        args.complexportal_predicted, "predicted", membrane_ids
    )
    corum = parse_corum(args.corum_human, membrane_ids)
    entities = cp_curated[0] + cp_predicted[0] + corum[0]
    components = cp_curated[1] + cp_predicted[1] + corum[1]
    binding_candidates = cp_curated[2] + cp_predicted[2] + corum[2]
    master, components, crosswalk, similarity_review = assign_canonical_ids(
        entities, components
    )
    source_to_target = {
        str(row["source_entity_id"]): str(row["complex_target_id"])
        for row in crosswalk
    }
    for row in binding_candidates:
        row["complex_target_id"] = source_to_target[str(row["source_entity_id"])]
    structures = build_structure_candidates(entities, args.pdbe_assembly)

    outputs = {
        "complex_target_master_v0_1.tsv": master,
        "complex_target_components_v0_1.tsv": components,
        "complex_source_crosswalk_v0_1.tsv": crosswalk,
        "complex_source_entities_v0_1.tsv": entities,
        "complex_identity_review_v0_1.tsv": similarity_review,
        "complex_binding_evidence_candidates_v0_1.tsv": binding_candidates,
        "complex_structure_evidence_candidates_v0_1.tsv": structures,
    }
    for filename, rows in outputs.items():
        write_tsv(args.output / filename, rows, fields_for(rows))

    master_ids = {str(row["complex_target_id"]) for row in master}
    component_ids = [str(row["component_assertion_id"]) for row in components]
    crosswalk_sources = [str(row["source_entity_id"]) for row in crosswalk]
    blocking_errors = {
        "duplicate_complex_target_ids": len(master_ids) != len(master),
        "duplicate_component_assertion_ids": len(set(component_ids)) != len(component_ids),
        "duplicate_crosswalk_source_entities": len(set(crosswalk_sources))
        != len(crosswalk_sources),
        "component_foreign_key_failures": sum(
            str(row["complex_target_id"]) not in master_ids for row in components
        ),
        "binding_candidate_foreign_key_failures": sum(
            str(row["complex_target_id"]) not in master_ids
            for row in binding_candidates
        ),
        "predicted_in_default_release": sum(
            int(row["default_release_inclusion"])
            for row in master
            if row["predicted_source_count"] and not row["curated_source_count"]
        ),
    }
    status = "PASS" if not any(bool(value) for value in blocking_errors.values()) else "FAIL"
    report = {
        "status": status,
        "generated_at_utc": now_utc(),
        "module_version": "complex_target_v0.1_v63_candidate",
        "v62_release_modified": False,
        "membrane_reference_count": len(membrane_ids),
        "counts": {
            "complex_targets": len(master),
            "default_curated_complex_targets": sum(
                int(row["default_release_inclusion"]) for row in master
            ),
            "predicted_only_complex_targets": sum(
                not int(row["curated_source_count"]) for row in master
            ),
            "component_assertions": len(components),
            "membrane_component_assertions": sum(
                int(row["membrane_component_flag"]) for row in components
            ),
            "source_crosswalks": len(crosswalk),
            "same_component_set_review_groups": len(similarity_review),
            "complex_binding_candidates": len(binding_candidates),
            "complex_structure_candidates": len(structures),
        },
        "source_counts": dict(Counter(str(row["source_database"]) for row in entities)),
        "policies": {
            "complexes_are_separate_entities": True,
            "single_protein_binding_evidence_promoted_to_complex": False,
            "component_set_similarity_auto_merged": False,
            "official_identity_xref_required_for_cross_source_merge": True,
            "required_optional_subunits_inferred": False,
            "predicted_complexes_default_release": False,
            "raw_drug_names_auto_mapped_to_compounds": False,
        },
        "blocking_errors": blocking_errors,
        "inputs": {
            str(path.name): {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in [
                args.protein_master,
                args.complexportal_curated,
                args.complexportal_predicted,
                args.corum_human,
                args.pdbe_assembly,
            ]
        },
        "outputs": {
            filename: {
                "rows": len(rows),
                "bytes": (args.output / filename).stat().st_size,
                "sha256": sha256(args.output / filename),
            }
            for filename, rows in outputs.items()
        },
    }
    report_path = args.output / "COMPLEX_TARGET_MODULE_V0_1_VALIDATION.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
