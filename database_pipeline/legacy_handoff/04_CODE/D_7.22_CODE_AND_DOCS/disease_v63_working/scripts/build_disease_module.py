#!/usr/bin/env python3
"""Build the V6.3 exact-only disease identity, hierarchy and annotation module."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


VERSION = "0.1.0"
TRUSTED_EQUIVALENCE_PREFIXES = {"OMIM", "ORPHANET", "EFO", "DOID"}
EVIDENCE_RANK = {"very_high": 3, "high": 2, "medium": 1, "low": 0}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--disease-relations", type=Path, required=True)
    p.add_argument("--protein-master", type=Path, required=True)
    p.add_argument("--mondo-json", type=Path, required=True)
    p.add_argument("--mondo-release-json", type=Path, required=True)
    p.add_argument("--doid-json", type=Path, required=True)
    p.add_argument("--uberon-obo", type=Path, required=True)
    p.add_argument("--opentargets-parquet", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--qa-dir", type=Path, required=True)
    return p.parse_args()


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_id(value: object, ontology_hint: str = "") -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    text = text.replace("MONDO_", "MONDO:").replace("EFO_", "EFO:")
    text = text.replace("DOID_", "DOID:").replace("Orphanet_", "ORPHANET:")
    if text.startswith("MIM_"):
        text = "OMIM:" + text.split("_", 1)[1]
    if text.upper().startswith("MIM:"):
        text = "OMIM:" + text.split(":", 1)[1]
    if ontology_hint.upper() == "OMIM" and re.fullmatch(r"\d+", text):
        text = f"OMIM:{text}"
    if ":" in text:
        prefix, suffix = text.split(":", 1)
        prefix = "ORPHANET" if prefix.upper() == "ORPHANET" else prefix.upper()
        return f"{prefix}:{suffix}"
    return text


def iri_to_curie(value: object) -> str:
    text = str(value or "")
    patterns = [
        (r"/MONDO_(\d+)$", "MONDO"),
        (r"/DOID_(\d+)$", "DOID"),
        (r"/UBERON_(\d+)$", "UBERON"),
        (r"/EFO_(\d+)$", "EFO"),
        (r"/Orphanet_(\d+)$", "ORPHANET"),
        (r"/omim/(\d+)$", "OMIM"),
        (r"/entry/(\d+)$", "OMIM"),
    ]
    for pattern, prefix in patterns:
        match = re.search(pattern, text, flags=re.I)
        if match:
            return f"{prefix}:{match.group(1)}"
    if re.match(r"^[A-Za-z][A-Za-z0-9_.-]*[:_]", text):
        return normalize_id(text)
    return ""


def write_tsv(path: Path, fields: list[str], rows) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        count = 0
        for row in rows:
            writer.writerow({key: "" if value is None else value for key, value in row.items()})
            count += 1
    return count


def truth(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def joined(values) -> str:
    return ";".join(sorted({str(value).strip() for value in values if str(value).strip() and str(value).strip().lower() != "nan"}))


def load_ontology_json(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    graph = payload["graphs"][0]
    nodes = {}
    for node in graph.get("nodes", []):
        curie = iri_to_curie(node.get("id"))
        if curie:
            nodes[curie] = node
    return graph, nodes


def node_obsolete(node: dict) -> bool:
    meta = node.get("meta") or {}
    return bool(meta.get("deprecated") or meta.get("obsolete"))


def node_definition(node: dict) -> str:
    definition = (node.get("meta") or {}).get("definition") or {}
    return str(definition.get("val") or "")


def node_synonyms(node: dict, predicate: str | None = None) -> list[str]:
    values = []
    for synonym in (node.get("meta") or {}).get("synonyms") or []:
        if predicate is None or synonym.get("pred") == predicate:
            if synonym.get("val"):
                values.append(str(synonym["val"]))
    return sorted(set(values))


def ontology_version(graph: dict) -> str:
    meta = graph.get("meta") or {}
    return str(meta.get("version") or "")


def build_parent_map(graph: dict, prefix: str) -> tuple[dict[str, set[str]], list[dict[str, str]]]:
    parents: dict[str, set[str]] = defaultdict(set)
    rows = []
    for edge in graph.get("edges", []):
        if edge.get("pred") != "is_a":
            continue
        child = iri_to_curie(edge.get("sub"))
        parent = iri_to_curie(edge.get("obj"))
        if child.startswith(prefix + ":") and parent.startswith(prefix + ":"):
            parents[child].add(parent)
            rows.append({"child_disease_id": child, "parent_disease_id": parent, "hierarchy_predicate": "is_a"})
    return parents, rows


def ancestor_distances(start: str, parents: dict[str, set[str]], max_depth: int = 40) -> dict[str, int]:
    result: dict[str, int] = {}
    queue = deque([(start, 0)])
    while queue:
        node, distance = queue.popleft()
        if distance >= max_depth:
            continue
        for parent in parents.get(node, set()):
            next_distance = distance + 1
            if parent not in result or next_distance < result[parent]:
                result[parent] = next_distance
                queue.append((parent, next_distance))
    return result


def is_descendant(node: str, root: str, parents: dict[str, set[str]]) -> bool:
    return node == root or root in ancestor_distances(node, parents)


def extract_mondo_mappings(mondo_nodes: dict[str, dict]):
    records: dict[tuple[str, str, str], dict[str, object]] = {}
    for mondo_id, node in mondo_nodes.items():
        if not mondo_id.startswith("MONDO:"):
            continue
        meta = node.get("meta") or {}
        for xref in meta.get("xrefs") or []:
            source_id = normalize_id(xref.get("val"))
            prefix = source_id.split(":", 1)[0] if ":" in source_id else ""
            if prefix not in TRUSTED_EQUIVALENCE_PREFIXES:
                continue
            key = (source_id, mondo_id, "mondo:xref_equivalence_proxy")
            records[key] = {
                "source_disease_id": source_id,
                "canonical_mondo_id": mondo_id,
                "mapping_predicate": "mondo:xref_equivalence_proxy",
                "mapping_relation": "exact_or_equivalent",
                "official_mapping_source": "Mondo meta.xrefs trusted equivalence prefixes",
                "identity_merge_candidate": 1,
            }
        for prop in meta.get("basicPropertyValues") or []:
            predicate = str(prop.get("pred") or "")
            if "skos/core#" not in predicate or not predicate.endswith("Match"):
                continue
            source_id = iri_to_curie(prop.get("val"))
            prefix = source_id.split(":", 1)[0] if ":" in source_id else ""
            if prefix not in TRUSTED_EQUIVALENCE_PREFIXES:
                continue
            short = predicate.rsplit("#", 1)[-1]
            exact = short == "exactMatch"
            key = (source_id, mondo_id, f"skos:{short}")
            records[key] = {
                "source_disease_id": source_id,
                "canonical_mondo_id": mondo_id,
                "mapping_predicate": f"skos:{short}",
                "mapping_relation": "exact_or_equivalent" if exact else short,
                "official_mapping_source": "Mondo basicPropertyValues",
                "identity_merge_candidate": int(exact),
            }
    exact_targets: dict[str, set[str]] = defaultdict(set)
    for row in records.values():
        if row["identity_merge_candidate"]:
            exact_targets[str(row["source_disease_id"])].add(str(row["canonical_mondo_id"]))
    rows = []
    for row in records.values():
        count = len(exact_targets.get(str(row["source_disease_id"]), set()))
        row = dict(row)
        row["exact_mapping_target_count"] = count
        row["automatic_identity_merge_allowed"] = int(bool(row["identity_merge_candidate"]) and count == 1)
        rows.append(row)
    return sorted(rows, key=lambda row: (str(row["source_disease_id"]), str(row["canonical_mondo_id"]), str(row["mapping_predicate"]))), exact_targets


def parse_uberon(path: Path):
    terms: dict[str, dict[str, object]] = {}
    current: dict[str, object] | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                if current and str(current.get("id", "")).startswith("UBERON:"):
                    terms[str(current["id"])] = current
                current = {"parents": set(), "part_of": set()}
            elif line.startswith("["):
                if current and str(current.get("id", "")).startswith("UBERON:"):
                    terms[str(current["id"])] = current
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:].strip()
                elif line.startswith("name: "):
                    current["name"] = line[6:].strip()
                elif line.startswith("is_a: "):
                    current["parents"].add(line[6:].split(" ! ", 1)[0].strip())
                elif line.startswith("relationship: part_of "):
                    current["part_of"].add(line.split()[2])
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
        if current and str(current.get("id", "")).startswith("UBERON:"):
            terms[str(current["id"])] = current
    return terms


def extract_anatomy_assertions(graph: dict, disease_prefix: str):
    assertions: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen = set()
    for axiom in graph.get("logicalDefinitionAxioms", []):
        disease_id = iri_to_curie(axiom.get("definedClassId"))
        if not disease_id.startswith(disease_prefix + ":"):
            continue
        for restriction in axiom.get("restrictions", []):
            anatomy_id = iri_to_curie(restriction.get("fillerId"))
            if not anatomy_id.startswith("UBERON:"):
                continue
            predicate = iri_to_curie(restriction.get("propertyId")) or str(restriction.get("propertyId") or "")
            key = (disease_id, anatomy_id, predicate, "logical_definition")
            if key not in seen:
                assertions[disease_id].append({"anatomy_id": anatomy_id, "predicate": predicate, "axiom_type": "logical_definition"})
                seen.add(key)
    for edge in graph.get("edges", []):
        disease_id = iri_to_curie(edge.get("sub"))
        anatomy_id = iri_to_curie(edge.get("obj"))
        if not disease_id.startswith(disease_prefix + ":") or not anatomy_id.startswith("UBERON:"):
            continue
        predicate = iri_to_curie(edge.get("pred")) or str(edge.get("pred") or "")
        key = (disease_id, anatomy_id, predicate, "edge")
        if key not in seen:
            assertions[disease_id].append({"anatomy_id": anatomy_id, "predicate": predicate, "axiom_type": "edge"})
            seen.add(key)
    return assertions


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.qa_dir.mkdir(parents=True, exist_ok=True)

    mondo_graph, all_mondo_nodes = load_ontology_json(args.mondo_json)
    doid_graph, all_doid_nodes = load_ontology_json(args.doid_json)
    mondo_nodes = {key: value for key, value in all_mondo_nodes.items() if key.startswith("MONDO:")}
    doid_nodes = {key: value for key, value in all_doid_nodes.items() if key.startswith("DOID:")}
    mondo_parents, hierarchy_rows = build_parent_map(mondo_graph, "MONDO")
    doid_parents, _ = build_parent_map(doid_graph, "DOID")
    mapping_rows, exact_targets = extract_mondo_mappings(mondo_nodes)

    relations = pd.read_csv(args.disease_relations, sep="\t", dtype=str, keep_default_na=False)
    proteins = pd.read_csv(args.protein_master, sep="\t", dtype=str, usecols=["target_uniprot_id"], keep_default_na=False)
    protein_ids = set(proteins["target_uniprot_id"])
    relations["normalized_source_disease_id"] = [normalize_id(value, hint) for value, hint in zip(relations["disease_id"], relations["disease_ontology"])]

    source_rows = []
    source_resolution: dict[str, dict[str, object]] = {}
    review_rows = []
    for normalized_id, frame in relations.groupby("normalized_source_disease_id", sort=True):
        raw_ids = sorted(set(frame["disease_id"]))
        names = sorted(set(name for name in frame["disease_name"] if name))
        ontologies = sorted(set(frame["disease_ontology"]))
        direct = normalized_id.startswith("MONDO:") and normalized_id in mondo_nodes
        candidates = {normalized_id} if direct else set(exact_targets.get(normalized_id, set()))
        if direct:
            canonical_id = normalized_id
            status = "direct_mondo_id"
            predicate = "identity"
        elif len(candidates) == 1:
            candidate = next(iter(candidates))
            if node_obsolete(mondo_nodes[candidate]):
                canonical_id = "SOURCE_ONLY:" + normalized_id.replace(":", "_")
                status = "unique_mapping_to_obsolete_mondo_review"
                predicate = "exact_but_obsolete"
            else:
                canonical_id = candidate
                status = "unique_official_exact_mapping"
                predicate = "official_exact_or_equivalent"
        elif len(candidates) > 1:
            canonical_id = "SOURCE_ONLY:" + normalized_id.replace(":", "_")
            status = "ambiguous_1_to_many_review"
            predicate = "exact_1_to_many_not_merged"
        else:
            canonical_id = "SOURCE_ONLY:" + normalized_id.replace(":", "_")
            status = "unmapped_source_only"
            predicate = "no_official_exact_mapping"
        non_disease = False
        if canonical_id.startswith("MONDO:") and not is_descendant(canonical_id, "MONDO:0000001", mondo_parents):
            non_disease = True
            status = "non_disease_mondo_entity_review"
        source_row = {
            "disease_source_entity_id": "SOURCE_DISEASE:" + normalized_id.replace(":", "_"),
            "normalized_source_disease_id": normalized_id,
            "raw_source_disease_ids": ";".join(raw_ids),
            "source_disease_names": ";".join(names),
            "source_ontologies": ";".join(ontologies),
            "source_databases": joined(frame["source_database"]),
            "source_relation_count": len(frame),
            "canonical_disease_id": canonical_id,
            "canonicalization_status": status,
            "mapping_predicate": predicate,
            "candidate_mondo_ids": ";".join(sorted(candidates)),
            "candidate_mondo_count": len(candidates),
            "name_based_mapping_used": 0,
        }
        source_rows.append(source_row)
        source_resolution[normalized_id] = source_row
        if "review" in status or status == "unmapped_source_only" or non_disease:
            review_rows.append({
                "review_id": "DMR-" + hashlib.sha1(normalized_id.encode()).hexdigest()[:14].upper(),
                "source_disease_id": normalized_id,
                "source_disease_name": names[0] if names else "",
                "review_type": status,
                "candidate_mondo_ids": ";".join(sorted(candidates)),
                "automatic_action": "retain_source_only_entity",
                "review_reason": "Name equality is never used; 1:N, obsolete, non-disease and unmapped records require authoritative review.",
            })

    active_canonical_ids = {str(row["canonical_disease_id"]) for row in source_rows}
    exact_by_mondo: dict[str, set[str]] = defaultdict(set)
    for row in mapping_rows:
        if row["automatic_identity_merge_allowed"]:
            exact_by_mondo[str(row["canonical_mondo_id"])].add(str(row["source_disease_id"]))

    entity_rows = []
    for mondo_id, node in sorted(mondo_nodes.items()):
        exact_ids = sorted(exact_by_mondo.get(mondo_id, set()))
        entity_rows.append({
            "canonical_disease_id": mondo_id,
            "preferred_name": node.get("lbl") or "",
            "entity_origin": "MONDO",
            "mondo_id": mondo_id,
            "is_obsolete": int(node_obsolete(node)),
            "definition": node_definition(node),
            "exact_synonyms": ";".join(node_synonyms(node, "hasExactSynonym")),
            "related_synonyms": ";".join(node_synonyms(node, "hasRelatedSynonym")),
            "official_exact_source_ids": ";".join(exact_ids),
            "official_exact_source_id_count": len(exact_ids),
            "in_protein_disease_relations": int(mondo_id in active_canonical_ids),
            "canonicalization_policy": "MONDO native entity",
        })
    for source in source_rows:
        canonical_id = str(source["canonical_disease_id"])
        if not canonical_id.startswith("SOURCE_ONLY:"):
            continue
        entity_rows.append({
            "canonical_disease_id": canonical_id,
            "preferred_name": str(source["source_disease_names"]).split(";", 1)[0],
            "entity_origin": "source_only",
            "mondo_id": "",
            "is_obsolete": 0,
            "definition": "",
            "exact_synonyms": "",
            "related_synonyms": "",
            "official_exact_source_ids": source["normalized_source_disease_id"],
            "official_exact_source_id_count": 0,
            "in_protein_disease_relations": 1,
            "canonicalization_policy": source["canonicalization_status"],
        })
    entity_ids = {str(row["canonical_disease_id"]) for row in entity_rows}

    source_fields = ["disease_source_entity_id", "normalized_source_disease_id", "raw_source_disease_ids", "source_disease_names", "source_ontologies", "source_databases", "source_relation_count", "canonical_disease_id", "canonicalization_status", "mapping_predicate", "candidate_mondo_ids", "candidate_mondo_count", "name_based_mapping_used"]
    mapping_fields = ["source_disease_id", "canonical_mondo_id", "mapping_predicate", "mapping_relation", "official_mapping_source", "identity_merge_candidate", "exact_mapping_target_count", "automatic_identity_merge_allowed"]
    entity_fields = ["canonical_disease_id", "preferred_name", "entity_origin", "mondo_id", "is_obsolete", "definition", "exact_synonyms", "related_synonyms", "official_exact_source_ids", "official_exact_source_id_count", "in_protein_disease_relations", "canonicalization_policy"]
    hierarchy_fields = ["child_disease_id", "parent_disease_id", "hierarchy_predicate", "source_ontology", "source_version"]
    for row in hierarchy_rows:
        row["source_ontology"] = "MONDO"
        row["source_version"] = ontology_version(mondo_graph)
    write_tsv(args.output_dir / "disease_source_entity_v0_1.tsv", source_fields, source_rows)
    write_tsv(args.output_dir / "disease_xref_mapping_v0_1.tsv", mapping_fields, mapping_rows)
    write_tsv(args.output_dir / "disease_entity_master_v0_1.tsv", entity_fields, entity_rows)
    write_tsv(args.output_dir / "disease_hierarchy_v0_1.tsv", hierarchy_fields, hierarchy_rows)

    source_evidence_rows = []
    aggregate: dict[tuple[str, str], dict[str, object]] = {}
    flag_fields = ["human_genetic_flag", "clinical_genetic_flag", "somatic_mutation_flag", "functional_flag", "animal_model_flag", "expression_flag", "literature_flag", "uniprot_disease_flag"]
    for row in relations.to_dict(orient="records"):
        source = source_resolution[row["normalized_source_disease_id"]]
        canonical_id = str(source["canonical_disease_id"])
        key = (row["target_uniprot_id"], canonical_id)
        relation_id = "PDR2-" + hashlib.sha1((key[0] + "|" + key[1]).encode()).hexdigest()[:18].upper()
        source_evidence_rows.append({
            "disease_relation_id_v2": relation_id,
            "source_disease_relation_id": row["disease_relation_id"],
            "target_uniprot_id": row["target_uniprot_id"],
            "canonical_disease_id": canonical_id,
            "source_disease_id": row["disease_id"],
            "normalized_source_disease_id": row["normalized_source_disease_id"],
            "mapping_predicate": source["mapping_predicate"],
            "canonicalization_status": source["canonicalization_status"],
            "source_database": row["source_database"],
            "source_record_id": row["source_record_id"],
            "source_url": row["source_url"],
            "evidence_level": row["disease_evidence_level"],
            "release_source": "V6.2 frozen",
        })
        if key not in aggregate:
            aggregate[key] = {
                "disease_relation_id_v2": relation_id,
                "target_uniprot_id": row["target_uniprot_id"],
                "approved_symbol": row["approved_symbol"],
                "canonical_disease_id": canonical_id,
                "canonical_disease_name": next((str(entity["preferred_name"]) for entity in entity_rows if entity["canonical_disease_id"] == canonical_id), row["disease_name"]),
                "source_disease_ids": set(), "source_disease_names": set(), "mapping_predicates": set(),
                "source_databases": set(), "source_record_ids": set(), "source_relation_ids": set(), "pubmed_ids": set(),
                "evidence_levels": set(), "ot_scores": [], "default_relation_inclusion": False,
                **{field: False for field in flag_fields},
            }
        item = aggregate[key]
        item["source_disease_ids"].add(row["disease_id"])
        item["source_disease_names"].add(row["disease_name"])
        item["mapping_predicates"].add(str(source["mapping_predicate"]))
        item["source_databases"].add(row["source_database"])
        item["source_record_ids"].add(row["source_record_id"])
        item["source_relation_ids"].add(row["disease_relation_id"])
        item["evidence_levels"].add(row["disease_evidence_level"])
        item["pubmed_ids"].update(value for value in row["pubmed_ids"].split(";") if value)
        if row["ot_overall_score"]:
            try:
                item["ot_scores"].append(float(row["ot_overall_score"]))
            except ValueError:
                pass
        item["default_relation_inclusion"] = item["default_relation_inclusion"] or truth(row["default_relation_inclusion"])
        for field in flag_fields:
            item[field] = item[field] or truth(row[field])

    relation_rows = []
    for item in aggregate.values():
        levels = [level for level in item["evidence_levels"] if level]
        best_level = max(levels, key=lambda level: EVIDENCE_RANK.get(level, -1)) if levels else ""
        relation_rows.append({
            "disease_relation_id_v2": item["disease_relation_id_v2"],
            "target_uniprot_id": item["target_uniprot_id"],
            "approved_symbol": item["approved_symbol"],
            "canonical_disease_id": item["canonical_disease_id"],
            "canonical_disease_name": item["canonical_disease_name"],
            "source_disease_ids": joined(item["source_disease_ids"]),
            "source_disease_names": joined(item["source_disease_names"]),
            "mapping_predicates": joined(item["mapping_predicates"]),
            "source_databases": joined(item["source_databases"]),
            "source_record_ids": joined(item["source_record_ids"]),
            "source_relation_ids": joined(item["source_relation_ids"]),
            "source_relation_count": len(item["source_relation_ids"]),
            "evidence_channels": joined(field.replace("_flag", "") for field in flag_fields if item[field]),
            "best_evidence_level": best_level,
            "max_ot_overall_score": max(item["ot_scores"]) if item["ot_scores"] else "",
            **{field: int(item[field]) for field in flag_fields},
            "pubmed_ids": joined(item["pubmed_ids"]),
            "default_relation_inclusion": int(item["default_relation_inclusion"]),
            "release_version": "V6.3-candidate",
        })
    relation_fields = ["disease_relation_id_v2", "target_uniprot_id", "approved_symbol", "canonical_disease_id", "canonical_disease_name", "source_disease_ids", "source_disease_names", "mapping_predicates", "source_databases", "source_record_ids", "source_relation_ids", "source_relation_count", "evidence_channels", "best_evidence_level", "max_ot_overall_score"] + flag_fields + ["pubmed_ids", "default_relation_inclusion", "release_version"]
    source_evidence_fields = ["disease_relation_id_v2", "source_disease_relation_id", "target_uniprot_id", "canonical_disease_id", "source_disease_id", "normalized_source_disease_id", "mapping_predicate", "canonicalization_status", "source_database", "source_record_id", "source_url", "evidence_level", "release_source"]
    write_tsv(args.output_dir / "protein_disease_relation_v2.tsv", relation_fields, sorted(relation_rows, key=lambda row: (str(row["target_uniprot_id"]), str(row["canonical_disease_id"]))))
    write_tsv(args.output_dir / "protein_disease_source_evidence_v2.tsv", source_evidence_fields, source_evidence_rows)

    ot = pd.read_parquet(args.opentargets_parquet)
    ot_by_id = {str(row.id): row for row in ot.itertuples(index=False)}
    therapeutic_rows = []
    therapy_coverage = []
    active_mondo_ids = sorted({str(row["canonical_disease_id"]) for row in relation_rows if str(row["canonical_disease_id"]).startswith("MONDO:")})
    for canonical_id in active_mondo_ids:
        ot_id = canonical_id.replace(":", "_")
        record = ot_by_id.get(ot_id)
        area_ids = [] if record is None else [str(value) for value in list(record.therapeuticAreas)]
        therapy_coverage.append({"canonical_disease_id": canonical_id, "opentargets_disease_id": ot_id if record is not None else "", "therapeutic_area_count": len(area_ids), "classification_status": "mapped_direct_ot_id" if record is not None else "no_direct_ot_26_06_entity"})
        for area_id in sorted(set(area_ids)):
            area = ot_by_id.get(area_id)
            therapeutic_rows.append({
                "canonical_disease_id": canonical_id,
                "opentargets_disease_id": ot_id,
                "therapeutic_area_id": area_id,
                "therapeutic_area_name": str(area.name) if area is not None else area_id,
                "membership_type": "multi_label_direct",
                "source_database": "Open Targets Platform",
                "source_version": "26.06",
                "keyword_classification_used": 0,
            })
    therapeutic_fields = ["canonical_disease_id", "opentargets_disease_id", "therapeutic_area_id", "therapeutic_area_name", "membership_type", "source_database", "source_version", "keyword_classification_used"]
    write_tsv(args.output_dir / "disease_therapeutic_area_v0_1.tsv", therapeutic_fields, therapeutic_rows)
    write_tsv(args.output_dir / "disease_therapeutic_area_coverage_v0_1.tsv", ["canonical_disease_id", "opentargets_disease_id", "therapeutic_area_count", "classification_status"], therapy_coverage)

    uberon = parse_uberon(args.uberon_obo)
    uberon_parents = {term_id: set(term.get("parents", set())) | set(term.get("part_of", set())) for term_id, term in uberon.items()}
    mondo_anatomy = extract_anatomy_assertions(mondo_graph, "MONDO")
    doid_anatomy = extract_anatomy_assertions(doid_graph, "DOID")
    doids_by_mondo = defaultdict(set)
    for row in mapping_rows:
        if row["automatic_identity_merge_allowed"] and str(row["source_disease_id"]).startswith("DOID:"):
            doids_by_mondo[str(row["canonical_mondo_id"])].add(str(row["source_disease_id"]))
    anatomy_rows_dict = {}
    for canonical_id in active_mondo_ids:
        direct_sources = [(canonical_id, 0, "MONDO")]
        direct_sources.extend((doid, 0, "DOID") for doid in sorted(doids_by_mondo.get(canonical_id, set())))
        for source_id, distance, ontology in direct_sources:
            assertion_map = mondo_anatomy if ontology == "MONDO" else doid_anatomy
            for assertion in assertion_map.get(source_id, []):
                anatomy_id = assertion["anatomy_id"]
                key = (canonical_id, anatomy_id, ontology, source_id, 0, assertion["predicate"])
                anatomy_rows_dict[key] = {
                    "canonical_disease_id": canonical_id, "anatomy_id": anatomy_id,
                    "anatomy_name": str(uberon.get(anatomy_id, {}).get("name", "")),
                    "assertion_status": "asserted", "source_disease_id": source_id,
                    "hierarchy_distance": 0, "mapping_predicate": assertion["predicate"],
                    "axiom_type": assertion["axiom_type"], "source_ontology": ontology,
                    "uberon_source": "Uberon basic", "review_flag": int(anatomy_id not in uberon),
                }
        for ancestor, distance in ancestor_distances(canonical_id, mondo_parents).items():
            for assertion in mondo_anatomy.get(ancestor, []):
                anatomy_id = assertion["anatomy_id"]
                key = (canonical_id, anatomy_id, "MONDO", ancestor, distance, assertion["predicate"])
                anatomy_rows_dict[key] = {
                    "canonical_disease_id": canonical_id, "anatomy_id": anatomy_id,
                    "anatomy_name": str(uberon.get(anatomy_id, {}).get("name", "")),
                    "assertion_status": "inferred_from_mondo_ancestor", "source_disease_id": ancestor,
                    "hierarchy_distance": distance, "mapping_predicate": assertion["predicate"],
                    "axiom_type": assertion["axiom_type"], "source_ontology": "MONDO",
                    "uberon_source": "Uberon basic", "review_flag": int(anatomy_id not in uberon),
                }
        for doid in sorted(doids_by_mondo.get(canonical_id, set())):
            for ancestor, distance in ancestor_distances(doid, doid_parents).items():
                for assertion in doid_anatomy.get(ancestor, []):
                    anatomy_id = assertion["anatomy_id"]
                    key = (canonical_id, anatomy_id, "DOID", ancestor, distance, assertion["predicate"])
                    anatomy_rows_dict[key] = {
                        "canonical_disease_id": canonical_id, "anatomy_id": anatomy_id,
                        "anatomy_name": str(uberon.get(anatomy_id, {}).get("name", "")),
                        "assertion_status": "inferred_from_doid_ancestor", "source_disease_id": ancestor,
                        "hierarchy_distance": distance, "mapping_predicate": assertion["predicate"],
                        "axiom_type": assertion["axiom_type"], "source_ontology": "DOID",
                        "uberon_source": "Uberon basic", "review_flag": int(anatomy_id not in uberon),
                    }
    anatomy_rows = list(anatomy_rows_dict.values())
    anatomy_fields = ["canonical_disease_id", "anatomy_id", "anatomy_name", "assertion_status", "source_disease_id", "hierarchy_distance", "mapping_predicate", "axiom_type", "source_ontology", "uberon_source", "review_flag"]
    write_tsv(args.output_dir / "disease_anatomy_v0_1.tsv", anatomy_fields, sorted(anatomy_rows, key=lambda row: (str(row["canonical_disease_id"]), int(row["hierarchy_distance"]), str(row["anatomy_id"]))))

    system_root = "UBERON:0000467"
    top_systems = {term_id for term_id, term in uberon.items() if system_root in set(term.get("parents", set()))}
    system_rows = []
    system_seen = set()
    for row in anatomy_rows:
        anatomy_id = str(row["anatomy_id"])
        anatomy_ancestors = ancestor_distances(anatomy_id, uberon_parents, max_depth=30)
        systems = ({anatomy_id} if anatomy_id in top_systems else set()) | (set(anatomy_ancestors) & top_systems)
        for system_id in systems:
            key = (row["canonical_disease_id"], system_id, row["assertion_status"], row["source_disease_id"])
            if key in system_seen:
                continue
            system_seen.add(key)
            system_rows.append({
                "canonical_disease_id": row["canonical_disease_id"],
                "anatomical_system_id": system_id,
                "anatomical_system_name": str(uberon.get(system_id, {}).get("name", "")),
                "source_anatomy_id": anatomy_id,
                "source_anatomy_name": row["anatomy_name"],
                "assertion_status": row["assertion_status"],
                "mapping_method": "Uberon is_a/part_of ancestry to child of anatomical system",
                "keyword_classification_used": 0,
            })
    system_fields = ["canonical_disease_id", "anatomical_system_id", "anatomical_system_name", "source_anatomy_id", "source_anatomy_name", "assertion_status", "mapping_method", "keyword_classification_used"]
    write_tsv(args.output_dir / "disease_anatomical_system_v0_1.tsv", system_fields, system_rows)

    review_fields = ["review_id", "source_disease_id", "source_disease_name", "review_type", "candidate_mondo_ids", "automatic_action", "review_reason"]
    write_tsv(args.output_dir / "disease_mapping_review_v0_1.tsv", review_fields, review_rows)

    release_meta = json.loads(args.mondo_release_json.read_text(encoding="utf-8"))
    source_registry = [
        {"source_name": "Mondo Disease Ontology", "version": ontology_version(mondo_graph), "release_label": release_meta.get("tag_name", ""), "license": "CC BY 4.0", "source_url": "https://purl.obolibrary.org/obo/mondo.json", "local_file": args.mondo_json.name, "bytes": args.mondo_json.stat().st_size, "sha256": hash_file(args.mondo_json)},
        {"source_name": "Human Disease Ontology", "version": ontology_version(doid_graph), "release_label": "", "license": "CC0 1.0", "source_url": "https://purl.obolibrary.org/obo/doid.json", "local_file": args.doid_json.name, "bytes": args.doid_json.stat().st_size, "sha256": hash_file(args.doid_json)},
        {"source_name": "Uberon basic", "version": "from frozen OBO header", "release_label": "", "license": "CC BY 3.0", "source_url": "https://purl.obolibrary.org/obo/uberon/uberon-basic.obo", "local_file": args.uberon_obo.name, "bytes": args.uberon_obo.stat().st_size, "sha256": hash_file(args.uberon_obo)},
        {"source_name": "Open Targets Platform disease/phenotype", "version": "26.06", "release_label": "26.06", "license": "Open Targets terms/data policy; retain attribution", "source_url": "https://platform.opentargets.org/downloads", "local_file": args.opentargets_parquet.name, "bytes": args.opentargets_parquet.stat().st_size, "sha256": hash_file(args.opentargets_parquet)},
        {"source_name": "MemPro protein-disease relations", "version": "V6.2 frozen", "release_label": "V6.2", "license": "project internal integration", "source_url": "local frozen release", "local_file": args.disease_relations.name, "bytes": args.disease_relations.stat().st_size, "sha256": hash_file(args.disease_relations)},
    ]
    source_registry_fields = ["source_name", "version", "release_label", "license", "source_url", "local_file", "bytes", "sha256"]
    write_tsv(args.output_dir / "disease_source_registry_v0_1.tsv", source_registry_fields, source_registry)

    relation_keys = [(row["target_uniprot_id"], row["canonical_disease_id"]) for row in relation_rows]
    source_ids_preserved = set(relations["disease_id"]) == {value for row in source_evidence_rows for value in [row["source_disease_id"]]}
    canonical_fk_ok = all(str(row["canonical_disease_id"]) in entity_ids for row in relation_rows)
    target_fk_ok = all(str(row["target_uniprot_id"]) in protein_ids for row in relation_rows)
    checks = {
        "v62_source_relation_rows_preserved": len(source_evidence_rows) == len(relations) == 10264,
        "source_disease_ids_preserved": source_ids_preserved,
        "protein_disease_key_unique": len(relation_keys) == len(set(relation_keys)),
        "protein_foreign_keys_complete": target_fk_ok,
        "canonical_disease_foreign_keys_complete": canonical_fk_ok,
        "no_name_based_automatic_mapping": all(int(row["name_based_mapping_used"]) == 0 for row in source_rows),
        "one_to_many_not_auto_merged": all(not (row["candidate_mondo_count"] > 1 and str(row["canonical_disease_id"]).startswith("MONDO:")) for row in source_rows),
        "non_exact_mappings_not_merged": all(not row["automatic_identity_merge_allowed"] for row in mapping_rows if row["mapping_relation"] != "exact_or_equivalent"),
        "source_only_records_retained": all(str(row["canonical_disease_id"]) in entity_ids for row in source_rows),
        "therapeutic_area_keyword_rules_absent": all(int(row["keyword_classification_used"]) == 0 for row in therapeutic_rows),
        "anatomy_keyword_rules_absent": all(int(row["keyword_classification_used"]) == 0 for row in system_rows),
        "review_queue_covers_ambiguous_and_unmapped": len(review_rows) == sum(1 for row in source_rows if "review" in str(row["canonicalization_status"]) or row["canonicalization_status"] == "unmapped_source_only"),
    }
    counts = {
        "v62_source_relation_rows": len(relations),
        "unique_source_diseases": len(source_rows),
        "mondo_entities": len(mondo_nodes),
        "source_only_entities": sum(1 for row in entity_rows if row["entity_origin"] == "source_only"),
        "disease_entities_total": len(entity_rows),
        "xref_mapping_rows": len(mapping_rows),
        "canonical_protein_disease_relations": len(relation_rows),
        "canonicalization_status": dict(Counter(str(row["canonicalization_status"]) for row in source_rows)),
        "review_rows": len(review_rows),
        "therapeutic_area_memberships": len(therapeutic_rows),
        "anatomy_mapping_rows": len(anatomy_rows),
        "anatomical_system_rows": len(system_rows),
        "anatomy_assertion_status": dict(Counter(str(row["assertion_status"]) for row in anatomy_rows)),
    }
    validation = {
        "module": "disease_identity_hierarchy_classification_anatomy",
        "module_version": VERSION,
        "candidate_release": "V6.3",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": counts,
        "policy": {
            "automatic_merge": "MONDO direct ID or official exact/equivalent mapping with exactly one MONDO target",
            "never_automatic": ["name equality", "broad match", "narrow match", "related/close match", "1:N mapping"],
            "source_preservation": "source disease ID, name, database, record ID and URL retained in source-evidence bridge",
            "classification": "Open Targets 26.06 therapeuticAreas; no disease-name keywords",
            "anatomy": "Mondo/DO logical axioms and Uberon ontology ancestry with asserted/inferred status",
        },
    }
    validation_path = args.qa_dir / "DISEASE_MODULE_V0_1_VALIDATION.json"
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")

    outputs = sorted(path for path in args.output_dir.iterdir() if path.is_file())
    manifest = {
        "module_version": VERSION,
        "candidate_release": "V6.3",
        "generated_at_utc": validation["generated_at_utc"],
        "outputs": [{"relative_path": path.name, "bytes": path.stat().st_size, "sha256": hash_file(path)} for path in outputs],
    }
    (args.qa_dir / "DISEASE_MODULE_V0_1_MANIFEST.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
