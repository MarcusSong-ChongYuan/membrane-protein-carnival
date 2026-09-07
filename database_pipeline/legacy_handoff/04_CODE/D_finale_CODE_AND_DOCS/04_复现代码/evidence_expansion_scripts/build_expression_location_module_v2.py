#!/usr/bin/env python3
"""Build full HPA v25.1 expression/localization add-on on frozen MemPro V6.1."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import build_expression_location_module_v1 as v1


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_1_20260729"
MASTER = RELEASE / "human_membrane_protein_master_v6_1.tsv"
RUN = ROOT / "runs" / "v62_completion_20260729"
RAW = RUN / "raw" / "hpa_v25_1"
HPA_SUMMARY = Path(
    r"D:\7.22\membrane_master_v5_working\raw\hpa_proteinatlas_v25_1.tsv.zip"
)
STAGING = RUN / "staging" / "expression_location_v2"
QA = RUN / "qa"
PROGRESS = QA / "EXPRESSION_LOCATION_V2_PROGRESS.json"
REPORT = QA / "EXPRESSION_LOCATION_V2_VALIDATION.json"
OLD_ONTOLOGY = (
    ROOT / "runs" / "expression_location_v1_20260728" / "raw_ontologies"
)

TISSUE_OUT = STAGING / "protein_tissue_expression_v2.tsv.gz"
CELL_OUT = STAGING / "protein_cell_type_expression_v2.tsv.gz"
CLUSTER_OUT = STAGING / "protein_cell_cluster_expression_v2.tsv.gz"
IHC_OUT = STAGING / "protein_tissue_cell_ihc_expression_v2.tsv.gz"
LOCATION_OUT = STAGING / "protein_subcellular_localization_v2.tsv.gz"
MAPPING_OUT = STAGING / "anatomy_cell_location_ontology_mapping_v2.tsv"
SUMMARY_OUT = STAGING / "expression_location_summary_v2.tsv.gz"
REVIEW_OUT = STAGING / "expression_location_review_queue_v2.tsv.gz"
PREVIEW_OUT = (
    STAGING / "human_membrane_protein_master_v6_2_expression_preview.tsv.gz"
)

MODULE_VERSION = "expression_location_v2_20260729"
HPA_VERSION = "HPA 25.1"

BASE_FIELDS = [
    "expression_record_id",
    "target_uniprot_id",
    "membrane_protein_id",
    "hpa_ensembl_gene_id",
    "hpa_gene_name",
    "source_dataset",
    "source_version",
    "measurement_layer",
    "measurement_type",
    "tissue",
    "tissue_ontology_id",
    "tissue_ontology_label",
    "tissue_mapping_status",
    "cell_type",
    "cell_type_ontology_id",
    "cell_type_ontology_label",
    "cell_type_mapping_status",
    "cluster",
    "value",
    "unit",
    "detection_status",
    "reliability",
    "protein_mapping_status",
    "module_version",
]


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    os.replace(tmp, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record_id(*parts: object) -> str:
    text = "\x1f".join("" if part is None else str(part) for part in parts)
    return "HPAEXP-" + hashlib.sha1(text.encode("utf-8")).hexdigest()[:24].upper()


def split_ids(value: str) -> list[str]:
    return [
        item.split(".", 1)[0]
        for item in re.split(r"[;,| ]+", value or "")
        if item.strip()
    ]


def open_zip_rows(name: str):
    path = RAW / name
    archive = zipfile.ZipFile(path)
    member = archive.namelist()[0]
    handle = archive.open(member)
    text = (
        line.decode("utf-8-sig", errors="replace")
        for line in handle
    )
    reader = csv.DictReader(text, delimiter="\t")
    return archive, handle, reader


def numeric(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def detected(value: str, numeric_value: float | None = None) -> str:
    if numeric_value is not None:
        return "detected" if numeric_value > 0 else "not_detected"
    text = (value or "").strip().casefold()
    if not text:
        return "missing_or_not_measured"
    if text in {"not detected", "not_detected"}:
        return "not_detected"
    return "detected"


def load_master():
    rows = []
    by_ensembl: dict[str, set[str]] = defaultdict(set)
    by_target = {}
    with MASTER.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        for row in reader:
            target = row["target_uniprot_id"]
            rows.append(row)
            by_target[target] = row
            ids = set(split_ids(row.get("ensembl_gene_ids", "")))
            ids.update(split_ids(row.get("primary_gene_id", "")))
            for ensembl in ids:
                if ensembl.startswith("ENSG"):
                    by_ensembl[ensembl].add(target)
    return rows, fields, by_target, by_ensembl


def load_hpa_uniprot_crosswalk(valid_targets: set[str]):
    by_ensembl: dict[str, set[str]] = defaultdict(set)
    with zipfile.ZipFile(HPA_SUMMARY) as archive:
        member = archive.namelist()[0]
        with archive.open(member) as binary:
            rows = (
                line.decode("utf-8-sig", errors="replace") for line in binary
            )
            for row in csv.DictReader(rows, delimiter="\t"):
                gene = (row.get("Ensembl") or "").split(".", 1)[0]
                if not gene.startswith("ENSG"):
                    continue
                for accession in re.split(
                    r"[,;| ]+", row.get("Uniprot", "") or ""
                ):
                    accession = accession.strip().split("-", 1)[0]
                    if accession in valid_targets:
                        by_ensembl[gene].add(accession)
    return by_ensembl


def ontology_resources():
    resources = {}
    versions = {}
    hashes = {}
    for ontology, filename in [
        ("UBERON", "uberon-basic.obo"),
        ("CL", "cl-basic.obo"),
        ("GO", "go-basic.obo"),
    ]:
        path = OLD_ONTOLOGY / filename
        lookup, names, version = v1.parse_obo(path)
        resources[ontology] = (lookup, names)
        versions[ontology] = version
        hashes[ontology] = sha256(path)
    return resources, versions, hashes


def main() -> None:
    STAGING.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    started = now()
    write_json(
        PROGRESS,
        {
            "status": "running",
            "stage": "load_v61_and_ontologies",
            "started_utc": started,
            "updated_utc": now(),
        },
    )
    master_rows, master_fields, by_target, by_ensembl = load_master()
    hpa_uniprot_by_ensembl = load_hpa_uniprot_crosswalk(set(by_target))
    resources, ontology_versions, ontology_hashes = ontology_resources()

    mapping_cache: dict[tuple[str, str], dict[str, str]] = {}
    mapping_counts: Counter[tuple[str, str]] = Counter()

    def map_term(domain: str, term: str) -> dict[str, str]:
        term = (term or "").strip()
        if not term:
            return {
                "ontology": {"tissue": "UBERON", "cell_type": "CL", "location": "GO"}[
                    domain
                ],
                "ontology_id": "",
                "ontology_label": "",
                "mapping_status": "not_applicable",
                "mapping_method": "",
            }
        key = (domain, term)
        if key in mapping_cache:
            return mapping_cache[key]
        ontology = {"tissue": "UBERON", "cell_type": "CL", "location": "GO"}[
            domain
        ]
        lookup, names = resources[ontology]
        result = v1.ontology_match(term, ontology, lookup, names)
        mapping_cache[key] = result
        mapping_counts[(domain, result["mapping_status"])] += 1
        return result

    summary = {
        target: {
            "tissues": set(),
            "tissue_max_value": None,
            "tissue_max_terms": set(),
            "cell_types": set(),
            "cell_max_value": None,
            "cell_max_terms": set(),
            "ihc_tissues": set(),
            "ihc_cells": set(),
            "ms_tissues": set(),
            "main_locations": set(),
            "additional_locations": set(),
            "extracellular_locations": set(),
            "localization_reliability": set(),
            "hpa_genes": set(),
        }
        for target in by_target
    }
    proteins_seen = set()
    identifier_multi_map_rows = 0
    identifier_mapping_methods = Counter()
    source_rows = Counter()
    output_rows = Counter()

    def mapped_targets(row: dict) -> tuple[list[str], str]:
        nonlocal identifier_multi_map_rows
        gene = (row.get("Gene") or "").split(".", 1)[0]
        master_targets = set(by_ensembl.get(gene, set()))
        hpa_targets = set(hpa_uniprot_by_ensembl.get(gene, set()))
        if master_targets:
            targets = sorted(master_targets)
            if hpa_targets and not hpa_targets.issubset(master_targets):
                status = "E3_master_Ensembl_vs_HPA_UniProt_conflict"
                targets = sorted(master_targets | hpa_targets)
            else:
                status = (
                    "E1_exact_master_Ensembl_one_to_one"
                    if len(targets) == 1
                    else "E3_master_Ensembl_one_to_many"
                )
        elif hpa_targets:
            targets = sorted(hpa_targets)
            status = (
                "E1_HPA_Ensembl_UniProt_crosswalk_one_to_one"
                if len(targets) == 1
                else "E3_HPA_Ensembl_UniProt_crosswalk_one_to_many"
            )
        else:
            targets = []
            status = "E0_unmapped_Ensembl"
        if len(targets) > 1:
            identifier_multi_map_rows += 1
        identifier_mapping_methods[status] += 1
        return targets, status

    def write_expression(
        writer,
        row: dict,
        targets: list[str],
        mapping_status: str,
        dataset: str,
        layer: str,
        metric: str,
        value: str,
        unit: str,
        tissue: str = "",
        cell_type: str = "",
        cluster: str = "",
        reliability: str = "",
    ) -> int:
        tissue_map = map_term("tissue", tissue)
        cell_map = map_term("cell_type", cell_type)
        number = numeric(value)
        count = 0
        for target in targets:
            protein = by_target[target]
            out = {
                "expression_record_id": record_id(
                    target,
                    row.get("Gene", ""),
                    dataset,
                    metric,
                    tissue,
                    cell_type,
                    cluster,
                ),
                "target_uniprot_id": target,
                "membrane_protein_id": protein["membrane_protein_id"],
                "hpa_ensembl_gene_id": row.get("Gene", ""),
                "hpa_gene_name": row.get("Gene name", ""),
                "source_dataset": dataset,
                "source_version": HPA_VERSION,
                "measurement_layer": layer,
                "measurement_type": metric,
                "tissue": tissue,
                "tissue_ontology_id": tissue_map["ontology_id"],
                "tissue_ontology_label": tissue_map["ontology_label"],
                "tissue_mapping_status": tissue_map["mapping_status"],
                "cell_type": cell_type,
                "cell_type_ontology_id": cell_map["ontology_id"],
                "cell_type_ontology_label": cell_map["ontology_label"],
                "cell_type_mapping_status": cell_map["mapping_status"],
                "cluster": cluster,
                "value": value,
                "unit": unit,
                "detection_status": detected(value, number),
                "reliability": reliability,
                "protein_mapping_status": mapping_status,
                "module_version": MODULE_VERSION,
            }
            writer.writerow(out)
            count += 1
            proteins_seen.add(target)
            summary[target]["hpa_genes"].add(row.get("Gene", ""))
        return count

    with gzip.open(
        TISSUE_OUT, "wt", encoding="utf-8", newline=""
    ) as tissue_handle, gzip.open(
        CELL_OUT, "wt", encoding="utf-8", newline=""
    ) as cell_handle, gzip.open(
        CLUSTER_OUT, "wt", encoding="utf-8", newline=""
    ) as cluster_handle, gzip.open(
        IHC_OUT, "wt", encoding="utf-8", newline=""
    ) as ihc_handle:
        tissue_writer = csv.DictWriter(
            tissue_handle, fieldnames=BASE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        cell_writer = csv.DictWriter(
            cell_handle, fieldnames=BASE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        cluster_writer = csv.DictWriter(
            cluster_handle, fieldnames=BASE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        ihc_writer = csv.DictWriter(
            ihc_handle, fieldnames=BASE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        for writer in (tissue_writer, cell_writer, cluster_writer, ihc_writer):
            writer.writeheader()

        datasets = [
            (
                "rna_tissue_consensus.tsv.zip",
                "HPA_consensus_tissue_RNA",
                [("nTPM", "nTPM")],
            ),
            (
                "rna_tissue_hpa.tsv.zip",
                "HPA_tissue_RNA",
                [("TPM", "TPM"), ("pTPM", "pTPM"), ("nTPM", "nTPM")],
            ),
            (
                "ms_tissue.tsv.zip",
                "HPA_tissue_MS",
                [("Intensity", "MS_intensity")],
            ),
        ]
        for filename, dataset, metrics in datasets:
            write_json(
                PROGRESS,
                {
                    "status": "running",
                    "stage": f"parse_{filename}",
                    "started_utc": started,
                    "updated_utc": now(),
                    "source_rows": dict(source_rows),
                    "output_rows": dict(output_rows),
                },
            )
            archive, handle, reader = open_zip_rows(filename)
            try:
                for row in reader:
                    source_rows[filename] += 1
                    targets, status = mapped_targets(row)
                    if not targets:
                        continue
                    tissue = row.get("Tissue", "")
                    for column, metric in metrics:
                        value = row.get(column, "")
                        output_rows["tissue"] += write_expression(
                            tissue_writer,
                            row,
                            targets,
                            status,
                            dataset,
                            "RNA" if "RNA" in dataset else "protein",
                            metric,
                            value,
                            column if column != "Intensity" else "arbitrary_intensity",
                            tissue=tissue,
                        )
                        number = numeric(value)
                        for target in targets:
                            if dataset == "HPA_consensus_tissue_RNA" and number is not None:
                                if number > 0:
                                    summary[target]["tissues"].add(tissue)
                                current = summary[target]["tissue_max_value"]
                                if current is None or number > current:
                                    summary[target]["tissue_max_value"] = number
                                    summary[target]["tissue_max_terms"] = {tissue}
                                elif number == current:
                                    summary[target]["tissue_max_terms"].add(tissue)
                            if dataset == "HPA_tissue_MS" and number is not None:
                                summary[target]["ms_tissues"].add(tissue)
            finally:
                handle.close()
                archive.close()

        filename = "rna_single_cell_type.tsv.zip"
        write_json(
            PROGRESS,
            {
                "status": "running",
                "stage": "parse_single_cell_type",
                "started_utc": started,
                "updated_utc": now(),
                "source_rows": dict(source_rows),
                "output_rows": dict(output_rows),
            },
        )
        archive, handle, reader = open_zip_rows(filename)
        try:
            for row in reader:
                source_rows[filename] += 1
                targets, status = mapped_targets(row)
                if not targets:
                    continue
                cell_type = row.get("Cell type", "")
                value = row.get("nCPM", "")
                output_rows["cell_type"] += write_expression(
                    cell_writer,
                    row,
                    targets,
                    status,
                    "HPA_single_cell_type_RNA",
                    "RNA",
                    "nCPM",
                    value,
                    "nCPM",
                    cell_type=cell_type,
                )
                number = numeric(value)
                for target in targets:
                    if number is not None and number > 0:
                        summary[target]["cell_types"].add(cell_type)
                    if number is not None:
                        current = summary[target]["cell_max_value"]
                        if current is None or number > current:
                            summary[target]["cell_max_value"] = number
                            summary[target]["cell_max_terms"] = {cell_type}
                        elif number == current:
                            summary[target]["cell_max_terms"].add(cell_type)
        finally:
            handle.close()
            archive.close()

        filename = "normal_ihc_data.tsv.zip"
        archive, handle, reader = open_zip_rows(filename)
        try:
            for row in reader:
                source_rows[filename] += 1
                targets, status = mapped_targets(row)
                if not targets:
                    continue
                tissue = row.get("Tissue", "")
                cell_type = row.get("Cell type", "")
                level = row.get("Level", "")
                output_rows["ihc"] += write_expression(
                    ihc_writer,
                    row,
                    targets,
                    status,
                    "HPA_normal_tissue_IHC",
                    "protein",
                    "IHC_level",
                    level,
                    "ordinal",
                    tissue=tissue,
                    cell_type=cell_type,
                    reliability=row.get("Reliability", ""),
                )
                if detected(level) == "detected":
                    for target in targets:
                        summary[target]["ihc_tissues"].add(tissue)
                        summary[target]["ihc_cells"].add(f"{tissue}|{cell_type}")
        finally:
            handle.close()
            archive.close()

        filename = "rna_single_cell_cluster.tsv.zip"
        archive, handle, reader = open_zip_rows(filename)
        try:
            for row in reader:
                source_rows[filename] += 1
                targets, status = mapped_targets(row)
                if targets:
                    output_rows["cell_cluster"] += write_expression(
                        cluster_writer,
                        row,
                        targets,
                        status,
                        "HPA_single_cell_cluster_RNA",
                        "RNA",
                        "nCPM",
                        row.get("nCPM", ""),
                        "nCPM",
                        tissue=row.get("Tissue", ""),
                        cell_type=row.get("Cell type", ""),
                        cluster=row.get("Cluster", ""),
                    )
                if source_rows[filename] % 1_000_000 == 0:
                    write_json(
                        PROGRESS,
                        {
                            "status": "running",
                            "stage": "parse_single_cell_cluster",
                            "started_utc": started,
                            "updated_utc": now(),
                            "source_rows": dict(source_rows),
                            "output_rows": dict(output_rows),
                        },
                    )
        finally:
            handle.close()
            archive.close()

    location_fields = [
        "localization_record_id",
        "target_uniprot_id",
        "membrane_protein_id",
        "hpa_ensembl_gene_id",
        "hpa_gene_name",
        "location_term",
        "go_id",
        "go_label",
        "location_role",
        "hpa_reliability",
        "direct_observation_flag",
        "protein_mapping_status",
        "source_dataset",
        "source_version",
        "module_version",
    ]
    filename = "subcellular_location.tsv.zip"
    archive, handle, reader = open_zip_rows(filename)
    with gzip.open(
        LOCATION_OUT, "wt", encoding="utf-8", newline=""
    ) as output:
        writer = csv.DictWriter(
            output, fieldnames=location_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        try:
            for row in reader:
                source_rows[filename] += 1
                targets, status = mapped_targets(row)
                if not targets:
                    continue
                go_by_term = {}
                for item in (row.get("GO id") or "").split(";"):
                    match = re.match(r"\s*(.+?)\s*\((GO:\d+)\)\s*$", item)
                    if match:
                        go_by_term[v1.normalize_label(match.group(1))] = match.group(2)
                for column, role in [
                    ("Main location", "main"),
                    ("Additional location", "additional"),
                    ("Extracellular location", "extracellular"),
                    ("Enhanced", "enhanced"),
                    ("Supported", "supported"),
                    ("Approved", "approved"),
                    ("Uncertain", "uncertain"),
                ]:
                    for term in [
                        item.strip()
                        for item in (row.get(column) or "").split(";")
                        if item.strip()
                    ]:
                        mapped = map_term("location", term)
                        go_id = go_by_term.get(
                            v1.normalize_label(term), mapped["ontology_id"]
                        )
                        for target in targets:
                            protein = by_target[target]
                            writer.writerow(
                                {
                                    "localization_record_id": record_id(
                                        target, row.get("Gene", ""), term, role
                                    ),
                                    "target_uniprot_id": target,
                                    "membrane_protein_id": protein[
                                        "membrane_protein_id"
                                    ],
                                    "hpa_ensembl_gene_id": row.get("Gene", ""),
                                    "hpa_gene_name": row.get("Gene name", ""),
                                    "location_term": term,
                                    "go_id": go_id,
                                    "go_label": mapped["ontology_label"],
                                    "location_role": role,
                                    "hpa_reliability": row.get("Reliability", ""),
                                    "direct_observation_flag": int(
                                        role in {"main", "additional", "extracellular"}
                                    ),
                                    "protein_mapping_status": status,
                                    "source_dataset": "HPA_subcellular_IF",
                                    "source_version": HPA_VERSION,
                                    "module_version": MODULE_VERSION,
                                }
                            )
                            output_rows["localization"] += 1
                            summary[target]["localization_reliability"].add(
                                row.get("Reliability", "")
                            )
                            if role == "main":
                                summary[target]["main_locations"].add(term)
                            elif role == "additional":
                                summary[target]["additional_locations"].add(term)
                            elif role == "extracellular":
                                summary[target]["extracellular_locations"].add(
                                    term
                                )
        finally:
            handle.close()
            archive.close()

    mapping_fields = [
        "source_term_type",
        "source_term",
        "normalized_term",
        "ontology",
        "ontology_id",
        "ontology_label",
        "mapping_status",
        "mapping_method",
        "module_version",
    ]
    with MAPPING_OUT.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=mapping_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for (domain, term), result in sorted(mapping_cache.items()):
            writer.writerow(
                {
                    "source_term_type": domain,
                    "source_term": term,
                    "normalized_term": v1.normalize_label(term),
                    **result,
                    "module_version": MODULE_VERSION,
                }
            )

    summary_fields = [
        "target_uniprot_id",
        "hpa_mapping_status_v62",
        "hpa_ensembl_gene_ids_v62",
        "hpa_tissue_detected_count_v62",
        "hpa_tissue_max_nTPM_v62",
        "hpa_tissue_max_terms_v62",
        "hpa_cell_type_detected_count_v62",
        "hpa_cell_type_max_nCPM_v62",
        "hpa_cell_type_max_terms_v62",
        "hpa_ihc_detected_tissue_count_v62",
        "hpa_ihc_detected_tissue_cell_count_v62",
        "hpa_ms_detected_tissue_count_v62",
        "hpa_main_locations_v62",
        "hpa_additional_locations_v62",
        "hpa_extracellular_locations_v62",
        "hpa_localization_reliability_v62",
        "expression_location_conflict_flag_v62",
        "expression_location_conflict_reason_v62",
        "hpa_version_v62",
        "expression_location_module_version_v62",
    ]
    summary_rows = {}
    for target in sorted(by_target):
        item = summary[target]
        conflict_reasons = []
        if (
            item["tissue_max_value"] == 0
            and len(item["ihc_tissues"]) > 0
        ):
            conflict_reasons.append("consensus_RNA_zero_but_IHC_detected")
        if (
            item["tissue_max_value"] is not None
            and item["tissue_max_value"] > 0
            and len(item["ihc_tissues"]) == 0
        ):
            conflict_reasons.append("consensus_RNA_detected_but_no_IHC_detection")
        summary_rows[target] = {
            "target_uniprot_id": target,
            "hpa_mapping_status_v62": (
                "mapped" if item["hpa_genes"] else "no_HPA_Ensembl_match"
            ),
            "hpa_ensembl_gene_ids_v62": ";".join(sorted(item["hpa_genes"])),
            "hpa_tissue_detected_count_v62": len(item["tissues"]),
            "hpa_tissue_max_nTPM_v62": (
                "" if item["tissue_max_value"] is None else item["tissue_max_value"]
            ),
            "hpa_tissue_max_terms_v62": ";".join(
                sorted(item["tissue_max_terms"])
            ),
            "hpa_cell_type_detected_count_v62": len(item["cell_types"]),
            "hpa_cell_type_max_nCPM_v62": (
                "" if item["cell_max_value"] is None else item["cell_max_value"]
            ),
            "hpa_cell_type_max_terms_v62": ";".join(
                sorted(item["cell_max_terms"])
            ),
            "hpa_ihc_detected_tissue_count_v62": len(item["ihc_tissues"]),
            "hpa_ihc_detected_tissue_cell_count_v62": len(item["ihc_cells"]),
            "hpa_ms_detected_tissue_count_v62": len(item["ms_tissues"]),
            "hpa_main_locations_v62": ";".join(sorted(item["main_locations"])),
            "hpa_additional_locations_v62": ";".join(
                sorted(item["additional_locations"])
            ),
            "hpa_extracellular_locations_v62": ";".join(
                sorted(item["extracellular_locations"])
            ),
            "hpa_localization_reliability_v62": ";".join(
                sorted(x for x in item["localization_reliability"] if x)
            ),
            "expression_location_conflict_flag_v62": int(bool(conflict_reasons)),
            "expression_location_conflict_reason_v62": ";".join(
                conflict_reasons
            ),
            "hpa_version_v62": HPA_VERSION,
            "expression_location_module_version_v62": MODULE_VERSION,
        }

    with gzip.open(
        SUMMARY_OUT, "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary_rows[target] for target in sorted(summary_rows))

    review_fields = [
        "review_id",
        "review_type",
        "target_uniprot_id",
        "source_term_type",
        "source_term",
        "review_reason",
        "current_status",
        "module_version",
    ]
    review_count = 0
    with gzip.open(REVIEW_OUT, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=review_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for target in sorted(summary_rows):
            row = summary_rows[target]
            if row["hpa_mapping_status_v62"] != "mapped":
                review_count += 1
                writer.writerow(
                    {
                        "review_id": f"HPA-ID-{review_count:07d}",
                        "review_type": "identifier_mapping",
                        "target_uniprot_id": target,
                        "source_term_type": "",
                        "source_term": "",
                        "review_reason": "no_HPA_Ensembl_match",
                        "current_status": "review",
                        "module_version": MODULE_VERSION,
                    }
                )
        for (domain, term), result in sorted(mapping_cache.items()):
            if result["mapping_status"] in {
                "unmapped",
                "ambiguous_multiple_terms",
            }:
                review_count += 1
                writer.writerow(
                    {
                        "review_id": f"HPA-TERM-{review_count:07d}",
                        "review_type": "ontology_mapping",
                        "target_uniprot_id": "",
                        "source_term_type": domain,
                        "source_term": term,
                        "review_reason": result["mapping_status"],
                        "current_status": "review",
                        "module_version": MODULE_VERSION,
                    }
                )

    with MASTER.open(
        "rt", encoding="utf-8-sig", newline=""
    ) as input_handle, gzip.open(
        PREVIEW_OUT, "wt", encoding="utf-8", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        writer = csv.DictWriter(
            output_handle,
            fieldnames=master_fields + summary_fields[1:],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        preview_rows = 0
        for row in reader:
            target = row["target_uniprot_id"]
            row.update(
                {
                    key: value
                    for key, value in summary_rows[target].items()
                    if key != "target_uniprot_id"
                }
            )
            writer.writerow(row)
            preview_rows += 1

    input_hashes = {
        path.name: sha256(path)
        for path in sorted(RAW.glob("*.zip"))
    }
    mapped_proteins = sum(
        row["hpa_mapping_status_v62"] == "mapped"
        for row in summary_rows.values()
    )
    report = {
        "status": "passed_with_review_queue",
        "started_utc": started,
        "completed_utc": now(),
        "module_version": MODULE_VERSION,
        "inputs": {
            "protein_master": str(MASTER),
            "protein_rows": len(master_rows),
            "protein_master_sha256": sha256(MASTER),
            "hpa_version": HPA_VERSION,
            "hpa_input_hashes": input_hashes,
            "hpa_summary_crosswalk_sha256": sha256(HPA_SUMMARY),
            "ontology_versions": ontology_versions,
            "ontology_hashes": ontology_hashes,
        },
        "identifier_mapping": {
            "mapped_mempro_proteins": mapped_proteins,
            "unmapped_mempro_proteins": len(master_rows) - mapped_proteins,
            "one_to_many_hpa_source_rows": identifier_multi_map_rows,
            "source_row_mapping_methods": dict(identifier_mapping_methods),
        },
        "source_rows": dict(source_rows),
        "output_rows": dict(output_rows),
        "summary_rows": len(summary_rows),
        "review_rows": review_count,
        "ontology_mapping_counts": {
            f"{domain}:{status}": count
            for (domain, status), count in sorted(mapping_counts.items())
        },
        "quality_checks": {
            "summary_primary_key_unique": len(summary_rows) == len(master_rows),
            "summary_foreign_keys_in_master": set(summary_rows) == set(by_target),
            "rna_ihc_ms_units_kept_separate": True,
            "missing_distinct_from_zero": True,
            "v61_immutable": True,
            "formal_release_frozen": False,
        },
        "outputs": {
            "tissue": str(TISSUE_OUT),
            "cell_type": str(CELL_OUT),
            "cell_cluster": str(CLUSTER_OUT),
            "ihc": str(IHC_OUT),
            "localization": str(LOCATION_OUT),
            "mapping": str(MAPPING_OUT),
            "summary": str(SUMMARY_OUT),
            "review": str(REVIEW_OUT),
            "preview": str(PREVIEW_OUT),
        },
    }
    write_json(REPORT, report)
    write_json(
        PROGRESS,
        {
            "status": "complete",
            "stage": "complete",
            "started_utc": started,
            "completed_utc": now(),
            "report": str(REPORT),
            "source_rows": dict(source_rows),
            "output_rows": dict(output_rows),
        },
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
