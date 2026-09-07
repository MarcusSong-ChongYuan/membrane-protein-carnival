#!/usr/bin/env python3
"""Build a normalized gene/canonical protein/isoform identity layer.

The script is deliberately path-agnostic. Frozen input paths and the output
directory are supplied on the command line so the workflow is portable.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator


VERSION = "0.1.0"
ISOFORM_RE = re.compile(r"^([A-Z0-9]+)-(\d+)$")
UNIPROT_HEADER_RE = re.compile(r"^>([^|]+)\|([^|]+)\|([^ ]+)(?:\s+(.*))?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-master", type=Path, required=True)
    parser.add_argument("--binding-evidence", type=Path, required=True)
    parser.add_argument("--binding-sites", type=Path, required=True)
    parser.add_argument("--disease-relations", type=Path, required=True)
    parser.add_argument("--complex-components", type=Path, required=True)
    parser.add_argument("--uniprot-fasta", type=Path, required=True)
    parser.add_argument("--uniprot-headers", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8-sig", newline="")


def split_ids(value: str) -> list[str]:
    if not value:
        return []
    return sorted({part.strip() for part in re.split(r"[;,|]", value) if part.strip()})


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        count = 0
        for row in rows:
            writer.writerow({key: "" if value is None else value for key, value in row.items()})
            count += 1
    return count


def fasta_records(path: Path) -> Iterator[tuple[str, str, str, str, str]]:
    header = None
    sequence: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    match = UNIPROT_HEADER_RE.match(header)
                    if not match:
                        raise ValueError(f"Unexpected FASTA header: {header}")
                    yield match.group(1), match.group(2), match.group(3), match.group(4) or "", "".join(sequence)
                header = line
                sequence = []
            elif line:
                sequence.append(line.strip())
        if header is not None:
            match = UNIPROT_HEADER_RE.match(header)
            if not match:
                raise ValueError(f"Unexpected FASTA header: {header}")
            yield match.group(1), match.group(2), match.group(3), match.group(4) or "", "".join(sequence)


def extract_uniprot_metadata(headers_path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for line in headers_path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip().lower()
        if key in {"x-uniprot-release", "x-uniprot-release-date", "x-api-deployment-date", "date"}:
            metadata[key] = value.strip()
    return metadata


def choose_gene_entity(row: dict[str, str]) -> tuple[str, str, list[str], list[str], list[str]]:
    hgnc = split_ids(row.get("hgnc_ids", ""))
    ensembl = split_ids(row.get("ensembl_gene_ids", ""))
    ncbi = split_ids(row.get("ncbi_gene_ids", ""))
    symbol = (row.get("approved_symbol") or "").strip()
    if len(hgnc) == 1:
        return hgnc[0], "exact_single_hgnc", hgnc, ensembl, ncbi
    if not hgnc and len(ensembl) == 1:
        return ensembl[0], "exact_single_ensembl", hgnc, ensembl, ncbi
    if not hgnc and not ensembl and len(ncbi) == 1:
        return f"NCBIGene:{ncbi[0]}", "exact_single_ncbi_gene", hgnc, ensembl, ncbi
    signature = "|".join([symbol, ";".join(hgnc), ";".join(ensembl), ";".join(ncbi)])
    if hgnc or ensembl or ncbi:
        return f"GENE_REVIEW:{sha256_text(signature)[:16]}", "multiple_or_conflicting_gene_ids_review", hgnc, ensembl, ncbi
    if symbol:
        return f"SYMBOL:{symbol}", "symbol_only_provisional", hgnc, ensembl, ncbi
    return f"GENE_UNRESOLVED:{sha256_text(row['target_uniprot_id'])[:16]}", "no_gene_identifier", hgnc, ensembl, ncbi


def build_gene_and_canonical(protein_path: Path):
    proteins: dict[str, dict[str, str]] = {}
    gene_accumulator: dict[str, dict[str, set[str] | Counter]] = {}
    links: list[dict[str, str]] = []
    canonical_rows: list[dict[str, str | int]] = []
    for row in read_tsv(protein_path):
        accession = row["target_uniprot_id"].strip()
        if accession in proteins:
            raise ValueError(f"Duplicate canonical UniProt accession: {accession}")
        proteins[accession] = row
        gene_id, mapping_status, hgnc, ensembl, ncbi = choose_gene_entity(row)
        if gene_id not in gene_accumulator:
            gene_accumulator[gene_id] = {
                "symbols": set(),
                "hgnc": set(),
                "ensembl": set(),
                "ncbi": set(),
                "mapping_status": Counter(),
                "proteins": set(),
            }
        record = gene_accumulator[gene_id]
        if row.get("approved_symbol"):
            record["symbols"].add(row["approved_symbol"].strip())  # type: ignore[union-attr]
        record["hgnc"].update(hgnc)  # type: ignore[union-attr]
        record["ensembl"].update(ensembl)  # type: ignore[union-attr]
        record["ncbi"].update(ncbi)  # type: ignore[union-attr]
        record["mapping_status"][mapping_status] += 1  # type: ignore[index]
        record["proteins"].add(accession)  # type: ignore[union-attr]
        protein_entity_id = f"UNIPROT:{accession}"
        links.append(
            {
                "canonical_protein_entity_id": protein_entity_id,
                "canonical_uniprot_accession": accession,
                "gene_entity_id": gene_id,
                "link_type": "canonical_protein_to_gene",
                "link_status": mapping_status,
                "same_gene_does_not_merge_proteins": "1",
            }
        )
        canonical_rows.append(
            {
                "canonical_protein_entity_id": protein_entity_id,
                "canonical_uniprot_accession": accession,
                "membrane_protein_id": row.get("membrane_protein_id", ""),
                "gene_entity_id": gene_id,
                "approved_symbol": row.get("approved_symbol", ""),
                "reviewed_status": row.get("reviewed", ""),
                "sequence_length": len(row.get("canonical_sequence", "")),
                "canonical_sequence": row.get("canonical_sequence", ""),
                "canonical_sequence_sha256": row.get("sequence_sha256_v53", "") or sha256_text(row.get("canonical_sequence", "")),
                "membrane_class": row.get("membrane_class_v52", ""),
                "membrane_evidence_level": row.get("evidence_level_v52", ""),
                "website_default": row.get("website_default_v52", ""),
                "source_release": "V6.2",
            }
        )
    gene_rows: list[dict[str, str | int]] = []
    identifier_rows: list[dict[str, str]] = []
    identifier_owners: dict[tuple[str, str], set[str]] = defaultdict(set)
    for gene_id, record in gene_accumulator.items():
        symbols = sorted(record["symbols"])  # type: ignore[arg-type]
        hgnc = sorted(record["hgnc"])  # type: ignore[arg-type]
        ensembl = sorted(record["ensembl"])  # type: ignore[arg-type]
        ncbi = sorted(record["ncbi"])  # type: ignore[arg-type]
        statuses = record["mapping_status"]  # type: ignore[assignment]
        status = statuses.most_common(1)[0][0]  # type: ignore[union-attr]
        conflict = int(len(symbols) > 1 or "review" in status)
        gene_rows.append(
            {
                "gene_entity_id": gene_id,
                "approved_symbol": ";".join(symbols),
                "hgnc_ids": ";".join(hgnc),
                "ensembl_gene_ids": ";".join(ensembl),
                "ncbi_gene_ids": ";".join(ncbi),
                "canonical_protein_count": len(record["proteins"]),  # type: ignore[arg-type]
                "gene_mapping_status": status,
                "gene_mapping_review_flag": conflict,
            }
        )
        for identifier_type, values in (("HGNC", hgnc), ("EnsemblGene", ensembl), ("NCBIGene", ncbi)):
            for identifier in values:
                identifier_rows.append(
                    {
                        "gene_entity_id": gene_id,
                        "identifier_type": identifier_type,
                        "identifier": identifier,
                        "identifier_status": "source_exact",
                    }
                )
                identifier_owners[(identifier_type, identifier)].add(gene_id)
    for row in identifier_rows:
        owners = identifier_owners[(row["identifier_type"], row["identifier"])]
        if len(owners) > 1:
            row["identifier_status"] = "cross_entity_collision_review"
    return proteins, gene_rows, identifier_rows, links, canonical_rows


def build_isoforms(
    fasta_path: Path,
    proteins: dict[str, dict[str, str]],
    metadata: dict[str, str],
):
    canonical_from_fasta: dict[str, str] = {}
    isoforms: dict[str, dict[str, object]] = {}
    fasta_total = 0
    relevant_fasta = 0
    for database, accession, entry_name, description, sequence in fasta_records(fasta_path):
        fasta_total += 1
        match = ISOFORM_RE.match(accession)
        base = match.group(1) if match else accession
        if base not in proteins:
            continue
        relevant_fasta += 1
        if not match:
            canonical_from_fasta[accession] = sequence
            continue
        if accession in isoforms:
            raise ValueError(f"Duplicate UniProt isoform in FASTA: {accession}")
        canonical_seq = proteins[base].get("canonical_sequence", "")
        isoforms[accession] = {
            "protein_isoform_entity_id": f"UNIPROT_ISOFORM:{accession}",
            "isoform_uniprot_accession": accession,
            "canonical_protein_entity_id": f"UNIPROT:{base}",
            "canonical_uniprot_accession": base,
            "isoform_number": match.group(2),
            "entry_name": entry_name,
            "database_section": database,
            "description": description,
            "sequence_length": len(sequence),
            "isoform_sequence": sequence,
            "sequence_sha256": sha256_text(sequence),
            "same_sequence_as_canonical": int(bool(canonical_seq) and sequence == canonical_seq),
            "source_database": "UniProtKB",
            "source_proteome": "UP000005640",
            "source_release": metadata.get("x-uniprot-release", ""),
            "source_release_date": metadata.get("x-uniprot-release-date", ""),
            "source_accession_explicit": "1",
        }
    canonical_comparison = Counter()
    canonical_review: list[dict[str, str]] = []
    for accession, row in proteins.items():
        source_sequence = canonical_from_fasta.get(accession)
        master_sequence = row.get("canonical_sequence", "")
        if source_sequence is None:
            status = "canonical_missing_from_reference_proteome_download"
        elif source_sequence == master_sequence:
            status = "exact_sequence_match"
        else:
            status = "sequence_mismatch_review"
        canonical_comparison[status] += 1
        if status != "exact_sequence_match":
            canonical_review.append(
                {
                    "canonical_uniprot_accession": accession,
                    "status": status,
                    "master_length": str(len(master_sequence)),
                    "download_length": "" if source_sequence is None else str(len(source_sequence)),
                    "master_sha256": sha256_text(master_sequence),
                    "download_sha256": "" if source_sequence is None else sha256_text(source_sequence),
                }
            )
    return list(isoforms.values()), isoforms, canonical_comparison, canonical_review, fasta_total, relevant_fasta


def resolve_uniprot(
    target_id: str,
    canonical_ids: set[str],
    isoform_ids: set[str],
) -> tuple[str, str, str, str, str]:
    target_id = target_id.strip()
    match = ISOFORM_RE.match(target_id)
    if match:
        base = match.group(1)
        if target_id in isoform_ids and base in canonical_ids:
            return "isoform", f"UNIPROT:{base}", f"UNIPROT_ISOFORM:{target_id}", "explicit_isoform_accession_exact", "0"
        return "unresolved_isoform_review", f"UNIPROT:{base}" if base in canonical_ids else "", "", "explicit_isoform_not_in_frozen_uniprot_set", "1"
    if target_id in canonical_ids:
        return "canonical_protein", f"UNIPROT:{target_id}", "", "explicit_or_legacy_normalized_canonical_accession", "0"
    if target_id:
        return "unresolved_protein_review", "", "", "target_accession_not_in_membrane_protein_master", "1"
    return "unresolved_target_review", "", "", "blank_target_identifier", "1"


def build_binding_resolution(
    input_path: Path,
    output_path: Path,
    protein_to_gene: dict[str, str],
    isoform_ids: set[str],
) -> Counter:
    counts = Counter()
    fieldnames = [
        "evidence_id", "source_database", "source_record_id", "source_target_identifier",
        "gene_entity_id", "canonical_protein_entity_id", "protein_isoform_entity_id",
        "target_resolution_level", "resolution_rule", "historical_isoform_specificity",
        "target_assignment_status_v62", "identity_review_flag",
    ]
    canonical_ids = set(protein_to_gene)
    with open_text(input_path) as source, open_text(output_path, "wt") as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(destination, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            target = row.get("target_uniprot_id", "")
            level, canonical_fk, isoform_fk, rule, review = resolve_uniprot(target, canonical_ids, isoform_ids)
            if row.get("target_assignment_status") in {"complex_context_ambiguous", "ambiguous_complex"}:
                level = "complex_context_ambiguous"
                rule = "v62_complex_assignment_requires_complex_entity_resolution"
                review = "1"
            writer.writerow({
                "evidence_id": row.get("evidence_id", ""),
                "source_database": row.get("source_database", ""),
                "source_record_id": row.get("source_record_id", ""),
                "source_target_identifier": target,
                "gene_entity_id": protein_to_gene.get(target.split("-", 1)[0], ""),
                "canonical_protein_entity_id": canonical_fk,
                "protein_isoform_entity_id": isoform_fk,
                "target_resolution_level": level,
                "resolution_rule": rule,
                "historical_isoform_specificity": "not_recoverable_from_v62_normalized_record" if not ISOFORM_RE.match(target) else "explicit_in_source_normalized_record",
                "target_assignment_status_v62": row.get("target_assignment_status", ""),
                "identity_review_flag": review,
            })
            counts[level] += 1
    return counts


def build_site_resolution(
    input_path: Path,
    output_path: Path,
    protein_to_gene: dict[str, str],
    isoform_ids: set[str],
) -> Counter:
    counts = Counter()
    fieldnames = [
        "binding_site_instance_id", "evidence_id", "source_database", "source_target_identifier",
        "gene_entity_id", "canonical_protein_entity_id", "protein_isoform_entity_id",
        "target_resolution_level", "resolution_rule", "historical_isoform_specificity", "identity_review_flag",
    ]
    canonical_ids = set(protein_to_gene)
    with open_text(input_path) as source, open_text(output_path, "wt") as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(destination, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            target = row.get("target_uniprot_id", "")
            level, canonical_fk, isoform_fk, rule, review = resolve_uniprot(target, canonical_ids, isoform_ids)
            writer.writerow({
                "binding_site_instance_id": row.get("binding_site_instance_id", ""),
                "evidence_id": row.get("evidence_id", ""),
                "source_database": row.get("source_database", ""),
                "source_target_identifier": target,
                "gene_entity_id": protein_to_gene.get(target.split("-", 1)[0], ""),
                "canonical_protein_entity_id": canonical_fk,
                "protein_isoform_entity_id": isoform_fk,
                "target_resolution_level": level,
                "resolution_rule": rule,
                "historical_isoform_specificity": "not_recoverable_from_v62_normalized_record" if not ISOFORM_RE.match(target) else "explicit_in_source_normalized_record",
                "identity_review_flag": review,
            })
            counts[level] += 1
    return counts


def build_disease_resolution(input_path: Path, output_path: Path, protein_to_gene: dict[str, str]) -> Counter:
    fieldnames = [
        "disease_relation_id", "source_database", "source_record_id", "source_target_identifier",
        "gene_entity_id", "canonical_protein_entity_id", "projected_canonical_protein_entity_id",
        "protein_isoform_entity_id", "target_resolution_level", "resolution_rule",
        "historical_isoform_specificity", "identity_review_flag",
    ]
    counts = Counter()
    rows = []
    for row in read_tsv(input_path):
        target = row.get("target_uniprot_id", "")
        gene_id = protein_to_gene.get(target, "")
        source = row.get("source_database", "")
        if source == "Open Targets":
            level = "gene_product_unspecified"
            canonical_fk = ""
            projected = f"UNIPROT:{target}" if target in protein_to_gene else ""
            rule = "open_targets_gene_level_association_no_isoform_assertion"
            historical = "gene_level_source"
        elif target in protein_to_gene:
            level = "canonical_protein"
            canonical_fk = f"UNIPROT:{target}"
            projected = canonical_fk
            rule = "uniprotkb_accession_level_disease_annotation"
            historical = "canonical_accession_no_isoform_assertion"
        else:
            level = "unresolved_target_review"
            canonical_fk = ""
            projected = ""
            rule = "target_not_in_membrane_protein_master"
            historical = "unresolved"
        rows.append({
            "disease_relation_id": row.get("disease_relation_id", ""),
            "source_database": source,
            "source_record_id": row.get("source_record_id", ""),
            "source_target_identifier": target,
            "gene_entity_id": gene_id,
            "canonical_protein_entity_id": canonical_fk,
            "projected_canonical_protein_entity_id": projected,
            "protein_isoform_entity_id": "",
            "target_resolution_level": level,
            "resolution_rule": rule,
            "historical_isoform_specificity": historical,
            "identity_review_flag": int(level == "unresolved_target_review"),
        })
        counts[level] += 1
    write_tsv(output_path, fieldnames, rows)
    return counts


def build_expression_resolution(output_path: Path, proteins: dict[str, dict[str, str]], protein_to_gene: dict[str, str]) -> Counter:
    fieldnames = [
        "gene_entity_id", "projected_canonical_protein_entity_id", "canonical_uniprot_accession",
        "hpa_mapping_status_v62", "source_target_level", "isoform_resolution_status",
        "projection_rule", "identity_review_flag",
    ]
    rows = []
    counts = Counter()
    for accession, row in proteins.items():
        hpa_status = row.get("hpa_mapping_status_v62", "")
        if hpa_status and "unmapped" not in hpa_status.lower() and "missing" not in hpa_status.lower():
            status = "gene_product_unspecified"
            review = "0"
        else:
            status = "hpa_gene_mapping_missing_or_unmapped"
            review = "1"
        rows.append({
            "gene_entity_id": protein_to_gene[accession],
            "projected_canonical_protein_entity_id": f"UNIPROT:{accession}",
            "canonical_uniprot_accession": accession,
            "hpa_mapping_status_v62": hpa_status,
            "source_target_level": "gene",
            "isoform_resolution_status": status,
            "projection_rule": "HPA_expression_is_gene_level_never_inferred_to_specific_isoform",
            "identity_review_flag": review,
        })
        counts[status] += 1
    write_tsv(output_path, fieldnames, rows)
    return counts


def build_complex_resolution(
    input_path: Path,
    output_path: Path,
    protein_to_gene: dict[str, str],
    isoform_ids: set[str],
) -> Counter:
    counts = Counter()
    canonical_ids = set(protein_to_gene)
    with open_text(input_path) as source, open_text(output_path, "wt") as destination:
        reader = csv.DictReader(source, delimiter="\t")
        extra = [
            "component_gene_entity_id_v0_3", "component_canonical_protein_entity_id_v0_3",
            "component_isoform_entity_id_v0_3", "component_reference_level_v0_3",
            "component_resolution_rule_v0_3", "component_identity_review_flag_v0_3",
        ]
        fieldnames = list(reader.fieldnames or []) + extra
        writer = csv.DictWriter(destination, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            isoform = (row.get("component_uniprot_isoform_id") or "").strip()
            canonical = (row.get("component_uniprot_id") or "").strip()
            if isoform:
                level, canonical_fk, isoform_fk, rule, review = resolve_uniprot(isoform, canonical_ids, isoform_ids)
            elif canonical:
                level, canonical_fk, isoform_fk, rule, review = resolve_uniprot(canonical, canonical_ids, isoform_ids)
                if level == "unresolved_protein_review":
                    level = "external_canonical_protein"
                    rule = "valid_protein_component_outside_membrane_protein_master"
                    review = "0"
            else:
                level = "non_protein_or_unresolved_component"
                canonical_fk = ""
                isoform_fk = ""
                rule = "no_uniprot_identifier_no_protein_identity_inference"
                review = "1" if row.get("component_entity_type") == "protein" else "0"
            base = canonical or (ISOFORM_RE.match(isoform).group(1) if ISOFORM_RE.match(isoform) else "")
            row.update({
                "component_gene_entity_id_v0_3": protein_to_gene.get(base, ""),
                "component_canonical_protein_entity_id_v0_3": canonical_fk,
                "component_isoform_entity_id_v0_3": isoform_fk,
                "component_reference_level_v0_3": level,
                "component_resolution_rule_v0_3": rule,
                "component_identity_review_flag_v0_3": review,
            })
            writer.writerow(row)
            counts[level] += 1
    return counts


def count_rows(path: Path) -> int:
    with open_text(path) as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> None:
    args = parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    qa_dir = out.parent / "qa"
    qa_dir.mkdir(parents=True, exist_ok=True)
    metadata = extract_uniprot_metadata(args.uniprot_headers)

    proteins, genes, identifiers, links, canonical_rows = build_gene_and_canonical(args.protein_master)
    protein_to_gene = {row["canonical_uniprot_accession"]: row["gene_entity_id"] for row in links}
    isoform_rows, isoforms, canonical_compare, canonical_review, fasta_total, relevant_fasta = build_isoforms(
        args.uniprot_fasta, proteins, metadata
    )

    write_tsv(out / "gene_entity_v0_1.tsv", [
        "gene_entity_id", "approved_symbol", "hgnc_ids", "ensembl_gene_ids", "ncbi_gene_ids",
        "canonical_protein_count", "gene_mapping_status", "gene_mapping_review_flag",
    ], sorted(genes, key=lambda row: str(row["gene_entity_id"])))
    write_tsv(out / "gene_identifier_v0_1.tsv", [
        "gene_entity_id", "identifier_type", "identifier", "identifier_status",
    ], sorted(identifiers, key=lambda row: (row["gene_entity_id"], row["identifier_type"], row["identifier"])))
    write_tsv(out / "canonical_protein_gene_link_v0_1.tsv", [
        "canonical_protein_entity_id", "canonical_uniprot_accession", "gene_entity_id", "link_type",
        "link_status", "same_gene_does_not_merge_proteins",
    ], sorted(links, key=lambda row: row["canonical_uniprot_accession"]))
    write_tsv(out / "canonical_protein_entity_v0_1.tsv.gz", [
        "canonical_protein_entity_id", "canonical_uniprot_accession", "membrane_protein_id", "gene_entity_id",
        "approved_symbol", "reviewed_status", "sequence_length", "canonical_sequence",
        "canonical_sequence_sha256", "membrane_class", "membrane_evidence_level", "website_default", "source_release",
    ], sorted(canonical_rows, key=lambda row: str(row["canonical_uniprot_accession"])))
    write_tsv(out / "protein_isoform_v0_1.tsv.gz", [
        "protein_isoform_entity_id", "isoform_uniprot_accession", "canonical_protein_entity_id",
        "canonical_uniprot_accession", "isoform_number", "entry_name", "database_section", "description",
        "sequence_length", "isoform_sequence", "sequence_sha256", "same_sequence_as_canonical",
        "source_database", "source_proteome", "source_release", "source_release_date", "source_accession_explicit",
    ], sorted(isoform_rows, key=lambda row: str(row["isoform_uniprot_accession"])))
    write_tsv(out / "canonical_sequence_review_queue_v0_1.tsv", [
        "canonical_uniprot_accession", "status", "master_length", "download_length",
        "master_sha256", "download_sha256",
    ], canonical_review)

    isoform_ids = set(isoforms)
    binding_counts = build_binding_resolution(
        args.binding_evidence, out / "binding_evidence_target_resolution_v0_1.tsv.gz", protein_to_gene, isoform_ids
    )
    site_counts = build_site_resolution(
        args.binding_sites, out / "binding_site_target_resolution_v0_1.tsv.gz", protein_to_gene, isoform_ids
    )
    disease_counts = build_disease_resolution(
        args.disease_relations, out / "disease_relation_target_resolution_v0_1.tsv", protein_to_gene
    )
    expression_counts = build_expression_resolution(
        out / "expression_target_resolution_v0_1.tsv", proteins, protein_to_gene
    )
    complex_counts = build_complex_resolution(
        args.complex_components, out / "complex_target_components_v0_3.tsv.gz", protein_to_gene, isoform_ids
    )

    policy_rows = [
        {"source_case": "explicit UniProt isoform accession", "target_level": "isoform", "automatic_action": "link exact protein_isoform foreign key", "prohibited_action": "do not collapse into canonical-only evidence"},
        {"source_case": "canonical UniProt accession only", "target_level": "canonical_protein", "automatic_action": "link canonical protein", "prohibited_action": "do not infer a specific isoform"},
        {"source_case": "gene-level record only", "target_level": "gene_product_unspecified", "automatic_action": "link gene; retain projected canonical record separately", "prohibited_action": "do not assign isoform"},
        {"source_case": "same gene has multiple canonical proteins", "target_level": "gene and protein remain separate", "automatic_action": "retain one row per canonical accession", "prohibited_action": "never merge proteins by gene name"},
        {"source_case": "complex component with explicit isoform", "target_level": "isoform component", "automatic_action": "link complex component to exact isoform", "prohibited_action": "do not infer isoform from gene or protein name"},
        {"source_case": "legacy V6.2 canonicalized evidence", "target_level": "canonical_protein", "automatic_action": "retain canonical link and mark isoform specificity unrecoverable", "prohibited_action": "do not back-infer historical isoform"},
    ]
    write_tsv(out / "target_resolution_policy_v0_1.tsv", ["source_case", "target_level", "automatic_action", "prohibited_action"], policy_rows)

    validation_checks = {
        "canonical_protein_count_unchanged": len(proteins) == 10997,
        "canonical_accessions_unique": len(proteins) == len(set(proteins)),
        "canonical_gene_link_complete": len(links) == len(proteins),
        "gene_fk_complete": all(row["gene_entity_id"] in {g["gene_entity_id"] for g in genes} for row in links),
        "isoform_ids_unique": len(isoform_ids) == len(isoform_rows),
        "isoform_canonical_fk_complete": all(row["canonical_uniprot_accession"] in proteins for row in isoform_rows),
        "binding_row_count_preserved": sum(binding_counts.values()) == count_rows(args.binding_evidence),
        "binding_site_row_count_preserved": sum(site_counts.values()) == count_rows(args.binding_sites),
        "disease_row_count_preserved": sum(disease_counts.values()) == count_rows(args.disease_relations),
        "expression_projection_count_complete": sum(expression_counts.values()) == len(proteins),
        "complex_component_row_count_preserved": sum(complex_counts.values()) == count_rows(args.complex_components),
        "no_isoform_inferred_for_binding": binding_counts.get("isoform", 0) == 0,
        "no_isoform_inferred_for_sites": site_counts.get("isoform", 0) == 0,
        "open_targets_is_gene_level": disease_counts.get("gene_product_unspecified", 0) == 5651,
    }
    validation = {
        "module": "gene_canonical_isoform_complex_identity_layer",
        "module_version": VERSION,
        "candidate_release": "V6.3",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(validation_checks.values()) else "FAIL",
        "checks": validation_checks,
        "counts": {
            "gene_entities": len(genes),
            "gene_identifier_rows": len(identifiers),
            "canonical_proteins": len(proteins),
            "canonical_proteins_with_multiple_members_per_gene": sum(1 for row in genes if int(row["canonical_protein_count"]) > 1),
            "uniprot_fasta_records_total": fasta_total,
            "uniprot_fasta_records_relevant_to_master": relevant_fasta,
            "protein_isoforms": len(isoform_rows),
            "canonical_sequence_comparison": dict(canonical_compare),
            "binding_resolution": dict(binding_counts),
            "binding_site_resolution": dict(site_counts),
            "disease_resolution": dict(disease_counts),
            "expression_resolution": dict(expression_counts),
            "complex_component_resolution": dict(complex_counts),
        },
        "source_metadata": metadata,
        "methodological_limit": "V6.2 canonicalized binding/site records contain no explicit isoform accession; historical isoform specificity cannot be reconstructed from the frozen normalized tables.",
    }
    validation_path = qa_dir / "IDENTITY_LAYER_V0_1_VALIDATION.json"
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")

    input_paths = [
        args.protein_master, args.binding_evidence, args.binding_sites, args.disease_relations,
        args.complex_components, args.uniprot_fasta, args.uniprot_headers,
    ]
    output_paths = sorted([path for path in out.iterdir() if path.is_file()]) + [validation_path]
    manifest = {
        "module_version": VERSION,
        "candidate_release": "V6.3",
        "generated_at_utc": validation["generated_at_utc"],
        "inputs": [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in input_paths],
        "outputs": [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in output_paths],
    }
    (qa_dir / "IDENTITY_LAYER_V0_1_MANIFEST.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
