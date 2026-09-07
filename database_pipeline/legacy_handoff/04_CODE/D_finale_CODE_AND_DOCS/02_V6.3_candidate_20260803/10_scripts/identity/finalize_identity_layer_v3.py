#!/usr/bin/env python3
"""Resolve direct external sequences and perform final relational QA."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode.replace("t", ""), encoding="utf-8-sig", newline="")


def rows(path: Path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write(path: Path, fields: list[str], records) -> int:
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        count = 0
        for record in records:
            writer.writerow(record)
            count += 1
    return count


def fasta(path: Path) -> dict[str, tuple[str, str, str, str]]:
    result, header, seq = {}, "", []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header:
                    parts = header.split("|", 2)
                    accession = parts[1]
                    tail = parts[2]
                    result[accession] = ("".join(seq), tail.split(" ", 1)[0], parts[0][1:], tail.split(" ", 1)[1] if " " in tail else "")
                header, seq = line, []
            elif line:
                seq.append(line.strip())
        if header:
            parts = header.split("|", 2)
            accession = parts[1]
            tail = parts[2]
            result[accession] = ("".join(seq), tail.split(" ", 1)[0], parts[0][1:], tail.split(" ", 1)[1] if " " in tail else "")
    return result


def hash_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--direct-fasta", type=Path, required=True)
    p.add_argument("--direct-manifest", type=Path, required=True)
    args = p.parse_args()
    out, qa = args.root / "output_v0_1", args.root / "qa"

    direct_sequences = fasta(args.direct_fasta)
    direct_manifest = json.loads(args.direct_manifest.read_text(encoding="utf-8"))
    manifest_by_request = {row["requested_accession"]: row for row in direct_manifest["records"]}
    external = list(rows(out / "external_complex_protein_entity_v0_1.tsv.gz"))
    review = []
    for row in external:
        accession = row["canonical_uniprot_accession"]
        if row["sequence_status"] == "resolved":
            row["resolution_method"] = "reference_proteome_bulk_fasta"
            row["resolved_accession"] = accession
            continue
        manifest_row = manifest_by_request.get(accession, {})
        returned = manifest_row.get("returned_accession", "")
        if returned in direct_sequences:
            sequence, entry, section, description = direct_sequences[returned]
            row.update({
                "entry_name": entry,
                "database_section": section,
                "description": description,
                "sequence_length": str(len(sequence)),
                "canonical_sequence": sequence,
                "sequence_sha256": hash_text(sequence),
                "sequence_status": "resolved_direct_nonreference",
                "resolution_method": "direct_uniprot_accession_api",
                "resolved_accession": returned,
            })
        else:
            row.update({
                "sequence_status": "unresolved_obsolete_or_deleted_review",
                "resolution_method": "direct_uniprot_accession_api_empty_response",
                "resolved_accession": "",
            })
            review.append({
                "external_protein_entity_id": row["external_protein_entity_id"],
                "source_accession": accession,
                "review_reason": "not_in_reference_proteome_and_direct_UniProt_FASTA_empty",
                "http_status": manifest_row.get("http_status", ""),
                "recommended_action": "audit source database version; map retired accession only with authoritative UniProt ID history",
            })
    external_fields = list(external[0].keys())
    write(out / "external_complex_protein_entity_v0_2.tsv.gz", external_fields, external)
    review_fields = ["external_protein_entity_id", "source_accession", "review_reason", "http_status", "recommended_action"]
    write(out / "external_complex_protein_review_v0_2.tsv", review_fields, review)

    gene_rows = list(rows(out / "gene_entity_v0_1.tsv"))
    canonical_rows = list(rows(out / "canonical_protein_entity_v0_1.tsv.gz"))
    isoform_rows = list(rows(out / "protein_isoform_v0_2.tsv.gz"))
    component_rows = list(rows(out / "complex_target_components_v0_4.tsv.gz"))
    gene_ids = {row["gene_entity_id"] for row in gene_rows}
    canonical_ids = {row["canonical_protein_entity_id"] for row in canonical_rows}
    external_ids = {row["external_protein_entity_id"] for row in external}
    isoform_ids = {row["protein_isoform_entity_id"] for row in isoform_rows}
    protein_union = canonical_ids | external_ids

    binding_counts, binding_gene_bad, binding_protein_bad = Counter(), 0, 0
    for row in rows(out / "binding_evidence_target_resolution_v0_1.tsv.gz"):
        binding_counts[row["target_resolution_level"]] += 1
        binding_gene_bad += bool(row["gene_entity_id"] and row["gene_entity_id"] not in gene_ids)
        binding_protein_bad += bool(row["canonical_protein_entity_id"] and row["canonical_protein_entity_id"] not in canonical_ids)
    site_counts = Counter(row["target_resolution_level"] for row in rows(out / "binding_site_target_resolution_v0_1.tsv.gz"))
    disease_counts = Counter(row["target_resolution_level"] for row in rows(out / "disease_relation_target_resolution_v0_1.tsv"))
    expression_counts = Counter(row["isoform_resolution_status"] for row in rows(out / "expression_target_resolution_v0_1.tsv"))
    component_counts = Counter(row["component_reference_level_v0_4"] for row in component_rows)
    search_counts = Counter(row["entity_type"] for row in rows(out / "website_entity_search_index_v0_1.tsv.gz"))

    component_fk_bad = 0
    for row in component_rows:
        protein_fk = row["component_canonical_protein_entity_id_v0_4"]
        isoform_fk = row["component_isoform_entity_id_v0_4"]
        component_fk_bad += bool(protein_fk and protein_fk not in protein_union)
        component_fk_bad += bool(isoform_fk and isoform_fk not in isoform_ids)

    checks = {
        "gene_primary_keys_unique": len(gene_rows) == len(gene_ids),
        "canonical_primary_keys_unique": len(canonical_rows) == len(canonical_ids) == 10997,
        "canonical_gene_foreign_keys_complete": all(row["gene_entity_id"] in gene_ids for row in canonical_rows),
        "isoform_primary_keys_unique": len(isoform_rows) == len(isoform_ids),
        "isoform_parent_foreign_keys_complete": all(row["canonical_protein_entity_id"] in protein_union for row in isoform_rows),
        "binding_foreign_keys_complete": binding_gene_bad == 0 and binding_protein_bad == 0,
        "binding_legacy_rows_complete": sum(binding_counts.values()) == 3003306,
        "binding_isoform_not_inferred": binding_counts.get("isoform", 0) == 0,
        "binding_site_rows_complete": sum(site_counts.values()) == 95598,
        "disease_rows_complete": sum(disease_counts.values()) == 10264,
        "open_targets_gene_product_unspecified": disease_counts.get("gene_product_unspecified", 0) == 5651,
        "expression_is_gene_level": expression_counts.get("gene_product_unspecified", 0) == 10997,
        "complex_component_rows_complete": len(component_rows) == 97993,
        "complex_component_foreign_keys_complete": component_fk_bad == 0,
        "all_explicit_complex_isoforms_resolved": component_counts.get("unresolved_isoform_review", 0) == 0,
        "all_external_sequences_have_terminal_status": all(row["sequence_status"] in {"resolved", "resolved_direct_nonreference", "unresolved_obsolete_or_deleted_review"} for row in external),
        "unresolved_external_records_queued": len(review) == direct_manifest["unresolved_count"],
        "website_supports_four_requested_entity_types": all(search_counts[key] > 0 for key in ("gene", "canonical_protein", "protein_isoform", "protein_complex")),
    }
    validation = {
        "module": "gene_canonical_isoform_complex_identity_layer",
        "module_version": "0.3.0",
        "candidate_release": "V6.3",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": {
            "gene_entities": len(gene_rows),
            "canonical_membrane_proteins": len(canonical_rows),
            "protein_isoforms": len(isoform_rows),
            "external_complex_canonical_proteins": len(external),
            "external_sequence_status": dict(Counter(row["sequence_status"] for row in external)),
            "binding_resolution": dict(binding_counts),
            "binding_site_resolution": dict(site_counts),
            "disease_resolution": dict(disease_counts),
            "expression_resolution": dict(expression_counts),
            "complex_component_resolution": dict(component_counts),
            "website_search_index": dict(search_counts),
        },
        "publication_interpretation": {
            "core_membrane_protein_count": "unchanged at 10,997; external complex components are context entities and are excluded from membrane-protein totals",
            "historical_isoform_limit": "V6.2 binding/site records contain canonical accessions only; isoform specificity was not inferred",
            "external_accession_review": f"{len(review)} source accessions returned no FASTA in UniProt 2026_02 and remain explicit review records, not silently deleted",
        },
        "supersedes": ["IDENTITY_LAYER_V0_1_VALIDATION.json", "IDENTITY_LAYER_V0_2_VALIDATION.json"],
    }
    validation_path = qa / "IDENTITY_LAYER_V0_3_VALIDATION.json"
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    deliverables = sorted(path for path in out.iterdir() if path.is_file()) + [validation_path]
    manifest = {
        "module_version": "0.3.0",
        "candidate_release": "V6.3",
        "generated_at_utc": validation["generated_at_utc"],
        "files": [{"path": str(path), "bytes": path.stat().st_size, "sha256": hash_file(path)} for path in deliverables],
    }
    (qa / "IDENTITY_LAYER_V0_3_MANIFEST.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
