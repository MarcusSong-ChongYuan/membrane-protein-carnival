#!/usr/bin/env python3
"""Extend the core identity layer to all protein components used by complexes."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ISOFORM_RE = re.compile(r"^([A-Z0-9]+)-(\d+)$")


def args_parser() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--protein-master", type=Path, required=True)
    p.add_argument("--gene-entities", type=Path, required=True)
    p.add_argument("--canonical-entities", type=Path, required=True)
    p.add_argument("--core-isoforms", type=Path, required=True)
    p.add_argument("--binding-evidence", type=Path, required=True)
    p.add_argument("--binding-resolution", type=Path, required=True)
    p.add_argument("--complex-components", type=Path, required=True)
    p.add_argument("--complex-master", type=Path, required=True)
    p.add_argument("--uniprot-fasta", type=Path, required=True)
    p.add_argument("--explicit-isoform-fasta", type=Path, required=True)
    p.add_argument("--explicit-isoform-manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--qa-dir", type=Path, required=True)
    return p.parse_args()


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8-sig", newline="")


def reader(path: Path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def row_count(path: Path) -> int:
    return sum(1 for _ in reader(path))


def sha256(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header = ""
    seq = []
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header:
                    parts = header.split("|", 2)
                    yield parts[0][1:], parts[1], parts[2].split(" ", 1)[0], parts[2].split(" ", 1)[1] if " " in parts[2] else "", "".join(seq)
                header, seq = line, []
            elif line:
                seq.append(line.strip())
        if header:
            parts = header.split("|", 2)
            yield parts[0][1:], parts[1], parts[2].split(" ", 1)[0], parts[2].split(" ", 1)[1] if " " in parts[2] else "", "".join(seq)


def write(path: Path, fields: list[str], rows) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(path, "wt") as handle:
        out = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        out.writeheader()
        for row in rows:
            out.writerow(row)
            count += 1
    return count


def main() -> None:
    args = args_parser()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.qa_dir.mkdir(parents=True, exist_ok=True)

    protein_rows = list(reader(args.protein_master))
    master_ids = {row["target_uniprot_id"] for row in protein_rows}
    protein_to_gene = {row["target_uniprot_id"]: row["gene_entity_id"] for row in reader(args.canonical_entities)}

    component_rows = list(reader(args.complex_components))
    component_bases = set()
    explicit_component_isoforms = set()
    for row in component_rows:
        canonical = (row.get("component_uniprot_id") or "").strip()
        isoform = (row.get("component_uniprot_isoform_id") or "").strip()
        if canonical:
            component_bases.add(canonical)
        match = ISOFORM_RE.match(isoform)
        if match:
            component_bases.add(match.group(1))
            explicit_component_isoforms.add(isoform)

    canonical_sequences: dict[str, tuple[str, str, str, str]] = {}
    bulk_isoforms: dict[str, dict[str, str]] = {}
    for section, accession, entry_name, description, sequence in fasta_records(args.uniprot_fasta):
        match = ISOFORM_RE.match(accession)
        base = match.group(1) if match else accession
        if base not in component_bases and base not in master_ids:
            continue
        if match:
            bulk_isoforms[accession] = {
                "protein_isoform_entity_id": f"UNIPROT_ISOFORM:{accession}",
                "isoform_uniprot_accession": accession,
                "canonical_protein_entity_id": f"UNIPROT:{base}",
                "canonical_uniprot_accession": base,
                "isoform_number": match.group(2),
                "entry_name": entry_name,
                "database_section": section,
                "description": description,
                "sequence_length": str(len(sequence)),
                "isoform_sequence": sequence,
                "sequence_sha256": sha256(sequence),
                "same_sequence_as_canonical": "",
                "source_database": "UniProtKB",
                "source_proteome": "UP000005640",
                "source_release": "2026_02",
                "source_release_date": "10-June-2026",
                "source_accession_explicit": "1",
                "entity_scope": "membrane_protein_master" if base in master_ids else "external_complex_component",
                "isoform_record_type": "bulk_reference_proteome_isoform",
            }
        else:
            canonical_sequences[accession] = (sequence, entry_name, section, description)

    combined_isoforms = {row["isoform_uniprot_accession"]: row for row in reader(args.core_isoforms)}
    for accession, row in bulk_isoforms.items():
        if accession not in combined_isoforms:
            combined_isoforms[accession] = row
    direct_isoform_count = 0
    for section, accession, entry_name, description, sequence in fasta_records(args.explicit_isoform_fasta):
        match = ISOFORM_RE.match(accession)
        if not match:
            continue
        base = match.group(1)
        direct_isoform_count += 1
        canonical_sequence = canonical_sequences.get(base, ("", "", "", ""))[0]
        row = combined_isoforms.get(accession)
        if row and row.get("sequence_sha256") and row["sequence_sha256"] != sha256(sequence):
            raise ValueError(f"UniProt sequence disagreement for {accession}")
        combined_isoforms[accession] = {
            "protein_isoform_entity_id": f"UNIPROT_ISOFORM:{accession}",
            "isoform_uniprot_accession": accession,
            "canonical_protein_entity_id": f"UNIPROT:{base}",
            "canonical_uniprot_accession": base,
            "isoform_number": match.group(2),
            "entry_name": entry_name,
            "database_section": section,
            "description": description,
            "sequence_length": str(len(sequence)),
            "isoform_sequence": sequence,
            "sequence_sha256": sha256(sequence),
            "same_sequence_as_canonical": str(int(bool(canonical_sequence) and sequence == canonical_sequence)),
            "source_database": "UniProtKB",
            "source_proteome": "UP000005640",
            "source_release": "2026_02",
            "source_release_date": "10-June-2026",
            "source_accession_explicit": "1",
            "entity_scope": "membrane_protein_master" if base in master_ids else "external_complex_component",
            "isoform_record_type": "bulk_and_direct_uniprot_isoform" if accession in bulk_isoforms else "direct_accession_api_isoform",
        }

    isoform_fields = [
        "protein_isoform_entity_id", "isoform_uniprot_accession", "canonical_protein_entity_id",
        "canonical_uniprot_accession", "isoform_number", "entry_name", "database_section", "description",
        "sequence_length", "isoform_sequence", "sequence_sha256", "same_sequence_as_canonical",
        "source_database", "source_proteome", "source_release", "source_release_date", "source_accession_explicit",
        "entity_scope", "isoform_record_type",
    ]
    write(args.output_dir / "protein_isoform_v0_2.tsv.gz", isoform_fields,
          (combined_isoforms[key] for key in sorted(combined_isoforms)))

    external_rows = []
    for accession in sorted(component_bases - master_ids):
        sequence, entry_name, section, description = canonical_sequences.get(accession, ("", "", "", ""))
        external_rows.append({
            "external_protein_entity_id": f"UNIPROT:{accession}",
            "canonical_uniprot_accession": accession,
            "entry_name": entry_name,
            "database_section": section,
            "description": description,
            "sequence_length": len(sequence),
            "canonical_sequence": sequence,
            "sequence_sha256": sha256(sequence) if sequence else "",
            "source_database": "UniProtKB",
            "source_proteome": "UP000005640",
            "source_release": "2026_02",
            "entity_scope": "external_complex_component_not_in_membrane_master",
            "sequence_status": "resolved" if sequence else "missing_review",
        })
    external_fields = [
        "external_protein_entity_id", "canonical_uniprot_accession", "entry_name", "database_section", "description",
        "sequence_length", "canonical_sequence", "sequence_sha256", "source_database", "source_proteome",
        "source_release", "entity_scope", "sequence_status",
    ]
    write(args.output_dir / "external_complex_protein_entity_v0_1.tsv.gz", external_fields, external_rows)
    external_ids = {row["canonical_uniprot_accession"] for row in external_rows}

    resolved_components = []
    component_counts = Counter()
    unresolved = []
    extra_fields = [
        "component_gene_entity_id_v0_4", "component_canonical_protein_entity_id_v0_4",
        "component_isoform_entity_id_v0_4", "component_reference_level_v0_4",
        "component_resolution_rule_v0_4", "component_identity_review_flag_v0_4",
    ]
    for row in component_rows:
        canonical = (row.get("component_uniprot_id") or "").strip()
        isoform = (row.get("component_uniprot_isoform_id") or "").strip()
        match = ISOFORM_RE.match(isoform)
        base = canonical or (match.group(1) if match else "")
        if isoform and isoform in combined_isoforms:
            level = "isoform" if base in master_ids else "external_isoform"
            canonical_fk = f"UNIPROT:{base}"
            isoform_fk = f"UNIPROT_ISOFORM:{isoform}"
            rule = "explicit_isoform_accession_verified_by_uniprot"
            review = "0"
        elif isoform:
            level, canonical_fk, isoform_fk = "unresolved_isoform_review", f"UNIPROT:{base}" if base else "", ""
            rule, review = "explicit_isoform_failed_uniprot_resolution", "1"
        elif canonical in master_ids:
            level, canonical_fk, isoform_fk = "canonical_protein", f"UNIPROT:{canonical}", ""
            rule, review = "canonical_membrane_protein_no_isoform_inference", "0"
        elif canonical in external_ids:
            level, canonical_fk, isoform_fk = "external_canonical_protein", f"UNIPROT:{canonical}", ""
            rule, review = "external_complex_protein_no_isoform_inference", "0"
        else:
            level, canonical_fk, isoform_fk = "non_protein_or_unresolved_component", "", ""
            rule = "no_resolvable_uniprot_identifier"
            review = "1" if row.get("component_entity_type") == "protein" else "0"
        row = dict(row)
        row.update({
            "component_gene_entity_id_v0_4": protein_to_gene.get(base, ""),
            "component_canonical_protein_entity_id_v0_4": canonical_fk,
            "component_isoform_entity_id_v0_4": isoform_fk,
            "component_reference_level_v0_4": level,
            "component_resolution_rule_v0_4": rule,
            "component_identity_review_flag_v0_4": review,
        })
        resolved_components.append(row)
        component_counts[level] += 1
        if review == "1":
            unresolved.append({
                "component_assertion_id": row.get("component_assertion_id", ""),
                "complex_target_id": row.get("complex_target_id", ""),
                "component_identifier": row.get("component_identifier", ""),
                "component_entity_type": row.get("component_entity_type", ""),
                "component_uniprot_id": canonical,
                "component_uniprot_isoform_id": isoform,
                "resolution_level": level,
                "review_reason": rule,
            })
    component_fields = list(component_rows[0].keys()) + extra_fields
    write(args.output_dir / "complex_target_components_v0_4.tsv.gz", component_fields, resolved_components)
    review_fields = [
        "component_assertion_id", "complex_target_id", "component_identifier", "component_entity_type",
        "component_uniprot_id", "component_uniprot_isoform_id", "resolution_level", "review_reason",
    ]
    write(args.output_dir / "complex_component_identity_review_v0_2.tsv", review_fields, unresolved)

    search_rows = []
    for row in reader(args.gene_entities):
        search_rows.append({"search_entity_id": row["gene_entity_id"], "entity_type": "gene", "primary_label": row.get("approved_symbol", ""), "aliases": ";".join(filter(None, [row.get("hgnc_ids", ""), row.get("ensembl_gene_ids", ""), row.get("ncbi_gene_ids", "")])), "parent_entity_id": "", "release_scope": "gene_layer"})
    for row in reader(args.canonical_entities):
        search_rows.append({"search_entity_id": row["canonical_protein_entity_id"], "entity_type": "canonical_protein", "primary_label": row.get("approved_symbol", "") or row["canonical_uniprot_accession"], "aliases": row["canonical_uniprot_accession"], "parent_entity_id": row.get("gene_entity_id", ""), "release_scope": "membrane_protein_master"})
    for row in combined_isoforms.values():
        search_rows.append({"search_entity_id": row["protein_isoform_entity_id"], "entity_type": "protein_isoform", "primary_label": row["isoform_uniprot_accession"], "aliases": row.get("entry_name", ""), "parent_entity_id": row["canonical_protein_entity_id"], "release_scope": row.get("entity_scope", "")})
    for row in reader(args.complex_master):
        search_rows.append({"search_entity_id": row["complex_target_id"], "entity_type": "protein_complex", "primary_label": row.get("complex_name", ""), "aliases": row.get("aliases", ""), "parent_entity_id": "", "release_scope": "complex_candidate_default" if row.get("default_release_inclusion") == "1" else "complex_audit"})
    search_fields = ["search_entity_id", "entity_type", "primary_label", "aliases", "parent_entity_id", "release_scope"]
    write(args.output_dir / "website_entity_search_index_v0_1.tsv.gz", search_fields, search_rows)

    explicit_manifest = json.loads(args.explicit_isoform_manifest.read_text(encoding="utf-8"))
    checks = {
        "canonical_membrane_protein_count_unchanged": len(master_ids) == 10997,
        "all_explicit_complex_isoforms_api_valid": explicit_manifest["failed_count"] == 0 and explicit_manifest["valid_fasta_count"] == len(explicit_component_isoforms),
        "all_explicit_complex_isoforms_resolved": all(i in combined_isoforms for i in explicit_component_isoforms),
        "complex_component_rows_preserved": len(resolved_components) == len(component_rows),
        "complex_explicit_isoform_unresolved_zero": component_counts["unresolved_isoform_review"] == 0,
        "external_canonical_sequences_complete": all(row["sequence_status"] == "resolved" for row in external_rows),
        "binding_bridge_rows_preserved": row_count(args.binding_evidence) == row_count(args.binding_resolution),
        "no_core_isoform_inferred_from_canonical_binding": all(row["target_resolution_level"] != "isoform" for row in reader(args.binding_resolution)),
    }
    counts = {
        "gene_entities": row_count(args.gene_entities),
        "canonical_membrane_proteins": len(master_ids),
        "protein_isoforms_total": len(combined_isoforms),
        "protein_isoforms_membrane_scope": sum(1 for row in combined_isoforms.values() if row.get("entity_scope") == "membrane_protein_master"),
        "protein_isoforms_external_complex_scope": sum(1 for row in combined_isoforms.values() if row.get("entity_scope") == "external_complex_component"),
        "explicit_complex_isoform_ids": len(explicit_component_isoforms),
        "external_complex_canonical_proteins": len(external_rows),
        "complex_component_rows": len(component_rows),
        "complex_component_resolution": dict(component_counts),
        "complex_identity_review_rows": len(unresolved),
        "website_search_entities": len(search_rows),
    }
    validation = {
        "module": "gene_canonical_isoform_complex_identity_layer",
        "module_version": "0.2.0",
        "candidate_release": "V6.3",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": counts,
        "supersedes": "IDENTITY_LAYER_V0_1_VALIDATION.json (physical-line counter corrected; external complex proteins added)",
        "methodological_limit": "Frozen V6.2 binding and binding-site records contain canonical accessions only, so historical isoform specificity is not recoverable and was not inferred.",
    }
    validation_path = args.qa_dir / "IDENTITY_LAYER_V0_2_VALIDATION.json"
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    outputs = sorted(path for path in args.output_dir.iterdir() if path.is_file())
    manifest = {
        "module_version": "0.2.0",
        "candidate_release": "V6.3",
        "generated_at_utc": validation["generated_at_utc"],
        "outputs": [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in outputs],
    }
    (args.qa_dir / "IDENTITY_LAYER_V0_2_OUTPUT_MANIFEST.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
