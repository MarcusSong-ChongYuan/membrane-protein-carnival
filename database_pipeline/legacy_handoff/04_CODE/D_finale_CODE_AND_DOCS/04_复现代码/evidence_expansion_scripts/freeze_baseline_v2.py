#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
BASE = Path(r"D:\7.22\small_molecule_v1_working\releases\release_small_molecule_v1_0")
PGD_BASE = Path(r"D:\7.22\membrane_master_v5_working\releases\release_v5_4_pgd_binding")

INPUTS = {
    "protein_master_v5_5": BASE / "human_membrane_audit_master_v5_5.tsv",
    "binding_evidence_v4_2": BASE / "binding_evidence_master_v4_2.tsv",
    "binding_sites_v4_2": BASE / "binding_site_instances_v4_2.tsv",
    "protein_compound_summary_v4_2": BASE / "protein_compound_summary_v4_2.tsv",
    "compound_master_v1_0": BASE / "small_molecule_master_v1_0.tsv",
    "compound_forms_v1_0": BASE / "compound_form_hierarchy_v1_0.tsv",
    "protein_gene_disease_v5_4": PGD_BASE / "protein_gene_disease_relations_v5_4.tsv",
    "small_molecule_validation_v1_0": BASE / "SMALL_MOLECULE_V1_VALIDATION_REPORT.json",
    "protein_binding_validation_v5_4": PGD_BASE / "V54_VALIDATION_REPORT.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def row_count(path: Path) -> int | None:
    if path.suffix.lower() != ".tsv":
        return None
    with path.open("rb") as handle:
        lines = sum(chunk.count(b"\n") for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""))
    return max(lines - 1, 0)


def build_protein_universe(path: Path, output: Path) -> dict:
    fields = [
        "target_uniprot_id",
        "approved_symbol",
        "protein_name",
        "primary_gene_id",
        "membrane_class_v52",
        "evidence_level_v52",
        "release_scope_v52",
        "website_default_v52",
        "functional_primary_class_v5",
        "functional_subclass_v5",
        "membrane_topology_v5",
        "ncbi_gene_ids",
        "ensembl_gene_ids",
    ]
    counters = {
        "rows": 0,
        "abc_rows": 0,
        "default_rows": 0,
        "class_counts": Counter(),
        "evidence_counts": Counter(),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with path.open("r", encoding="utf-8-sig", newline="") as source, output.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            counters["rows"] += 1
            membrane_class = row.get("membrane_class_v52", "")
            evidence = row.get("evidence_level_v52", "")
            counters["class_counts"][membrane_class] += 1
            counters["evidence_counts"][evidence] += 1
            if membrane_class not in {"A", "B", "C"}:
                continue
            counters["abc_rows"] += 1
            if row.get("website_default_v52") == "1":
                counters["default_rows"] += 1
            writer.writerow({field: row.get(field, "") for field in fields})
    return {
        "rows": counters["rows"],
        "abc_rows": counters["abc_rows"],
        "default_rows": counters["default_rows"],
        "class_counts": dict(counters["class_counts"]),
        "evidence_counts": dict(counters["evidence_counts"]),
    }


def main() -> None:
    missing = [str(path) for path in INPUTS.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing baseline inputs: " + "; ".join(missing))

    manifest_rows = []
    for key, path in INPUTS.items():
        manifest_rows.append(
            {
                "input_id": key,
                "absolute_path": str(path),
                "size_bytes": path.stat().st_size,
                "row_count": row_count(path),
                "sha256": sha256(path),
            }
        )

    report_dir = ROOT / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    with (report_dir / "baseline_file_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["input_id", "absolute_path", "size_bytes", "row_count", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    protein_stats = build_protein_universe(
        INPUTS["protein_master_v5_5"], ROOT / "intermediate" / "protein_universe_abc_v2.tsv"
    )
    with INPUTS["small_molecule_validation_v1_0"].open("r", encoding="utf-8") as handle:
        compound_validation = json.load(handle)
    with INPUTS["protein_binding_validation_v5_4"].open("r", encoding="utf-8") as handle:
        binding_validation = json.load(handle)

    report = {
        "freeze_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": "frozen",
        "baseline_inputs": {key: str(value) for key, value in INPUTS.items()},
        "protein_universe": protein_stats,
        "canonical_counts": {
            "binding_evidence_rows_v4_2": compound_validation["binding_evidence_rows_v4_2"],
            "binding_site_rows_v4_2": compound_validation["binding_site_rows_v4_2"],
            "protein_compound_pair_rows_v4_2": compound_validation[
                "protein_compound_pair_rows_v4_2"
            ],
            "parent_master_rows": compound_validation["parent_master_rows"],
            "parent_compounds_with_default_binding": compound_validation["counters"][
                "parent_compounds_with_default_binding"
            ],
            "proteins_with_canonical_core_compounds": compound_validation["counters"][
                "proteins_with_canonical_core_compounds"
            ],
            "disease_relation_rows": binding_validation["disease_relation_rows"],
            "proteins_with_supported_disease": binding_validation[
                "proteins_with_supported_disease"
            ],
        },
        "validation_inputs_passed": {
            "small_molecule_v1_0": compound_validation.get("status") == "passed",
            "protein_binding_v5_4": binding_validation.get("status") == "passed",
        },
    }
    with (report_dir / "BASELINE_FREEZE_V2.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)

    with (ROOT / "config" / "baseline_paths.json").open("w", encoding="utf-8") as handle:
        json.dump({key: str(value) for key, value in INPUTS.items()}, handle, indent=2)

    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
