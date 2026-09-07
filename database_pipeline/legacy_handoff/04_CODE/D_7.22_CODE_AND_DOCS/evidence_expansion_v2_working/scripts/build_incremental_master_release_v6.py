#!/usr/bin/env python3
"""Build the four linked MemProDB V6 master tables from validated evidence."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import pathlib
import shutil
import sys
from collections import defaultdict


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "incremental_20260727"
OUT = RUN / "release_candidate_v6_0"
QA_DIR = RUN / "qa"
CONFIG = ROOT / "config" / "baseline_paths.json"

PAIR_SUMMARY = OUT / "protein_compound_summary_v6_0.tsv.gz"
RELEASE_VERSION = "v6.0"
RELEASE_DATE = "2026-07-27"
TIER_RANK = {"BE1": 1, "BE2": 2, "BE3": 3}


def split_values(value: str) -> list[str]:
    import re

    return [
        part.strip()
        for part in re.split(r"[;|,]", value or "")
        if part.strip()
    ]


def numeric(value: str) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def integer(value: str) -> int:
    try:
        return int(float(str(value).strip() or 0))
    except (TypeError, ValueError):
        return 0


def new_summary() -> dict[str, object]:
    return {
        "pair_count": 0,
        "best_rank": 99,
        "best_tier_counts": defaultdict(int),
        "be1_evidence": 0,
        "be2_evidence": 0,
        "be3_evidence": 0,
        "evidence_count": 0,
        "sources": set(),
        "best_value": None,
    }


def update_summary(summary: dict[str, object], row: dict[str, str]) -> None:
    summary["pair_count"] += 1
    tier = row.get("best_binding_evidence_level", "")
    if tier in TIER_RANK:
        summary["best_rank"] = min(summary["best_rank"], TIER_RANK[tier])
        summary["best_tier_counts"][tier] += 1
    summary["be1_evidence"] += integer(row.get("BE1_evidence_count", ""))
    summary["be2_evidence"] += integer(row.get("BE2_evidence_count", ""))
    summary["be3_evidence"] += integer(row.get("BE3_evidence_count", ""))
    summary["evidence_count"] += integer(
        row.get("binding_evidence_count", "")
    )
    summary["sources"].update(split_values(row.get("independent_sources", "")))
    value = numeric(row.get("best_standard_value_nM", ""))
    if value is not None:
        if summary["best_value"] is None:
            summary["best_value"] = value
        else:
            summary["best_value"] = min(summary["best_value"], value)


def best_tier(summary: dict[str, object]) -> str:
    rank = summary["best_rank"]
    return next(
        (tier for tier, tier_rank in TIER_RANK.items() if tier_rank == rank),
        "",
    )


def main() -> int:
    evidence_qa_path = QA_DIR / "INCREMENTAL_EVIDENCE_BUILD_V6_QA.json"
    if not evidence_qa_path.exists():
        raise FileNotFoundError(evidence_qa_path)
    evidence_qa = json.loads(evidence_qa_path.read_text(encoding="utf-8"))
    if evidence_qa.get("status") != "passed":
        raise RuntimeError("Evidence build QA did not pass")
    if not PAIR_SUMMARY.exists():
        raise FileNotFoundError(PAIR_SUMMARY)

    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    protein_baseline = pathlib.Path(config["protein_master_v5_5"])
    compound_baseline = pathlib.Path(config["compound_master_v1_0"])
    form_baseline = pathlib.Path(config["compound_forms_v1_0"])
    disease_baseline = pathlib.Path(config["protein_gene_disease_v5_4"])

    protein_summary: dict[str, dict[str, object]] = defaultdict(new_summary)
    compound_summary: dict[str, dict[str, object]] = defaultdict(new_summary)
    pair_rows = 0
    with gzip.open(
        PAIR_SUMMARY, "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pair_rows += 1
            update_summary(protein_summary[row["target_uniprot_id"]], row)
            update_summary(compound_summary[row["compound_internal_id"]], row)

    protein_output = OUT / "human_membrane_protein_master_v6_0.tsv"
    protein_fields_appended = [
        "canonical_core_small_molecule_count_v60",
        "canonical_BE1_compound_count_v60",
        "canonical_BE2_compound_count_v60",
        "canonical_BE3_compound_count_v60",
        "canonical_best_binding_evidence_level_v60",
        "canonical_binding_evidence_count_v60",
        "canonical_binding_source_count_v60",
        "canonical_binding_sources_v60",
        "canonical_best_standard_value_nM_v60",
        "compound_master_version_v60",
        "incremental_evidence_run_v60",
        "v60_release_date",
    ]
    protein_rows = 0
    with protein_baseline.open(
        "r", encoding="utf-8-sig", newline=""
    ) as source, protein_output.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or []) + protein_fields_appended
        writer = csv.DictWriter(
            destination,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            summary = protein_summary.get(
                row.get("target_uniprot_id", ""), new_summary()
            )
            tier_counts = summary["best_tier_counts"]
            row.update(
                {
                    "canonical_core_small_molecule_count_v60": summary[
                        "pair_count"
                    ],
                    "canonical_BE1_compound_count_v60": tier_counts["BE1"],
                    "canonical_BE2_compound_count_v60": tier_counts["BE2"],
                    "canonical_BE3_compound_count_v60": tier_counts["BE3"],
                    "canonical_best_binding_evidence_level_v60": best_tier(
                        summary
                    ),
                    "canonical_binding_evidence_count_v60": summary[
                        "evidence_count"
                    ],
                    "canonical_binding_source_count_v60": len(
                        summary["sources"]
                    ),
                    "canonical_binding_sources_v60": ";".join(
                        sorted(summary["sources"])
                    ),
                    "canonical_best_standard_value_nM_v60": (
                        summary["best_value"]
                        if summary["best_value"] is not None
                        else ""
                    ),
                    "compound_master_version_v60": "v1.1",
                    "incremental_evidence_run_v60": "incremental_20260727",
                    "v60_release_date": RELEASE_DATE,
                }
            )
            writer.writerow(row)
            protein_rows += 1

    compound_output = OUT / "small_molecule_master_v1_1.tsv"
    compound_rows = 0
    with compound_baseline.open(
        "r", encoding="utf-8-sig", newline=""
    ) as source, compound_output.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if "incremental_evidence_run_v11" not in fields:
            fields.append("incremental_evidence_run_v11")
        writer = csv.DictWriter(
            destination,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            summary = compound_summary.get(
                row["compound_internal_id"], new_summary()
            )
            row.update(
                {
                    "protein_target_count": summary["pair_count"],
                    "best_binding_evidence_level": best_tier(summary),
                    "BE1_evidence_count": summary["be1_evidence"],
                    "BE2_evidence_count": summary["be2_evidence"],
                    "BE3_evidence_count": summary["be3_evidence"],
                    "binding_evidence_count": summary["evidence_count"],
                    "binding_source_count": len(summary["sources"]),
                    "binding_sources": ";".join(sorted(summary["sources"])),
                    "best_standard_value_nM": (
                        summary["best_value"]
                        if summary["best_value"] is not None
                        else ""
                    ),
                    "release_version": "small-molecule-v1.1",
                    "release_date": RELEASE_DATE,
                    "incremental_evidence_run_v11": "incremental_20260727",
                }
            )
            writer.writerow(row)
            compound_rows += 1

    form_output = OUT / "compound_form_hierarchy_v1_1.tsv"
    form_rows = 0
    with form_baseline.open(
        "r", encoding="utf-8-sig", newline=""
    ) as source, form_output.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or []) + ["release_version_v11"]
        writer = csv.DictWriter(
            destination,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            row["release_version_v11"] = "small-molecule-v1.1"
            writer.writerow(row)
            form_rows += 1

    disease_output = OUT / "protein_gene_disease_relations_v6_0.tsv"
    disease_rows = 0
    with disease_baseline.open(
        "r", encoding="utf-8-sig", newline=""
    ) as source, disease_output.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or []) + [
            "release_version_v60",
            "incremental_change_v60",
        ]
        writer = csv.DictWriter(
            destination,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            row.update(
                {
                    "release_version_v60": RELEASE_VERSION,
                    "incremental_change_v60": (
                        "relation_content_unchanged_from_v5_4"
                    ),
                }
            )
            writer.writerow(row)
            disease_rows += 1

    report = {
        "status": "passed",
        "run_type": "incremental_supplement_not_rebuild",
        "release_version": RELEASE_VERSION,
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "counts": {
            "protein_master_rows": protein_rows,
            "protein_compound_pair_rows": pair_rows,
            "compound_master_rows": compound_rows,
            "compound_form_rows": form_rows,
            "protein_gene_disease_rows": disease_rows,
            "proteins_with_compounds": sum(
                1
                for summary in protein_summary.values()
                if summary["pair_count"] > 0
            ),
            "compounds_with_proteins": sum(
                1
                for summary in compound_summary.values()
                if summary["pair_count"] > 0
            ),
        },
        "outputs": {
            "protein_master": str(protein_output),
            "protein_gene_disease": str(disease_output),
            "small_molecule_master": str(compound_output),
            "compound_form_hierarchy": str(form_output),
            "protein_compound_summary": str(PAIR_SUMMARY),
        },
    }
    qa_path = QA_DIR / "INCREMENTAL_MASTER_BUILD_V6_QA.json"
    temporary = qa_path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(qa_path)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
