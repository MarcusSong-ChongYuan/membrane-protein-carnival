#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
NORMALIZED = ROOT / "intermediate" / "gpcrdb_normalized_evidence_v2.tsv"
DEFAULT = ROOT / "intermediate" / "gpcrdb_default_new_evidence_v2.tsv"
PAIR_SUMMARY = ROOT / "intermediate" / "gpcrdb_default_pair_summary_v2.tsv"
REPORT = ROOT / "reports" / "GPCRDB_NORMALIZATION_V2_QA.json"


def main() -> None:
    errors = Counter()
    normalized_rows = 0
    evidence_ids = set()
    mapping_counts = Counter()
    unresolved_by_source = Counter()
    unresolved_by_ligand_type = Counter()

    with NORMALIZED.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            normalized_rows += 1
            evidence_id = row["source_evidence_id"]
            if evidence_id in evidence_ids:
                errors["duplicate_source_evidence_id"] += 1
            evidence_ids.add(evidence_id)
            mapping_counts[row["compound_mapping_status"]] += 1
            if row["compound_mapping_status"] == "unresolved":
                unresolved_by_source[row["source_database"]] += 1
                unresolved_by_ligand_type[row["source_ligand_type"]] += 1

    pair_rows: dict[tuple[str, str], dict] = {}
    source_evidence_rows = Counter()
    source_pairs: dict[str, set[tuple[str, str]]] = defaultdict(set)
    source_novel_pairs: dict[str, set[tuple[str, str]]] = defaultdict(set)
    source_strengthened_pairs: dict[str, set[tuple[str, str]]] = defaultdict(set)
    source_proteins: dict[str, set[str]] = defaultdict(set)
    source_compounds: dict[str, set[str]] = defaultdict(set)
    evidence_tiers = Counter()
    activity_types = Counter()
    default_rows = 0

    with DEFAULT.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            default_rows += 1
            source = row["source_database"]
            target = row["target_uniprot_id"]
            compound = row["compound_internal_id"]
            pair = (target, compound)
            existing_pair = row["existing_pair_overlap"] == "1"
            if row["default_release_inclusion"] != "1":
                errors["default_file_contains_nondefault_row"] += 1
            if row["compound_mapping_status"] != "mapped_core":
                errors["default_row_not_mapped_core"] += 1
            if row["evidence_tier"] not in {"BE2", "BE3"}:
                errors["default_row_bad_evidence_tier"] += 1
            if source == "ChEMBL":
                errors["default_row_chembl_proxy"] += 1
            try:
                if float(row["standard_value_nM"]) <= 0:
                    errors["default_row_nonpositive_standard_value"] += 1
            except ValueError:
                errors["default_row_missing_standard_value"] += 1

            source_evidence_rows[source] += 1
            source_pairs[source].add(pair)
            source_proteins[source].add(target)
            source_compounds[source].add(compound)
            if existing_pair:
                source_strengthened_pairs[source].add(pair)
            else:
                source_novel_pairs[source].add(pair)
            evidence_tiers[row["evidence_tier"]] += 1
            activity_types[row["activity_type"]] += 1

            summary = pair_rows.setdefault(
                pair,
                {
                    "target_uniprot_id": target,
                    "approved_symbol": row["approved_symbol"],
                    "compound_internal_id": compound,
                    "compound_name_examples": set(),
                    "source_databases": set(),
                    "evidence_tiers": set(),
                    "activity_types": set(),
                    "evidence_count": 0,
                    "best_standard_value_nM": None,
                    "existing_pair_overlap": int(existing_pair),
                },
            )
            summary["compound_name_examples"].add(row["compound_name"])
            summary["source_databases"].add(source)
            summary["evidence_tiers"].add(row["evidence_tier"])
            summary["activity_types"].add(row["activity_type"])
            summary["evidence_count"] += 1
            value = float(row["standard_value_nM"])
            current = summary["best_standard_value_nM"]
            summary["best_standard_value_nM"] = (
                value if current is None else min(current, value)
            )
            summary["existing_pair_overlap"] = max(
                summary["existing_pair_overlap"], int(existing_pair)
            )

    fields = [
        "target_uniprot_id",
        "approved_symbol",
        "compound_internal_id",
        "compound_name_examples",
        "source_databases",
        "best_evidence_tier",
        "activity_types",
        "evidence_count",
        "best_standard_value_nM",
        "existing_pair_overlap",
        "pair_action",
    ]
    with PAIR_SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for pair in sorted(pair_rows):
            summary = pair_rows[pair]
            row = dict(summary)
            row["compound_name_examples"] = ";".join(
                sorted(summary["compound_name_examples"])[:5]
            )
            row["source_databases"] = ";".join(sorted(summary["source_databases"]))
            row["best_evidence_tier"] = (
                "BE2" if "BE2" in summary["evidence_tiers"] else "BE3"
            )
            row["activity_types"] = ";".join(sorted(summary["activity_types"]))
            row["pair_action"] = (
                "strengthen_existing_pair"
                if summary["existing_pair_overlap"]
                else "add_new_pair"
            )
            writer.writerow({field: row.get(field, "") for field in fields})

    source_summary = {}
    for source in sorted(source_evidence_rows):
        source_summary[source] = {
            "evidence_rows": source_evidence_rows[source],
            "unique_pairs": len(source_pairs[source]),
            "novel_pairs": len(source_novel_pairs[source]),
            "existing_pairs_strengthened": len(source_strengthened_pairs[source]),
            "proteins": len(source_proteins[source]),
            "compounds": len(source_compounds[source]),
        }

    report = {
        "status": "passed" if not errors else "failed",
        "normalized_rows": normalized_rows,
        "unique_source_evidence_ids": len(evidence_ids),
        "default_candidate_rows": default_rows,
        "default_unique_pairs": len(pair_rows),
        "novel_protein_compound_pairs": sum(
            row["existing_pair_overlap"] == 0 for row in pair_rows.values()
        ),
        "existing_pairs_strengthened": sum(
            row["existing_pair_overlap"] == 1 for row in pair_rows.values()
        ),
        "default_unique_proteins": len(
            {row["target_uniprot_id"] for row in pair_rows.values()}
        ),
        "default_unique_compounds": len(
            {row["compound_internal_id"] for row in pair_rows.values()}
        ),
        "source_summary": source_summary,
        "evidence_tier_rows": dict(evidence_tiers),
        "activity_type_rows": dict(activity_types),
        "compound_mapping_counts_all_normalized": dict(mapping_counts),
        "unresolved_rows_by_source": dict(unresolved_by_source),
        "unresolved_rows_by_ligand_type": dict(unresolved_by_ligand_type),
        "qa_error_counts": dict(errors),
        "pair_summary_file": str(PAIR_SUMMARY),
    }
    with REPORT.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
