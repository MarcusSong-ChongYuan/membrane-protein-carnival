#!/usr/bin/env python3
"""Build an early pair-level delta from completed GPCRdb/PDBe/PDBbind staging."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import pathlib
import sqlite3
import sys
from collections import Counter


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "runs" / "incremental_20260727"
OUT_DIR = RUN_ROOT / "merge_candidate"
INDEX_DB = ROOT / "intermediate" / "integration_index_v2.sqlite"
GPCR_PAIRS = ROOT / "intermediate" / "gpcrdb_default_pair_summary_v2.tsv"
PDBE = RUN_ROOT / "staging" / "pdbe" / "pdbe_normalized_evidence_v3.tsv.gz"
PDBBIND = (
    RUN_ROOT / "staging" / "pdbbind" / "pdbbind_normalized_evidence_v3.tsv.gz"
)

TIER_RANK = {"BE1": 1, "BE2": 2, "BE3": 3, "MECH": 4, "nondefault": 9}


def split_values(value: str) -> list[str]:
    import re

    return [
        part.strip()
        for part in re.split(r"[;,|]", value or "")
        if part.strip()
    ]


def existing_pairs() -> dict[tuple[str, str], dict[str, object]]:
    con = sqlite3.connect(INDEX_DB)
    result = {}
    for target, compound, tier, count, sources in con.execute(
        """
        SELECT target_uniprot_id, compound_internal_id, best_evidence_tier,
               evidence_count, independent_sources
        FROM existing_pair
        """
    ):
        result[(target, compound)] = {
            "best_evidence_tier": tier or "",
            "evidence_count": count or 0,
            "independent_sources": sources or "",
        }
    con.close()
    return result


def new_aggregate():
    return {
        "approved_symbol": "",
        "compound_names": set(),
        "sources": set(),
        "evidence_tiers": set(),
        "activity_types": set(),
        "evidence_ids": set(),
        "pdb_by_source": {},
        "standard_values_nM": [],
    }


def add_evidence(
    pairs: dict[tuple[str, str], dict[str, object]],
    target: str,
    compound: str,
    symbol: str,
    name: str,
    source: str,
    tier: str,
    activity_types: list[str],
    evidence_ids: list[str],
    pdb_ids: list[str],
    standard_values: list[float],
) -> None:
    key = (target, compound)
    item = pairs.setdefault(key, new_aggregate())
    item["approved_symbol"] = item["approved_symbol"] or symbol
    if name:
        item["compound_names"].add(name)
    item["sources"].add(source)
    if tier:
        item["evidence_tiers"].add(tier)
    item["activity_types"].update(activity_types)
    item["evidence_ids"].update(evidence_ids)
    item["pdb_by_source"].setdefault(source, set()).update(pdb_ids)
    item["standard_values_nM"].extend(standard_values)


def numeric(value: str) -> list[float]:
    try:
        return [float(value)] if str(value).strip() else []
    except ValueError:
        return []


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    baseline = existing_pairs()
    pairs: dict[tuple[str, str], dict[str, object]] = {}
    counts: Counter[str] = Counter()

    with GPCR_PAIRS.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            add_evidence(
                pairs,
                row["target_uniprot_id"],
                row["compound_internal_id"],
                row.get("approved_symbol", ""),
                (row.get("compound_name_examples") or "").split(";")[0],
                "GPCRdb_curated_sources",
                row.get("best_evidence_tier", ""),
                split_values(row.get("activity_types", "")),
                [f"GPCRPAIR:{row['target_uniprot_id']}:{row['compound_internal_id']}"],
                [],
                numeric(row.get("best_standard_value_nM", "")),
            )
            counts["gpcrdb_pair_rows"] += 1

    for path, source_label in [(PDBE, "PDBe"), (PDBBIND, "PDBbind")]:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("default_release_inclusion") != "1":
                    continue
                target = row.get("target_uniprot_id", "")
                compound = row.get("compound_internal_id", "")
                if not target or not compound:
                    counts[f"{source_label}_default_rows_missing_key"] += 1
                    continue
                add_evidence(
                    pairs,
                    target,
                    compound,
                    row.get("approved_symbol", ""),
                    row.get("compound_name", ""),
                    source_label,
                    row.get("evidence_tier", ""),
                    split_values(row.get("activity_type", "")),
                    [row.get("source_evidence_id", "")],
                    split_values(row.get("pdb_ids", "")),
                    numeric(row.get("standard_value_nM", "")),
                )
                counts[f"{source_label}_default_evidence_rows"] += 1

    output_path = OUT_DIR / "pre_pubchem_pair_delta_v3.tsv.gz"
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    fields = [
        "target_uniprot_id",
        "approved_symbol",
        "compound_internal_id",
        "compound_name_examples",
        "new_source_databases",
        "new_source_count",
        "new_evidence_count",
        "new_best_evidence_tier",
        "new_activity_types",
        "new_pdb_ids",
        "pdbe_pdb_count",
        "pdbbind_pdb_count",
        "shared_pdbe_pdbbind_structure_count",
        "best_new_standard_value_nM",
        "existing_pair_overlap",
        "existing_best_evidence_tier",
        "existing_evidence_count",
        "existing_sources",
        "pair_action",
    ]
    structural_overlap_path = OUT_DIR / "pdbe_pdbbind_structure_overlap_v3.tsv.gz"
    structural_tmp = structural_overlap_path.with_suffix(
        structural_overlap_path.suffix + ".tmp"
    )
    with gzip.open(
        temporary, "wt", encoding="utf-8", newline="", compresslevel=6
    ) as output_handle, gzip.open(
        structural_tmp, "wt", encoding="utf-8", newline="", compresslevel=6
    ) as overlap_handle:
        writer = csv.DictWriter(
            output_handle,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )
        overlap_writer = csv.DictWriter(
            overlap_handle,
            delimiter="\t",
            fieldnames=[
                "target_uniprot_id",
                "compound_internal_id",
                "pdb_id",
                "overlap_type",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        overlap_writer.writeheader()
        for (target, compound), item in sorted(pairs.items()):
            existing = baseline.get((target, compound))
            action = "strengthen_existing_pair" if existing else "add_new_pair"
            counts[action] += 1
            pdbe_pdbs = item["pdb_by_source"].get("PDBe", set())
            pdbbind_pdbs = item["pdb_by_source"].get("PDBbind", set())
            shared = pdbe_pdbs & pdbbind_pdbs
            all_pdbs = set().union(*item["pdb_by_source"].values())
            for pdb in sorted(shared):
                overlap_writer.writerow(
                    {
                        "target_uniprot_id": target,
                        "compound_internal_id": compound,
                        "pdb_id": pdb,
                        "overlap_type": "same_target_compound_pdb_in_pdbe_and_pdbbind",
                    }
                )
                counts["shared_pdbe_pdbbind_structures"] += 1
            best_tier = min(
                item["evidence_tiers"],
                key=lambda value: TIER_RANK.get(value, 99),
                default="",
            )
            values = item["standard_values_nM"]
            writer.writerow(
                {
                    "target_uniprot_id": target,
                    "approved_symbol": item["approved_symbol"],
                    "compound_internal_id": compound,
                    "compound_name_examples": ";".join(
                        sorted(item["compound_names"])[:5]
                    ),
                    "new_source_databases": ";".join(sorted(item["sources"])),
                    "new_source_count": len(item["sources"]),
                    "new_evidence_count": len(item["evidence_ids"]),
                    "new_best_evidence_tier": best_tier,
                    "new_activity_types": ";".join(
                        sorted(item["activity_types"])
                    ),
                    "new_pdb_ids": ";".join(sorted(all_pdbs)),
                    "pdbe_pdb_count": len(pdbe_pdbs),
                    "pdbbind_pdb_count": len(pdbbind_pdbs),
                    "shared_pdbe_pdbbind_structure_count": len(shared),
                    "best_new_standard_value_nM": min(values) if values else "",
                    "existing_pair_overlap": int(existing is not None),
                    "existing_best_evidence_tier": (
                        existing["best_evidence_tier"] if existing else ""
                    ),
                    "existing_evidence_count": (
                        existing["evidence_count"] if existing else 0
                    ),
                    "existing_sources": (
                        existing["independent_sources"] if existing else ""
                    ),
                    "pair_action": action,
                }
            )
            counts["output_pair_rows"] += 1

    temporary.replace(output_path)
    structural_tmp.replace(structural_overlap_path)
    report = {
        "status": "passed",
        "run_type": "incremental_supplement_not_rebuild",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "scope": "GPCRdb + PDBe + PDBbind only; PubChem and BRENDA pending",
        "counts": dict(sorted(counts.items())),
        "outputs": {
            "pair_delta": str(output_path),
            "pdbe_pdbbind_structure_overlap": str(structural_overlap_path),
        },
    }
    qa_path = RUN_ROOT / "qa" / "PRE_PUBCHEM_INCREMENTAL_DELTA_V3_QA.json"
    tmp = qa_path.with_suffix(".json.tmp")
    tmp.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    tmp.replace(qa_path)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
