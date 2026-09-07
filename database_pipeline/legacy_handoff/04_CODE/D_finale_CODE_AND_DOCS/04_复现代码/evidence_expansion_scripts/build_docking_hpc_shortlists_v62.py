#!/usr/bin/env python3
"""Build auditable HPC docking shortlists from the frozen MemPro V6.2 release."""

from __future__ import annotations

import csv
import gzip
import hashlib
import heapq
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


RELEASE = Path(
    r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730"
)
BASE_RUN = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\docking_selection_v1_v62_20260730"
)
RUN = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\docking_selection_v2_v62_20260730"
)
STAGING = RUN / "staging"
QA = RUN / "qa"

RECOMMENDED = BASE_RUN / "staging" / "docking_recommended_pairs_v1_v62.tsv.gz"
NEGATIVE_SUMMARY = RELEASE / "protein_compound_negative_summary_v6_2.tsv.gz"

FULL_NONCONFLICT = STAGING / "docking_full_nonconflict_v2_v62.tsv.gz"
CONFLICT_REVIEW = STAGING / "docking_conflict_review_v2_v62.tsv.gz"
STANDARD = STAGING / "docking_hpc_standard_v2_v62.tsv.gz"
PILOT = STAGING / "docking_hpc_pilot_v2_v62.tsv"
TARGET_SUMMARY = STAGING / "docking_hpc_target_summary_v2_v62.tsv"
VALIDATION = QA / "DOCKING_HPC_SHORTLIST_V2_V62_VALIDATION.json"
DICTIONARY = RUN / "DOCKING_HPC_SHORTLIST_V2_V62_DATA_DICTIONARY.md"

STANDARD_CAP = {"R0": 3, "D1": 20, "D2": 10, "D3": 5}
PILOT_CAP = {"R0": 1, "D1": 5, "D2": 3, "D3": 2}
TIER_BASE = {"R0": 100.0, "D1": 90.0, "D2": 70.0, "D3": 50.0}


def text(value: object | None) -> str:
    return "" if value is None else str(value).strip()


def integer(value: str | None) -> int:
    try:
        return int(float(text(value)))
    except (TypeError, ValueError):
        return 0


def number(value: str | None) -> float | None:
    try:
        value = text(value)
        return float(value) if value else None
    except (TypeError, ValueError):
        return None


def truth(value: str | None) -> bool:
    return text(value).lower() in {"1", "true", "yes", "y"}


def first_pdb(row: dict[str, str]) -> str:
    fields = (
        "pair_site_pdb_ids",
        "opm_pdb_ids",
        "pdbtm_pdb_ids",
        "all_pdb_ids",
    )
    for field in fields:
        raw = text(row.get(field))
        if not raw:
            continue
        for token in raw.replace(",", ";").split(";"):
            token = text(token).upper()
            if len(token) == 4 and token.isalnum():
                return token
    return ""


def ranking_score(row: dict[str, str]) -> float:
    tier = text(row["docking_tier"])
    score = TIER_BASE[tier]
    sources = integer(row["independent_source_count"])
    be1 = integer(row["BE1_evidence_count"])
    be2 = integer(row["BE2_evidence_count"])
    potency = number(row["best_standard_value_nM"])
    mw = number(row["molecular_weight"])
    tpsa = number(row["tpsa"])
    rotatable = number(row["rotatable_bond_count"])
    charge = number(row["formal_charge"])

    score += min(sources, 5) * 3
    score += min(be1, 5) * 4
    score += min(be2, 5) * 2
    score += 5 if truth(row["is_approved_drug"]) else 0
    score += 4 if truth(row["is_clinical_candidate"]) else 0
    score += 3 if truth(row["is_chemical_probe"]) else 0
    score += 2 if truth(row["is_endogenous_ligand"]) else 0
    if potency is not None:
        if potency <= 10:
            score += 15
        elif potency <= 100:
            score += 10
        elif potency <= 1000:
            score += 5
    if mw is not None and 150 <= mw <= 550:
        score += 2
    if tpsa is not None and tpsa <= 140:
        score += 1
    if rotatable is not None and rotatable <= 12:
        score += 1
    if charge is not None and abs(charge) <= 1:
        score += 1
    return score


def load_conflicts() -> dict[tuple[str, str], int]:
    conflicts: dict[tuple[str, str], int] = {}
    with gzip.open(NEGATIVE_SUMMARY, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if truth(row["positive_negative_conflict_flag"]):
                conflicts[
                    (text(row["target_uniprot_id"]), text(row["compound_internal_id"]))
                ] = integer(row["negative_evidence_count"])
    return conflicts


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build() -> dict:
    STAGING.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    conflicts = load_conflicts()

    heaps: dict[tuple[str, str], list[tuple[float, int, dict[str, str]]]] = defaultdict(list)
    full_counts: Counter[str] = Counter()
    conflict_counts: Counter[str] = Counter()
    sequence = 0

    with gzip.open(
        RECOMMENDED, "rt", encoding="utf-8-sig", newline=""
    ) as source, gzip.open(
        FULL_NONCONFLICT, "wt", encoding="utf-8", newline=""
    ) as full_handle, gzip.open(
        CONFLICT_REVIEW, "wt", encoding="utf-8", newline=""
    ) as conflict_handle:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or []) + [
            "positive_negative_conflict_flag_v62",
            "negative_evidence_count_v62",
            "ranking_score_v2",
            "candidate_receptor_pdb_id_v2",
            "receptor_preparation_status_v2",
            "ligand_preparation_status_v2",
            "box_definition_status_v2",
        ]
        full_writer = csv.DictWriter(
            full_handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        conflict_writer = csv.DictWriter(
            conflict_handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        full_writer.writeheader()
        conflict_writer.writeheader()

        for row in reader:
            tier = text(row["docking_tier"])
            target = text(row["target_uniprot_id"])
            compound = text(row["compound_internal_id"])
            negative_count = conflicts.get((target, compound), 0)
            has_conflict = negative_count > 0
            score = ranking_score(row)
            candidate_pdb = first_pdb(row)
            row.update(
                {
                    "positive_negative_conflict_flag_v62": int(has_conflict),
                    "negative_evidence_count_v62": negative_count,
                    "ranking_score_v2": round(score, 2),
                    "candidate_receptor_pdb_id_v2": candidate_pdb,
                    "receptor_preparation_status_v2": (
                        "review_pair_specific_chain_state"
                        if tier == "R0"
                        else "select_chain_state_and_biological_assembly"
                    ),
                    "ligand_preparation_status_v2": (
                        "enumerate_protonation_tautomer_and_stereochemistry"
                    ),
                    "box_definition_status_v2": (
                        "derive_from_cocrystal_ligand"
                        if tier == "R0"
                        else "transfer_validated_site_or_define_pocket"
                    ),
                }
            )
            if has_conflict:
                conflict_writer.writerow(row)
                conflict_counts[tier] += 1
                continue

            full_writer.writerow(row)
            full_counts[tier] += 1
            cap = STANDARD_CAP[tier]
            key = (target, tier)
            sequence += 1
            item = (score, sequence, dict(row))
            if len(heaps[key]) < cap:
                heapq.heappush(heaps[key], item)
            elif score > heaps[key][0][0]:
                heapq.heapreplace(heaps[key], item)

    selected: list[dict[str, str]] = []
    for (target, tier), heap in heaps.items():
        ranked = sorted(heap, key=lambda item: (item[0], item[1]), reverse=True)
        for rank, (_, _, row) in enumerate(ranked, start=1):
            row["within_target_tier_rank_v2"] = rank
            row["hpc_selection_set_v2"] = "standard"
            selected.append(row)

    tier_order = {"R0": 0, "D1": 1, "D2": 2, "D3": 3}
    selected.sort(
        key=lambda row: (
            tier_order[text(row["docking_tier"])],
            text(row["target_uniprot_id"]),
            integer(row["within_target_tier_rank_v2"]),
        )
    )
    output_fields = list(selected[0].keys()) if selected else []
    if "within_target_tier_rank_v2" not in output_fields:
        output_fields.append("within_target_tier_rank_v2")
    if "hpc_selection_set_v2" not in output_fields:
        output_fields.append("hpc_selection_set_v2")

    standard_counts: Counter[str] = Counter()
    pilot_counts: Counter[str] = Counter()
    target_counts: dict[str, Counter[str]] = defaultdict(Counter)
    target_meta: dict[str, tuple[str, str]] = {}
    standard_keys: set[tuple[str, str]] = set()
    pilot_keys: set[tuple[str, str]] = set()
    standard_duplicate_keys = 0
    pilot_duplicate_keys = 0
    missing_candidate_pdb = 0

    with gzip.open(
        STANDARD, "wt", encoding="utf-8", newline=""
    ) as standard_handle, PILOT.open(
        "wt", encoding="utf-8", newline=""
    ) as pilot_handle:
        standard_writer = csv.DictWriter(
            standard_handle,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        pilot_writer = csv.DictWriter(
            pilot_handle,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        standard_writer.writeheader()
        pilot_writer.writeheader()
        for row in selected:
            tier = text(row["docking_tier"])
            target = text(row["target_uniprot_id"])
            pair = (target, text(row["compound_internal_id"]))
            if pair in standard_keys:
                standard_duplicate_keys += 1
            standard_keys.add(pair)
            if not text(row["candidate_receptor_pdb_id_v2"]):
                missing_candidate_pdb += 1
            standard_writer.writerow(row)
            standard_counts[tier] += 1
            target_counts[target][tier] += 1
            target_meta[target] = (
                text(row["approved_symbol"]),
                text(row["protein_name"]),
            )
            if integer(row["within_target_tier_rank_v2"]) <= PILOT_CAP[tier]:
                pilot_row = dict(row)
                pilot_row["hpc_selection_set_v2"] = "pilot"
                if pair in pilot_keys:
                    pilot_duplicate_keys += 1
                pilot_keys.add(pair)
                pilot_writer.writerow(pilot_row)
                pilot_counts[tier] += 1

    with TARGET_SUMMARY.open("wt", encoding="utf-8", newline="") as handle:
        fields = [
            "target_uniprot_id",
            "approved_symbol",
            "protein_name",
            "R0_standard_count",
            "D1_standard_count",
            "D2_standard_count",
            "D3_standard_count",
            "standard_pair_count",
            "pilot_pair_count",
        ]
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for target in sorted(target_counts):
            counts = target_counts[target]
            standard_total = sum(counts[tier] for tier in TIER_BASE)
            pilot_total = sum(
                min(counts[tier], PILOT_CAP[tier]) for tier in TIER_BASE
            )
            symbol, name = target_meta[target]
            writer.writerow(
                {
                    "target_uniprot_id": target,
                    "approved_symbol": symbol,
                    "protein_name": name,
                    "R0_standard_count": counts["R0"],
                    "D1_standard_count": counts["D1"],
                    "D2_standard_count": counts["D2"],
                    "D3_standard_count": counts["D3"],
                    "standard_pair_count": standard_total,
                    "pilot_pair_count": pilot_total,
                }
            )

    standard_cap_violations = sum(
        1
        for counts in target_counts.values()
        for tier, cap in STANDARD_CAP.items()
        if counts[tier] > cap
    )
    blocking_errors = {
        key: value
        for key, value in {
            "standard_duplicate_pair_keys": standard_duplicate_keys,
            "pilot_duplicate_pair_keys": pilot_duplicate_keys,
            "pilot_pairs_not_in_standard": len(pilot_keys - standard_keys),
            "standard_cap_violations": standard_cap_violations,
            "missing_candidate_receptor_pdb": missing_candidate_pdb,
        }.items()
        if value
    }
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if not blocking_errors else "FAIL",
        "release_basis": "MemPro V6.2 frozen release",
        "selection_policy": "docking_hpc_shortlist_v2",
        "full_recommended_input_rows": sum(full_counts.values())
        + sum(conflict_counts.values()),
        "full_nonconflict_rows": sum(full_counts.values()),
        "conflict_review_rows": sum(conflict_counts.values()),
        "full_nonconflict_tier_counts": dict(full_counts),
        "conflict_tier_counts": dict(conflict_counts),
        "standard_caps_per_target": STANDARD_CAP,
        "pilot_caps_per_target": PILOT_CAP,
        "standard_counts": dict(standard_counts),
        "standard_total": sum(standard_counts.values()),
        "pilot_counts": dict(pilot_counts),
        "pilot_total": sum(pilot_counts.values()),
        "standard_target_count": len(target_counts),
        "standard_unique_pair_count": len(standard_keys),
        "pilot_unique_pair_count": len(pilot_keys),
        "blocking_errors": blocking_errors,
        "outputs": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in (
                FULL_NONCONFLICT,
                CONFLICT_REVIEW,
                STANDARD,
                PILOT,
                TARGET_SUMMARY,
            )
        },
    }
    VALIDATION.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    DICTIONARY.write_text(
        f"""# Docking HPC shortlist V2 (MemPro V6.2)

The source release is immutable. Positive-negative conflicts are excluded from
default HPC lists and retained in a separate review table.

## Workload choices

| Set | Rule | Rows |
|---|---|---:|
| Pilot | Per target: R0 1, D1 5, D2 3, D3 2 | {sum(pilot_counts.values()):,} |
| Standard | Per target: R0 3, D1 20, D2 10, D3 5 | {sum(standard_counts.values()):,} |
| Full non-conflict | All R0-D3 pairs without V6.2 positive-negative conflict | {sum(full_counts.values()):,} |

R0 is a redocking/protocol-validation set. D1-D3 are prospective or mechanistic
sets. A candidate receptor PDB is provided for triage only; receptor state,
chain, biological assembly, microstates, and the docking box still require
preparation and validation.

## Ranking

Ranking is deterministic within target and tier. It favors independent sources,
BE1/BE2 support, stronger quantitative potency, approved/clinical/probe/
endogenous status, and tractable physicochemical properties. Ranking is not a
binding-affinity prediction.
""",
        encoding="utf-8",
    )
    return report


if __name__ == "__main__":
    print(json.dumps(build(), ensure_ascii=False, indent=2))
