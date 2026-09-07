#!/usr/bin/env python3
"""Rebase and harden the subunit/biological-assembly module on MemPro V6.1."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_1_20260729"
MASTER = RELEASE / "human_membrane_protein_master_v6_1.tsv"
OLD = (
    ROOT / "runs" / "subunit_assembly_v1_20260728" / "staging"
)
RUN = ROOT / "runs" / "v62_completion_20260729"
STAGING = RUN / "staging" / "subunit_assembly_v2"
QA = RUN / "qa"
PROGRESS = QA / "SUBUNIT_ASSEMBLY_V2_PROGRESS.json"
REPORT = QA / "SUBUNIT_ASSEMBLY_V2_VALIDATION.json"

UNIPROT_V1 = OLD / "uniprot_subunit_annotations_v1.tsv.gz"
ASSEMBLY_V1 = OLD / "membrane_protein_biological_assembly_v1.tsv.gz"
UNIPROT_OUT = STAGING / "uniprot_subunit_annotations_v2.tsv.gz"
ASSEMBLY_OUT = STAGING / "protein_biological_assembly_v2.tsv.gz"
COMPONENT_OUT = STAGING / "biological_assembly_components_v2.tsv.gz"
SUMMARY_OUT = STAGING / "protein_subunit_summary_v2.tsv.gz"
REVIEW_OUT = STAGING / "subunit_conflict_review_v2.tsv.gz"
PREVIEW_OUT = (
    STAGING / "human_membrane_protein_master_v6_2_subunit_preview.tsv.gz"
)
MODULE_VERSION = "subunit_assembly_v2_20260729"

NUMBER_WORDS = {
    "dimer": 2,
    "trimer": 3,
    "tetramer": 4,
    "pentamer": 5,
    "hexamer": 6,
    "heptamer": 7,
    "octamer": 8,
    "nonamer": 9,
    "decamer": 10,
    "dodecamer": 12,
}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    os.replace(tmp, path)


def split_values(value: str) -> list[str]:
    return [
        item.strip()
        for item in re.split(r"[;|]+", value or "")
        if item.strip()
    ]


def parse_comment(comment: str) -> dict:
    """Extract only explicit stoichiometric statements from a UniProt comment."""
    text = (comment or "").casefold()
    states: set[str] = set()
    counts: set[int] = set()
    reasons: set[str] = set()

    if re.search(r"\bmonomer(?:ic)?\b", text):
        states.add("monomer")
        counts.add(1)
    if re.search(r"\bhomo[- ]?and/or[- ]?heterodimer(?:ic)?\b", text):
        states.update({"homodimer", "heterodimer"})
        counts.add(2)
    for stem, count in NUMBER_WORDS.items():
        if re.search(rf"\bhomo[- ]?{stem}(?:ic)?\b", text):
            states.add("homo" + stem)
            counts.add(count)
        if re.search(rf"\bhetero[- ]?{stem}(?:ic)?\b", text):
            states.add("hetero" + stem)
            counts.add(count)
        if re.search(
            rf"(?<!homo)(?<!hetero)(?<!homo-)(?<!hetero-)\b{stem}(?:ic)?\b",
            text,
        ):
            states.add(stem + "_unspecified")
            counts.add(count)

    if re.search(
        r"two identical heavy chains and two identical light chains", text
    ):
        states.add("heterotetramer_complex")
        counts.add(4)
        reasons.add("immunoglobulin_two_heavy_two_light")
    if re.search(r"\b(?:homo)?oligomer(?:ic)?\b", text):
        if not any(state.startswith("homo") for state in states):
            states.add("homooligomer_count_unknown")
    if re.search(r"\bheterooligomer(?:ic)?\b", text):
        states.add("heterooligomer_count_unknown")
    if (
        re.search(r"\b(?:heteromeric|multimeric) complex\b", text)
        and not states
    ):
        states.add("heteromeric_complex_count_unknown")
    if "interacts with" in text and not states:
        reasons.add("interaction_only_no_stoichiometry")
    if comment and not states:
        reasons.add("annotation_without_explicit_stoichiometry")

    return {
        "states": sorted(states),
        "counts": sorted(counts),
        "parser_notes": sorted(reasons),
    }


def load_master():
    rows = []
    by_target = {}
    pdb_by_target = {}
    with MASTER.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        for row in reader:
            target = row["target_uniprot_id"]
            rows.append(row)
            by_target[target] = row
            pdb_by_target[target] = {
                item.strip().lower()
                for item in re.split(r"[;,| ]+", row.get("pdb_ids", ""))
                if item.strip()
            }
    return rows, fields, by_target, pdb_by_target


def main() -> None:
    STAGING.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    started = now()
    write_json(
        PROGRESS,
        {
            "status": "running",
            "stage": "load_v61_and_cached_annotations",
            "started_utc": started,
            "updated_utc": now(),
        },
    )
    master_rows, master_fields, by_target, pdb_by_target = load_master()

    uniprot = {}
    with gzip.open(
        UNIPROT_V1, "rt", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            uniprot[row["target_uniprot_id"]] = row

    parsed_uniprot = {}
    uniprot_fields = [
        "target_uniprot_id",
        "fetch_status",
        "subunit_states_v2",
        "subunit_count_values_v2",
        "direct_experimental_flag",
        "evidence_grade",
        "pubmed_ids",
        "eco_codes",
        "parser_notes_v2",
        "subunit_comment",
        "retrieval_utc",
        "module_version",
    ]
    with gzip.open(
        UNIPROT_OUT, "wt", encoding="utf-8", newline=""
    ) as output:
        writer = csv.DictWriter(
            output,
            fieldnames=uniprot_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for target in sorted(by_target):
            old = uniprot.get(
                target,
                {
                    "fetch_status": "missing_cache",
                    "direct_experimental_flag": "0",
                    "pubmed_ids": "",
                    "eco_codes": "",
                    "subunit_comment": "",
                    "retrieval_utc": "",
                },
            )
            parsed = parse_comment(old.get("subunit_comment", ""))
            direct = old.get("direct_experimental_flag") == "1"
            if direct and parsed["states"]:
                grade = "SA1_direct_experimental"
            elif parsed["states"]:
                grade = "SA2_curated_annotation"
            else:
                grade = "SA0_unresolved"
            row = {
                "target_uniprot_id": target,
                "fetch_status": old.get("fetch_status", ""),
                "subunit_states_v2": ";".join(parsed["states"]),
                "subunit_count_values_v2": ";".join(
                    map(str, parsed["counts"])
                ),
                "direct_experimental_flag": int(direct),
                "evidence_grade": grade,
                "pubmed_ids": old.get("pubmed_ids", ""),
                "eco_codes": old.get("eco_codes", ""),
                "parser_notes_v2": ";".join(parsed["parser_notes"]),
                "subunit_comment": old.get("subunit_comment", ""),
                "retrieval_utc": old.get("retrieval_utc", ""),
                "module_version": MODULE_VERSION,
            }
            parsed_uniprot[target] = row
            writer.writerow(row)

    assembly_fields = [
        "target_uniprot_id",
        "pdb_id",
        "biological_assembly_id",
        "assembly_details",
        "assembly_assignment_class",
        "structure_experimental_flag",
        "assembly_direct_experimental_confirmation_flag",
        "assembly_complex_state",
        "assembly_total_protein_subunit_count",
        "distinct_protein_entity_count",
        "target_copy_count",
        "target_state_from_assembly",
        "target_state_inference_status",
        "subunit_composition",
        "assembly_composition",
        "polymeric_count",
        "molecular_weight_kda",
        "source_database",
        "source_endpoint",
        "module_version",
    ]
    component_fields = [
        "assembly_component_id",
        "target_uniprot_id",
        "pdb_id",
        "biological_assembly_id",
        "component_index",
        "component_name",
        "copy_count",
        "target_component_status",
        "module_version",
    ]
    assemblies_by_target: dict[str, list[dict]] = defaultdict(list)
    assembly_count = 0
    component_count = 0
    stale_pdb_rows = 0
    with gzip.open(
        ASSEMBLY_V1, "rt", encoding="utf-8", newline=""
    ) as input_handle, gzip.open(
        ASSEMBLY_OUT, "wt", encoding="utf-8", newline=""
    ) as assembly_handle, gzip.open(
        COMPONENT_OUT, "wt", encoding="utf-8", newline=""
    ) as component_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        assembly_writer = csv.DictWriter(
            assembly_handle,
            fieldnames=assembly_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        component_writer = csv.DictWriter(
            component_handle,
            fieldnames=component_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        assembly_writer.writeheader()
        component_writer.writeheader()
        for old in reader:
            target = old["target_uniprot_id"]
            if target not in by_target:
                stale_pdb_rows += 1
                continue
            pdb_id = old["pdb_id"].lower()
            if pdb_id not in pdb_by_target[target]:
                stale_pdb_rows += 1
                continue
            distinct = int(old.get("distinct_protein_entity_count") or 0)
            total = int(old.get("protein_subunit_count") or 0)
            assignment = old.get("assembly_assignment_class", "unspecified")
            if distinct == 1 and total > 0:
                target_copy = total
                target_state = old.get("oligomeric_state", "")
                inference = (
                    "author_defined_single_protein_entity"
                    if assignment == "author_defined"
                    else "software_or_unspecified_single_protein_entity"
                )
            elif distinct > 1:
                target_copy = ""
                target_state = "heteromeric_complex_target_copy_unknown"
                inference = "component_to_uniprot_not_resolved"
            else:
                target_copy = ""
                target_state = ""
                inference = "no_polypeptide_entity"
            row = {
                "target_uniprot_id": target,
                "pdb_id": pdb_id,
                "biological_assembly_id": old.get(
                    "biological_assembly_id", ""
                ),
                "assembly_details": old.get("assembly_details", ""),
                "assembly_assignment_class": assignment,
                "structure_experimental_flag": old.get(
                    "structure_experimental_flag", "1"
                ),
                "assembly_direct_experimental_confirmation_flag": old.get(
                    "assembly_direct_experimental_confirmation_flag", "0"
                ),
                "assembly_complex_state": old.get("oligomeric_state", ""),
                "assembly_total_protein_subunit_count": total or "",
                "distinct_protein_entity_count": distinct,
                "target_copy_count": target_copy,
                "target_state_from_assembly": target_state,
                "target_state_inference_status": inference,
                "subunit_composition": old.get("subunit_composition", ""),
                "assembly_composition": old.get("assembly_composition", ""),
                "polymeric_count": old.get("polymeric_count", ""),
                "molecular_weight_kda": old.get("molecular_weight_kda", ""),
                "source_database": old.get("source_database", "PDBe"),
                "source_endpoint": old.get("source_endpoint", ""),
                "module_version": MODULE_VERSION,
            }
            assembly_writer.writerow(row)
            assemblies_by_target[target].append(row)
            assembly_count += 1

            components = [
                item.strip()
                for item in old.get("subunit_composition", "").split(";")
                if item.strip()
            ]
            for index, component in enumerate(components, start=1):
                match = re.match(r"(.+?)\s+x(\d+)\s*$", component)
                name = match.group(1).strip() if match else component
                copies = match.group(2) if match else ""
                component_count += 1
                component_writer.writerow(
                    {
                        "assembly_component_id": (
                            f"ASMCOMP-{hashlib.sha1((target+'|'+pdb_id+'|'+old.get('biological_assembly_id','')+'|'+str(index)).encode()).hexdigest()[:20].upper()}"
                        ),
                        "target_uniprot_id": target,
                        "pdb_id": pdb_id,
                        "biological_assembly_id": old.get(
                            "biological_assembly_id", ""
                        ),
                        "component_index": index,
                        "component_name": name,
                        "copy_count": copies,
                        "target_component_status": (
                            "target_only_component"
                            if distinct == 1
                            else "unresolved_component_to_uniprot"
                        ),
                        "module_version": MODULE_VERSION,
                    }
                )

    summary_fields = [
        "membrane_protein_id",
        "target_uniprot_id",
        "approved_symbol",
        "protein_name",
        "oligomeric_state_consensus_v62",
        "oligomeric_state_values_v62",
        "subunit_count_values_v62",
        "subunit_count_min_v62",
        "subunit_count_max_v62",
        "subunit_composition_summary_v62",
        "biological_assembly_ids_v62",
        "pdbe_pdb_count_v62",
        "pdbe_assembly_count_v62",
        "experimental_subunit_evidence_flag_v62",
        "subunit_evidence_grade_v62",
        "subunit_conflict_flag_v62",
        "subunit_conflict_reason_v62",
        "large_assembly_review_flag_v62",
        "uniprot_subunit_comment_v62",
        "uniprot_pubmed_ids_v62",
        "subunit_module_version_v62",
    ]
    summary_rows = {}
    review_rows = []
    state_distribution = Counter()
    grade_distribution = Counter()
    for target in sorted(by_target):
        protein = by_target[target]
        uni = parsed_uniprot[target]
        uni_states = set(split_values(uni["subunit_states_v2"]))
        uni_counts = {
            int(value)
            for value in split_values(uni["subunit_count_values_v2"])
            if value.isdigit()
        }
        assemblies = assemblies_by_target.get(target, [])
        author_target_states = {
            row["target_state_from_assembly"]
            for row in assemblies
            if row["assembly_assignment_class"] == "author_defined"
            and row["target_state_inference_status"]
            == "author_defined_single_protein_entity"
            and row["target_state_from_assembly"]
        }
        other_target_states = {
            row["target_state_from_assembly"]
            for row in assemblies
            if row["target_state_inference_status"]
            == "software_or_unspecified_single_protein_entity"
            and row["target_state_from_assembly"]
        }
        assembly_counts = {
            int(row["target_copy_count"])
            for row in assemblies
            if str(row["target_copy_count"]).isdigit()
        }
        all_states = uni_states | author_target_states | other_target_states
        all_counts = uni_counts | assembly_counts
        if uni_states:
            consensus_states = uni_states
            grade = uni["evidence_grade"]
        elif author_target_states:
            consensus_states = author_target_states
            grade = "SA2_author_defined_structure"
        elif other_target_states:
            consensus_states = other_target_states
            grade = "SA3_software_inferred_assembly"
        else:
            consensus_states = set()
            grade = "SA0_unresolved"
        consensus = (
            next(iter(consensus_states))
            if len(consensus_states) == 1
            else "multiple_or_state_dependent"
            if len(consensus_states) > 1
            else "unknown"
        )

        conflict_reasons = []
        if uni_states and author_target_states and uni_states.isdisjoint(
            author_target_states
        ):
            conflict_reasons.append("UniProt_vs_author_defined_PDBe")
        if len(author_target_states | other_target_states) > 1:
            conflict_reasons.append("multiple_PDBe_target_states")
        if uni_states and len(uni_states) > 1:
            conflict_reasons.append("multiple_UniProt_states_or_contexts")
        large = any(
            int(row["assembly_total_protein_subunit_count"] or 0) >= 24
            for row in assemblies
        )
        composition_values = sorted(
            {
                row["subunit_composition"]
                for row in assemblies
                if row["subunit_composition"]
            }
        )
        assembly_ids = sorted(
            {
                f"{row['pdb_id']}:{row['biological_assembly_id']}"
                for row in assemblies
            }
        )
        pdb_ids = {row["pdb_id"] for row in assemblies}
        summary = {
            "membrane_protein_id": protein["membrane_protein_id"],
            "target_uniprot_id": target,
            "approved_symbol": protein.get("approved_symbol", ""),
            "protein_name": protein.get("protein_name", ""),
            "oligomeric_state_consensus_v62": consensus,
            "oligomeric_state_values_v62": ";".join(sorted(all_states)),
            "subunit_count_values_v62": ";".join(
                map(str, sorted(all_counts))
            ),
            "subunit_count_min_v62": min(all_counts) if all_counts else "",
            "subunit_count_max_v62": max(all_counts) if all_counts else "",
            "subunit_composition_summary_v62": " || ".join(
                composition_values[:10]
            ),
            "biological_assembly_ids_v62": ";".join(assembly_ids),
            "pdbe_pdb_count_v62": len(pdb_ids),
            "pdbe_assembly_count_v62": len(assemblies),
            "experimental_subunit_evidence_flag_v62": uni[
                "direct_experimental_flag"
            ],
            "subunit_evidence_grade_v62": grade,
            "subunit_conflict_flag_v62": int(bool(conflict_reasons)),
            "subunit_conflict_reason_v62": ";".join(conflict_reasons),
            "large_assembly_review_flag_v62": int(large),
            "uniprot_subunit_comment_v62": uni["subunit_comment"],
            "uniprot_pubmed_ids_v62": uni["pubmed_ids"],
            "subunit_module_version_v62": MODULE_VERSION,
        }
        summary_rows[target] = summary
        state_distribution[consensus] += 1
        grade_distribution[grade] += 1
        review_reasons = list(conflict_reasons)
        if large:
            review_reasons.append("large_assembly_total_ge_24")
        if "immunoglobulin" in uni["parser_notes_v2"]:
            review_reasons.append("immunoglobulin_complex_context")
        if "immunoglobulin_two_heavy_two_light" in uni["parser_notes_v2"]:
            review_reasons.append("immunoglobulin_complex_context")
        if review_reasons:
            review_rows.append(
                {
                    "review_id": f"SUBUNIT-REV-{len(review_rows)+1:07d}",
                    "target_uniprot_id": target,
                    "approved_symbol": protein.get("approved_symbol", ""),
                    "review_reasons": ";".join(sorted(set(review_reasons))),
                    "oligomeric_state_consensus_v62": consensus,
                    "oligomeric_state_values_v62": ";".join(
                        sorted(all_states)
                    ),
                    "uniprot_subunit_comment": uni["subunit_comment"],
                    "pdbe_pdb_count": len(pdb_ids),
                    "pdbe_assembly_count": len(assemblies),
                    "review_status": "pending_prioritized_review",
                    "module_version": MODULE_VERSION,
                }
            )

    with gzip.open(
        SUMMARY_OUT, "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=summary_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows[target] for target in sorted(summary_rows))

    review_fields = [
        "review_id",
        "target_uniprot_id",
        "approved_symbol",
        "review_reasons",
        "oligomeric_state_consensus_v62",
        "oligomeric_state_values_v62",
        "uniprot_subunit_comment",
        "pdbe_pdb_count",
        "pdbe_assembly_count",
        "review_status",
        "module_version",
    ]
    with gzip.open(
        REVIEW_OUT, "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=review_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(review_rows)

    master_add_fields = summary_fields[4:]
    with MASTER.open(
        "rt", encoding="utf-8-sig", newline=""
    ) as input_handle, gzip.open(
        PREVIEW_OUT, "wt", encoding="utf-8", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        writer = csv.DictWriter(
            output_handle,
            fieldnames=master_fields + master_add_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            target = row["target_uniprot_id"]
            row.update(
                {
                    field: summary_rows[target][field]
                    for field in master_add_fields
                }
            )
            writer.writerow(row)

    expected_pdb_pairs = {
        (target, pdb)
        for target, pdbs in pdb_by_target.items()
        for pdb in pdbs
    }
    observed_pdb_pairs = {
        (target, row["pdb_id"])
        for target, rows in assemblies_by_target.items()
        for row in rows
    }
    report = {
        "status": "passed_with_review_queue",
        "started_utc": started,
        "completed_utc": now(),
        "module_version": MODULE_VERSION,
        "inputs": {
            "protein_master": str(MASTER),
            "protein_rows": len(master_rows),
            "protein_master_sha256": sha256(MASTER),
            "uniprot_cache": str(UNIPROT_V1),
            "uniprot_cache_sha256": sha256(UNIPROT_V1),
            "assembly_cache_export": str(ASSEMBLY_V1),
            "assembly_cache_export_sha256": sha256(ASSEMBLY_V1),
        },
        "outputs": {
            "uniprot_rows": len(parsed_uniprot),
            "assembly_rows": assembly_count,
            "component_rows": component_count,
            "summary_rows": len(summary_rows),
            "review_rows": len(review_rows),
        },
        "state_distribution": dict(state_distribution),
        "evidence_grade_distribution": dict(grade_distribution),
        "cache_reconciliation": {
            "expected_v61_protein_pdb_pairs": len(expected_pdb_pairs),
            "observed_pairs_with_assembly": len(observed_pdb_pairs),
            "pairs_without_assembly_export": len(
                expected_pdb_pairs - observed_pdb_pairs
            ),
            "stale_v1_assembly_rows_excluded": stale_pdb_rows,
        },
        "quality_checks": {
            "summary_primary_key_unique": len(summary_rows)
            == len(master_rows),
            "summary_foreign_keys_in_master": set(summary_rows)
            == set(by_target),
            "pdb_asymmetric_unit_not_used_as_native_state": True,
            "heteromer_target_copy_not_inferred_without_component_mapping": True,
            "pdb_experimental_structure_not_equal_direct_stoichiometry": True,
            "large_complexes_review_flagged": True,
            "v61_immutable": True,
            "formal_release_frozen": False,
        },
        "output_files": {
            "uniprot": str(UNIPROT_OUT),
            "assembly": str(ASSEMBLY_OUT),
            "components": str(COMPONENT_OUT),
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
        },
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
