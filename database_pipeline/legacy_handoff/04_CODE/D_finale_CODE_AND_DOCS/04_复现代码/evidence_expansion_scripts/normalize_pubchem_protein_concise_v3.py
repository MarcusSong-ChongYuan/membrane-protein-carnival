#!/usr/bin/env python3
"""Watch and normalize exact-UniProt PubChem concise files incrementally.

Each completed protein is written atomically to its own gzip TSV.  This allows
normalization to run beside the downloader without reading partial files or
rewriting already completed targets.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import pathlib
import re
import shutil
import sqlite3
import sys
import time
from collections import Counter, defaultdict


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "runs" / "incremental_20260727"
OUT_DIR = RUN_ROOT / "staging" / "pubchem" / "by_target"
QA_DIR = RUN_ROOT / "staging" / "pubchem" / "target_qa"
CHECKPOINT = ROOT / "raw" / "pubchem_202607" / "pubchem_protein_concise_v2.jsonl"
RAW_DIR = ROOT / "raw" / "pubchem_202607" / "protein_concise"
PROTEIN_UNIVERSE = ROOT / "intermediate" / "protein_universe_abc_v2.tsv"
INDEX_DB = ROOT / "intermediate" / "integration_index_v2.sqlite"
SCHEMA_PATH = ROOT / "config" / "normalized_source_schema.tsv"
PROGRESS = RUN_ROOT / "staging" / "pubchem" / "normalization_progress_v3.jsonl"
SUMMARY = RUN_ROOT / "qa" / "PUBCHEM_NORMALIZATION_V3_PROGRESS.json"
MIN_FREE_BYTES = 25 * 1024**3

EXTRA_FIELDS = [
    "pubchem_aid",
    "pubchem_sid",
    "pubchem_cid",
    "pubchem_target_gene_id",
    "pubchem_assay_type",
    "pubchem_activity_name",
    "pubchem_activity_value_uM",
    "release_bucket",
]

DIRECT_BINDING_TYPES = {
    "KD",
    "KI",
    "KB",
    "KA",
    "KE",
}
QUANTITATIVE_ACTIVITY_TYPES = {
    "KD",
    "KI",
    "KB",
    "KA",
    "KE",
    "IC50",
    "EC50",
    "AC50",
    "DC50",
    "ED50",
    "GI50",
    "ID50",
    "KM",
    "POTENCY",
    "INHIBITION",
}


def stable_id(prefix: str, *parts: object) -> str:
    material = "\x1f".join("" if part is None else str(part) for part in parts)
    return prefix + hashlib.sha256(material.encode("utf-8")).hexdigest()[:24].upper()


def load_schema() -> list[str]:
    with SCHEMA_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        return [row["column_name"] for row in csv.DictReader(handle, delimiter="\t")]


def load_proteins() -> dict[str, str]:
    result = {}
    with PROTEIN_UNIVERSE.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["target_uniprot_id"]] = row.get("approved_symbol", "")
    return result


def load_cid_mappings():
    con = sqlite3.connect(INDEX_DB)
    candidates: dict[str, dict[tuple[str, str], dict[str, str]]] = defaultdict(dict)
    for row in con.execute(
        """
        SELECT id_value, compound_internal_id, compound_form_id, mapping_status,
               scope_status, preferred_source_name
        FROM compound_xref
        WHERE id_type='PubChem CID'
        """
    ):
        cid, internal, form, status, scope, name = row
        if not internal:
            continue
        candidates[str(cid)][(str(internal), str(form or ""))] = {
            "compound_internal_id": str(internal),
            "compound_form_id": str(form or ""),
            "compound_mapping_status": str(status or "mapped_core"),
            "compound_scope_status": str(scope or ""),
            "compound_name": str(name or ""),
        }
    con.close()
    unique = {}
    ambiguous = set()
    for cid, mappings in candidates.items():
        parents = {key[0] for key in mappings}
        if len(parents) == 1:
            unique[cid] = sorted(
                mappings.values(),
                key=lambda row: (
                    not bool(row["compound_form_id"]),
                    row["compound_form_id"],
                ),
            )[0]
        else:
            ambiguous.add(cid)
    return unique, ambiguous


def load_download_state() -> dict[str, dict[str, object]]:
    result = {}
    if not CHECKPOINT.exists():
        return result
    with CHECKPOINT.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if row.get("status") in {"ok", "no_data"}:
                result[str(row["target_uniprot_id"])] = row
    return result


def load_completed() -> set[str]:
    completed = set()
    if not PROGRESS.exists():
        return completed
    with PROGRESS.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if row.get("status") == "ok":
                completed.add(str(row["target_uniprot_id"]))
    return completed


def clean_activity_name(value: str) -> str:
    return re.sub(r"\s+", " ", (value or "").strip())


def classify_row(
    outcome: str,
    assay_type: str,
    activity_name: str,
    activity_value: str,
    mapping_status: str,
) -> tuple[str, str, str, list[str]]:
    outcome_upper = outcome.upper()
    assay_upper = assay_type.upper()
    activity_upper = activity_name.upper()
    numeric = False
    try:
        float(activity_value)
        numeric = bool(activity_value.strip())
    except (TypeError, ValueError):
        numeric = False

    if outcome_upper == "INACTIVE":
        return "negative", "negative", "pubchem_inactive", []
    if outcome_upper == "INCONCLUSIVE":
        return (
            "nondefault",
            "nondefault",
            "pubchem_inconclusive",
            ["inconclusive assay outcome"],
        )

    reasons = []
    positive_signal = outcome_upper == "ACTIVE" or (
        numeric and activity_upper in QUANTITATIVE_ACTIVITY_TYPES
    )
    if not positive_signal:
        reasons.append("neither active nor recognized quantitative result")
    if assay_upper == "SCREENING":
        reasons.append("primary screening assay without confirmatory designation")
    if mapping_status not in {"mapped_core", "mapped_extended"}:
        reasons.append("compound identity is not mapped to frozen V1.0")

    if activity_upper in DIRECT_BINDING_TYPES and numeric:
        tier = "BE2"
        evidence_type = "pubchem_direct_quantitative_binding"
    elif numeric and activity_upper in QUANTITATIVE_ACTIVITY_TYPES:
        tier = "BE3"
        evidence_type = "pubchem_target_specific_quantitative_activity"
    else:
        tier = "BE3"
        evidence_type = "pubchem_target_specific_activity"
    bucket = "default_candidate" if not reasons else "nondefault"
    return bucket, tier, evidence_type, reasons


def process_target(
    target: str,
    approved_symbol: str,
    cid_map: dict[str, dict[str, str]],
    ambiguous_cids: set[str],
    fields: list[str],
) -> dict[str, object]:
    source = RAW_DIR / f"{target}.csv.gz"
    destination = OUT_DIR / f"{target}.tsv.gz"
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    qa_path = QA_DIR / f"{target}.json"
    counts: Counter[str] = Counter()
    aids = set()
    cids = set()
    mapped_compounds = set()
    seen_evidence_ids = set()

    with gzip.open(source, "rt", encoding="utf-8-sig", newline="") as input_handle:
        reader = csv.DictReader(input_handle)
        required = {
            "AID",
            "SID",
            "CID",
            "Activity Outcome",
            "Activity Name",
            "Activity Value [uM]",
            "Assay Type",
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Unexpected PubChem header in {source}")
        with gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=6
        ) as output_handle:
            writer = csv.DictWriter(
                output_handle,
                delimiter="\t",
                fieldnames=fields,
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()
            for raw in reader:
                aid = (raw.get("AID") or "").strip()
                sid = (raw.get("SID") or "").strip()
                cid = (raw.get("CID") or "").strip()
                outcome = (raw.get("Activity Outcome") or "").strip()
                assay_type = (raw.get("Assay Type") or "").strip()
                activity_name = clean_activity_name(raw.get("Activity Name") or "")
                qualifier = (raw.get("Activity Qualifier") or "").strip()
                value = (raw.get("Activity Value [uM]") or "").strip()
                mapping = cid_map.get(cid)
                if cid in ambiguous_cids:
                    mapping_status = "ambiguous_pubchem_cid"
                elif mapping:
                    mapping_status = mapping["compound_mapping_status"]
                else:
                    mapping_status = "unresolved"
                bucket, tier, evidence_type, reasons = classify_row(
                    outcome, assay_type, activity_name, value, mapping_status
                )
                standard_value = ""
                if value:
                    try:
                        standard_value = str(float(value) * 1000.0)
                    except ValueError:
                        reasons.append("quantitative value could not be parsed")
                evidence_id = stable_id(
                    "PUBCHEM-",
                    target,
                    aid,
                    sid,
                    cid,
                    outcome,
                    activity_name,
                    qualifier,
                    value,
                    assay_type,
                    raw.get("Assay Name") or "",
                )
                if evidence_id in seen_evidence_ids:
                    counts["exact_duplicate_rows_skipped"] += 1
                    continue
                seen_evidence_ids.add(evidence_id)
                internal = mapping["compound_internal_id"] if mapping else ""
                form = mapping["compound_form_id"] if mapping else ""
                row = {
                    "source_evidence_id": evidence_id,
                    "source_database": "PubChem BioAssay",
                    "source_version": "current_2026-07",
                    "source_record_id": f"AID:{aid};SID:{sid};CID:{cid}",
                    "source_url": f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{aid}",
                    "retrieval_date": "2026-07-27",
                    "target_uniprot_id": target,
                    "target_mapping_status": "single_human_uniprot",
                    "approved_symbol": approved_symbol,
                    "compound_source_id": f"PubChem CID:{cid}",
                    "compound_name": mapping["compound_name"] if mapping else "",
                    "compound_internal_id": internal,
                    "compound_form_id": form,
                    "compound_mapping_status": mapping_status,
                    "evidence_tier": tier,
                    "evidence_type": evidence_type,
                    "activity_type": activity_name,
                    "activity_relation": qualifier,
                    "activity_value": value,
                    "activity_unit": "uM" if value else "",
                    "standard_value_nM": standard_value,
                    "activity_outcome": outcome.lower(),
                    "assay_or_mechanism": (raw.get("Assay Name") or "").strip(),
                    "pdb_ids": "",
                    "binding_site_residues": "",
                    "pubmed_ids": (raw.get("PubMed ID") or "").strip(),
                    "doi": "",
                    "record_qc_status": "ok" if not reasons else "review",
                    "default_release_inclusion": str(
                        int(bucket == "default_candidate")
                    ),
                    "exclusion_reason": "; ".join(reasons),
                    "pubchem_aid": aid,
                    "pubchem_sid": sid,
                    "pubchem_cid": cid,
                    "pubchem_target_gene_id": (
                        raw.get("Target GeneID") or ""
                    ).strip(),
                    "pubchem_assay_type": assay_type,
                    "pubchem_activity_name": activity_name,
                    "pubchem_activity_value_uM": value,
                    "release_bucket": bucket,
                }
                writer.writerow(row)
                counts["normalized_rows"] += 1
                counts[f"bucket_{bucket}"] += 1
                counts[f"tier_{tier}"] += 1
                counts[f"mapping_{mapping_status}"] += 1
                aids.add(aid)
                cids.add(cid)
                if internal:
                    mapped_compounds.add(internal)

    temporary.replace(destination)
    qa = {
        "status": "ok",
        "target_uniprot_id": target,
        "source_path": str(source),
        "output_path": str(destination),
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "counts": dict(sorted(counts.items())),
        "unique_aids": len(aids),
        "unique_cids": len(cids),
        "unique_mapped_compounds": len(mapped_compounds),
        "compressed_output_bytes": destination.stat().st_size,
    }
    qa_tmp = qa_path.with_suffix(".json.tmp")
    qa_tmp.write_text(
        json.dumps(qa, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    qa_tmp.replace(qa_path)
    return qa


def write_summary(
    download_state: dict[str, dict[str, object]],
    completed: set[str],
    run_counts: Counter[str],
) -> None:
    payload = {
        "status": (
            "passed"
            if len(download_state) == 10556
            and all(
                row.get("status") in {"ok", "no_data"}
                for row in download_state.values()
            )
            and all(
                target in completed
                for target, row in download_state.items()
                if row.get("status") == "ok"
            )
            else "in_progress"
        ),
        "run_type": "incremental_supplement_not_rebuild",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "download_completed_targets": len(download_state),
        "download_ok_targets": sum(
            row.get("status") == "ok" for row in download_state.values()
        ),
        "download_no_data_targets": sum(
            row.get("status") == "no_data" for row in download_state.values()
        ),
        "normalized_ok_targets": len(completed),
        "run_counts": dict(sorted(run_counts.items())),
        "output_directory": str(OUT_DIR),
        "progress_checkpoint": str(PROGRESS),
    }
    temporary = SUMMARY.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    temporary.replace(SUMMARY)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    fields = load_schema() + EXTRA_FIELDS
    proteins = load_proteins()
    cid_map, ambiguous_cids = load_cid_mappings()
    completed = load_completed()
    run_counts: Counter[str] = Counter()

    while True:
        if shutil.disk_usage(ROOT.drive + "\\").free < MIN_FREE_BYTES:
            run_counts["stopped_low_disk"] += 1
            write_summary(load_download_state(), completed, run_counts)
            return 2
        download_state = load_download_state()
        available = sorted(
            target
            for target, row in download_state.items()
            if row.get("status") == "ok"
            and target not in completed
            and (RAW_DIR / f"{target}.csv.gz").exists()
        )
        for target in available:
            try:
                qa = process_target(
                    target,
                    proteins.get(target, ""),
                    cid_map,
                    ambiguous_cids,
                    fields,
                )
            except Exception as exc:
                run_counts["target_errors"] += 1
                error_path = QA_DIR / f"{target}.error.json"
                error_path.write_text(
                    json.dumps(
                        {
                            "status": "error",
                            "target_uniprot_id": target,
                            "error": repr(exc),
                            "generated_at_utc": dt.datetime.now(
                                dt.timezone.utc
                            ).isoformat(),
                        },
                        ensure_ascii=False,
                        indent=2,
                    )
                    + "\n",
                    encoding="utf-8",
                )
                continue
            with PROGRESS.open("a", encoding="utf-8") as handle:
                handle.write(
                    json.dumps(
                        {
                            "status": "ok",
                            "target_uniprot_id": target,
                            "normalized_rows": qa["counts"].get(
                                "normalized_rows", 0
                            ),
                            "completed_at_utc": qa["generated_at_utc"],
                        },
                        ensure_ascii=False,
                    )
                    + "\n"
                )
            completed.add(target)
            run_counts["targets_normalized_this_process"] += 1
            run_counts["rows_normalized_this_process"] += int(
                qa["counts"].get("normalized_rows", 0)
            )
            write_summary(download_state, completed, run_counts)

        write_summary(download_state, completed, run_counts)
        fully_downloaded = (
            len(download_state) == 10556
            and set(row.get("status") for row in download_state.values())
            <= {"ok", "no_data"}
        )
        all_ok_normalized = all(
            target in completed
            for target, row in download_state.items()
            if row.get("status") == "ok"
        )
        if not args.watch or (fully_downloaded and all_ok_normalized):
            break
        time.sleep(max(10, args.poll_seconds))
    return 0


if __name__ == "__main__":
    sys.exit(main())
