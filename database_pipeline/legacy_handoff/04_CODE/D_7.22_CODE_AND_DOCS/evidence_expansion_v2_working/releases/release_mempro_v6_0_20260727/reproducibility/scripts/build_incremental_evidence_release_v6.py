#!/usr/bin/env python3
"""Build the high-confidence incremental evidence layer for MemProDB V6.

The published V5.5/V5.4/V4.2/V1.0 baselines are read-only.  This builder:
* retains uniquely mapped positive evidence in the binding master;
* keeps inactive PubChem results in a separate negative-evidence table;
* keeps unresolved, inconclusive, artifact and unsuitable-ligand records in a
  review table;
* flags exact quantitative source mirrors and excludes them from default
  pair-level counting without deleting their provenance;
* rebuilds binding-site and protein-compound pair summaries.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import pathlib
import re
import sqlite3
import sys
from collections import Counter, defaultdict


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "incremental_20260727"
BASELINE_CONFIG = ROOT / "config" / "baseline_paths.json"
OUT = RUN / "release_candidate_v6_0"
QA_DIR = RUN / "qa"
DB_PATH = RUN / "merge_candidate" / "incremental_merge_v6.sqlite"

GPCR = ROOT / "intermediate" / "gpcrdb_default_new_evidence_v2.tsv"
PDBe = RUN / "staging" / "pdbe" / "pdbe_normalized_evidence_v3.tsv.gz"
PDBBIND = RUN / "staging" / "pdbbind" / "pdbbind_normalized_evidence_v3.tsv.gz"
BRENDA = RUN / "staging" / "brenda" / "brenda_normalized_evidence_v3.tsv.gz"
PUBCHEM_DIR = RUN / "staging" / "pubchem" / "by_target"

RELEASE_VERSION = "v6.0"
RELEASE_DATE = "2026-07-27"
TIER_RANK = {"BE1": 1, "BE2": 2, "BE3": 3}
VALID_MAPPING = {"mapped_core", "mapped_extended"}


def open_text(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def split_values(value: str) -> list[str]:
    return [
        part.strip()
        for part in re.split(r"[;|,]", value or "")
        if part.strip()
    ]


def numeric(value: str) -> float | None:
    try:
        result = float(str(value).strip())
        return result if result >= 0 else None
    except (TypeError, ValueError):
        return None


def integer(value: str) -> int:
    try:
        return int(float(str(value).strip() or 0))
    except (TypeError, ValueError):
        return 0


def quantitative_key(
    row: dict[str, str],
) -> tuple[str, str, str, str, str] | None:
    value = numeric(row.get("standard_value_nM", ""))
    target = (row.get("target_uniprot_id") or "").strip()
    compound = (row.get("compound_internal_id") or "").strip()
    activity = (row.get("activity_type") or "").strip().upper()
    assay = " ".join(
        (row.get("assay_or_mechanism") or "").lower().split()
    )
    context = assay[:240] or "|".join(
        sorted(split_values(row.get("pubmed_ids", "")))
    )
    if (
        value is None
        or value <= 0
        or not target
        or not compound
        or not activity
        or not context
    ):
        return None
    rounded = f"{value:.8g}"
    return target, compound, activity, rounded, context


def stable_id(prefix: str, *parts: str) -> str:
    payload = "\x1f".join(parts).encode("utf-8")
    return f"{prefix}_{hashlib.sha1(payload).hexdigest()[:20].upper()}"


def atomic_gzip_writer(path: pathlib.Path, fieldnames: list[str]):
    temporary = path.with_name(path.name + ".tmp")
    handle = gzip.open(
        temporary,
        "wt",
        encoding="utf-8",
        newline="",
        compresslevel=6,
    )
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fieldnames,
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    return temporary, handle, writer


def source_id_type(value: str) -> str:
    value = value or ""
    if value.startswith("PubChem CID:"):
        return "PubChem CID"
    if value.startswith("PDBCCD:"):
        return "PDB CCD ID"
    if value.startswith("PDBBIND:"):
        return "PDBbind complex ID"
    if value.startswith("GPCRLIG-"):
        return "GPCRdb ligand record"
    if value.startswith("BRENDA-NAME-"):
        return "BRENDA normalized name"
    return "source compound ID"


def directness(tier: str) -> str:
    return {
        "BE1": "direct_structural",
        "BE2": "direct_quantitative",
        "BE3": "functional_or_pharmacological",
    }.get(tier, "supporting")


def classify_new_row(row: dict[str, str], source_group: str) -> tuple[str, str]:
    outcome = (row.get("activity_outcome") or "").strip().lower()
    tier = (row.get("evidence_tier") or "").strip()
    mapping = (row.get("compound_mapping_status") or "").strip()
    target = (row.get("target_uniprot_id") or "").strip()
    compound = (row.get("compound_internal_id") or "").strip()

    if outcome == "inactive" or tier == "negative":
        return "negative", "reported inactive"
    if not target:
        return "review", "missing protein identity"
    if mapping not in VALID_MAPPING or not compound:
        return "review", "compound identity unresolved"
    if outcome == "inconclusive" or tier == "nondefault":
        return "review", "inconclusive activity outcome"
    if source_group == "PDBe" and (row.get("artifact_class") or "").strip():
        return "review", "structural artifact or nondefault component"
    if source_group == "PDBbind":
        flags = [
            row.get("incomplete_ligand_flag") == "1",
            row.get("peptide_or_oligomer_flag") == "1",
        ]
        if any(flags):
            return "review", "incomplete or non-small-molecule ligand"
    return "binding", ""


def binding_row(
    row: dict[str, str],
    baseline_fields: list[str],
    protein: dict[str, dict[str, str]],
    compound: dict[str, dict[str, str]],
    source_group: str,
    duplicate_status: str,
    effective_default: str,
) -> dict[str, str]:
    target = row.get("target_uniprot_id", "")
    internal = row.get("compound_internal_id", "")
    pinfo = protein[target]
    cinfo = compound[internal]
    result = {field: "" for field in baseline_fields}
    result.update(
        {
            "evidence_id": row.get("source_evidence_id", ""),
            "target_uniprot_id": target,
            "approved_symbol": row.get("approved_symbol", "")
            or pinfo.get("approved_symbol", ""),
            "protein_release_tier": pinfo.get("evidence_level_v52", ""),
            "protein_website_default": pinfo.get("website_default_v52", ""),
            "compound_id": row.get("compound_source_id", ""),
            "compound_id_type": source_id_type(
                row.get("compound_source_id", "")
            ),
            "compound_name": row.get("compound_name", "")
            or cinfo.get("preferred_name", ""),
            "small_molecule_scope_status": cinfo.get(
                "compound_scope_status", ""
            ),
            "is_core_small_molecule": str(
                int(cinfo.get("compound_scope_status") == "core")
            ),
            "evidence_tier": row.get("evidence_tier", ""),
            "evidence_type": row.get("evidence_type", ""),
            "evidence_directness": directness(row.get("evidence_tier", "")),
            "target_assignment_status": row.get("target_mapping_status", ""),
            "relationship_context_id": row.get("source_record_id", ""),
            "activity_type": row.get("activity_type", ""),
            "activity_relation": row.get("activity_relation", ""),
            "activity_value": row.get("activity_value", ""),
            "activity_unit": row.get("activity_unit", ""),
            "standard_value_nM": row.get("standard_value_nM", ""),
            "assay_or_mechanism": row.get("assay_or_mechanism", ""),
            "pdb_ids": row.get("pdb_ids", ""),
            "binding_site_residues": row.get("binding_site_residues", ""),
            "ligand_het_id": row.get("ligand_het_id", ""),
            "pubmed_ids": row.get("pubmed_ids", ""),
            "doi": row.get("doi", ""),
            "source_database": row.get("source_database", ""),
            "source_version": row.get("source_version", ""),
            "source_record_id": row.get("source_record_id", ""),
            "source_url": row.get("source_url", ""),
            "record_qc_status": row.get("record_qc_status", ""),
            "default_release_inclusion": effective_default,
            "originating_dataset": source_group,
            "compound_internal_id": internal,
            "compound_form_id": row.get("compound_form_id", ""),
            "compound_mapping_status": row.get(
                "compound_mapping_status", ""
            ),
            "compound_scope_status_v1": cinfo.get(
                "compound_scope_status", ""
            ),
            "default_release_inclusion_v4_2": effective_default,
            "compound_master_version": "v1.0",
            "activity_outcome_v60": row.get("activity_outcome", ""),
            "incremental_source_v60": source_group,
            "duplicate_status_v60": duplicate_status,
            "identity_status_v60": "canonical_compound_and_protein",
            "conflict_status_v60": "context_preserved",
            "release_version_v60": RELEASE_VERSION,
        }
    )
    return result


def candidate_row(
    row: dict[str, str], source_group: str, bucket: str, reason: str
) -> dict[str, str]:
    return {
        "source_evidence_id": row.get("source_evidence_id", ""),
        "target_uniprot_id": row.get("target_uniprot_id", ""),
        "approved_symbol": row.get("approved_symbol", ""),
        "compound_source_id": row.get("compound_source_id", ""),
        "compound_name": row.get("compound_name", ""),
        "compound_internal_id": row.get("compound_internal_id", ""),
        "compound_form_id": row.get("compound_form_id", ""),
        "compound_mapping_status": row.get("compound_mapping_status", ""),
        "evidence_tier": row.get("evidence_tier", ""),
        "evidence_type": row.get("evidence_type", ""),
        "activity_type": row.get("activity_type", ""),
        "activity_relation": row.get("activity_relation", ""),
        "activity_value": row.get("activity_value", ""),
        "activity_unit": row.get("activity_unit", ""),
        "standard_value_nM": row.get("standard_value_nM", ""),
        "activity_outcome": row.get("activity_outcome", ""),
        "assay_or_mechanism": row.get("assay_or_mechanism", ""),
        "pdb_ids": row.get("pdb_ids", ""),
        "binding_site_residues": row.get("binding_site_residues", ""),
        "source_database": row.get("source_database", ""),
        "source_version": row.get("source_version", ""),
        "source_record_id": row.get("source_record_id", ""),
        "source_url": row.get("source_url", ""),
        "record_qc_status": row.get("record_qc_status", ""),
        "original_default_release_inclusion": row.get(
            "default_release_inclusion", ""
        ),
        "original_exclusion_reason": row.get("exclusion_reason", ""),
        "review_bucket": bucket,
        "review_reason": reason,
        "incremental_source": source_group,
        "release_version": RELEASE_VERSION,
    }


def negative_row(row: dict[str, str], source_group: str) -> dict[str, str]:
    result = candidate_row(
        row, source_group, "negative_evidence", "reported inactive"
    )
    result.update(
        {
            "pubchem_aid": row.get("pubchem_aid", ""),
            "pubchem_sid": row.get("pubchem_sid", ""),
            "pubchem_cid": row.get("pubchem_cid", ""),
        }
    )
    return result


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    config = json.loads(BASELINE_CONFIG.read_text(encoding="utf-8"))
    protein_path = pathlib.Path(config["protein_master_v5_5"])
    compound_path = pathlib.Path(config["compound_master_v1_0"])
    form_path = pathlib.Path(config["compound_forms_v1_0"])
    baseline_binding = pathlib.Path(config["binding_evidence_v4_2"])
    baseline_sites = pathlib.Path(config["binding_sites_v4_2"])
    baseline_pairs = pathlib.Path(config["protein_compound_summary_v4_2"])

    protein: dict[str, dict[str, str]] = {}
    with protein_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target = row.get("target_uniprot_id", "")
            if row.get("membrane_class_v52") in {"A", "B", "C"}:
                protein[target] = row

    compound: dict[str, dict[str, str]] = {}
    with compound_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            compound[row["compound_internal_id"]] = {
                "preferred_name": row.get("preferred_name", ""),
                "compound_scope_status": row.get(
                    "compound_scope_status", ""
                ),
            }

    valid_forms: set[str] = set()
    with form_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            valid_forms.add(row["compound_form_id"])

    if DB_PATH.exists():
        DB_PATH.unlink()
    con = sqlite3.connect(DB_PATH)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA temp_store=FILE")
    con.executescript(
        """
        CREATE TABLE baseline_pair(
            target TEXT NOT NULL,
            compound TEXT NOT NULL,
            approved_symbol TEXT,
            preferred_name TEXT,
            best_tier TEXT,
            be1 INTEGER,
            be2 INTEGER,
            be3 INTEGER,
            evidence_count INTEGER,
            sources TEXT,
            best_value REAL,
            form_count INTEGER,
            PRIMARY KEY(target, compound)
        );
        CREATE TABLE new_pair_raw(
            target TEXT NOT NULL,
            compound TEXT NOT NULL,
            approved_symbol TEXT,
            tier TEXT,
            tier_rank INTEGER,
            source TEXT,
            value REAL,
            form_id TEXT,
            evidence_id TEXT NOT NULL
        );
        CREATE TABLE quant_obs(
            target TEXT NOT NULL,
            compound TEXT NOT NULL,
            activity_type TEXT NOT NULL,
            value REAL NOT NULL,
            source TEXT NOT NULL
        );
        """
    )

    with baseline_pairs.open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        batch = []
        for row in csv.DictReader(handle, delimiter="\t"):
            batch.append(
                (
                    row["target_uniprot_id"],
                    row["compound_internal_id"],
                    row.get("approved_symbol", ""),
                    row.get("preferred_name", ""),
                    row.get("best_binding_evidence_level", ""),
                    integer(row.get("BE1_evidence_count", "")),
                    integer(row.get("BE2_evidence_count", "")),
                    integer(row.get("BE3_evidence_count", "")),
                    integer(row.get("binding_evidence_count", "")),
                    row.get("independent_sources", ""),
                    numeric(row.get("best_standard_value_nM", "")),
                    integer(row.get("compound_form_count", "")),
                )
            )
            if len(batch) >= 5000:
                con.executemany(
                    "INSERT INTO baseline_pair VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                batch.clear()
        if batch:
            con.executemany(
                "INSERT INTO baseline_pair VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                batch,
            )
    con.commit()

    with baseline_binding.open(
        "r", encoding="utf-8-sig", newline=""
    ) as input_handle:
        baseline_reader = csv.DictReader(input_handle, delimiter="\t")
        baseline_fields = list(baseline_reader.fieldnames or [])

    appended_binding_fields = [
        "activity_outcome_v60",
        "incremental_source_v60",
        "duplicate_status_v60",
        "identity_status_v60",
        "conflict_status_v60",
        "release_version_v60",
    ]
    binding_fields = baseline_fields + appended_binding_fields
    candidate_fields = [
        "source_evidence_id",
        "target_uniprot_id",
        "approved_symbol",
        "compound_source_id",
        "compound_name",
        "compound_internal_id",
        "compound_form_id",
        "compound_mapping_status",
        "evidence_tier",
        "evidence_type",
        "activity_type",
        "activity_relation",
        "activity_value",
        "activity_unit",
        "standard_value_nM",
        "activity_outcome",
        "assay_or_mechanism",
        "pdb_ids",
        "binding_site_residues",
        "source_database",
        "source_version",
        "source_record_id",
        "source_url",
        "record_qc_status",
        "original_default_release_inclusion",
        "original_exclusion_reason",
        "review_bucket",
        "review_reason",
        "incremental_source",
        "release_version",
    ]
    negative_fields = candidate_fields + [
        "pubchem_aid",
        "pubchem_sid",
        "pubchem_cid",
    ]
    conflict_fields = [
        "target_uniprot_id",
        "compound_internal_id",
        "pubchem_cid",
        "pubchem_aid",
        "active_row_count",
        "inactive_row_count",
        "inconclusive_row_count",
        "conflict_type",
        "handling",
    ]

    binding_path = OUT / "binding_evidence_master_v6_0.tsv.gz"
    candidate_path = OUT / "binding_evidence_review_queue_v6_0.tsv.gz"
    negative_path = OUT / "negative_binding_evidence_v1_0.tsv.gz"
    context_conflict_path = OUT / "pubchem_context_conflicts_v1_0.tsv.gz"
    binding_tmp, binding_handle, binding_writer = atomic_gzip_writer(
        binding_path, binding_fields
    )
    candidate_tmp, candidate_handle, candidate_writer = atomic_gzip_writer(
        candidate_path, candidate_fields
    )
    negative_tmp, negative_handle, negative_writer = atomic_gzip_writer(
        negative_path, negative_fields
    )
    conflict_tmp, conflict_handle, conflict_writer = atomic_gzip_writer(
        context_conflict_path, conflict_fields
    )

    counts: Counter[str] = Counter()
    source_counts: dict[str, Counter[str]] = defaultdict(Counter)
    source_default_targets: dict[str, set[str]] = defaultdict(set)
    source_default_compounds: dict[str, set[str]] = defaultdict(set)
    baseline_quant: set[tuple[str, str, str, str, str]] = set()
    new_quant_seen: dict[tuple[str, str, str, str, str], str] = {}
    new_evidence_ids: set[str] = set()
    site_map: dict[tuple[str, ...], dict[str, object]] = {}
    pair_batch: list[tuple[object, ...]] = []
    quant_batch: list[tuple[object, ...]] = []

    with baseline_binding.open(
        "r", encoding="utf-8-sig", newline=""
    ) as input_handle:
        for row in csv.DictReader(input_handle, delimiter="\t"):
            output = dict(row)
            output.update(
                {
                    "activity_outcome_v60": "positive_or_observed",
                    "incremental_source_v60": "frozen_v4_2_baseline",
                    "duplicate_status_v60": "baseline_record",
                    "identity_status_v60": "baseline_canonical",
                    "conflict_status_v60": "baseline_policy_preserved",
                    "release_version_v60": RELEASE_VERSION,
                }
            )
            binding_writer.writerow(output)
            counts["baseline_binding_rows"] += 1
            key = quantitative_key(row)
            if key:
                baseline_quant.add(key)
                if (
                    row.get("default_release_inclusion_v4_2") == "1"
                    and numeric(row.get("standard_value_nM", "")) not in (None, 0)
                ):
                    quant_batch.append(
                        (
                            row.get("target_uniprot_id", ""),
                            row.get("compound_internal_id", ""),
                            (row.get("activity_type") or "").upper(),
                            numeric(row.get("standard_value_nM", "")),
                            row.get("source_database", ""),
                        )
                    )
                    if len(quant_batch) >= 5000:
                        con.executemany(
                            "INSERT INTO quant_obs VALUES(?,?,?,?,?)",
                            quant_batch,
                        )
                        quant_batch.clear()
    if quant_batch:
        con.executemany("INSERT INTO quant_obs VALUES(?,?,?,?,?)", quant_batch)
        quant_batch.clear()
    con.commit()

    def flush_pair_batches() -> None:
        nonlocal pair_batch, quant_batch
        if pair_batch:
            con.executemany(
                "INSERT INTO new_pair_raw VALUES(?,?,?,?,?,?,?,?,?)",
                pair_batch,
            )
            pair_batch.clear()
        if quant_batch:
            con.executemany(
                "INSERT INTO quant_obs VALUES(?,?,?,?,?)", quant_batch
            )
            quant_batch.clear()

    def process_row(row: dict[str, str], source_group: str) -> None:
        source_stat = source_counts[source_group]
        source_stat["input_rows"] += 1
        bucket, reason = classify_new_row(row, source_group)
        target = (row.get("target_uniprot_id") or "").strip()
        internal = (row.get("compound_internal_id") or "").strip()
        form = (row.get("compound_form_id") or "").strip()

        if bucket == "negative":
            negative_writer.writerow(negative_row(row, source_group))
            counts["negative_rows"] += 1
            source_stat["negative_rows"] += 1
            return
        if target not in protein:
            bucket, reason = "review", "protein outside canonical A/B/C universe"
        if internal and internal not in compound:
            bucket, reason = "review", "compound ID absent from canonical master"
        if form and form not in valid_forms:
            bucket, reason = "review", "compound form ID absent from form hierarchy"
        evidence_id = (row.get("source_evidence_id") or "").strip()
        if bucket == "binding" and evidence_id in new_evidence_ids:
            bucket, reason = "review", "duplicate source evidence ID"

        if bucket != "binding":
            candidate_writer.writerow(
                candidate_row(row, source_group, bucket, reason)
            )
            counts["review_rows"] += 1
            source_stat["review_rows"] += 1
            return

        new_evidence_ids.add(evidence_id)
        qkey = quantitative_key(row)
        duplicate_status = "unique_or_complementary"
        if qkey and qkey in baseline_quant:
            duplicate_status = "baseline_quantitative_mirror_exact"
        elif qkey and qkey in new_quant_seen:
            duplicate_status = "incremental_quantitative_mirror_exact"
        if qkey and qkey not in new_quant_seen:
            new_quant_seen[qkey] = source_group

        effective_default = row.get("default_release_inclusion", "0")
        if duplicate_status != "unique_or_complementary":
            effective_default = "0"
            source_stat["exact_mirror_rows"] += 1
            counts["exact_mirror_rows"] += 1

        binding_writer.writerow(
            binding_row(
                row,
                baseline_fields,
                protein,
                compound,
                source_group,
                duplicate_status,
                effective_default,
            )
        )
        counts["new_binding_rows"] += 1
        source_stat["binding_rows"] += 1
        if effective_default == "1":
            counts["new_default_binding_rows"] += 1
            source_stat["default_binding_rows"] += 1
            source_default_targets[source_group].add(target)
            source_default_compounds[source_group].add(internal)
            tier = row.get("evidence_tier", "")
            value = numeric(row.get("standard_value_nM", ""))
            pair_batch.append(
                (
                    target,
                    internal,
                    row.get("approved_symbol", "")
                    or protein[target].get("approved_symbol", ""),
                    tier,
                    TIER_RANK.get(tier, 99),
                    row.get("source_database", "") or source_group,
                    value,
                    form,
                    evidence_id,
                )
            )
            if qkey and value not in (None, 0):
                quant_batch.append(
                    (
                        target,
                        internal,
                        (row.get("activity_type") or "").upper(),
                        value,
                        row.get("source_database", "") or source_group,
                    )
                )
            if row.get("binding_site_residues") and row.get("pdb_ids"):
                for pdb_id in split_values(row.get("pdb_ids", "")):
                    coordinate = (
                        row.get("residue_index_type")
                        or row.get("pocket_residue_index_type")
                        or "unspecified"
                    )
                    chains = row.get("pdb_chain_ids", "")
                    key = (
                        target,
                        internal,
                        pdb_id.lower(),
                        coordinate,
                        chains,
                        source_group,
                    )
                    item = site_map.setdefault(
                        key,
                        {
                            "evidence_ids": set(),
                            "residues": set(),
                            "compound_source_ids": set(),
                            "form_ids": set(),
                            "mapping_status": row.get(
                                "compound_mapping_status", ""
                            ),
                        },
                    )
                    item["evidence_ids"].add(evidence_id)
                    item["residues"].update(
                        split_values(row.get("binding_site_residues", ""))
                    )
                    if row.get("compound_source_id"):
                        item["compound_source_ids"].add(
                            row["compound_source_id"]
                        )
                    if form:
                        item["form_ids"].add(form)
        if len(pair_batch) >= 5000 or len(quant_batch) >= 5000:
            flush_pair_batches()

    source_specs = [
        (GPCR, "GPCRdb"),
        (PDBe, "PDBe"),
        (PDBBIND, "PDBbind"),
        (BRENDA, "BRENDA"),
    ]
    for path, source_group in source_specs:
        with open_text(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                process_row(row, source_group)
        flush_pair_batches()
        con.commit()
        print(
            f"{source_group}: {dict(source_counts[source_group])}",
            flush=True,
        )

    pubchem_files = sorted(PUBCHEM_DIR.glob("*.tsv.gz"))
    for file_number, path in enumerate(pubchem_files, 1):
        context: dict[tuple[str, str, str], list[object]] = {}
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                aid = row.get("pubchem_aid", "")
                cid = row.get("pubchem_cid", "")
                internal = row.get("compound_internal_id", "")
                context_key = (aid, cid, internal)
                state = context.setdefault(context_key, [0, 0, 0])
                outcome = (row.get("activity_outcome") or "").lower()
                if outcome == "active":
                    state[0] += 1
                elif outcome == "inactive":
                    state[1] += 1
                elif outcome == "inconclusive":
                    state[2] += 1
                process_row(row, "PubChem BioAssay")
        target = path.name.removesuffix(".tsv.gz")
        for (aid, cid, internal), state in context.items():
            if state[0] and state[1]:
                conflict_writer.writerow(
                    {
                        "target_uniprot_id": target,
                        "compound_internal_id": internal,
                        "pubchem_cid": cid,
                        "pubchem_aid": aid,
                        "active_row_count": state[0],
                        "inactive_row_count": state[1],
                        "inconclusive_row_count": state[2],
                        "conflict_type": "same_assay_active_inactive",
                        "handling": (
                            "preserve both outcomes outside default pair "
                            "count; require assay-context review"
                        ),
                    }
                )
                counts["pubchem_exact_context_conflicts"] += 1
        if file_number % 250 == 0:
            flush_pair_batches()
            con.commit()
            print(
                f"PubChem files={file_number}/{len(pubchem_files)} "
                f"input={source_counts['PubChem BioAssay']['input_rows']} "
                f"binding={source_counts['PubChem BioAssay']['binding_rows']} "
                f"negative={source_counts['PubChem BioAssay']['negative_rows']}",
                flush=True,
            )

    flush_pair_batches()
    con.commit()
    binding_handle.close()
    candidate_handle.close()
    negative_handle.close()
    conflict_handle.close()
    binding_tmp.replace(binding_path)
    candidate_tmp.replace(candidate_path)
    negative_tmp.replace(negative_path)
    conflict_tmp.replace(context_conflict_path)

    con.executescript(
        """
        CREATE INDEX idx_new_pair_key
        ON new_pair_raw(target, compound);
        CREATE INDEX idx_quant_key
        ON quant_obs(target, compound, activity_type);
        CREATE TABLE new_pair_agg AS
        SELECT target, compound,
               MAX(approved_symbol) AS approved_symbol,
               MIN(tier_rank) AS best_rank,
               SUM(CASE WHEN tier='BE1' THEN 1 ELSE 0 END) AS be1,
               SUM(CASE WHEN tier='BE2' THEN 1 ELSE 0 END) AS be2,
               SUM(CASE WHEN tier='BE3' THEN 1 ELSE 0 END) AS be3,
               COUNT(*) AS evidence_count,
               GROUP_CONCAT(DISTINCT source) AS sources,
               MIN(value) AS best_value,
               COUNT(DISTINCT CASE WHEN form_id<>'' THEN form_id END) AS form_count
        FROM new_pair_raw
        GROUP BY target, compound;
        CREATE UNIQUE INDEX idx_new_pair_agg_key
        ON new_pair_agg(target, compound);
        """
    )
    con.commit()

    pair_fields = [
        "target_uniprot_id",
        "approved_symbol",
        "compound_internal_id",
        "preferred_name",
        "best_binding_evidence_level",
        "BE1_evidence_count",
        "BE2_evidence_count",
        "BE3_evidence_count",
        "binding_evidence_count",
        "independent_source_count",
        "independent_sources",
        "best_standard_value_nM",
        "compound_form_count",
        "compound_master_version",
        "new_evidence_count_v60",
        "new_sources_v60",
        "pair_action_v60",
        "release_version_v60",
    ]
    pair_path = OUT / "protein_compound_summary_v6_0.tsv.gz"
    pair_tmp, pair_handle, pair_writer = atomic_gzip_writer(
        pair_path, pair_fields
    )
    pair_query = """
        WITH all_keys AS (
            SELECT target, compound FROM baseline_pair
            UNION
            SELECT target, compound FROM new_pair_agg
        )
        SELECT k.target, k.compound,
               b.approved_symbol, b.preferred_name, b.best_tier,
               b.be1, b.be2, b.be3, b.evidence_count, b.sources,
               b.best_value, b.form_count,
               n.approved_symbol, n.best_rank, n.be1, n.be2, n.be3,
               n.evidence_count, n.sources, n.best_value, n.form_count
        FROM all_keys k
        LEFT JOIN baseline_pair b
          ON b.target=k.target AND b.compound=k.compound
        LEFT JOIN new_pair_agg n
          ON n.target=k.target AND n.compound=k.compound
        ORDER BY k.target, k.compound
    """
    for row in con.execute(pair_query):
        (
            target,
            internal,
            b_symbol,
            b_name,
            b_tier,
            b_be1,
            b_be2,
            b_be3,
            b_count,
            b_sources,
            b_value,
            b_forms,
            n_symbol,
            n_rank,
            n_be1,
            n_be2,
            n_be3,
            n_count,
            n_sources,
            n_value,
            n_forms,
        ) = row
        ranks = [
            TIER_RANK[tier]
            for tier in [b_tier]
            if tier in TIER_RANK
        ]
        if n_rank and n_rank in TIER_RANK.values():
            ranks.append(n_rank)
        best_rank = min(ranks) if ranks else None
        best_tier = (
            next(
                (tier for tier, rank in TIER_RANK.items() if rank == best_rank),
                "",
            )
            if best_rank
            else ""
        )
        sources = sorted(
            set(split_values(b_sources or ""))
            | set(split_values(n_sources or ""))
        )
        values = [
            value
            for value in [b_value, n_value]
            if value is not None
        ]
        pair_writer.writerow(
            {
                "target_uniprot_id": target,
                "approved_symbol": b_symbol or n_symbol or "",
                "compound_internal_id": internal,
                "preferred_name": compound.get(internal, {}).get(
                    "preferred_name", b_name or ""
                ),
                "best_binding_evidence_level": best_tier,
                "BE1_evidence_count": (b_be1 or 0) + (n_be1 or 0),
                "BE2_evidence_count": (b_be2 or 0) + (n_be2 or 0),
                "BE3_evidence_count": (b_be3 or 0) + (n_be3 or 0),
                "binding_evidence_count": (b_count or 0) + (n_count or 0),
                "independent_source_count": len(sources),
                "independent_sources": ";".join(sources),
                "best_standard_value_nM": min(values) if values else "",
                "compound_form_count": max(b_forms or 0, n_forms or 0),
                "compound_master_version": "v1.0",
                "new_evidence_count_v60": n_count or 0,
                "new_sources_v60": ";".join(
                    sorted(set(split_values(n_sources or "")))
                ),
                "pair_action_v60": (
                    "strengthen_existing_pair"
                    if b_count is not None and n_count is not None
                    else "add_new_pair"
                    if n_count is not None
                    else "baseline_unchanged"
                ),
                "release_version_v60": RELEASE_VERSION,
            }
        )
        counts["pair_summary_rows"] += 1
        if b_count is None:
            counts["new_pair_rows"] += 1
        elif n_count is not None:
            counts["strengthened_pair_rows"] += 1
    pair_handle.close()
    pair_tmp.replace(pair_path)

    with baseline_sites.open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        site_reader = csv.DictReader(handle, delimiter="\t")
        site_fields = list(site_reader.fieldnames or []) + [
            "residue_index_type_v60",
            "pdb_chain_ids_v60",
            "site_merge_status_v60",
            "release_version_v60",
        ]
    site_path = OUT / "binding_site_instances_v6_0.tsv.gz"
    site_tmp, site_handle, site_writer = atomic_gzip_writer(
        site_path, site_fields
    )
    with baseline_sites.open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row.update(
                {
                    "residue_index_type_v60": "",
                    "pdb_chain_ids_v60": "",
                    "site_merge_status_v60": "baseline_site",
                    "release_version_v60": RELEASE_VERSION,
                }
            )
            site_writer.writerow(row)
            counts["baseline_site_rows"] += 1
    for key, item in sorted(site_map.items()):
        target, internal, pdb_id, coordinate, chains, source_group = key
        evidence_ids = sorted(item["evidence_ids"])
        residues = sorted(item["residues"])
        source_ids = sorted(item["compound_source_ids"])
        forms = sorted(item["form_ids"])
        site_writer.writerow(
            {
                "binding_site_instance_id": stable_id(
                    "BSI6", *key, ";".join(residues)
                ),
                "evidence_id": ";".join(evidence_ids),
                "target_uniprot_id": target,
                "compound_id": ";".join(source_ids),
                "site_type": (
                    "experimental_structure_residue_contact"
                    if source_group == "PDBe"
                    else "experimental_complex_pocket"
                ),
                "pdb_ids": pdb_id,
                "residue_or_site_description": ";".join(residues),
                "site_compound_specificity": "compound_specific",
                "source_database": source_group,
                "record_qc_status": "ok",
                "compound_internal_id": internal,
                "compound_form_id": ";".join(forms),
                "compound_mapping_status": item["mapping_status"],
                "compound_master_version": "v1.0",
                "residue_index_type_v60": coordinate,
                "pdb_chain_ids_v60": chains,
                "site_merge_status_v60": (
                    "merged_same_source_structure_context"
                    if len(evidence_ids) > 1
                    else "single_structure_context"
                ),
                "release_version_v60": RELEASE_VERSION,
            }
        )
        counts["new_site_rows"] += 1
    site_handle.close()
    site_tmp.replace(site_path)

    heterogeneity_fields = [
        "target_uniprot_id",
        "compound_internal_id",
        "activity_type",
        "minimum_standard_value_nM",
        "maximum_standard_value_nM",
        "max_min_ratio",
        "observation_count",
        "source_count",
        "sources",
        "interpretation",
    ]
    heterogeneity_path = OUT / "quantitative_heterogeneity_audit_v1_0.tsv.gz"
    h_tmp, h_handle, h_writer = atomic_gzip_writer(
        heterogeneity_path, heterogeneity_fields
    )
    for row in con.execute(
        """
        SELECT target, compound, activity_type, MIN(value), MAX(value),
               COUNT(*), COUNT(DISTINCT source),
               GROUP_CONCAT(DISTINCT source)
        FROM quant_obs
        WHERE value>0 AND activity_type<>''
        GROUP BY target, compound, activity_type
        HAVING COUNT(DISTINCT source)>=2 AND MAX(value)/MIN(value)>=100
        ORDER BY MAX(value)/MIN(value) DESC
        """
    ):
        target, internal, activity, minimum, maximum, n, ns, sources = row
        h_writer.writerow(
            {
                "target_uniprot_id": target,
                "compound_internal_id": internal,
                "activity_type": activity,
                "minimum_standard_value_nM": minimum,
                "maximum_standard_value_nM": maximum,
                "max_min_ratio": maximum / minimum,
                "observation_count": n,
                "source_count": ns,
                "sources": sources,
                "interpretation": (
                    "contextual assay heterogeneity; preserve source-specific "
                    "values and do not average"
                ),
            }
        )
        counts["quantitative_heterogeneity_rows"] += 1
    h_handle.close()
    h_tmp.replace(heterogeneity_path)

    coverage_path = OUT / "incremental_source_coverage_v6_0.tsv"
    with coverage_path.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = [
            "source_group",
            "input_rows",
            "binding_rows",
            "default_binding_rows",
            "negative_rows",
            "review_rows",
            "exact_mirror_rows",
            "default_unique_targets",
            "default_unique_compounds",
        ]
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        for source_group in [
            "GPCRdb",
            "PDBe",
            "PDBbind",
            "BRENDA",
            "PubChem BioAssay",
        ]:
            stat = source_counts[source_group]
            writer.writerow(
                {
                    "source_group": source_group,
                    **{field: stat.get(field, 0) for field in fields[1:7]},
                    "default_unique_targets": len(
                        source_default_targets[source_group]
                    ),
                    "default_unique_compounds": len(
                        source_default_compounds[source_group]
                    ),
                }
            )

    input_total = sum(
        stat.get("input_rows", 0) for stat in source_counts.values()
    )
    reconciled = (
        counts["new_binding_rows"]
        + counts["negative_rows"]
        + counts["review_rows"]
    )
    blocking: list[str] = []
    if input_total != reconciled:
        blocking.append(
            f"source row reconciliation failed: input={input_total}, "
            f"assigned={reconciled}"
        )
    if counts["new_default_binding_rows"] <= 0:
        blocking.append("no new default binding evidence rows")

    report = {
        "status": "passed" if not blocking else "failed",
        "run_type": "incremental_supplement_not_rebuild",
        "release_version": RELEASE_VERSION,
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "baseline_policy": "read_only_no_overwrite",
        "counts": dict(sorted(counts.items())),
        "source_counts": {
            source: dict(sorted(stat.items()))
            for source, stat in sorted(source_counts.items())
        },
        "row_reconciliation": {
            "input_source_rows": input_total,
            "assigned_source_rows": reconciled,
            "difference": input_total - reconciled,
        },
        "identity_policy": (
            "only canonical A/B/C proteins and uniquely mapped V1.0 compounds "
            "enter the binding master; unresolved identities remain in review"
        ),
        "duplicate_policy": (
            "exact target-compound-activity-value mirrors are retained for "
            "provenance but excluded from default pair counting"
        ),
        "conflict_policy": (
            "negative results are separated; same-assay active/inactive "
            "contexts and cross-source quantitative heterogeneity are audited "
            "without averaging"
        ),
        "outputs": {
            "binding_evidence": str(binding_path),
            "negative_evidence": str(negative_path),
            "review_queue": str(candidate_path),
            "binding_sites": str(site_path),
            "protein_compound_summary": str(pair_path),
            "pubchem_context_conflicts": str(context_conflict_path),
            "quantitative_heterogeneity": str(heterogeneity_path),
            "source_coverage": str(coverage_path),
            "merge_database": str(DB_PATH),
        },
        "blocking_findings": blocking,
    }
    qa_path = QA_DIR / "INCREMENTAL_EVIDENCE_BUILD_V6_QA.json"
    qa_tmp = qa_path.with_suffix(".json.tmp")
    qa_tmp.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    qa_tmp.replace(qa_path)
    con.close()
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not blocking else 1


if __name__ == "__main__":
    sys.exit(main())
