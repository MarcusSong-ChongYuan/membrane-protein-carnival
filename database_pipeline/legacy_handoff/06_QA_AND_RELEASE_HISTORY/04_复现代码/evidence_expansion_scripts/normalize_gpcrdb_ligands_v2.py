#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import html
import json
import math
import re
import sqlite3
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RAW_DIR = ROOT / "raw" / "gpcrdb_202607" / "ligands"
CHECKPOINT = ROOT / "raw" / "gpcrdb_202607" / "gpcrdb_ligand_fetch_v2.jsonl"
INDEX = ROOT / "intermediate" / "integration_index_v2.sqlite"
COMPOUND_MASTER = Path(
    r"D:\7.22\small_molecule_v1_working\releases\release_small_molecule_v1_0"
    r"\small_molecule_master_v1_0.tsv"
)
COMPOUND_FORMS = Path(
    r"D:\7.22\small_molecule_v1_working\releases\release_small_molecule_v1_0"
    r"\compound_form_hierarchy_v1_0.tsv"
)
BASELINE_EVIDENCE = Path(
    r"D:\7.22\small_molecule_v1_working\releases\release_small_molecule_v1_0"
    r"\binding_evidence_master_v4_2.tsv"
)
OUTPUT = ROOT / "intermediate" / "gpcrdb_normalized_evidence_v2.tsv"
DEFAULT_CANDIDATES = ROOT / "intermediate" / "gpcrdb_default_new_evidence_v2.tsv"
UNRESOLVED = ROOT / "intermediate" / "gpcrdb_unresolved_compounds_v2.tsv"
REPORT = ROOT / "reports" / "GPCRDB_NORMALIZATION_V2_REPORT.json"
LOG = ROOT / "logs" / "gpcrdb_normalization_v2.log"

SOURCE_CANONICAL = {
    "ChEMBL": "ChEMBL",
    "Drug Central": "DrugCentral",
    "Guide to Pharmacology": "IUPHAR/BPS Guide to PHARMACOLOGY",
    "PDSP KiDatabase": "PDSP KiDatabase",
}
DIRECT_TYPES = {"Ki", "Kd"}
FUNCTIONAL_TYPES = {
    "IC50",
    "EC50",
    "AC50",
    "ED50",
    "Potency",
    "Kb",
    "KB",
    "A2",
}
PUBMED_RE = re.compile(r"(?:pubmed(?:\.ncbi\.nlm\.nih\.gov)?/|pubmed/)(\d+)", re.I)
TAG_RE = re.compile(r"<[^>]+>")

NORMALIZED_FIELDS = [
    "source_evidence_id",
    "source_database",
    "source_version",
    "source_record_id",
    "source_url",
    "retrieval_date",
    "target_uniprot_id",
    "target_mapping_status",
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
    "pubmed_ids",
    "doi",
    "record_qc_status",
    "default_release_inclusion",
    "exclusion_reason",
    "gpcrdb_entry_name",
    "gpcrdb_protein_name",
    "source_ligand_type",
    "source_smiles",
    "exact_standard_smiles",
    "exact_inchikey",
    "parent_standard_smiles",
    "parent_inchikey",
    "source_p_activity_type",
    "source_activity_range_p",
    "source_assay_type",
    "existing_pair_overlap",
    "existing_same_source_relation",
]


def clean_text(value: object) -> str:
    text = "" if value is None else str(value)
    return TAG_RE.sub("", html.unescape(text)).strip()


def float_or_none(value: object) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def stable_hash(*parts: object, prefix: str) -> str:
    text = "\x1f".join("" if part is None else str(part) for part in parts)
    return prefix + hashlib.sha256(text.encode("utf-8")).hexdigest()[:20].upper()


def log(message: str) -> None:
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def standardize_smiles(smiles: str) -> dict[str, str]:
    result = {
        "parse_status": "no_structure",
        "exact_standard_smiles": "",
        "exact_inchikey": "",
        "parent_standard_smiles": "",
        "parent_inchikey": "",
    }
    if not smiles:
        return result
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result["parse_status"] = "parse_failed"
        return result
    try:
        exact = rdMolStandardize.Cleanup(mol)
    except Exception:
        exact = mol
    try:
        parent = rdMolStandardize.FragmentParent(exact)
        if parent is None or parent.GetNumAtoms() == 0:
            parent = exact
    except Exception:
        parent = exact
    try:
        parent = rdMolStandardize.Uncharger().uncharge(parent)
    except Exception:
        pass
    try:
        exact_smiles = Chem.MolToSmiles(exact, canonical=True, isomericSmiles=True)
        parent_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        exact_inchi = Chem.MolToInchi(exact)
        parent_inchi = Chem.MolToInchi(parent)
        exact_key = Chem.InchiToInchiKey(exact_inchi) if exact_inchi else ""
        parent_key = Chem.InchiToInchiKey(parent_inchi) if parent_inchi else ""
    except Exception:
        result["parse_status"] = "standardization_failed"
        return result
    result.update(
        {
            "parse_status": "parsed",
            "exact_standard_smiles": exact_smiles,
            "exact_inchikey": exact_key,
            "parent_standard_smiles": parent_smiles,
            "parent_inchikey": parent_key,
        }
    )
    return result


def unique_map(mapping: dict[str, set[tuple]]) -> dict[str, tuple]:
    return {key: next(iter(values)) for key, values in mapping.items() if len(values) == 1}


def load_compound_maps() -> tuple[dict[str, tuple], dict[str, tuple]]:
    parent_sets: dict[str, set[tuple]] = defaultdict(set)
    with COMPOUND_MASTER.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("standard_inchikey", "").strip()
            if key:
                parent_sets[key].add(
                    (
                        row["compound_internal_id"].strip(),
                        row.get("compound_scope_status", "").strip(),
                    )
                )
    form_sets: dict[str, set[tuple]] = defaultdict(set)
    with COMPOUND_FORMS.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("exact_inchikey", "").strip()
            if key:
                form_sets[key].add(
                    (
                        row["compound_internal_id"].strip(),
                        row["compound_form_id"].strip(),
                    )
                )
    return unique_map(parent_sets), unique_map(form_sets)


def load_index() -> tuple[dict[str, str], set[tuple[str, str]]]:
    connection = sqlite3.connect(INDEX)
    symbols = dict(
        connection.execute("SELECT target_uniprot_id, approved_symbol FROM protein")
    )
    pairs = set(
        connection.execute(
            "SELECT target_uniprot_id, compound_internal_id FROM existing_pair"
        )
    )
    connection.close()
    return symbols, pairs


def load_existing_source_relations() -> set[tuple[str, str, str]]:
    relations: set[tuple[str, str, str]] = set()
    with BASELINE_EVIDENCE.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target = row.get("target_uniprot_id", "").strip()
            compound = row.get("compound_internal_id", "").strip()
            for source in row.get("source_database", "").split(";"):
                source = source.strip()
                if target and compound and source:
                    relations.add((target, compound, source))
    return relations


def load_target_map() -> dict[str, dict]:
    targets = {}
    with CHECKPOINT.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if record.get("status") == "ok":
                targets[record["entry_name"]] = record
    return targets


def normalize_activity(row: dict) -> tuple[str, float | None, float | None]:
    raw_type = clean_text(row.get("Value type"))
    activity_type = raw_type[1:] if raw_type.startswith("p") else raw_type
    p_value = float_or_none(row.get("Activity value (P)"))
    standard = float_or_none(row.get("Activity value (Standard)"))
    if standard is None or standard <= 0:
        standard = 10 ** (9 - p_value) if p_value is not None else None
    return activity_type, p_value, standard


def extract_reference(value: object) -> tuple[str, str]:
    text = clean_text(value)
    if not text or text.lower() == "not available":
        return "", ""
    pubmed = PUBMED_RE.search(text)
    if pubmed:
        return pubmed.group(1), ""
    if "doi.org/" in text.lower():
        return "", text.split("doi.org/", 1)[-1].strip()
    if text.lower().startswith("10."):
        return "", text
    return "", ""


def main() -> None:
    parent_map, form_map = load_compound_maps()
    symbols, existing_pairs = load_index()
    existing_sources = load_existing_source_relations()
    target_map = load_target_map()
    structure_cache: dict[str, dict[str, str]] = {}
    seen_evidence: set[str] = set()
    counters = Counter()
    source_counts = Counter()
    default_source_counts = Counter()
    mapping_counts = Counter()
    unresolved_summary: dict[str, dict] = {}

    with (
        OUTPUT.open("w", encoding="utf-8", newline="") as output_handle,
        DEFAULT_CANDIDATES.open("w", encoding="utf-8", newline="") as default_handle,
    ):
        writer = csv.DictWriter(
            output_handle,
            fieldnames=NORMALIZED_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        default_writer = csv.DictWriter(
            default_handle,
            fieldnames=NORMALIZED_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        default_writer.writeheader()

        for entry_name, target_record in sorted(target_map.items()):
            raw_path = RAW_DIR / f"{entry_name}.json.gz"
            if not raw_path.exists():
                counters["missing_raw_file"] += 1
                continue
            with gzip.open(raw_path, "rt", encoding="utf-8") as handle:
                rows = json.load(handle)
            target = target_record["target_uniprot_id"]
            for row in rows:
                counters["raw_rows"] += 1
                raw_source = clean_text(row.get("Source"))
                source = SOURCE_CANONICAL.get(raw_source, raw_source)
                ligand_name = clean_text(row.get("Ligand name"))
                ligand_type = clean_text(row.get("Ligand type")).lower()
                smiles = clean_text(row.get("Smiles"))
                activity_type, p_value, standard_nm = normalize_activity(row)
                pubmed_ids, doi = extract_reference(row.get("DOI"))
                assay_description = clean_text(row.get("Assay description"))
                assay_type = clean_text(row.get("Assay type"))
                activity_range = clean_text(row.get("Activity ranges (P)"))

                fingerprint = stable_hash(
                    source,
                    target,
                    ligand_name,
                    smiles,
                    clean_text(row.get("Value type")),
                    p_value,
                    standard_nm,
                    activity_range,
                    assay_type,
                    assay_description,
                    pubmed_ids,
                    doi,
                    prefix="GPEV-",
                )
                if fingerprint in seen_evidence:
                    counters["exact_duplicates_dropped"] += 1
                    continue
                seen_evidence.add(fingerprint)

                skip_structure = ligand_type in {"protein", "peptide"} or len(smiles) > 2000
                if skip_structure:
                    structure = {
                        "parse_status": "excluded_non_small_molecule_structure",
                        "exact_standard_smiles": "",
                        "exact_inchikey": "",
                        "parent_standard_smiles": "",
                        "parent_inchikey": "",
                    }
                else:
                    if smiles not in structure_cache:
                        structure_cache[smiles] = standardize_smiles(smiles)
                    structure = structure_cache[smiles]
                compound_id = ""
                compound_form_id = ""
                compound_scope = ""
                mapping_status = "unresolved"
                if structure["exact_inchikey"] in form_map:
                    compound_id, compound_form_id = form_map[structure["exact_inchikey"]]
                    parent_hit = parent_map.get(structure["parent_inchikey"])
                    if parent_hit and parent_hit[0] == compound_id:
                        compound_scope = parent_hit[1]
                    mapping_status = (
                        "mapped_core" if compound_scope == "core" else "mapped_extended"
                    )
                elif structure["parent_inchikey"] in parent_map:
                    compound_id, compound_scope = parent_map[structure["parent_inchikey"]]
                    mapping_status = (
                        "mapped_core" if compound_scope == "core" else "mapped_extended"
                    )

                evidence_tier = (
                    "BE2"
                    if activity_type in DIRECT_TYPES
                    else "BE3"
                    if activity_type in FUNCTIONAL_TYPES
                    else "nondefault"
                )
                evidence_type = (
                    "direct_quantitative_binding"
                    if evidence_tier == "BE2"
                    else "target_specific_functional_pharmacology"
                    if evidence_tier == "BE3"
                    else "unclassified_gpcr_activity"
                )
                existing_pair = int((target, compound_id) in existing_pairs)
                same_source = int(
                    bool(compound_id)
                    and (target, compound_id, source) in existing_sources
                )
                exclusions = []
                is_small_molecule = ligand_type == "small molecule" or (
                    ligand_type in {"", "na"} and mapping_status == "mapped_core"
                )
                if not is_small_molecule:
                    exclusions.append("not_confirmed_small_molecule")
                if mapping_status != "mapped_core":
                    exclusions.append(mapping_status)
                if evidence_tier not in {"BE2", "BE3"}:
                    exclusions.append("unclassified_activity_type")
                if p_value is None or standard_nm is None:
                    exclusions.append("missing_quantitative_value")
                if source == "ChEMBL":
                    exclusions.append("chembl_proxy_duplicate_risk")
                if same_source:
                    exclusions.append("existing_same_source_relation")
                default_inclusion = int(not exclusions)
                record_qc = (
                    "excluded"
                    if "not_confirmed_small_molecule" in exclusions
                    else "review"
                    if structure["parse_status"] != "parsed"
                    or mapping_status == "unresolved"
                    else "ok"
                )
                source_record_id = stable_hash(
                    source,
                    target,
                    ligand_name,
                    smiles,
                    activity_type,
                    p_value,
                    standard_nm,
                    assay_description,
                    prefix="GPCRDB-",
                )
                normalized = {
                    "source_evidence_id": fingerprint,
                    "source_database": source,
                    "source_version": "GPCRdb retrieval 2026-07-27",
                    "source_record_id": source_record_id,
                    "source_url": f"https://gpcrdb.org/services/ligands/{entry_name}/",
                    "retrieval_date": date.today().isoformat(),
                    "target_uniprot_id": target,
                    "target_mapping_status": "single_human_uniprot",
                    "approved_symbol": symbols.get(target, ""),
                    "compound_source_id": stable_hash(
                        source, smiles, prefix="GPCRLIG-"
                    ),
                    "compound_name": ligand_name,
                    "compound_internal_id": compound_id,
                    "compound_form_id": compound_form_id,
                    "compound_mapping_status": mapping_status,
                    "evidence_tier": evidence_tier,
                    "evidence_type": evidence_type,
                    "activity_type": activity_type,
                    "activity_relation": "=",
                    "activity_value": "" if p_value is None else p_value,
                    "activity_unit": "pActivity",
                    "standard_value_nM": "" if standard_nm is None else standard_nm,
                    "activity_outcome": "active",
                    "assay_or_mechanism": assay_description,
                    "pdb_ids": "",
                    "binding_site_residues": "",
                    "pubmed_ids": pubmed_ids,
                    "doi": doi,
                    "record_qc_status": record_qc,
                    "default_release_inclusion": default_inclusion,
                    "exclusion_reason": ";".join(exclusions),
                    "gpcrdb_entry_name": entry_name,
                    "gpcrdb_protein_name": clean_text(row.get("Protein name")),
                    "source_ligand_type": ligand_type,
                    "source_smiles": smiles,
                    "exact_standard_smiles": structure["exact_standard_smiles"],
                    "exact_inchikey": structure["exact_inchikey"],
                    "parent_standard_smiles": structure["parent_standard_smiles"],
                    "parent_inchikey": structure["parent_inchikey"],
                    "source_p_activity_type": clean_text(row.get("Value type")),
                    "source_activity_range_p": activity_range,
                    "source_assay_type": assay_type,
                    "existing_pair_overlap": existing_pair,
                    "existing_same_source_relation": same_source,
                }
                writer.writerow(normalized)
                if default_inclusion:
                    default_writer.writerow(normalized)

                counters["normalized_unique_rows"] += 1
                counters["default_candidate_rows"] += default_inclusion
                counters["existing_pair_overlap_rows"] += existing_pair
                counters["existing_same_source_rows"] += same_source
                source_counts[source] += 1
                default_source_counts[source] += default_inclusion
                mapping_counts[mapping_status] += 1

                if mapping_status == "unresolved" and structure["parse_status"] == "parsed":
                    key = structure["parent_inchikey"] or structure["parent_standard_smiles"]
                    summary = unresolved_summary.setdefault(
                        key,
                        {
                            "parent_inchikey": structure["parent_inchikey"],
                            "parent_standard_smiles": structure[
                                "parent_standard_smiles"
                            ],
                            "example_name": ligand_name,
                            "source_databases": set(),
                            "target_uniprot_ids": set(),
                            "evidence_count": 0,
                            "best_evidence_tier": evidence_tier,
                        },
                    )
                    summary["source_databases"].add(source)
                    summary["target_uniprot_ids"].add(target)
                    summary["evidence_count"] += 1
                    if evidence_tier == "BE2":
                        summary["best_evidence_tier"] = "BE2"

                if counters["raw_rows"] % 5000 == 0:
                    output_handle.flush()
                    default_handle.flush()
                    log(
                        f"raw_rows={counters['raw_rows']} "
                        f"normalized={counters['normalized_unique_rows']} "
                        f"default={counters['default_candidate_rows']} "
                        f"structures={len(structure_cache)}"
                    )

    unresolved_fields = [
        "parent_inchikey",
        "parent_standard_smiles",
        "example_name",
        "source_databases",
        "target_count",
        "target_uniprot_ids",
        "evidence_count",
        "best_evidence_tier",
    ]
    with UNRESOLVED.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=unresolved_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for _, summary in sorted(unresolved_summary.items()):
            row = dict(summary)
            row["source_databases"] = ";".join(sorted(summary["source_databases"]))
            row["target_count"] = len(summary["target_uniprot_ids"])
            row["target_uniprot_ids"] = ";".join(sorted(summary["target_uniprot_ids"]))
            writer.writerow(row)

    report = {
        "status": "passed",
        "retrieval_date": date.today().isoformat(),
        "input_target_files": len(target_map),
        "counts": dict(counters),
        "source_normalized_rows": dict(source_counts),
        "source_default_candidate_rows": dict(default_source_counts),
        "compound_mapping_counts": dict(mapping_counts),
        "unique_structures_processed": len(structure_cache),
        "unresolved_parent_structures": len(unresolved_summary),
        "output_files": {
            "all_normalized": str(OUTPUT),
            "default_new_candidates": str(DEFAULT_CANDIDATES),
            "unresolved_compounds": str(UNRESOLVED),
        },
        "policy_notes": [
            "Names alone never trigger compound merges.",
            "Only unique exact-form or parent InChIKey matches map automatically.",
            "ChEMBL proxy rows are retained for audit but excluded from default additions.",
            "Existing target-compound-source relations are retained for audit but excluded from default additions.",
        ],
    }
    with REPORT.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
