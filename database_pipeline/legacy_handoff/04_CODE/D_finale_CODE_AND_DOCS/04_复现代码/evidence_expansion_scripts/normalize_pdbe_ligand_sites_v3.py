#!/usr/bin/env python3
"""Normalize PDBe UniProt ligand-site responses for incremental supplementation.

Raw PDBe files and the frozen releases are read-only.  Source-specific gzip TSV
files are written atomically under the isolated run staging directory.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import pathlib
import re
import sys
from collections import Counter, defaultdict


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "runs" / "incremental_20260727"
OUT_DIR = RUN_ROOT / "staging" / "pdbe"
RAW_DIR = ROOT / "raw" / "pdbe_202607" / "ligand_sites"
PROTEIN_UNIVERSE = ROOT / "intermediate" / "protein_universe_abc_v2.tsv"
BASELINE_PATHS = ROOT / "config" / "baseline_paths.json"

# High-confidence crystallization/solvent artifacts.  They are preserved in an
# audit table but are not candidates for the default binding release.
SOLVENT_BUFFER_HET_IDS = {
    "HOH",
    "DOD",
    "EDO",
    "GOL",
    "PEG",
    "PG4",
    "PGE",
    "1PE",
    "MPD",
    "TRS",
    "MES",
    "HEP",
    "EOH",
    "IPA",
    "DMS",
    "ACY",
}
COMMON_ION_HET_IDS = {
    "NA",
    "K",
    "CL",
    "BR",
    "IOD",
    "SO4",
    "PO4",
    "NO3",
    "NH4",
    "MG",
    "CA",
}

EXTRA_FIELDS = [
    "ligand_het_id",
    "ligand_atom_count",
    "artifact_class",
    "pdb_entity_id",
    "pdb_chain_ids",
    "residue_index_type",
    "residue_qc_mismatch_count",
]


def stable_id(prefix: str, *parts: object) -> str:
    material = "\x1f".join("" if part is None else str(part) for part in parts)
    return prefix + hashlib.sha256(material.encode("utf-8")).hexdigest()[:24].upper()


def atomic_gzip_writer(path: pathlib.Path, fields: list[str]):
    temporary = path.with_suffix(path.suffix + ".tmp")
    handle = gzip.open(temporary, "wt", encoding="utf-8", newline="")
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fields,
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    return temporary, handle, writer


def load_schema() -> list[str]:
    with (ROOT / "config" / "normalized_source_schema.tsv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        return [row["column_name"] for row in csv.DictReader(handle, delimiter="\t")]


def load_proteins() -> dict[str, dict[str, str]]:
    result = {}
    with PROTEIN_UNIVERSE.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["target_uniprot_id"]] = row
    return result


def load_het_crosswalk() -> tuple[dict[str, dict[str, str]], set[str]]:
    """Use only unambiguous HET mappings already present in frozen V4.2."""
    paths = json.loads(BASELINE_PATHS.read_text(encoding="utf-8"))
    evidence_path = pathlib.Path(paths["binding_evidence_v4_2"])
    candidates: dict[str, dict[tuple[str, str], dict[str, str]]] = defaultdict(dict)
    with evidence_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            het = (row.get("ligand_het_id") or "").strip().upper()
            internal = (row.get("compound_internal_id") or "").strip()
            if not het or not internal:
                continue
            form = (row.get("compound_form_id") or "").strip()
            candidates[het][(internal, form)] = {
                "compound_internal_id": internal,
                "compound_form_id": form,
                "compound_mapping_status": (
                    row.get("compound_mapping_status") or "mapped_core"
                ),
                "compound_scope_status": (
                    row.get("compound_scope_status_v1") or ""
                ),
            }
    unambiguous: dict[str, dict[str, str]] = {}
    ambiguous: set[str] = set()
    for het, mappings in candidates.items():
        internal_ids = {key[0] for key in mappings}
        if len(internal_ids) == 1:
            # Prefer an exact form if present; otherwise retain parent mapping.
            chosen = sorted(
                mappings.values(),
                key=lambda row: (not bool(row["compound_form_id"]), row["compound_form_id"]),
            )[0]
            unambiguous[het] = chosen
        else:
            ambiguous.add(het)
    return unambiguous, ambiguous


def classify_artifact(
    het: str, additional: dict[str, object], atom_count: float | None
) -> tuple[str, str]:
    if bool(additional.get("isSolvent")):
        return "solvent_pdbe", "PDBe isSolvent=true"
    if het in SOLVENT_BUFFER_HET_IDS:
        return "solvent_or_buffer_curated", "curated solvent/buffer component"
    if het in COMMON_ION_HET_IDS:
        return "common_ion_nondefault", "common ion/buffer component"
    if atom_count is not None and atom_count < 3:
        return "very_small_component_nondefault", "fewer than three atoms"
    return "", ""


def split_chain_ids(value: object) -> list[str]:
    if value is None:
        return [""]
    values = [part.strip() for part in re.split(r"[,;]", str(value)) if part.strip()]
    return values or [""]


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    schema = load_schema()
    fields = schema + EXTRA_FIELDS
    proteins = load_proteins()
    het_map, ambiguous_hets = load_het_crosswalk()

    evidence_path = OUT_DIR / "pdbe_normalized_evidence_v3.tsv.gz"
    artifact_path = OUT_DIR / "pdbe_excluded_artifacts_v3.tsv.gz"
    ev_tmp, ev_handle, ev_writer = atomic_gzip_writer(evidence_path, fields)
    ar_tmp, ar_handle, ar_writer = atomic_gzip_writer(artifact_path, fields)

    counts: Counter[str] = Counter()
    ligands: set[str] = set()
    pdb_ids: set[str] = set()
    target_ids: set[str] = set()
    mapped_compounds: set[str] = set()
    target_ligand_pairs: set[tuple[str, str]] = set()
    errors: list[dict[str, str]] = []

    try:
        for path in sorted(RAW_DIR.glob("*.json.gz")):
            target = path.name.removesuffix(".json.gz")
            if target not in proteins:
                counts["files_target_not_in_abc_universe"] += 1
                continue
            try:
                with gzip.open(path, "rt", encoding="utf-8") as handle:
                    payload = json.load(handle)
            except Exception as exc:  # retain a bounded audit of malformed source files
                counts["files_parse_error"] += 1
                if len(errors) < 100:
                    errors.append({"path": str(path), "error": repr(exc)})
                continue

            record = payload.get(target)
            if not isinstance(record, dict):
                counts["files_missing_target_record"] += 1
                continue
            sequence = str(record.get("sequence") or "")
            target_row = proteins[target]
            data = record.get("data") or []
            counts["files_parsed"] += 1
            target_ids.add(target)

            for ligand in data:
                if not isinstance(ligand, dict):
                    counts["malformed_ligand_records"] += 1
                    continue
                het = str(ligand.get("accession") or "").strip().upper()
                name = str(ligand.get("name") or "").strip()
                additional = ligand.get("additionalData") or {}
                try:
                    atom_count = (
                        float(additional.get("numAtoms"))
                        if additional.get("numAtoms") is not None
                        else None
                    )
                except (TypeError, ValueError):
                    atom_count = None
                artifact_class, artifact_reason = classify_artifact(
                    het, additional, atom_count
                )
                mapping = het_map.get(het, {})
                compound_internal_id = mapping.get("compound_internal_id", "")
                compound_form_id = mapping.get("compound_form_id", "")
                if het in ambiguous_hets:
                    mapping_status = "ambiguous_existing_het_mapping"
                elif compound_internal_id:
                    mapping_status = mapping.get("compound_mapping_status", "mapped_core")
                else:
                    mapping_status = "unresolved"

                # Aggregate contacts separately for each PDB/entity/chain.
                grouped: dict[
                    tuple[str, str, str], list[tuple[int, str, str]]
                ] = defaultdict(list)
                for residue in ligand.get("residues") or []:
                    if not isinstance(residue, dict):
                        continue
                    try:
                        position = int(residue.get("startIndex"))
                    except (TypeError, ValueError):
                        counts["residues_without_numeric_position"] += 1
                        continue
                    amino_acid = str(residue.get("startCode") or "").upper()
                    index_type = str(residue.get("indexType") or "")
                    interactions = residue.get("interactingPDBEntries") or []
                    for interaction in interactions:
                        pdb_id = str(interaction.get("pdbId") or "").lower()
                        entity = str(interaction.get("entityId") or "")
                        for chain in split_chain_ids(interaction.get("chainIds")):
                            grouped[(pdb_id, entity, chain)].append(
                                (position, amino_acid, index_type)
                            )

                for (pdb_id, entity, chain), residues in sorted(grouped.items()):
                    unique_residues = sorted(set(residues), key=lambda x: (x[0], x[1], x[2]))
                    residue_text = ";".join(
                        f"{aa}{position}" if aa else str(position)
                        for position, aa, _ in unique_residues
                    )
                    index_types = ";".join(sorted({item[2] for item in unique_residues}))
                    mismatch_count = 0
                    if index_types == "UNIPROT" and sequence:
                        for position, amino_acid, _ in unique_residues:
                            if (
                                position < 1
                                or position > len(sequence)
                                or (
                                    amino_acid
                                    and sequence[position - 1].upper()
                                    != {
                                        "ALA": "A",
                                        "ARG": "R",
                                        "ASN": "N",
                                        "ASP": "D",
                                        "CYS": "C",
                                        "GLN": "Q",
                                        "GLU": "E",
                                        "GLY": "G",
                                        "HIS": "H",
                                        "ILE": "I",
                                        "LEU": "L",
                                        "LYS": "K",
                                        "MET": "M",
                                        "PHE": "F",
                                        "PRO": "P",
                                        "SER": "S",
                                        "THR": "T",
                                        "TRP": "W",
                                        "TYR": "Y",
                                        "VAL": "V",
                                    }.get(amino_acid, sequence[position - 1].upper())
                                )
                            ):
                                mismatch_count += 1

                    exclusion_reasons = []
                    if artifact_reason:
                        exclusion_reasons.append(artifact_reason)
                    if not compound_internal_id:
                        exclusion_reasons.append("compound not mapped to frozen V1.0")
                    elif mapping_status not in {"mapped_core", "mapped_extended"}:
                        exclusion_reasons.append(
                            "frozen compound mapping is not release-eligible"
                        )
                    if mismatch_count:
                        exclusion_reasons.append("UniProt residue code mismatch")
                    if not residue_text:
                        exclusion_reasons.append("no residue-level contacts")
                    default = int(not exclusion_reasons)
                    qc_status = "ok" if not mismatch_count else "review"
                    source_record_id = f"{target}:{pdb_id}:{het}:{entity}:{chain}"
                    evidence_id = stable_id("PDBe-", source_record_id, residue_text)
                    row = {
                        "source_evidence_id": evidence_id,
                        "source_database": "PDBe",
                        "source_version": "current_2026-07",
                        "source_record_id": source_record_id,
                        "source_url": f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}",
                        "retrieval_date": "2026-07-27",
                        "target_uniprot_id": target,
                        "target_mapping_status": "single_human_uniprot",
                        "approved_symbol": target_row.get("approved_symbol", ""),
                        "compound_source_id": f"PDBCCD:{het}",
                        "compound_name": name,
                        "compound_internal_id": compound_internal_id,
                        "compound_form_id": compound_form_id,
                        "compound_mapping_status": mapping_status,
                        "evidence_tier": "BE1",
                        "evidence_type": "experimental_structure_residue_contact",
                        "activity_type": "structure_contact",
                        "activity_relation": "",
                        "activity_value": "",
                        "activity_unit": "",
                        "standard_value_nM": "",
                        "activity_outcome": "observed",
                        "assay_or_mechanism": "PDBe ligand binding site",
                        "pdb_ids": pdb_id,
                        "binding_site_residues": residue_text,
                        "pubmed_ids": "",
                        "doi": "",
                        "record_qc_status": qc_status,
                        "default_release_inclusion": str(default),
                        "exclusion_reason": "; ".join(exclusion_reasons),
                        "ligand_het_id": het,
                        "ligand_atom_count": "" if atom_count is None else atom_count,
                        "artifact_class": artifact_class,
                        "pdb_entity_id": entity,
                        "pdb_chain_ids": chain,
                        "residue_index_type": index_types,
                        "residue_qc_mismatch_count": mismatch_count,
                    }
                    ev_writer.writerow(row)
                    counts["evidence_rows"] += 1
                    counts[f"mapping_{mapping_status}"] += 1
                    counts[f"artifact_{artifact_class or 'none'}"] += 1
                    counts[f"default_{default}"] += 1
                    counts["residue_contacts"] += len(unique_residues)
                    if artifact_class:
                        ar_writer.writerow(row)
                        counts["artifact_audit_rows"] += 1
                    ligands.add(het)
                    pdb_ids.add(pdb_id)
                    if compound_internal_id:
                        mapped_compounds.add(compound_internal_id)
                        target_ligand_pairs.add((target, compound_internal_id))
    finally:
        ev_handle.close()
        ar_handle.close()

    ev_tmp.replace(evidence_path)
    ar_tmp.replace(artifact_path)

    report = {
        "status": "passed" if not errors else "passed_with_parse_errors",
        "run_type": "incremental_supplement_not_rebuild",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "raw_directory": str(RAW_DIR),
        "baseline_policy": "read_only",
        "counts": dict(sorted(counts.items())),
        "unique_targets": len(target_ids),
        "unique_ligand_het_ids": len(ligands),
        "unique_pdb_ids": len(pdb_ids),
        "unique_mapped_compounds": len(mapped_compounds),
        "unique_mapped_target_compound_pairs": len(target_ligand_pairs),
        "unambiguous_frozen_het_crosswalk_size": len(het_map),
        "ambiguous_frozen_het_ids": len(ambiguous_hets),
        "parse_errors_sample": errors,
        "outputs": {
            "normalized_evidence": str(evidence_path),
            "excluded_artifacts": str(artifact_path),
        },
    }
    report_path = RUN_ROOT / "qa" / "PDBE_NORMALIZATION_V3_QA.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    temporary.replace(report_path)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
