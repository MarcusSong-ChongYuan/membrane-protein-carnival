#!/usr/bin/env python3
"""Normalize the selected PDBbind membrane subset and make a docking manifest."""

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
OUT_DIR = RUN_ROOT / "staging" / "pdbbind"
MANIFEST = ROOT / "intermediate" / "pdbbind_2020r1_membrane_subset_manifest.tsv"
BASELINE_PATHS = ROOT / "config" / "baseline_paths.json"
SCHEMA_PATH = ROOT / "config" / "normalized_source_schema.tsv"

EXTRA_FIELDS = [
    "ligand_het_id",
    "pdbbind_resolution",
    "pdbbind_release_year",
    "pdbbind_ligand_annotation",
    "pdbbind_index_comment",
    "pocket_residue_count",
    "pocket_residue_index_type",
    "ligand_sdf_sha256",
    "ligand_atom_count",
    "ligand_heavy_atom_count",
    "covalent_complex_flag",
    "incomplete_ligand_flag",
    "peptide_or_oligomer_flag",
    "compound_mapping_method",
]

AFFINITY_RE = re.compile(
    r"^(?P<type>Kd|Ki|IC50)(?P<relation><=|>=|=|<|>|~)"
    r"(?P<value>[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?\d+)?)"
    r"(?P<unit>pM|nM|uM|µM|mM)$",
    re.I,
)
UNIT_TO_NM = {"PM": 0.001, "NM": 1.0, "UM": 1000.0, "ΜM": 1000.0, "MM": 1_000_000.0}
PEPTIDE_RE = re.compile(r"^\d+-mer$", re.I)
HET_RE = re.compile(r"^[A-Za-z0-9]{1,5}$")


def stable_id(prefix: str, *parts: object) -> str:
    material = "\x1f".join("" if part is None else str(part) for part in parts)
    return prefix + hashlib.sha256(material.encode("utf-8")).hexdigest()[:24].upper()


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_schema() -> list[str]:
    with SCHEMA_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        return [row["column_name"] for row in csv.DictReader(handle, delimiter="\t")]


def split_values(value: str) -> list[str]:
    return [
        part.strip()
        for part in re.split(r"[;,|]", value or "")
        if part.strip()
    ]


def load_frozen_structure_crosswalks():
    """Build unambiguous PDB/target and HET mappings from frozen V4.2."""
    paths = json.loads(BASELINE_PATHS.read_text(encoding="utf-8"))
    evidence_path = pathlib.Path(paths["binding_evidence_v4_2"])
    by_pdb_target: dict[tuple[str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    by_pdb: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    by_het: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    with evidence_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            compound = (row.get("compound_internal_id") or "").strip()
            if not compound:
                continue
            form = (row.get("compound_form_id") or "").strip()
            key = compound + "\x1f" + form
            mapping = {
                "compound_internal_id": compound,
                "compound_form_id": form,
                "compound_mapping_status": (
                    row.get("compound_mapping_status") or "mapped_core"
                ),
                "compound_scope_status": (
                    row.get("compound_scope_status_v1") or ""
                ),
                "compound_name": (row.get("compound_name") or "").strip(),
            }
            target = (row.get("target_uniprot_id") or "").strip()
            for pdb in split_values(row.get("pdb_ids") or ""):
                pdb = pdb.lower()
                by_pdb_target[(pdb, target)][key] = mapping
                by_pdb[pdb][key] = mapping
            het = (row.get("ligand_het_id") or "").strip().upper()
            if het:
                by_het[het][key] = mapping

    def unique(table):
        result = {}
        ambiguous = set()
        for key, mappings in table.items():
            internal_ids = {row["compound_internal_id"] for row in mappings.values()}
            if len(internal_ids) == 1:
                result[key] = sorted(
                    mappings.values(),
                    key=lambda row: (
                        not bool(row["compound_form_id"]),
                        row["compound_form_id"],
                    ),
                )[0]
            else:
                ambiguous.add(key)
        return result, ambiguous

    return (*unique(by_pdb_target), *unique(by_pdb), *unique(by_het))


def parse_affinity(raw: str):
    match = AFFINITY_RE.match((raw or "").strip())
    if not match:
        return None
    unit = match.group("unit").upper().replace("µ", "U")
    value = float(match.group("value"))
    return {
        "activity_type": match.group("type"),
        "activity_relation": match.group("relation"),
        "activity_value": match.group("value"),
        "activity_unit": match.group("unit"),
        "standard_value_nM": value * UNIT_TO_NM[unit],
    }


def ligand_annotation_info(annotation: str, comment: str):
    annotation = (annotation or "").strip()
    peptide = bool(PEPTIDE_RE.fullmatch(annotation))
    het = (
        annotation.upper()
        if annotation and HET_RE.fullmatch(annotation) and not peptide
        else ""
    )
    name = annotation
    match = re.search(r"ligand is (.+?)(?:[.;]|$)", comment or "", flags=re.I)
    if match:
        name = match.group(1).strip()
    return het, name, peptide


def parse_sdf_counts(path: pathlib.Path) -> tuple[int, int]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = [handle.readline() for _ in range(4)]
        if len(lines) < 4:
            return 0, 0
        counts = lines[3]
        try:
            atom_count = int(counts[0:3])
        except ValueError:
            parts = counts.split()
            atom_count = int(parts[0]) if parts else 0
        heavy = 0
        for _ in range(atom_count):
            line = handle.readline()
            if not line:
                break
            element = line[31:34].strip()
            if not element:
                parts = line.split()
                element = parts[3] if len(parts) > 3 else ""
            if element.upper() not in {"H", "D", "T"}:
                heavy += 1
    return atom_count, heavy


def parse_pocket(path: pathlib.Path) -> list[tuple[str, str, str, str]]:
    residues: set[tuple[str, str, str, str]] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            residue = line[17:20].strip()
            chain = line[21:22].strip()
            number = line[22:26].strip()
            insertion = line[26:27].strip()
            residues.add((chain, residue, number, insertion))
    def key(item):
        chain, residue, number, insertion = item
        try:
            numeric = int(number)
        except ValueError:
            numeric = 10**9
        return chain, numeric, insertion, residue
    return sorted(residues, key=key)


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


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = load_schema() + EXTRA_FIELDS
    (
        pdb_target_map,
        ambiguous_pdb_target,
        pdb_map,
        ambiguous_pdb,
        het_map,
        ambiguous_het,
    ) = load_frozen_structure_crosswalks()

    evidence_path = OUT_DIR / "pdbbind_normalized_evidence_v3.tsv.gz"
    docking_path = OUT_DIR / "pdbbind_docking_structure_manifest_v3.tsv.gz"
    ev_tmp, ev_handle, ev_writer = atomic_gzip_writer(evidence_path, fields)
    docking_fields = [
        "pdb_id",
        "target_uniprot_ids",
        "approved_symbols",
        "compound_internal_id",
        "compound_form_id",
        "compound_mapping_status",
        "compound_mapping_method",
        "ligand_het_id",
        "ligand_name",
        "activity_type",
        "activity_relation",
        "activity_value",
        "activity_unit",
        "standard_value_nM",
        "resolution",
        "release_year",
        "protein_pdb_path",
        "pocket_pdb_path",
        "ligand_sdf_path",
        "ligand_mol2_path",
        "ligand_sdf_sha256",
        "ligand_atom_count",
        "ligand_heavy_atom_count",
        "pocket_residue_count",
        "pocket_residues_pdb_numbering",
        "covalent_complex_flag",
        "incomplete_ligand_flag",
        "peptide_or_oligomer_flag",
        "redocking_candidate",
        "redocking_exclusion_reason",
    ]
    dk_tmp, dk_handle, dk_writer = atomic_gzip_writer(docking_path, docking_fields)

    counts: Counter[str] = Counter()
    targets: set[str] = set()
    compounds: set[str] = set()
    target_compound_pairs: set[tuple[str, str]] = set()
    pdb_ids: set[str] = set()
    failures: list[dict[str, str]] = []

    try:
        with MANIFEST.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                counts["manifest_rows"] += 1
                pdb = row["pdb_id"].lower()
                directory = pathlib.Path(row["structure_directory"])
                protein_path = directory / f"{pdb}_protein.pdb"
                pocket_path = directory / f"{pdb}_pocket.pdb"
                sdf_path = directory / f"{pdb}_ligand.sdf"
                mol2_path = directory / f"{pdb}_ligand.mol2"
                required = [protein_path, pocket_path, sdf_path, mol2_path]
                if not all(path.exists() for path in required):
                    counts["missing_structure_file_rows"] += 1
                    if len(failures) < 100:
                        failures.append(
                            {
                                "pdb_id": pdb,
                                "error": "one or more required structure files missing",
                            }
                        )
                    continue
                affinity = parse_affinity(row["binding_data_raw"])
                if affinity is None:
                    counts["unparsed_affinity_rows"] += 1
                    continue

                targets_for_pdb = split_values(row["target_uniprot_ids"])
                symbols = split_values(row["approved_symbols"])
                comment = row["index_comment"]
                het, ligand_name, peptide = ligand_annotation_info(
                    row["ligand_annotation"], comment
                )
                incomplete = "incomplete ligand" in comment.lower()
                covalent = "covalent complex" in comment.lower()
                atom_count, heavy_count = parse_sdf_counts(sdf_path)
                ligand_hash = sha256(sdf_path)
                pocket_residues = parse_pocket(pocket_path)
                pocket_text = ";".join(
                    f"{chain}:{residue}{number}{insertion}"
                    for chain, residue, number, insertion in pocket_residues
                )

                # Prefer the frozen mapping for the same PDB and target.  Use a
                # PDB-level mapping only if unique across all frozen evidence,
                # then fall back to an unambiguous HET crosswalk.
                mapping = None
                mapping_method = ""
                target_specific_mappings = {
                    pdb_target_map[(pdb, target)]["compound_internal_id"]: (
                        pdb_target_map[(pdb, target)]
                    )
                    for target in targets_for_pdb
                    if (pdb, target) in pdb_target_map
                }
                if len(target_specific_mappings) == 1:
                    mapping = next(iter(target_specific_mappings.values()))
                    mapping_method = "frozen_same_pdb_target"
                elif pdb in pdb_map:
                    mapping = pdb_map[pdb]
                    mapping_method = "frozen_same_pdb"
                elif het and het in het_map:
                    mapping = het_map[het]
                    mapping_method = "frozen_unambiguous_het"
                elif len(target_specific_mappings) > 1:
                    mapping_method = "ambiguous_same_pdb_target"
                elif pdb in ambiguous_pdb:
                    mapping_method = "ambiguous_same_pdb"
                elif het and het in ambiguous_het:
                    mapping_method = "ambiguous_het"
                else:
                    mapping_method = "unresolved_structure_identity"

                compound_internal = (
                    mapping["compound_internal_id"] if mapping else ""
                )
                compound_form = mapping["compound_form_id"] if mapping else ""
                mapping_status = (
                    mapping["compound_mapping_status"] if mapping else "unresolved"
                )
                frozen_name = mapping["compound_name"] if mapping else ""
                if not ligand_name or PEPTIDE_RE.fullmatch(ligand_name):
                    ligand_name = frozen_name or ligand_name

                common_exclusions = []
                if len(targets_for_pdb) != 1:
                    common_exclusions.append(
                        "PDBbind complex does not resolve to one membrane protein"
                    )
                if not compound_internal:
                    common_exclusions.append("compound identity unresolved")
                elif mapping_status not in {"mapped_core", "mapped_extended"}:
                    common_exclusions.append(
                        "frozen compound mapping is not release-eligible"
                    )
                if peptide:
                    common_exclusions.append("peptide or oligomer ligand")
                if incomplete:
                    common_exclusions.append("incomplete ligand structure")

                for target in targets_for_pdb:
                    source_record_id = f"{pdb}:{target}"
                    evidence_id = stable_id(
                        "PDBBIND-",
                        source_record_id,
                        row["binding_data_raw"],
                        ligand_hash,
                    )
                    row_output = {
                        "source_evidence_id": evidence_id,
                        "source_database": "PDBbind",
                        "source_version": "2020.R1_reprocessed_2024",
                        "source_record_id": source_record_id,
                        "source_url": f"https://www.pdbbind-plus.org.cn/",
                        "retrieval_date": "2026-07-27",
                        "target_uniprot_id": target,
                        "target_mapping_status": (
                            "single_human_uniprot"
                            if len(targets_for_pdb) == 1
                            else "ambiguous_complex"
                        ),
                        "approved_symbol": ";".join(symbols),
                        "compound_source_id": (
                            f"PDBCCD:{het}" if het else f"PDBBIND:{pdb}"
                        ),
                        "compound_name": ligand_name,
                        "compound_internal_id": compound_internal,
                        "compound_form_id": compound_form,
                        "compound_mapping_status": mapping_status,
                        "evidence_tier": "BE1",
                        "evidence_type": "experimental_complex_with_affinity",
                        "activity_type": affinity["activity_type"],
                        "activity_relation": affinity["activity_relation"],
                        "activity_value": affinity["activity_value"],
                        "activity_unit": affinity["activity_unit"],
                        "standard_value_nM": affinity["standard_value_nM"],
                        "activity_outcome": "observed",
                        "assay_or_mechanism": comment,
                        "pdb_ids": pdb,
                        "binding_site_residues": pocket_text,
                        "pubmed_ids": "",
                        "doi": "",
                        "record_qc_status": (
                            "ok" if not common_exclusions else "review"
                        ),
                        "default_release_inclusion": str(
                            int(not common_exclusions)
                        ),
                        "exclusion_reason": "; ".join(common_exclusions),
                        "ligand_het_id": het,
                        "pdbbind_resolution": row["resolution"],
                        "pdbbind_release_year": row["release_year"],
                        "pdbbind_ligand_annotation": row["ligand_annotation"],
                        "pdbbind_index_comment": comment,
                        "pocket_residue_count": len(pocket_residues),
                        "pocket_residue_index_type": "PDB",
                        "ligand_sdf_sha256": ligand_hash,
                        "ligand_atom_count": atom_count,
                        "ligand_heavy_atom_count": heavy_count,
                        "covalent_complex_flag": int(covalent),
                        "incomplete_ligand_flag": int(incomplete),
                        "peptide_or_oligomer_flag": int(peptide),
                        "compound_mapping_method": mapping_method,
                    }
                    ev_writer.writerow(row_output)
                    counts["evidence_rows"] += 1
                    counts[
                        f"default_{row_output['default_release_inclusion']}"
                    ] += 1
                    counts[f"mapping_method_{mapping_method}"] += 1
                    counts[f"activity_type_{affinity['activity_type']}"] += 1
                    targets.add(target)
                    if compound_internal:
                        compounds.add(compound_internal)
                        target_compound_pairs.add((target, compound_internal))

                redocking_exclusions = list(common_exclusions)
                if covalent:
                    redocking_exclusions.append(
                        "covalent complex requires a dedicated covalent workflow"
                    )
                if heavy_count < 3:
                    redocking_exclusions.append("ligand has fewer than three heavy atoms")
                dk_writer.writerow(
                    {
                        "pdb_id": pdb,
                        "target_uniprot_ids": ";".join(targets_for_pdb),
                        "approved_symbols": ";".join(symbols),
                        "compound_internal_id": compound_internal,
                        "compound_form_id": compound_form,
                        "compound_mapping_status": mapping_status,
                        "compound_mapping_method": mapping_method,
                        "ligand_het_id": het,
                        "ligand_name": ligand_name,
                        **affinity,
                        "resolution": row["resolution"],
                        "release_year": row["release_year"],
                        "protein_pdb_path": str(protein_path),
                        "pocket_pdb_path": str(pocket_path),
                        "ligand_sdf_path": str(sdf_path),
                        "ligand_mol2_path": str(mol2_path),
                        "ligand_sdf_sha256": ligand_hash,
                        "ligand_atom_count": atom_count,
                        "ligand_heavy_atom_count": heavy_count,
                        "pocket_residue_count": len(pocket_residues),
                        "pocket_residues_pdb_numbering": pocket_text,
                        "covalent_complex_flag": int(covalent),
                        "incomplete_ligand_flag": int(incomplete),
                        "peptide_or_oligomer_flag": int(peptide),
                        "redocking_candidate": int(not redocking_exclusions),
                        "redocking_exclusion_reason": "; ".join(
                            redocking_exclusions
                        ),
                    }
                )
                counts["docking_manifest_rows"] += 1
                counts[
                    "redocking_candidate"
                    if not redocking_exclusions
                    else "redocking_noncandidate"
                ] += 1
                pdb_ids.add(pdb)
    finally:
        ev_handle.close()
        dk_handle.close()

    ev_tmp.replace(evidence_path)
    dk_tmp.replace(docking_path)
    report = {
        "status": "passed" if not failures else "passed_with_missing_files",
        "run_type": "incremental_supplement_not_rebuild",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "baseline_policy": "read_only",
        "source_manifest": str(MANIFEST),
        "source_version": "2020.R1 complexes reprocessed with PDBbind v2024 workflow",
        "counts": dict(sorted(counts.items())),
        "unique_pdb_ids": len(pdb_ids),
        "unique_targets": len(targets),
        "unique_mapped_compounds": len(compounds),
        "unique_mapped_target_compound_pairs": len(target_compound_pairs),
        "failures_sample": failures,
        "outputs": {
            "normalized_evidence": str(evidence_path),
            "docking_structure_manifest": str(docking_path),
        },
    }
    qa_path = RUN_ROOT / "qa" / "PDBBIND_NORMALIZATION_V3_QA.json"
    temporary = qa_path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    temporary.replace(qa_path)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
