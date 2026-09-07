#!/usr/bin/env python3
"""Standardize local PDBbind ligand SDF files and compare them with V6.0."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
DEFAULT_RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
DEFAULT_RUN = ROOT / "runs" / "incremental_v61_20260727"
DEFAULT_SDF_ROOT = ROOT / "raw" / "pdbbind_2020r1" / "membrane_subset"
DEFAULT_ANNOTATIONS = (
    ROOT
    / "runs"
    / "incremental_20260727"
    / "staging"
    / "pdbbind"
    / "pdbbind_normalized_evidence_v3.tsv.gz"
)


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def load_master(path: Path):
    full: dict[str, list[tuple[str, str]]] = defaultdict(list)
    connectivity: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            compound_id = row["compound_internal_id"]
            name = row["preferred_name"]
            key = row["standard_inchikey"].strip().upper()
            conn = row["connectivity_key"].strip().upper()
            if key:
                full[key].append((compound_id, name))
            if conn:
                connectivity[conn].append((compound_id, name, key))
    return full, connectivity


def load_annotations(path: Path) -> dict[str, dict]:
    result: dict[str, dict] = {}
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            complex_id = row["source_record_id"].split(":", 1)[0].lower()
            item = result.setdefault(
                complex_id,
                {
                    "compound_name": row.get("compound_name", ""),
                    "ligand_het_id": row.get("ligand_het_id", ""),
                    "covalent": False,
                    "incomplete": False,
                    "peptide": False,
                    "review": False,
                    "target_ids": set(),
                },
            )
            item["covalent"] |= row.get("covalent_complex_flag", "") == "1"
            item["incomplete"] |= row.get("incomplete_ligand_flag", "") == "1"
            item["peptide"] |= row.get("peptide_or_oligomer_flag", "") == "1"
            item["review"] |= row.get("record_qc_status", "") == "review"
            if row.get("target_uniprot_id", ""):
                item["target_ids"].add(row["target_uniprot_id"])
    return result


def parse_sdf(path: Path):
    error = ""
    mol = None
    try:
        supplier = Chem.SDMolSupplier(
            str(path), removeHs=False, sanitize=True, strictParsing=False
        )
        mol = next((candidate for candidate in supplier if candidate is not None), None)
    except Exception as exc:
        error = f"sanitized_read:{type(exc).__name__}:{exc}"
    if mol is None:
        try:
            block = path.read_text(encoding="utf-8", errors="replace")
            mol = Chem.MolFromMolBlock(
                block, sanitize=False, removeHs=False, strictParsing=False
            )
            if mol is not None:
                Chem.SanitizeMol(mol)
        except Exception as exc:
            error = f"fallback_read:{type(exc).__name__}:{exc}"
            mol = None
    if mol is None:
        return None, error or "RDKit returned no molecule"
    try:
        mol = Chem.RemoveHs(mol, sanitize=True)
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        connectivity_smiles = Chem.MolToSmiles(
            mol, canonical=True, isomericSmiles=False
        )
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi) if inchi else ""
        fragments = len(Chem.GetMolFrags(mol))
        formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        return (
            {
                "standard_smiles": smiles,
                "connectivity_smiles": connectivity_smiles,
                "standard_inchi": inchi,
                "standard_inchikey": inchikey,
                "connectivity_key": inchikey[:14] if inchikey else "",
                "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
                "molecular_weight": round(Descriptors.MolWt(mol), 6),
                "formal_charge": formal_charge,
                "atom_count": mol.GetNumAtoms(),
                "heavy_atom_count": mol.GetNumHeavyAtoms(),
                "fragment_count": fragments,
            },
            error,
        )
    except Exception as exc:
        return None, f"standardization:{type(exc).__name__}:{exc}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--sdf-root", type=Path, default=DEFAULT_SDF_ROOT)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    args = parser.parse_args()

    output_dir = args.run / "staging" / "pdbbind_identity"
    qa_dir = args.run / "qa"
    output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir.mkdir(parents=True, exist_ok=True)
    progress = qa_dir / "PDBBIND_IDENTITY_V61_PROGRESS.json"
    started = now()

    full_map, connectivity_map = load_master(
        args.release / "small_molecule_master_v1_1.tsv"
    )
    annotations = load_annotations(args.annotations)
    sdf_paths = sorted(args.sdf_root.rglob("*_ligand.sdf"))

    fields = [
        "complex_id",
        "sdf_path",
        "compound_name",
        "ligand_het_id",
        "target_uniprot_ids",
        "standard_smiles",
        "connectivity_smiles",
        "standard_inchi",
        "standard_inchikey",
        "connectivity_key",
        "molecular_formula",
        "molecular_weight",
        "formal_charge",
        "atom_count",
        "heavy_atom_count",
        "fragment_count",
        "covalent_complex_flag",
        "incomplete_ligand_flag",
        "peptide_or_oligomer_flag",
        "identity_resolution",
        "matched_compound_internal_ids",
        "matched_compound_names",
        "structure_eligibility",
        "standardization_error",
    ]
    output = output_dir / "pdbbind_ligand_identity_v61.tsv.gz"
    status_counts: dict[str, int] = defaultdict(int)
    new_clusters: dict[str, dict] = {}

    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for index, path in enumerate(sdf_paths, 1):
            complex_id = path.name.removesuffix("_ligand.sdf").lower()
            annotation = annotations.get(complex_id, {})
            structure, error = parse_sdf(path)
            row = {
                "complex_id": complex_id,
                "sdf_path": str(path),
                "compound_name": annotation.get("compound_name", ""),
                "ligand_het_id": annotation.get("ligand_het_id", ""),
                "target_uniprot_ids": ";".join(
                    sorted(annotation.get("target_ids", set()))
                ),
                "covalent_complex_flag": int(annotation.get("covalent", False)),
                "incomplete_ligand_flag": int(annotation.get("incomplete", False)),
                "peptide_or_oligomer_flag": int(annotation.get("peptide", False)),
                "standardization_error": error,
            }
            if structure is None:
                row.update(
                    {
                        "identity_resolution": "unresolved_invalid_structure",
                        "structure_eligibility": "manual_review",
                    }
                )
            else:
                row.update(structure)
                key = structure["standard_inchikey"].upper()
                conn_key = structure["connectivity_key"].upper()
                exact = full_map.get(key, [])
                connected = connectivity_map.get(conn_key, [])
                if exact:
                    row["identity_resolution"] = "exact_full_inchikey"
                    matches = exact
                elif connected:
                    row["identity_resolution"] = "connectivity_only_form_or_stereo"
                    matches = [(x[0], x[1]) for x in connected]
                else:
                    row["identity_resolution"] = "new_structure_candidate"
                    matches = []
                    if key:
                        cluster = new_clusters.setdefault(
                            key,
                            {
                                "standard_inchikey": key,
                                "connectivity_key": conn_key,
                                "standard_smiles": structure["standard_smiles"],
                                "standard_inchi": structure["standard_inchi"],
                                "molecular_formula": structure["molecular_formula"],
                                "molecular_weight": structure["molecular_weight"],
                                "complex_ids": [],
                            },
                        )
                        cluster["complex_ids"].append(complex_id)
                row["matched_compound_internal_ids"] = ";".join(
                    sorted({item[0] for item in matches})
                )
                row["matched_compound_names"] = ";".join(
                    sorted({item[1] for item in matches})
                )
                if annotation.get("incomplete") or annotation.get("peptide"):
                    row["structure_eligibility"] = "exclude_or_manual_review"
                elif structure["heavy_atom_count"] < 2:
                    row["structure_eligibility"] = "exclude_atomic_ion"
                elif structure["fragment_count"] > 1:
                    row["structure_eligibility"] = "form_hierarchy_review"
                elif annotation.get("covalent"):
                    row["structure_eligibility"] = "covalent_separate_layer"
                else:
                    row["structure_eligibility"] = "eligible_small_molecule"

            status_counts[row["identity_resolution"]] += 1
            writer.writerow({field: row.get(field, "") for field in fields})
            if index % 500 == 0:
                write_json(
                    progress,
                    {
                        "status": "running",
                        "started_utc": started,
                        "updated_utc": now(),
                        "sdf_total": len(sdf_paths),
                        "sdf_processed": index,
                        "identity_resolution_counts": status_counts,
                    },
                )

    cluster_output = output_dir / "pdbbind_new_structure_clusters_v61.tsv.gz"
    with gzip.open(cluster_output, "wt", encoding="utf-8", newline="") as handle:
        cluster_fields = [
            "standard_inchikey",
            "connectivity_key",
            "standard_smiles",
            "standard_inchi",
            "molecular_formula",
            "molecular_weight",
            "complex_count",
            "complex_ids",
        ]
        writer = csv.DictWriter(
            handle, fieldnames=cluster_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(new_clusters):
            item = new_clusters[key]
            writer.writerow(
                {
                    **item,
                    "complex_count": len(item["complex_ids"]),
                    "complex_ids": ";".join(sorted(item["complex_ids"])),
                }
            )

    summary = {
        "status": "complete",
        "started_utc": started,
        "completed_utc": now(),
        "sdf_total": len(sdf_paths),
        "identity_resolution_counts": status_counts,
        "new_structure_cluster_count": len(new_clusters),
        "output": str(output),
        "new_cluster_output": str(cluster_output),
        "rdkit_version": Chem.rdBase.rdkitVersion,
    }
    write_json(progress, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
