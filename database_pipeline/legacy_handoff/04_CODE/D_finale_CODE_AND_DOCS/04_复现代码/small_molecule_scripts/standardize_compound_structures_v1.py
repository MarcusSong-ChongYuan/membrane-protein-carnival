from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize


RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")
csv.field_size_limit(100_000_000)

ROOT = Path(r"D:\7.22")
RELEASE = ROOT / "membrane_master_v5_working" / "releases" / "release_v5_4_pgd_binding"
WORK = ROOT / "small_molecule_v1_working"
RAW = WORK / "raw"
INTERMEDIATE = WORK / "intermediate"
CHUNKS = INTERMEDIATE / "standardized_chunks"
REPORTS = WORK / "reports"
CHUNKS.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

INDEX = RELEASE / "small_molecule_index_v4_1.tsv"
CHEMBL = RAW / "chembl37" / "chembl37_selected_chemreps.tsv"
PUBCHEM = RAW / "pubchem" / "pubchem_missing_properties_v1.jsonl"
REPORT = REPORTS / "STRUCTURE_STANDARDIZATION_PROGRESS.json"

CHUNK_SIZE = 5_000
AMIDE = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
METALS = {
    3,
    4,
    11,
    12,
    13,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    84,
    87,
    88,
    89,
    90,
    91,
    92,
}

OUTPUT_FIELDS = [
    "source_record_id",
    "input_row_number",
    "source_compound_id",
    "source_compound_id_type",
    "source_compound_name",
    "source_scope_status",
    "source_scope_class",
    "source_is_core_small_molecule",
    "source_databases",
    "source_versions",
    "source_record_qc_status",
    "pubchem_cids",
    "chembl_ids",
    "structure_input_source",
    "input_smiles",
    "input_inchikey",
    "structure_parse_status",
    "exact_standard_smiles",
    "exact_standard_inchi",
    "exact_standard_inchikey",
    "exact_identity_key",
    "parent_standard_smiles",
    "parent_standard_inchi",
    "parent_standard_inchikey",
    "parent_identity_key",
    "exact_connectivity_key",
    "parent_connectivity_key",
    "fragment_count",
    "form_type",
    "exact_molecular_formula",
    "exact_molecular_weight",
    "exact_formal_charge",
    "parent_molecular_formula",
    "parent_molecular_weight",
    "parent_exact_mass",
    "parent_xlogp",
    "parent_tpsa",
    "parent_hbond_donor_count",
    "parent_hbond_acceptor_count",
    "parent_rotatable_bond_count",
    "parent_formal_charge",
    "parent_heavy_atom_count",
    "parent_ring_count",
    "parent_max_ring_size",
    "parent_amide_bond_count",
    "computed_structural_class",
    "contains_metal",
    "standardization_actions",
    "structure_qc_note",
]
_WORKER_CHEMBL: dict[str, dict] = {}
_WORKER_PUBCHEM: dict[str, dict] = {}


def normalize_chembl_id(value: str) -> str:
    value = (value or "").strip().upper()
    if value.startswith("CHEMBL:"):
        value = value.split(":", 1)[1]
    return value if value.startswith("CHEMBL") else ""


def split_ids(value: str) -> list[str]:
    return [
        item.strip()
        for item in (value or "").replace("|", ";").replace(",", ";").split(";")
        if item.strip()
    ]


def load_chembl() -> dict[str, dict]:
    result = {}
    with CHEMBL.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[row["chembl_id"].upper()] = row
    return result


def load_pubchem() -> dict[str, dict]:
    result = {}
    with PUBCHEM.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                row = json.loads(line)
                result[str(row["CID"])] = row
    return result


def initialize_worker() -> None:
    global _WORKER_CHEMBL, _WORKER_PUBCHEM
    RDLogger.DisableLog("rdApp.warning")
    RDLogger.DisableLog("rdApp.error")
    _WORKER_CHEMBL = load_chembl()
    _WORKER_PUBCHEM = load_pubchem()


def standardize_task(item: tuple[int, dict]) -> dict:
    row_number, row = item
    return standardize_row(
        row,
        row_number,
        _WORKER_CHEMBL,
        _WORKER_PUBCHEM,
    )


def identity_key(inchikey: str, smiles: str) -> str:
    if inchikey:
        return f"IK:{inchikey}"
    if smiles:
        return "SMI:" + hashlib.sha256(smiles.encode("utf-8")).hexdigest()[:32]
    return ""


def inchi_and_key(mol: Chem.Mol) -> tuple[str, str]:
    try:
        inchi = Chem.MolToInchi(mol)
        key = Chem.InchiToInchiKey(inchi) if inchi else ""
        return inchi, key
    except Exception:
        return "", ""


def molecule_summary(mol: Chem.Mol) -> dict:
    atom_rings = mol.GetRingInfo().AtomRings()
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    oxygen_count = sum(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    nitrogen_count = sum(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    metal_count = sum(atom.GetAtomicNum() in METALS for atom in mol.GetAtoms())
    amide_count = len(mol.GetSubstructMatches(AMIDE)) if AMIDE else 0
    max_ring = max((len(ring) for ring in atom_rings), default=0)
    mw = Descriptors.MolWt(mol)
    if metal_count and carbon_count:
        structural_class = "organometallic_or_metal_containing"
    elif metal_count or carbon_count == 0:
        structural_class = "inorganic_chemical_entity"
    elif max_ring >= 12:
        structural_class = "macrocyclic_compound"
    elif amide_count >= 4 and mw >= 400:
        structural_class = "peptidic_or_peptidomimetic"
    elif (
        carbon_count >= 4
        and oxygen_count >= 3
        and oxygen_count / max(carbon_count, 1) >= 0.55
        and len(atom_rings) >= 1
    ):
        structural_class = "carbohydrate_like"
    elif carbon_count >= 16 and oxygen_count <= 6 and nitrogen_count <= 4:
        structural_class = "lipid_like"
    elif len(atom_rings) >= 3 and 17 <= carbon_count <= 35:
        structural_class = "polycyclic_or_steroid_like"
    elif len(atom_rings) >= 1:
        structural_class = "cyclic_organic_compound"
    else:
        structural_class = "acyclic_organic_compound"
    return {
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(mw, 6),
        "exact_mass": round(Descriptors.ExactMolWt(mol), 6),
        "xlogp": round(Crippen.MolLogP(mol), 6),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 6),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable": int(Lipinski.NumRotatableBonds(mol)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "ring_count": int(len(atom_rings)),
        "max_ring": int(max_ring),
        "amide_count": int(amide_count),
        "structural_class": structural_class,
        "contains_metal": int(metal_count > 0),
    }


def select_structure(
    row: dict, chembl: dict[str, dict], pubchem: dict[str, dict]
) -> tuple[str, str]:
    original = row.get("canonical_smiles", "").strip()
    if original:
        return original, "v4.1_index"

    compound_id = row.get("compound_id", "").strip()
    if row.get("compound_id_type") == "PubChem CID" and compound_id in pubchem:
        record = pubchem[compound_id]
        smiles = str(record.get("SMILES") or record.get("ConnectivitySMILES") or "")
        if smiles:
            return smiles, "PubChem_PUG_REST_2026-07-26"

    chembl_candidates = []
    normalized = normalize_chembl_id(compound_id)
    if normalized:
        chembl_candidates.append(normalized)
    for item in split_ids(row.get("chembl_ids", "")):
        normalized = normalize_chembl_id(item)
        if normalized:
            chembl_candidates.append(normalized)
    for chembl_id in chembl_candidates:
        record = chembl.get(chembl_id)
        if record and record.get("canonical_smiles", "").strip():
            return record["canonical_smiles"].strip(), "ChEMBL_37_chemreps"
    return "", "missing"


def standardize_row(
    row: dict, row_number: int, chembl: dict[str, dict], pubchem: dict[str, dict]
) -> dict:
    compound_id = row.get("compound_id", "").strip()
    compound_type = row.get("compound_id_type", "").strip()
    source_record_id = "SRC-" + hashlib.sha256(
        f"{compound_type}|{compound_id}".encode("utf-8")
    ).hexdigest()[:20].upper()
    result = {
        "source_record_id": source_record_id,
        "input_row_number": row_number,
        "source_compound_id": compound_id,
        "source_compound_id_type": compound_type,
        "source_compound_name": row.get("compound_name", "").strip(),
        "source_scope_status": row.get("small_molecule_scope_status", "").strip(),
        "source_scope_class": row.get("small_molecule_scope_class", "").strip(),
        "source_is_core_small_molecule": row.get("is_core_small_molecule", "").strip(),
        "source_databases": row.get("source_databases", "").strip(),
        "source_versions": row.get("source_versions", "").strip(),
        "source_record_qc_status": row.get("record_qc_status", "").strip(),
        "pubchem_cids": row.get("pubchem_cids", "").strip(),
        "chembl_ids": row.get("chembl_ids", "").strip(),
    }
    smiles, structure_source = select_structure(row, chembl, pubchem)
    result["structure_input_source"] = structure_source
    result["input_smiles"] = smiles
    result["input_inchikey"] = row.get("inchikey", "").strip()
    if not smiles:
        result.update(
            {
                "structure_parse_status": "no_structure",
                "structure_qc_note": "No parseable structure supplied by current sources",
            }
        )
        return result

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result.update(
            {
                "structure_parse_status": "parse_failed",
                "structure_qc_note": "RDKit could not parse the supplied SMILES",
            }
        )
        return result

    actions = []
    try:
        cleaned = rdMolStandardize.Cleanup(mol)
        if Chem.MolToSmiles(mol, isomericSmiles=True) != Chem.MolToSmiles(
            cleaned, isomericSmiles=True
        ):
            actions.append("cleanup_normalized")
    except Exception as exc:
        cleaned = mol
        actions.append("cleanup_failed_fallback_to_parsed")
        result["structure_qc_note"] = f"Cleanup failed: {type(exc).__name__}"

    try:
        parent = rdMolStandardize.FragmentParent(cleaned)
    except Exception:
        parent = cleaned
        actions.append("fragment_parent_failed")
    if parent is None or parent.GetNumAtoms() == 0:
        parent = cleaned
        actions.append("empty_fragment_parent_fallback")

    try:
        uncharged_parent = rdMolStandardize.Uncharger().uncharge(parent)
        if Chem.MolToSmiles(parent, isomericSmiles=True) != Chem.MolToSmiles(
            uncharged_parent, isomericSmiles=True
        ):
            actions.append("parent_charge_normalized")
        parent = uncharged_parent
    except Exception:
        actions.append("parent_uncharging_failed")

    for label, current in (("exact", cleaned), ("parent", parent)):
        try:
            current.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(current)
        except Exception:
            actions.append(f"{label}_ring_initialization_failed")

    exact_smiles = Chem.MolToSmiles(cleaned, canonical=True, isomericSmiles=True)
    parent_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
    exact_inchi, exact_inchikey = inchi_and_key(cleaned)
    parent_inchi, parent_inchikey = inchi_and_key(parent)
    fragments = len(Chem.GetMolFrags(cleaned))
    if fragments > 1 and exact_smiles != parent_smiles:
        form_type = "salt_or_multicomponent_form"
        actions.append("fragments_removed_for_parent")
    elif exact_smiles != parent_smiles:
        form_type = "charge_or_standardized_form"
    else:
        form_type = "parent_form"

    exact_summary = molecule_summary(cleaned)
    parent_summary = molecule_summary(parent)
    qc_notes = []
    input_key = result["input_inchikey"]
    if input_key and exact_inchikey and input_key != exact_inchikey:
        qc_notes.append("input_inchikey_differs_from_computed_standard")
    if not exact_inchikey:
        qc_notes.append("standard_inchikey_unavailable")
    if result.get("structure_qc_note"):
        qc_notes.append(result["structure_qc_note"])

    result.update(
        {
            "structure_parse_status": "parsed",
            "exact_standard_smiles": exact_smiles,
            "exact_standard_inchi": exact_inchi,
            "exact_standard_inchikey": exact_inchikey,
            "exact_identity_key": identity_key(exact_inchikey, exact_smiles),
            "parent_standard_smiles": parent_smiles,
            "parent_standard_inchi": parent_inchi,
            "parent_standard_inchikey": parent_inchikey,
            "parent_identity_key": identity_key(parent_inchikey, parent_smiles),
            "exact_connectivity_key": exact_inchikey[:14] if exact_inchikey else "",
            "parent_connectivity_key": parent_inchikey[:14] if parent_inchikey else "",
            "fragment_count": fragments,
            "form_type": form_type,
            "exact_molecular_formula": exact_summary["molecular_formula"],
            "exact_molecular_weight": exact_summary["molecular_weight"],
            "exact_formal_charge": exact_summary["formal_charge"],
            "parent_molecular_formula": parent_summary["molecular_formula"],
            "parent_molecular_weight": parent_summary["molecular_weight"],
            "parent_exact_mass": parent_summary["exact_mass"],
            "parent_xlogp": parent_summary["xlogp"],
            "parent_tpsa": parent_summary["tpsa"],
            "parent_hbond_donor_count": parent_summary["hbd"],
            "parent_hbond_acceptor_count": parent_summary["hba"],
            "parent_rotatable_bond_count": parent_summary["rotatable"],
            "parent_formal_charge": parent_summary["formal_charge"],
            "parent_heavy_atom_count": parent_summary["heavy_atoms"],
            "parent_ring_count": parent_summary["ring_count"],
            "parent_max_ring_size": parent_summary["max_ring"],
            "parent_amide_bond_count": parent_summary["amide_count"],
            "computed_structural_class": parent_summary["structural_class"],
            "contains_metal": parent_summary["contains_metal"],
            "standardization_actions": ";".join(actions) if actions else "none",
            "structure_qc_note": ";".join(qc_notes),
        }
    )
    return result


def completed_chunks() -> set[int]:
    result = set()
    for path in CHUNKS.glob("standardized_chunk_*.tsv"):
        try:
            result.add(int(path.stem.rsplit("_", 1)[1]))
        except ValueError:
            continue
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-chunks", type=int, default=3)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    with INDEX.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    total_chunks = (len(rows) + CHUNK_SIZE - 1) // CHUNK_SIZE
    done = completed_chunks()
    processed_now = 0
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=initialize_worker,
    ) as executor:
        for chunk_index in range(total_chunks):
            if chunk_index in done:
                continue
            start = chunk_index * CHUNK_SIZE
            end = min(len(rows), start + CHUNK_SIZE)
            output_path = CHUNKS / f"standardized_chunk_{chunk_index:04d}.tsv"
            temp_path = output_path.with_suffix(".tsv.tmp")
            tasks = list(enumerate(rows[start:end], start=start + 1))
            results = executor.map(standardize_task, tasks, chunksize=50)
            with temp_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=OUTPUT_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="ignore",
                )
                writer.writeheader()
                writer.writerows(results)
            os.replace(temp_path, output_path)
            processed_now += 1
            done.add(chunk_index)
            print(
                json.dumps(
                    {
                        "completed_chunk": chunk_index,
                        "completed_chunks": len(done),
                        "total_chunks": total_chunks,
                        "rows_in_chunk": end - start,
                    }
                ),
                flush=True,
            )
            if processed_now >= args.max_chunks:
                break

    report = {
        "input_rows": len(rows),
        "chunk_size": CHUNK_SIZE,
        "total_chunks": total_chunks,
        "completed_chunks": len(done),
        "completed": len(done) == total_chunks,
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "standardization_policy": "MemProDB-compound-standardization-v1",
    }
    REPORT.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
