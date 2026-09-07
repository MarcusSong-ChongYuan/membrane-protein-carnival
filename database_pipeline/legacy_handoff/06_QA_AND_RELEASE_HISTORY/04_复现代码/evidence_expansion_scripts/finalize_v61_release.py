#!/usr/bin/env python3
"""Build, validate, and freeze the formal MemPro V6.1 release.

V6.0 is immutable. V6.1 promotes only structure-resolved records from the
V6.0 review queue. Full InChIKey identity is required for automatic merging;
connectivity-only, ambiguous, generic, and covalent-special cases remain in
review. Every V6.0 review row receives one and only one V6.1 disposition.
"""

from __future__ import annotations

import csv
import datetime as dt
import gc
import gzip
import hashlib
import heapq
import json
import math
import os
import re
import shutil
import sqlite3
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize


RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
BASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
RUN = ROOT / "runs" / "incremental_v61_20260727"
CANDIDATE = RUN / "release_candidate_v6_1"
FINAL = ROOT / "releases" / "release_mempro_v6_1_20260729"
QA = RUN / "qa"
DB = RUN / "merge_candidate" / "v61_formal_build.sqlite"
PROGRESS = QA / "V61_FORMAL_BUILD_PROGRESS.json"

MASTER0 = BASE / "small_molecule_master_v1_1.tsv"
FORM0 = BASE / "compound_form_hierarchy_v1_1.tsv"
EVIDENCE0 = BASE / "binding_evidence_master_v6_0.tsv.gz"
REVIEW0 = BASE / "binding_evidence_review_queue_v6_0.tsv.gz"
NEGATIVE0 = BASE / "negative_binding_evidence_v1_0.tsv.gz"
PAIR0 = BASE / "protein_compound_summary_v6_0.tsv.gz"
SITE0 = BASE / "binding_site_instances_v6_0.tsv.gz"
PROTEIN0 = BASE / "human_membrane_protein_master_v6_0.tsv"
DISEASE0 = BASE / "protein_gene_disease_relations_v6_0.tsv"

PUBCHEM = [
    RUN / "staging" / "pubchem_identity" / f"pubchem_identity_crosswalk_{tag}_v61.tsv.gz"
    for tag in ("p1", "p2", "p3")
]
PDBE = RUN / "staging" / "pdbe_identity" / "pdbe_ccd_identity_v61.tsv.gz"
PDBBIND = (
    RUN / "staging" / "pdbbind_identity" / "pdbbind_ligand_identity_v61.tsv.gz"
)
BRENDA = RUN / "staging" / "brenda_identity" / "brenda_name_identity_v61.tsv.gz"
GATE = QA / "V61_FINALIZATION_GATE.json"

RELEASE_VERSION = "v6.1"
RELEASE_DATE = "2026-07-29"
PUBCHEM_RE = re.compile(r"(?:PUBCHEM:|CID[: ]?)(\d+)", re.IGNORECASE)
TIER_RANK = {"BE1": 1, "BE2": 2, "BE3": 3}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def integer(value: object) -> int:
    try:
        return int(float(clean(value)))
    except (TypeError, ValueError):
        return 0


def number(value: object) -> float | None:
    try:
        value = clean(value)
        return float(value) if value else None
    except (TypeError, ValueError):
        return None


def split_values(value: object) -> list[str]:
    return [
        part.strip()
        for part in clean(value).replace("|", ";").replace(",", ";").split(";")
        if part.strip()
    ]


def normalize_name(value: str) -> str:
    return " ".join(value.casefold().split())


def stable_id(prefix: str, *parts: str, length: int = 20) -> str:
    token = hashlib.sha256("\x1f".join(parts).encode("utf-8")).hexdigest()[:length]
    return f"{prefix}_{token.upper()}"


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def update_progress(stage: str, status: str = "running", **details: object) -> None:
    write_json(
        PROGRESS,
        {
            "status": status,
            "stage": stage,
            "updated_utc": now(),
            **details,
        },
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("rt", encoding="utf-8-sig", newline="")


@contextmanager
def atomic_text(path: Path, fields: list[str], gzipped: bool = False):
    temporary = path.with_suffix(path.suffix + ".tmp")
    opener = gzip.open if gzipped else open
    with opener(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        yield writer
    temporary.replace(path)


def require_inputs() -> None:
    required = [
        MASTER0,
        FORM0,
        EVIDENCE0,
        REVIEW0,
        NEGATIVE0,
        SITE0,
        PROTEIN0,
        DISEASE0,
        PDBE,
        PDBBIND,
        BRENDA,
        GATE,
        *PUBCHEM,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing V6.1 inputs: " + "; ".join(missing))
    gate = json.loads(GATE.read_text(encoding="utf-8"))
    if gate.get("status") != "ready_for_final_refresh":
        raise RuntimeError(f"V6.1 finalization gate is not ready: {gate.get('status')}")
    if FINAL.exists():
        raise FileExistsError(
            f"Formal V6.1 release already exists and will not be overwritten: {FINAL}"
        )


def connect_db() -> sqlite3.Connection:
    if DB.exists():
        DB.unlink()
    DB.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(DB)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA temp_store=FILE")
    con.executescript(
        """
        CREATE TABLE identity (
          source_type TEXT NOT NULL,
          source_id TEXT NOT NULL,
          source_database TEXT NOT NULL,
          exact_key TEXT,
          connectivity_key TEXT,
          smiles TEXT,
          inchi TEXT,
          formula TEXT,
          molecular_weight REAL,
          formal_charge INTEGER,
          preferred_name TEXT,
          identity_resolution TEXT,
          source_action TEXT,
          safe_status TEXT NOT NULL,
          existing_compound_id TEXT,
          compound_internal_id TEXT,
          compound_form_id TEXT,
          notes TEXT,
          PRIMARY KEY(source_type, source_id)
        );
        CREATE TABLE candidate_source (
          exact_key TEXT NOT NULL,
          source_database TEXT NOT NULL,
          source_id TEXT NOT NULL,
          PRIMARY KEY(exact_key, source_database, source_id)
        );
        CREATE TABLE candidate_structure (
          exact_key TEXT PRIMARY KEY,
          connectivity_key TEXT,
          smiles TEXT,
          inchi TEXT,
          formula TEXT,
          molecular_weight REAL,
          formal_charge INTEGER,
          preferred_name TEXT
        );
        CREATE TABLE disposition (
          source_evidence_id TEXT PRIMARY KEY,
          action TEXT NOT NULL,
          source_type TEXT,
          source_id TEXT,
          identity_safe_status TEXT,
          identity_resolution TEXT,
          exact_key TEXT,
          existing_compound_id TEXT,
          reason TEXT
        );
        CREATE TABLE referenced_candidate (
          exact_key TEXT PRIMARY KEY
        );
        CREATE TABLE canonical_form (
          exact_key TEXT PRIMARY KEY,
          parent_key TEXT NOT NULL,
          compound_internal_id TEXT NOT NULL,
          compound_form_id TEXT NOT NULL,
          parent_relation TEXT NOT NULL,
          exact_smiles TEXT,
          exact_inchi TEXT,
          exact_formula TEXT,
          exact_mw REAL,
          exact_charge INTEGER,
          parent_smiles TEXT,
          parent_inchi TEXT,
          parent_formula TEXT,
          parent_mw REAL,
          parent_charge INTEGER,
          connectivity_key TEXT,
          fragment_count INTEGER,
          stereocenter_count INTEGER,
          unassigned_stereocenter_count INTEGER,
          standardization_status TEXT,
          standardization_notes TEXT
        );
        CREATE TABLE pair_raw (
          target TEXT NOT NULL,
          compound TEXT NOT NULL,
          approved_symbol TEXT,
          tier TEXT,
          tier_rank INTEGER,
          source TEXT,
          value_nm REAL,
          form_id TEXT,
          activity_type TEXT,
          origin TEXT NOT NULL
        );
        CREATE TABLE quant_raw (
          target TEXT NOT NULL,
          compound TEXT NOT NULL,
          activity_type TEXT,
          value_nm REAL,
          source TEXT
        );
        CREATE TABLE site_new (
          site_key TEXT PRIMARY KEY,
          evidence_id TEXT,
          target TEXT,
          compound TEXT,
          form_id TEXT,
          pdb_ids TEXT,
          residues TEXT,
          source TEXT,
          source_compound_id TEXT
        );
        """
    )
    return con


def load_baseline_identity() -> tuple[dict, dict, dict, int, int]:
    master_by_id: dict[str, dict[str, str]] = {}
    master_by_key: dict[str, str] = {}
    connectivity: dict[str, list[str]] = defaultdict(list)
    max_compound = 0
    with MASTER0.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            compound_id = row["compound_internal_id"]
            master_by_id[compound_id] = row
            key = row["standard_inchikey"].upper()
            if key:
                if key in master_by_key and master_by_key[key] != compound_id:
                    raise RuntimeError(f"Duplicate baseline full InChIKey: {key}")
                master_by_key[key] = compound_id
            conn_key = row["connectivity_key"].upper()
            if conn_key:
                connectivity[conn_key].append(compound_id)
            match = re.search(r"(\d+)$", compound_id)
            if match:
                max_compound = max(max_compound, int(match.group(1)))
    form_by_key: dict[str, str] = {}
    max_form = 0
    with FORM0.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["exact_inchikey"].upper()
            if key and key not in form_by_key:
                form_by_key[key] = row["compound_form_id"]
            match = re.search(r"(\d+)$", row["compound_form_id"])
            if match:
                max_form = max(max_form, int(match.group(1)))
    return master_by_id, master_by_key, form_by_key, max_compound, max_form


def add_identity(
    con: sqlite3.Connection,
    *,
    source_type: str,
    source_id: str,
    source_database: str,
    exact_key: str = "",
    connectivity_key: str = "",
    smiles: str = "",
    inchi: str = "",
    formula: str = "",
    molecular_weight: float | None = None,
    formal_charge: int | None = None,
    preferred_name: str = "",
    identity_resolution: str = "",
    source_action: str = "",
    safe_status: str,
    existing_compound_id: str = "",
    notes: str = "",
) -> None:
    exact_key = exact_key.upper()
    con.execute(
        """
        INSERT OR REPLACE INTO identity
        (source_type,source_id,source_database,exact_key,connectivity_key,
         smiles,inchi,formula,molecular_weight,formal_charge,preferred_name,
         identity_resolution,source_action,safe_status,existing_compound_id,
         compound_internal_id,compound_form_id,notes)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
        (
            source_type,
            source_id,
            source_database,
            exact_key,
            connectivity_key.upper(),
            smiles,
            inchi,
            formula,
            molecular_weight,
            formal_charge,
            preferred_name,
            identity_resolution,
            source_action,
            safe_status,
            existing_compound_id,
            existing_compound_id,
            "",
            notes,
        ),
    )
    if safe_status == "exact_new" and exact_key and smiles:
        con.execute(
            """
            INSERT OR IGNORE INTO candidate_structure
            (exact_key,connectivity_key,smiles,inchi,formula,molecular_weight,
             formal_charge,preferred_name)
            VALUES (?,?,?,?,?,?,?,?)
            """,
            (
                exact_key,
                connectivity_key.upper(),
                smiles,
                inchi,
                formula,
                molecular_weight,
                formal_charge,
                preferred_name,
            ),
        )
        con.execute(
            "INSERT OR IGNORE INTO candidate_source VALUES (?,?,?)",
            (exact_key, source_database, source_id),
        )


def load_identity_sources(con: sqlite3.Connection) -> dict:
    counts: Counter[str] = Counter()
    for path in PUBCHEM:
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                existing = split_values(row["matched_compound_internal_ids"])
                conflict = row["cid_structure_conflict"] == "1"
                if (
                    row["release_action"] == "map_existing_exact"
                    and len(existing) == 1
                    and not conflict
                ):
                    safe = "exact_existing"
                elif (
                    row["release_action"] == "create_candidate_after_global_dedup"
                    and row["inchikey"]
                    and row["smiles"]
                    and not conflict
                ):
                    safe = "exact_new"
                else:
                    safe = "identity_review"
                add_identity(
                    con,
                    source_type="PUBCHEM_CID",
                    source_id=row["cid"],
                    source_database="PubChem BioAssay",
                    exact_key=row["inchikey"],
                    connectivity_key=row["connectivity_key"],
                    smiles=row["smiles"],
                    inchi=row["inchi"],
                    formula=row["molecular_formula"],
                    molecular_weight=number(row["molecular_weight"]),
                    formal_charge=integer(row["formal_charge"]),
                    preferred_name=f"PubChem CID {row['cid']}",
                    identity_resolution=row["identity_resolution"],
                    source_action=row["release_action"],
                    safe_status=safe,
                    existing_compound_id=existing[0] if len(existing) == 1 else "",
                    notes=(
                        "CID structure conflict"
                        if conflict
                        else row["fetch_error_message"]
                    ),
                )
                counts[f"PubChem:{safe}"] += 1
        con.commit()

    with gzip.open(PDBE, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            existing = split_values(row["matched_compound_internal_ids"])
            eligible = row["release_eligibility"] == "eligible_small_molecule"
            if eligible and row["identity_resolution"] == "exact_full_inchikey" and len(existing) == 1:
                safe = "exact_existing"
            elif eligible and row["identity_resolution"] == "new_structure_candidate" and row["inchikey"] and row["smiles_stereo"]:
                safe = "exact_new"
            elif row["release_eligibility"].startswith("exclude"):
                safe = "excluded_identity"
            else:
                safe = "identity_review"
            add_identity(
                con,
                source_type="PDB_CCD",
                source_id=row["component_id"].upper(),
                source_database="PDBe",
                exact_key=row["inchikey"],
                connectivity_key=row["connectivity_key"],
                smiles=row["smiles_stereo"] or row["smiles"],
                inchi=row["inchi"],
                formula=row["formula"],
                molecular_weight=number(row["formula_weight"]),
                formal_charge=integer(row["formal_charge"]),
                preferred_name=row["preferred_name"] or row["ccd_name"],
                identity_resolution=row["identity_resolution"],
                source_action=row["release_eligibility"],
                safe_status=safe,
                existing_compound_id=existing[0] if len(existing) == 1 else "",
                notes=row["error_message"],
            )
            counts[f"PDBe:{safe}"] += 1
    con.commit()

    with gzip.open(PDBBIND, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            existing = split_values(row["matched_compound_internal_ids"])
            eligible = row["structure_eligibility"] == "eligible_small_molecule"
            if eligible and row["identity_resolution"] == "exact_full_inchikey" and len(existing) == 1:
                safe = "exact_existing"
            elif eligible and row["identity_resolution"] == "new_structure_candidate" and row["standard_inchikey"] and row["standard_smiles"]:
                safe = "exact_new"
            elif row["structure_eligibility"] == "covalent_separate_layer":
                safe = "covalent_review"
            elif row["structure_eligibility"].startswith("exclude"):
                safe = "excluded_identity"
            else:
                safe = "identity_review"
            add_identity(
                con,
                source_type="PDBBIND_COMPLEX",
                source_id=row["complex_id"].lower(),
                source_database="PDBbind",
                exact_key=row["standard_inchikey"],
                connectivity_key=row["connectivity_key"],
                smiles=row["standard_smiles"],
                inchi=row["standard_inchi"],
                formula=row["molecular_formula"],
                molecular_weight=number(row["molecular_weight"]),
                formal_charge=integer(row["formal_charge"]),
                preferred_name=row["compound_name"],
                identity_resolution=row["identity_resolution"],
                source_action=row["structure_eligibility"],
                safe_status=safe,
                existing_compound_id=existing[0] if len(existing) == 1 else "",
                notes=row["standardization_error"],
            )
            counts[f"PDBbind:{safe}"] += 1
    con.commit()

    with gzip.open(BRENDA, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            status = row["resolution_status"]
            existing = split_values(row["matched_compound_internal_ids"])
            keys = split_values(row["pubchem_result_inchikeys"])
            smiles = split_values(row["structure_smiles"])
            if status in {
                "verified_existing_exact_structure",
                "resolved_unique_name_to_existing_exact_structure",
            } and len(existing) == 1:
                safe = "exact_existing"
            elif (
                status == "new_structure_candidate_from_unique_name"
                and len(keys) == 1
                and len(smiles) == 1
            ):
                safe = "exact_new"
            elif status == "excluded_generic_or_common":
                safe = "excluded_identity"
            else:
                safe = "identity_review"
            add_identity(
                con,
                source_type="BRENDA_NAME",
                source_id=row["normalized_name"],
                source_database="BRENDA",
                exact_key=keys[0] if len(keys) == 1 else "",
                connectivity_key=(keys[0][:14] if len(keys) == 1 else ""),
                smiles=smiles[0] if len(smiles) == 1 else "",
                preferred_name=row["preferred_name"],
                identity_resolution=status,
                source_action=status,
                safe_status=safe,
                existing_compound_id=existing[0] if len(existing) == 1 else "",
                notes=row["message"],
            )
            counts[f"BRENDA:{safe}"] += 1
    con.commit()
    con.execute("CREATE INDEX idx_identity_exact ON identity(exact_key)")
    con.execute("CREATE INDEX idx_identity_compound ON identity(existing_compound_id)")
    con.commit()
    return dict(counts)


def source_lookup(row: dict[str, str]) -> tuple[str, str]:
    source = row["source_database"]
    if source.startswith("PubChem"):
        match = PUBCHEM_RE.search(row["compound_source_id"])
        return "PUBCHEM_CID", match.group(1) if match else ""
    if source == "PDBe":
        return "PDB_CCD", row["compound_source_id"].split(":", 1)[-1].upper()
    if source == "PDBbind":
        complex_id = (
            row["source_record_id"].split(":", 1)[0].lower()
            or row["pdb_ids"].split(";", 1)[0].lower()
        )
        return "PDBBIND_COMPLEX", complex_id
    if source == "BRENDA":
        return "BRENDA_NAME", normalize_name(row["compound_name"])
    return "", ""


def classify_dispositions(con: sqlite3.Connection) -> dict:
    counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    total = 0
    with gzip.open(REVIEW0, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            total += 1
            source_type, source_id = source_lookup(row)
            identity = None
            if source_type and source_id:
                identity = con.execute(
                    """
                    SELECT safe_status,identity_resolution,exact_key,
                           existing_compound_id,source_action,notes
                    FROM identity WHERE source_type=? AND source_id=?
                    """,
                    (source_type, source_id),
                ).fetchone()
            outcome = row["activity_outcome"].casefold()
            if identity is None:
                action, reason = "review", "no_identity_lookup"
                safe = resolution = exact_key = existing = ""
            else:
                safe, resolution, exact_key, existing, source_action, notes = identity
                if safe == "excluded_identity":
                    action, reason = "excluded", source_action
                elif safe in {"identity_review", "covalent_review"}:
                    action, reason = "review", safe
                elif outcome == "inactive":
                    action, reason = "negative", "reported_inactive"
                elif outcome == "inconclusive":
                    action, reason = "review", "inconclusive_activity"
                elif outcome in {"active", "observed", "positive_or_observed"}:
                    action, reason = "positive", "structure_resolved_positive"
                else:
                    action, reason = "review", "unspecified_activity_outcome"
            con.execute(
                """
                INSERT OR REPLACE INTO disposition
                (source_evidence_id,action,source_type,source_id,
                 identity_safe_status,identity_resolution,exact_key,
                 existing_compound_id,reason)
                VALUES (?,?,?,?,?,?,?,?,?)
                """,
                (
                    row["source_evidence_id"],
                    action,
                    source_type,
                    source_id,
                    safe,
                    resolution,
                    exact_key,
                    existing,
                    reason,
                ),
            )
            if action in {"positive", "negative"} and safe == "exact_new" and exact_key:
                con.execute(
                    "INSERT OR IGNORE INTO referenced_candidate VALUES (?)",
                    (exact_key,),
                )
            counts[action] += 1
            source_counts[f"{row['source_database']}:{action}"] += 1
            if total % 100000 == 0:
                con.commit()
                update_progress(
                    "classify_review_queue",
                    rows_processed=total,
                    disposition_counts=dict(counts),
                )
    con.commit()
    return {
        "review_rows_total": total,
        "disposition_counts": dict(counts),
        "source_disposition_counts": dict(source_counts),
    }


def mol_identity(mol: Chem.Mol) -> tuple[str, str, str]:
    smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    inchi = Chem.MolToInchi(mol)
    key = Chem.InchiToInchiKey(inchi) if inchi else ""
    return smiles, inchi, key


def mol_properties(mol: Chem.Mol) -> dict[str, object]:
    formula = rdMolDescriptors.CalcMolFormula(mol)
    return {
        "formula": formula,
        "mw": float(Descriptors.MolWt(mol)),
        "exact_mass": float(Descriptors.ExactMolWt(mol)),
        "xlogp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable": int(Lipinski.NumRotatableBonds(mol)),
        "charge": int(Chem.GetFormalCharge(mol)),
        "heavy": int(mol.GetNumHeavyAtoms()),
        "rings": int(rdMolDescriptors.CalcNumRings(mol)),
        "amide": int(rdMolDescriptors.CalcNumAmideBonds(mol)),
    }


def canonicalize_candidate_worker(payload: tuple) -> tuple[str, dict | None, str, str]:
    """Canonicalize one exact structure without touching shared state."""
    key, conn_key, smiles, inchi, formula, mw, charge, preferred = payload
    if smiles is None:
        return key, None, "missing_candidate_structure", "missing candidate structure"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("RDKit returned no molecule")
        cleaned = rdMolStandardize.Cleanup(mol)
        exact_smiles, exact_inchi, computed_key = mol_identity(cleaned)
        if computed_key != key:
            raise ValueError(f"full_inchikey_mismatch:{computed_key or 'missing'}")
        fragments = len(Chem.GetMolFrags(cleaned))
        fragment_parent = rdMolStandardize.FragmentParent(cleaned)
        charge_parent = rdMolStandardize.ChargeParent(fragment_parent)
        parent_smiles, parent_inchi, parent_key = mol_identity(charge_parent)
        if not parent_key:
            raise ValueError("missing_parent_inchikey")
        centers = Chem.FindMolChiralCenters(
            cleaned, includeUnassigned=True, includeCIP=True
        )
        parent_props = mol_properties(charge_parent)
        exact_props = mol_properties(cleaned)
        if fragments > 1 and parent_key != key:
            relation = "salt_or_multicomponent_form"
        elif parent_key != key:
            relation = "charged_or_protonation_form"
        elif centers:
            relation = "stereochemically_defined_parent"
        else:
            relation = "parent_form"
        record = {
            "exact_key": key,
            "parent_key": parent_key,
            "relation": relation,
            "exact_smiles": exact_smiles,
            "exact_inchi": exact_inchi,
            "exact_props": exact_props,
            "parent_smiles": parent_smiles,
            "parent_inchi": parent_inchi,
            "parent_props": parent_props,
            "connectivity_key": parent_key[:14],
            "fragment_count": fragments,
            "stereocenter_count": len(centers),
            "unassigned_stereocenter_count": sum(
                1 for _atom, label in centers if label == "?"
            ),
            "preferred_name": preferred,
        }
        return key, record, "", ""
    except Exception as exc:
        return key, None, type(exc).__name__, str(exc)


def canonicalize_candidates(
    con: sqlite3.Connection,
    master_by_key: dict[str, str],
    form_by_key: dict[str, str],
    max_compound: int,
    max_form: int,
) -> dict:
    referenced_total = con.execute(
        "SELECT COUNT(*) FROM referenced_candidate"
    ).fetchone()[0]
    records: dict[str, dict] = {}
    failures: Counter[str] = Counter()
    derived_parent_structures: dict[str, dict] = {}
    candidate_rows = con.execute(
        """
        SELECT r.exact_key,c.connectivity_key,c.smiles,c.inchi,c.formula,
               c.molecular_weight,c.formal_charge,c.preferred_name
        FROM referenced_candidate r
        LEFT JOIN candidate_structure c ON c.exact_key=r.exact_key
        ORDER BY r.exact_key
        """
    )
    worker_count = min(8, max(2, (os.cpu_count() or 2) // 2))
    with ProcessPoolExecutor(max_workers=worker_count) as pool:
        results = pool.map(
            canonicalize_candidate_worker,
            (tuple(row) for row in candidate_rows),
            chunksize=200,
        )
        for index, (key, record, error_type, error_message) in enumerate(results, 1):
            if record is not None:
                records[key] = record
                parent_key = record["parent_key"]
                derived_parent_structures[parent_key] = {
                    "smiles": record["parent_smiles"],
                    "inchi": record["parent_inchi"],
                    "properties": record["parent_props"],
                    "preferred_name": record["preferred_name"],
                }
            else:
                failures[error_type or "unknown_error"] += 1
                con.execute(
                    "UPDATE disposition SET action='review',reason=? WHERE exact_key=?",
                    (f"canonicalization_failed:{error_message}", key),
                )
                con.execute(
                    "DELETE FROM referenced_candidate WHERE exact_key=?", (key,)
                )
            if index % 10000 == 0:
                con.commit()
                update_progress(
                    "canonicalize_parent_form",
                    candidates_processed=index,
                    candidates_total=referenced_total,
                    failures=sum(failures.values()),
                    workers=worker_count,
                )
    con.commit()

    new_parent_keys = sorted(
        {
            record["parent_key"]
            for record in records.values()
            if record["parent_key"] not in master_by_key
        }
    )
    parent_ids = dict(master_by_key)
    for offset, key in enumerate(new_parent_keys, 1):
        parent_ids[key] = f"HMPD-CMPD-{max_compound + offset:07d}"
    new_form_keys = sorted(key for key in records if key not in form_by_key)
    form_ids = dict(form_by_key)
    for offset, key in enumerate(new_form_keys, 1):
        form_ids[key] = f"HMPD-FORM-{max_form + offset:07d}"

    for key, record in records.items():
        compound_id = parent_ids[record["parent_key"]]
        form_id = form_ids[key]
        props = record["exact_props"]
        parent_props = record["parent_props"]
        con.execute(
            """
            INSERT OR REPLACE INTO canonical_form
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """,
            (
                key,
                record["parent_key"],
                compound_id,
                form_id,
                record["relation"],
                record["exact_smiles"],
                record["exact_inchi"],
                props["formula"],
                props["mw"],
                props["charge"],
                record["parent_smiles"],
                record["parent_inchi"],
                parent_props["formula"],
                parent_props["mw"],
                parent_props["charge"],
                record["connectivity_key"],
                record["fragment_count"],
                record["stereocenter_count"],
                record["unassigned_stereocenter_count"],
                "ok",
                "",
            ),
        )
        con.execute(
            """
            UPDATE identity SET compound_internal_id=?,compound_form_id=?
            WHERE exact_key=? AND safe_status='exact_new'
            """,
            (compound_id, form_id, key),
        )
    con.commit()
    return {
        "referenced_new_exact_structures": referenced_total,
        "canonicalized_exact_structures": len(records),
        "canonicalization_failures": dict(failures),
        "new_parent_compounds": len(new_parent_keys),
        "new_exact_forms": len(new_form_keys),
        "parent_ids": parent_ids,
        "derived_parent_structures": derived_parent_structures,
        "records": records,
    }


def bool_text(value: bool) -> str:
    return "1" if value else "0"


def structural_class(props: dict[str, object], mol: Chem.Mol) -> str:
    atoms = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if not atoms.intersection({"C"}):
        return "inorganic"
    if atoms.intersection(
        {"Fe", "Zn", "Mg", "Mn", "Co", "Ni", "Cu", "Ca", "Na", "K", "Li"}
    ):
        return "metal_containing"
    if int(props["rings"]) >= 2:
        return "polycyclic_organic"
    if int(props["rings"]) == 1:
        return "cyclic_organic"
    return "acyclic_organic"


def build_master_and_forms(
    con: sqlite3.Connection,
    master_by_id: dict[str, dict[str, str]],
    canonical: dict,
    baseline_form_ids: set[str] | None = None,
) -> tuple[dict[str, dict[str, str]], dict]:
    CANDIDATE.mkdir(parents=True, exist_ok=True)
    master_fields: list[str]
    with MASTER0.open("rt", encoding="utf-8-sig", newline="") as handle:
        master_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    for field in (
        "identity_refresh_v12",
        "new_source_record_count_v12",
        "release_version_v12",
        "release_date_v12",
    ):
        if field not in master_fields:
            master_fields.append(field)

    updates: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"sources": set(), "pubchem": set()}
    )
    for source_db, source_id, compound_id in con.execute(
        """
        SELECT source_database,source_id,compound_internal_id
        FROM identity
        WHERE safe_status='exact_existing' AND compound_internal_id<>''
        """
    ):
        updates[compound_id]["sources"].add(source_db)
        if source_db == "PubChem BioAssay":
            updates[compound_id]["pubchem"].add(source_id)

    con.executescript(
        """
        DROP TABLE IF EXISTS compound_source_agg;
        CREATE TABLE compound_source_agg AS
        SELECT cf.compound_internal_id,
               GROUP_CONCAT(DISTINCT cs.source_database) AS sources,
               COUNT(*) AS source_record_count,
               GROUP_CONCAT(
                 DISTINCT CASE WHEN cs.source_database='PubChem BioAssay'
                               THEN cs.source_id END
               ) AS pubchem_cids
        FROM candidate_source cs
        JOIN canonical_form cf ON cf.exact_key=cs.exact_key
        GROUP BY cf.compound_internal_id;
        CREATE UNIQUE INDEX idx_compound_source_agg_id
        ON compound_source_agg(compound_internal_id);
        DROP TABLE IF EXISTS compound_form_count;
        CREATE TABLE compound_form_count AS
        SELECT compound_internal_id,COUNT(*) AS form_count
        FROM canonical_form GROUP BY compound_internal_id;
        CREATE UNIQUE INDEX idx_compound_form_count_id
        ON compound_form_count(compound_internal_id);
        """
    )
    con.commit()
    parent_ids = canonical["parent_ids"]
    derived = canonical["derived_parent_structures"]
    new_master_count = 0
    new_master_total = len(derived)
    aggregate_rows = iter(
        con.execute(
            """
            SELECT s.compound_internal_id,COALESCE(s.sources,''),
                   s.source_record_count,COALESCE(s.pubchem_cids,''),
                   COALESCE(f.form_count,0)
            FROM compound_source_agg s
            LEFT JOIN compound_form_count f USING(compound_internal_id)
            ORDER BY s.compound_internal_id
            """
        )
    )
    aggregate_row = next(aggregate_rows, None)
    all_master: dict[str, dict[str, str]] = {}
    master_path = CANDIDATE / "small_molecule_master_v1_2.tsv"
    with atomic_text(master_path, master_fields) as writer:
        with MASTER0.open("rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                compound_id = row["compound_internal_id"]
                update = updates.get(compound_id)
                if update:
                    sources = sorted(
                        set(split_values(row["source_databases"]))
                        | update["sources"]
                    )
                    pubchem = sorted(
                        set(split_values(row["pubchem_cids"]))
                        | update["pubchem"],
                        key=lambda value: int(value),
                    )
                    row["source_databases"] = ";".join(sources)
                    row["source_database_count"] = str(len(sources))
                    row["pubchem_cids"] = ";".join(pubchem)
                    row["identity_refresh_v12"] = "existing_identity_enriched"
                    row["new_source_record_count_v12"] = str(
                        len(update["sources"]) + len(update["pubchem"])
                    )
                else:
                    row["identity_refresh_v12"] = "baseline_preserved"
                    row["new_source_record_count_v12"] = "0"
                row["release_version_v12"] = "small-molecule-v1.2"
                row["release_date_v12"] = RELEASE_DATE
                all_master[compound_id] = {
                    "preferred_name": row["preferred_name"],
                    "compound_scope_status": row["compound_scope_status"],
                    "standard_inchikey": row["standard_inchikey"],
                }
                writer.writerow({field: row.get(field, "") for field in master_fields})
        for parent_key, compound_id in sorted(parent_ids.items()):
            if compound_id in master_by_id:
                continue
            structure = derived.pop(parent_key)
            mol = Chem.MolFromSmiles(structure["smiles"])
            props = structure["properties"]
            while aggregate_row is not None and aggregate_row[0] < compound_id:
                aggregate_row = next(aggregate_rows, None)
            if aggregate_row is not None and aggregate_row[0] == compound_id:
                (
                    _aggregate_id,
                    source_text,
                    source_record_count,
                    pubchem_text,
                    form_count,
                ) = aggregate_row
                aggregate_row = next(aggregate_rows, None)
            else:
                source_text, source_record_count, pubchem_text, form_count = (
                    "",
                    0,
                    "",
                    0,
                )
            sources = sorted(split_values(source_text))
            pubchem = sorted(
                split_values(pubchem_text),
                key=lambda value: int(value),
            )
            preferred = clean(structure["preferred_name"]) or f"Structure {parent_key}"
            mw = float(props["mw"])
            heavy = int(props["heavy"])
            charge = int(props["charge"])
            rotatable = int(props["rotatable"])
            scope = (
                "core"
                if 50 <= mw <= 1000 and 2 <= heavy <= 100 and abs(charge) <= 4
                else "extended"
            )
            row = {field: "" for field in master_fields}
            row.update(
                {
                    "compound_internal_id": compound_id,
                    "preferred_name": preferred,
                    "preferred_name_source": "V6.1 resolved source",
                    "compound_scope_status": scope,
                    "compound_scope_class": "structure_resolved_small_molecule",
                    "identity_confidence": (
                        "high" if len(sources) >= 2 else "medium"
                    ),
                    "standard_smiles": structure["smiles"],
                    "standard_inchi": structure["inchi"],
                    "standard_inchikey": parent_key,
                    "identity_key": f"IK:{parent_key}",
                    "connectivity_key": parent_key[:14],
                    "molecular_formula": props["formula"],
                    "molecular_weight": f"{mw:.6g}",
                    "exact_mass": f"{float(props['exact_mass']):.6g}",
                    "xlogp": f"{float(props['xlogp']):.6g}",
                    "tpsa": f"{float(props['tpsa']):.6g}",
                    "hbond_donor_count": props["hbd"],
                    "hbond_acceptor_count": props["hba"],
                    "rotatable_bond_count": rotatable,
                    "formal_charge": charge,
                    "heavy_atom_count": heavy,
                    "ring_count": props["rings"],
                    "amide_bond_count": props["amide"],
                    "computed_structural_class": structural_class(props, mol),
                    "molecule_types": "small_molecule",
                    "form_count": form_count,
                    "source_record_count": source_record_count,
                    "source_database_count": len(sources),
                    "source_databases": ";".join(sources),
                    "pubchem_cids": ";".join(pubchem),
                    "record_qc_status": "ok",
                    "qc_notes": "",
                    "release_version": "small-molecule-v1.2",
                    "release_date": RELEASE_DATE,
                    "standardization_policy": "V6.1 full-InChIKey; RDKit cleanup/fragment/charge parent",
                    "identity_refresh_v12": "new_canonical_parent",
                    "new_source_record_count_v12": source_record_count,
                    "release_version_v12": "small-molecule-v1.2",
                    "release_date_v12": RELEASE_DATE,
                }
            )
            row = {key: clean(value) for key, value in row.items()}
            all_master[compound_id] = {
                "preferred_name": row["preferred_name"],
                "compound_scope_status": row["compound_scope_status"],
                "standard_inchikey": row["standard_inchikey"],
            }
            writer.writerow({field: row.get(field, "") for field in master_fields})
            new_master_count += 1
            if new_master_count % 50000 == 0:
                update_progress(
                    "build_master_and_forms",
                    new_master_rows_written=new_master_count,
                    new_master_rows_total=new_master_total,
                )

    with FORM0.open("rt", encoding="utf-8-sig", newline="") as handle:
        form_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    for field in (
        "parent_inchikey_v12",
        "parent_relation_v12",
        "standardization_status_v12",
        "release_version_v12",
    ):
        if field not in form_fields:
            form_fields.append(field)
    form_path = CANDIDATE / "compound_form_hierarchy_v1_2.tsv"
    appended_forms = 0
    baseline_form_ids = baseline_form_ids or set()
    con.executescript(
        """
        DROP TABLE IF EXISTS candidate_source_agg;
        CREATE TABLE candidate_source_agg AS
        SELECT exact_key,
               GROUP_CONCAT(DISTINCT source_database) AS sources,
               GROUP_CONCAT(source_database || ':' || source_id) AS source_ids,
               COUNT(*) AS source_record_count
        FROM candidate_source
        GROUP BY exact_key;
        CREATE UNIQUE INDEX idx_candidate_source_agg_key
        ON candidate_source_agg(exact_key);
        """
    )
    con.commit()
    with atomic_text(form_path, form_fields) as writer:
        with FORM0.open("rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row.update(
                    {
                        "parent_inchikey_v12": all_master[
                            row["compound_internal_id"]
                        ]["standard_inchikey"],
                        "parent_relation_v12": "baseline_form",
                        "standardization_status_v12": "baseline_preserved",
                        "release_version_v12": "small-molecule-v1.2",
                    }
                )
                writer.writerow({field: row.get(field, "") for field in form_fields})
        for row in con.execute(
            """
            SELECT exact_key,parent_key,compound_internal_id,compound_form_id,
                   parent_relation,exact_smiles,exact_inchi,exact_formula,
                   exact_mw,exact_charge,standardization_status,
                   standardization_notes,COALESCE(a.sources,''),
                   COALESCE(a.source_ids,''),COALESCE(a.source_record_count,0)
            FROM canonical_form c
            LEFT JOIN candidate_source_agg a USING(exact_key)
            ORDER BY compound_form_id
            """
        ):
            (
                exact_key,
                parent_key,
                compound_id,
                form_id,
                relation,
                smiles,
                inchi,
                formula,
                mw,
                charge,
                status,
                notes,
                source_text,
                source_id_text,
                source_record_count,
            ) = row
            if form_id in baseline_form_ids:
                continue
            sources = sorted(split_values(source_text))
            source_ids = sorted(split_values(source_id_text))
            writer.writerow(
                {
                    "compound_form_id": form_id,
                    "compound_internal_id": compound_id,
                    "form_type": relation,
                    "exact_smiles": smiles,
                    "exact_inchi": inchi,
                    "exact_inchikey": exact_key,
                    "exact_identity_key": f"IK:{exact_key}",
                    "molecular_formula": formula,
                    "molecular_weight": mw,
                    "formal_charge": charge,
                    "source_record_count": source_record_count,
                    "source_databases": ";".join(sources),
                    "source_compound_ids": ";".join(source_ids),
                    "record_qc_status": "ok",
                    "qc_notes": notes,
                    "release_version_v11": "small-molecule-v1.2",
                    "parent_inchikey_v12": parent_key,
                    "parent_relation_v12": relation,
                    "standardization_status_v12": status,
                    "release_version_v12": "small-molecule-v1.2",
                }
            )
            appended_forms += 1
    return all_master, {
        "master_rows": len(all_master),
        "new_master_rows": new_master_count,
        "new_form_rows": appended_forms,
    }


def evidence_directness(tier: str) -> str:
    return {"BE1": "direct_structural", "BE2": "direct", "BE3": "functional"}.get(
        tier, ""
    )


def source_id_type(source: str) -> str:
    return {
        "PubChem BioAssay": "PubChem CID",
        "PDBe": "PDB CCD",
        "PDBbind": "PDBbind complex",
        "BRENDA": "BRENDA compound name",
    }.get(source, "source identifier")


def build_evidence_layers(
    con: sqlite3.Connection,
    all_master: dict[str, dict[str, str]],
    form_by_key: dict[str, str],
) -> dict:
    with gzip.open(EVIDENCE0, "rt", encoding="utf-8-sig", newline="") as handle:
        evidence_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    for field in (
        "source_evidence_id_v61",
        "evidence_action_v61",
        "identity_resolution_v61",
        "release_version_v61",
    ):
        if field not in evidence_fields:
            evidence_fields.append(field)

    with gzip.open(NEGATIVE0, "rt", encoding="utf-8-sig", newline="") as handle:
        negative_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    for field in ("identity_resolution_v61", "release_action_v61", "release_version_v61"):
        if field not in negative_fields:
            negative_fields.append(field)

    with gzip.open(REVIEW0, "rt", encoding="utf-8-sig", newline="") as handle:
        review_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    for field in (
        "identity_resolution_v61",
        "resolved_compound_internal_id_v61",
        "resolved_compound_form_id_v61",
        "release_disposition_v61",
        "release_reason_v61",
        "release_version_v61",
    ):
        if field not in review_fields:
            review_fields.append(field)

    conflict_fields = [
        "source_evidence_id",
        "target_uniprot_id",
        "source_database",
        "source_record_id",
        "compound_source_id",
        "compound_name",
        "identity_resolution",
        "identity_safe_status",
        "release_disposition",
        "reason",
    ]
    duplicate_fields = [
        "source_evidence_id",
        "target_uniprot_id",
        "compound_internal_id",
        "source_database",
        "source_record_id",
        "duplicate_of_evidence_id",
        "duplicate_key_sha256",
    ]
    sample_fields = [
        "sample_rank",
        "source_evidence_id",
        "target_uniprot_id",
        "approved_symbol",
        "compound_internal_id",
        "compound_form_id",
        "compound_name",
        "source_database",
        "source_record_id",
        "source_url",
        "evidence_tier",
        "activity_type",
        "activity_value",
        "activity_unit",
        "standard_value_nM",
        "verification_status",
        "verification_notes",
    ]

    evidence_path = CANDIDATE / "binding_evidence_master_v6_1.tsv.gz"
    negative_path = CANDIDATE / "negative_binding_evidence_v1_1.tsv.gz"
    review_path = CANDIDATE / "binding_evidence_review_queue_v6_1.tsv.gz"
    conflict_path = CANDIDATE / "identity_conflict_audit_v1_0.tsv.gz"
    duplicate_path = CANDIDATE / "duplicate_evidence_audit_v1_0.tsv.gz"

    counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    new_evidence_digests: dict[str, str] = {}
    new_negative_digests: set[str] = set()
    sample_heap: list[tuple[int, dict[str, str]]] = []
    sample_size = 500

    with (
        atomic_text(evidence_path, evidence_fields, gzipped=True) as evidence_writer,
        atomic_text(negative_path, negative_fields, gzipped=True) as negative_writer,
        atomic_text(review_path, review_fields, gzipped=True) as review_writer,
        atomic_text(conflict_path, conflict_fields, gzipped=True) as conflict_writer,
        atomic_text(duplicate_path, duplicate_fields, gzipped=True) as duplicate_writer,
    ):
        with gzip.open(EVIDENCE0, "rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row.update(
                    {
                        "source_evidence_id_v61": "",
                        "evidence_action_v61": "baseline_preserved",
                        "identity_resolution_v61": "baseline_canonical",
                        "release_version_v61": RELEASE_VERSION,
                    }
                )
                evidence_writer.writerow(
                    {field: row.get(field, "") for field in evidence_fields}
                )
                counts["baseline_positive"] += 1
                if row.get("default_release_inclusion") == "1":
                    tier = row["evidence_tier"]
                    value = number(row["standard_value_nM"])
                    con.execute(
                        "INSERT INTO pair_raw VALUES (?,?,?,?,?,?,?,?,?,?)",
                        (
                            row["target_uniprot_id"],
                            row["compound_internal_id"],
                            row["approved_symbol"],
                            tier,
                            TIER_RANK.get(tier, 9),
                            row["source_database"],
                            value,
                            row["compound_form_id"],
                            row["activity_type"],
                            "baseline",
                        ),
                    )
                    if value is not None:
                        con.execute(
                            "INSERT INTO quant_raw VALUES (?,?,?,?,?)",
                            (
                                row["target_uniprot_id"],
                                row["compound_internal_id"],
                                row["activity_type"],
                                value,
                                row["source_database"],
                            ),
                        )

        with gzip.open(NEGATIVE0, "rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row.update(
                    {
                        "identity_resolution_v61": "baseline_canonical",
                        "release_action_v61": "baseline_negative_preserved",
                        "release_version_v61": RELEASE_VERSION,
                    }
                )
                negative_writer.writerow(
                    {field: row.get(field, "") for field in negative_fields}
                )
                counts["baseline_negative"] += 1

        con.commit()
        with gzip.open(REVIEW0, "rt", encoding="utf-8-sig", newline="") as handle:
            for row_number, row in enumerate(
                csv.DictReader(handle, delimiter="\t"), 1
            ):
                disposition = con.execute(
                    """
                    SELECT action,source_type,source_id,identity_safe_status,
                           identity_resolution,exact_key,existing_compound_id,reason
                    FROM disposition WHERE source_evidence_id=?
                    """,
                    (row["source_evidence_id"],),
                ).fetchone()
                if disposition is None:
                    raise RuntimeError(
                        f"Missing disposition: {row['source_evidence_id']}"
                    )
                (
                    action,
                    source_type,
                    source_id,
                    safe,
                    resolution,
                    exact_key,
                    existing,
                    reason,
                ) = disposition
                compound_id = existing
                form_id = ""
                if safe == "exact_new" and exact_key:
                    mapping = con.execute(
                        """
                        SELECT compound_internal_id,compound_form_id
                        FROM canonical_form WHERE exact_key=?
                        """,
                        (exact_key,),
                    ).fetchone()
                    if mapping:
                        compound_id, form_id = mapping
                    else:
                        action, reason = "review", "canonical_form_not_available"
                elif safe == "exact_existing" and exact_key:
                    compound_id = existing
                    form_id = form_by_key.get(exact_key, "")
                if action in {"positive", "negative"} and (
                    not compound_id or compound_id not in all_master
                ):
                    action, reason = "review", "compound_foreign_key_unavailable"

                digest_parts = [
                    row["target_uniprot_id"],
                    compound_id,
                    row["source_database"],
                    row["source_record_id"],
                    row["evidence_type"],
                    row["activity_type"],
                    row["activity_relation"],
                    row["activity_value"],
                    row["activity_unit"],
                    row["pdb_ids"],
                    row["binding_site_residues"],
                    row["activity_outcome"],
                ]
                digest = hashlib.sha256(
                    "\x1f".join(digest_parts).encode("utf-8")
                ).hexdigest()

                if action == "positive":
                    evidence_id = stable_id(
                        "BE61",
                        row["source_evidence_id"],
                        row["target_uniprot_id"],
                        compound_id,
                    )
                    if digest in new_evidence_digests:
                        action = "duplicate"
                        duplicate_writer.writerow(
                            {
                                "source_evidence_id": row["source_evidence_id"],
                                "target_uniprot_id": row["target_uniprot_id"],
                                "compound_internal_id": compound_id,
                                "source_database": row["source_database"],
                                "source_record_id": row["source_record_id"],
                                "duplicate_of_evidence_id": new_evidence_digests[
                                    digest
                                ],
                                "duplicate_key_sha256": digest,
                            }
                        )
                        counts["duplicate"] += 1
                    else:
                        new_evidence_digests[digest] = evidence_id
                        master = all_master[compound_id]
                        output = {
                            "evidence_id": evidence_id,
                            "target_uniprot_id": row["target_uniprot_id"],
                            "approved_symbol": row["approved_symbol"],
                            "protein_release_tier": "",
                            "protein_website_default": "1",
                            "compound_id": row["compound_source_id"],
                            "compound_id_type": source_id_type(
                                row["source_database"]
                            ),
                            "compound_name": row["compound_name"],
                            "small_molecule_scope_status": master[
                                "compound_scope_status"
                            ],
                            "is_core_small_molecule": bool_text(
                                master["compound_scope_status"] == "core"
                            ),
                            "evidence_tier": row["evidence_tier"],
                            "evidence_type": row["evidence_type"],
                            "evidence_directness": evidence_directness(
                                row["evidence_tier"]
                            ),
                            "target_assignment_status": "single_uniprot_mapping",
                            "relationship_context_id": row["source_record_id"],
                            "activity_type": row["activity_type"],
                            "activity_relation": row["activity_relation"],
                            "activity_value": row["activity_value"],
                            "activity_unit": row["activity_unit"],
                            "standard_value_nM": row["standard_value_nM"],
                            "assay_or_mechanism": row["assay_or_mechanism"],
                            "pdb_ids": row["pdb_ids"],
                            "binding_site_residues": row[
                                "binding_site_residues"
                            ],
                            "ligand_het_id": (
                                row["compound_source_id"].split(":", 1)[-1]
                                if row["source_database"] == "PDBe"
                                else ""
                            ),
                            "pubmed_ids": "",
                            "doi": "",
                            "source_database": row["source_database"],
                            "source_version": row["source_version"],
                            "source_record_id": row["source_record_id"],
                            "source_url": row["source_url"],
                            "record_qc_status": row["record_qc_status"],
                            "default_release_inclusion": "1",
                            "originating_dataset": row["incremental_source"],
                            "compound_internal_id": compound_id,
                            "compound_form_id": form_id,
                            "compound_mapping_status": (
                                "mapped_core"
                                if master["compound_scope_status"] == "core"
                                else "mapped_extended"
                            ),
                            "compound_scope_status_v1": master[
                                "compound_scope_status"
                            ],
                            "default_release_inclusion_v4_2": "0",
                            "compound_master_version": "v1.2",
                            "activity_outcome_v60": row["activity_outcome"],
                            "incremental_source_v60": row["incremental_source"],
                            "duplicate_status_v60": "new_unique_v61",
                            "identity_status_v60": resolution,
                            "conflict_status_v60": "no_blocking_conflict",
                            "release_version_v60": "v6.1",
                            "source_evidence_id_v61": row["source_evidence_id"],
                            "evidence_action_v61": "promoted_from_v60_review",
                            "identity_resolution_v61": resolution,
                            "release_version_v61": RELEASE_VERSION,
                        }
                        evidence_writer.writerow(
                            {
                                field: output.get(field, "")
                                for field in evidence_fields
                            }
                        )
                        counts["new_positive"] += 1
                        source_counts[f"{row['source_database']}:positive"] += 1
                        tier = row["evidence_tier"]
                        value = number(row["standard_value_nM"])
                        con.execute(
                            "INSERT INTO pair_raw VALUES (?,?,?,?,?,?,?,?,?,?)",
                            (
                                row["target_uniprot_id"],
                                compound_id,
                                row["approved_symbol"],
                                tier,
                                TIER_RANK.get(tier, 9),
                                row["source_database"],
                                value,
                                form_id,
                                row["activity_type"],
                                "v61",
                            ),
                        )
                        if value is not None:
                            con.execute(
                                "INSERT INTO quant_raw VALUES (?,?,?,?,?)",
                                (
                                    row["target_uniprot_id"],
                                    compound_id,
                                    row["activity_type"],
                                    value,
                                    row["source_database"],
                                ),
                            )
                        if (
                            tier == "BE1"
                            and row["pdb_ids"]
                            and row["binding_site_residues"]
                        ):
                            site_key = hashlib.sha256(
                                "\x1f".join(
                                    [
                                        row["target_uniprot_id"],
                                        compound_id,
                                        row["pdb_ids"],
                                        row["binding_site_residues"],
                                        row["source_database"],
                                    ]
                                ).encode("utf-8")
                            ).hexdigest()
                            con.execute(
                                """
                                INSERT OR IGNORE INTO site_new
                                VALUES (?,?,?,?,?,?,?,?,?)
                                """,
                                (
                                    site_key,
                                    evidence_id,
                                    row["target_uniprot_id"],
                                    compound_id,
                                    form_id,
                                    row["pdb_ids"],
                                    row["binding_site_residues"],
                                    row["source_database"],
                                    row["compound_source_id"],
                                ),
                            )
                        sample_hash = int(
                            hashlib.sha256(
                                row["source_evidence_id"].encode("utf-8")
                            ).hexdigest(),
                            16,
                        )
                        sample_row = {
                            "source_evidence_id": row["source_evidence_id"],
                            "target_uniprot_id": row["target_uniprot_id"],
                            "approved_symbol": row["approved_symbol"],
                            "compound_internal_id": compound_id,
                            "compound_form_id": form_id,
                            "compound_name": row["compound_name"],
                            "source_database": row["source_database"],
                            "source_record_id": row["source_record_id"],
                            "source_url": row["source_url"],
                            "evidence_tier": row["evidence_tier"],
                            "activity_type": row["activity_type"],
                            "activity_value": row["activity_value"],
                            "activity_unit": row["activity_unit"],
                            "standard_value_nM": row["standard_value_nM"],
                            "verification_status": "pending_manual_review",
                            "verification_notes": "",
                        }
                        if len(sample_heap) < sample_size:
                            heapq.heappush(sample_heap, (-sample_hash, sample_row))
                        elif sample_hash < -sample_heap[0][0]:
                            heapq.heapreplace(
                                sample_heap, (-sample_hash, sample_row)
                            )

                if action == "negative":
                    if digest in new_negative_digests:
                        action = "duplicate"
                        counts["duplicate_negative"] += 1
                    else:
                        new_negative_digests.add(digest)
                        output = dict(row)
                        output.update(
                            {
                                "compound_internal_id": compound_id,
                                "compound_form_id": form_id,
                                "compound_mapping_status": (
                                    "mapped_core"
                                    if all_master[compound_id][
                                        "compound_scope_status"
                                    ]
                                    == "core"
                                    else "mapped_extended"
                                ),
                                "review_bucket": "negative_evidence",
                                "review_reason": "reported inactive; identity resolved in V6.1",
                                "release_version": RELEASE_VERSION,
                                "identity_resolution_v61": resolution,
                                "release_action_v61": "negative_identity_resolved",
                                "release_version_v61": RELEASE_VERSION,
                            }
                        )
                        negative_writer.writerow(
                            {
                                field: output.get(field, "")
                                for field in negative_fields
                            }
                        )
                        counts["new_negative"] += 1
                        source_counts[f"{row['source_database']}:negative"] += 1

                if action not in {"positive", "negative"}:
                    output = dict(row)
                    output.update(
                        {
                            "identity_resolution_v61": resolution,
                            "resolved_compound_internal_id_v61": compound_id,
                            "resolved_compound_form_id_v61": form_id,
                            "release_disposition_v61": action,
                            "release_reason_v61": reason,
                            "release_version_v61": RELEASE_VERSION,
                            "release_version": RELEASE_VERSION,
                        }
                    )
                    review_writer.writerow(
                        {field: output.get(field, "") for field in review_fields}
                    )
                    counts[f"retained_{action}"] += 1
                    source_counts[f"{row['source_database']}:{action}"] += 1
                    if action in {"review", "excluded"}:
                        conflict_writer.writerow(
                            {
                                "source_evidence_id": row["source_evidence_id"],
                                "target_uniprot_id": row["target_uniprot_id"],
                                "source_database": row["source_database"],
                                "source_record_id": row["source_record_id"],
                                "compound_source_id": row["compound_source_id"],
                                "compound_name": row["compound_name"],
                                "identity_resolution": resolution,
                                "identity_safe_status": safe,
                                "release_disposition": action,
                                "reason": reason,
                            }
                        )

                if row_number % 100000 == 0:
                    con.commit()
                    update_progress(
                        "rebuild_evidence_layers",
                        review_rows_processed=row_number,
                        evidence_counts=dict(counts),
                    )
        con.commit()

    sample_path = CANDIDATE / "manual_validation_sample_pending_v6_1.tsv"
    sample_rows = [
        item[1]
        for item in sorted(sample_heap, key=lambda item: -item[0])
    ]
    with atomic_text(sample_path, sample_fields) as writer:
        for rank, row in enumerate(sample_rows, 1):
            row["sample_rank"] = rank
            writer.writerow(row)
    return {
        "counts": dict(counts),
        "source_counts": dict(source_counts),
        "manual_validation_sample_rows": len(sample_rows),
    }


def build_pair_site_protein_outputs(
    con: sqlite3.Connection, all_master: dict[str, dict[str, str]]
) -> dict:
    con.executescript(
        """
        CREATE INDEX idx_pair_raw_key ON pair_raw(target,compound);
        CREATE INDEX idx_quant_raw_key
          ON quant_raw(target,compound,activity_type,value_nm);
        CREATE TABLE pair_agg AS
        SELECT target,compound,
               MAX(approved_symbol) AS approved_symbol,
               MIN(tier_rank) AS best_rank,
               SUM(CASE WHEN tier='BE1' THEN 1 ELSE 0 END) AS be1,
               SUM(CASE WHEN tier='BE2' THEN 1 ELSE 0 END) AS be2,
               SUM(CASE WHEN tier='BE3' THEN 1 ELSE 0 END) AS be3,
               COUNT(*) AS evidence_count,
               COUNT(DISTINCT source) AS source_count,
               GROUP_CONCAT(DISTINCT source) AS sources,
               MIN(value_nm) AS best_value,
               COUNT(DISTINCT CASE WHEN form_id<>'' THEN form_id END) AS form_count,
               SUM(CASE WHEN origin='v61' THEN 1 ELSE 0 END) AS new_count,
               GROUP_CONCAT(DISTINCT CASE WHEN origin='v61' THEN source END)
                 AS new_sources
        FROM pair_raw
        GROUP BY target,compound;
        CREATE UNIQUE INDEX idx_pair_agg_key ON pair_agg(target,compound);
        CREATE TABLE target_agg AS
        SELECT target,
               COUNT(DISTINCT compound) AS compound_count,
               COUNT(DISTINCT CASE WHEN best_rank=1 THEN compound END) AS be1_compounds,
               COUNT(DISTINCT CASE WHEN best_rank=2 THEN compound END) AS be2_compounds,
               COUNT(DISTINCT CASE WHEN best_rank=3 THEN compound END) AS be3_compounds,
               SUM(evidence_count) AS evidence_count,
               COUNT(DISTINCT source_count) AS source_count_placeholder
        FROM pair_agg GROUP BY target;
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
        "new_evidence_count_v61",
        "new_sources_v61",
        "pair_action_v61",
        "release_version_v61",
    ]
    pair_path = CANDIDATE / "protein_compound_summary_v6_1.tsv.gz"
    pair_rows = 0
    with atomic_text(pair_path, pair_fields, gzipped=True) as writer:
        for row in con.execute(
            """
            SELECT target,compound,approved_symbol,best_rank,be1,be2,be3,
                   evidence_count,source_count,sources,best_value,form_count,
                   new_count,new_sources
            FROM pair_agg ORDER BY target,compound
            """
        ):
            (
                target,
                compound_id,
                symbol,
                rank,
                be1,
                be2,
                be3,
                evidence_count,
                source_count,
                sources,
                best_value,
                form_count,
                new_count,
                new_sources,
            ) = row
            writer.writerow(
                {
                    "target_uniprot_id": target,
                    "approved_symbol": symbol,
                    "compound_internal_id": compound_id,
                    "preferred_name": all_master[compound_id]["preferred_name"],
                    "best_binding_evidence_level": {
                        1: "BE1",
                        2: "BE2",
                        3: "BE3",
                    }.get(rank, ""),
                    "BE1_evidence_count": be1,
                    "BE2_evidence_count": be2,
                    "BE3_evidence_count": be3,
                    "binding_evidence_count": evidence_count,
                    "independent_source_count": source_count,
                    "independent_sources": ";".join(
                        sorted(split_values(sources))
                    ),
                    "best_standard_value_nM": (
                        "" if best_value is None else best_value
                    ),
                    "compound_form_count": form_count,
                    "compound_master_version": "v1.2",
                    "new_evidence_count_v61": new_count,
                    "new_sources_v61": ";".join(
                        sorted(split_values(new_sources))
                    ),
                    "pair_action_v61": (
                        "new_or_strengthened_from_identity_refresh"
                        if new_count
                        else "baseline_unchanged"
                    ),
                    "release_version_v61": RELEASE_VERSION,
                }
            )
            pair_rows += 1

    with gzip.open(SITE0, "rt", encoding="utf-8-sig", newline="") as handle:
        site_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    for field in ("site_action_v61", "release_version_v61"):
        if field not in site_fields:
            site_fields.append(field)
    site_path = CANDIDATE / "binding_site_instances_v6_1.tsv.gz"
    site_rows = 0
    with atomic_text(site_path, site_fields, gzipped=True) as writer:
        with gzip.open(SITE0, "rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row.update(
                    {
                        "site_action_v61": "baseline_preserved",
                        "release_version_v61": RELEASE_VERSION,
                    }
                )
                writer.writerow({field: row.get(field, "") for field in site_fields})
                site_rows += 1
        for row in con.execute(
            """
            SELECT site_key,evidence_id,target,compound,form_id,pdb_ids,
                   residues,source,source_compound_id
            FROM site_new ORDER BY target,compound,pdb_ids,residues
            """
        ):
            (
                site_key,
                evidence_id,
                target,
                compound,
                form_id,
                pdb_ids,
                residues,
                source,
                source_compound_id,
            ) = row
            writer.writerow(
                {
                    "binding_site_instance_id": f"BSI61_{site_key[:20].upper()}",
                    "evidence_id": evidence_id,
                    "target_uniprot_id": target,
                    "compound_id": source_compound_id,
                    "site_type": "experimental_structure_residue_contact",
                    "pdb_ids": pdb_ids,
                    "residue_or_site_description": residues,
                    "site_compound_specificity": "compound_specific",
                    "source_database": source,
                    "record_qc_status": "ok",
                    "compound_internal_id": compound,
                    "compound_form_id": form_id,
                    "compound_mapping_status": "mapped_v61",
                    "compound_master_version": "v1.2",
                    "residue_index_type_v60": "",
                    "pdb_chain_ids_v60": "",
                    "site_merge_status_v60": "new_v61_site",
                    "release_version_v60": "v6.1",
                    "site_action_v61": "new_identity_resolved_site",
                    "release_version_v61": RELEASE_VERSION,
                }
            )
            site_rows += 1

    with PROTEIN0.open("rt", encoding="utf-8-sig", newline="") as handle:
        protein_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    for field in (
        "canonical_core_small_molecule_count_v61",
        "canonical_BE1_compound_count_v61",
        "canonical_BE2_compound_count_v61",
        "canonical_BE3_compound_count_v61",
        "canonical_best_binding_evidence_level_v61",
        "canonical_binding_evidence_count_v61",
        "canonical_binding_source_count_v61",
        "canonical_binding_sources_v61",
        "canonical_best_standard_value_nM_v61",
        "compound_master_version_v61",
        "v61_release_date",
    ):
        if field not in protein_fields:
            protein_fields.append(field)
    protein_path = CANDIDATE / "human_membrane_protein_master_v6_1.tsv"
    protein_rows = 0
    with atomic_text(protein_path, protein_fields) as writer:
        with PROTEIN0.open("rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                target = row["target_uniprot_id"]
                agg = con.execute(
                    """
                    SELECT COUNT(*),
                           SUM(CASE WHEN best_rank=1 THEN 1 ELSE 0 END),
                           SUM(CASE WHEN best_rank=2 THEN 1 ELSE 0 END),
                           SUM(CASE WHEN best_rank=3 THEN 1 ELSE 0 END),
                           MIN(best_rank),SUM(evidence_count),MIN(best_value)
                    FROM pair_agg WHERE target=?
                    """,
                    (target,),
                ).fetchone()
                source_rows = con.execute(
                    """
                    SELECT DISTINCT source FROM pair_raw WHERE target=?
                    """,
                    (target,),
                ).fetchall()
                sources = sorted(value[0] for value in source_rows)
                row.update(
                    {
                        "canonical_core_small_molecule_count_v61": agg[0] or 0,
                        "canonical_BE1_compound_count_v61": agg[1] or 0,
                        "canonical_BE2_compound_count_v61": agg[2] or 0,
                        "canonical_BE3_compound_count_v61": agg[3] or 0,
                        "canonical_best_binding_evidence_level_v61": {
                            1: "BE1",
                            2: "BE2",
                            3: "BE3",
                        }.get(agg[4], ""),
                        "canonical_binding_evidence_count_v61": agg[5] or 0,
                        "canonical_binding_source_count_v61": len(sources),
                        "canonical_binding_sources_v61": ";".join(sources),
                        "canonical_best_standard_value_nM_v61": (
                            "" if agg[6] is None else agg[6]
                        ),
                        "compound_master_version_v61": "v1.2",
                        "v61_release_date": RELEASE_DATE,
                    }
                )
                writer.writerow(
                    {field: row.get(field, "") for field in protein_fields}
                )
                protein_rows += 1

    with DISEASE0.open("rt", encoding="utf-8-sig", newline="") as handle:
        disease_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    if "release_version_v61" not in disease_fields:
        disease_fields.append("release_version_v61")
    disease_path = CANDIDATE / "protein_gene_disease_relations_v6_1.tsv"
    disease_rows = 0
    with atomic_text(disease_path, disease_fields) as writer:
        with DISEASE0.open("rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row["release_version_v61"] = RELEASE_VERSION
                writer.writerow(
                    {field: row.get(field, "") for field in disease_fields}
                )
                disease_rows += 1

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
    heterogeneity_path = CANDIDATE / "quantitative_heterogeneity_audit_v1_1.tsv.gz"
    heterogeneity_rows = 0
    with atomic_text(
        heterogeneity_path, heterogeneity_fields, gzipped=True
    ) as writer:
        for row in con.execute(
            """
            SELECT target,compound,activity_type,MIN(value_nm),MAX(value_nm),
                   COUNT(*),COUNT(DISTINCT source),
                   GROUP_CONCAT(DISTINCT source)
            FROM quant_raw
            WHERE value_nm>0 AND activity_type<>''
            GROUP BY target,compound,activity_type
            HAVING COUNT(*)>=2 AND MAX(value_nm)/MIN(value_nm)>=100
            ORDER BY MAX(value_nm)/MIN(value_nm) DESC
            """
        ):
            target, compound, activity, minimum, maximum, n, ns, sources = row
            writer.writerow(
                {
                    "target_uniprot_id": target,
                    "compound_internal_id": compound,
                    "activity_type": activity,
                    "minimum_standard_value_nM": minimum,
                    "maximum_standard_value_nM": maximum,
                    "max_min_ratio": maximum / minimum,
                    "observation_count": n,
                    "source_count": ns,
                    "sources": sources,
                    "interpretation": "preserve source-specific values; do not average",
                }
            )
            heterogeneity_rows += 1
    return {
        "pair_rows": pair_rows,
        "site_rows": site_rows,
        "protein_rows": protein_rows,
        "disease_rows": disease_rows,
        "heterogeneity_rows": heterogeneity_rows,
    }


def count_rows(path: Path) -> int:
    with open_text(path) as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def duplicate_and_fk_checks(
    con: sqlite3.Connection,
    all_master: dict[str, dict[str, str]],
    reconciliation: dict,
) -> dict:
    protein_ids = set()
    with (CANDIDATE / "human_membrane_protein_master_v6_1.tsv").open(
        "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        protein_rows = list(csv.DictReader(handle, delimiter="\t"))
        protein_ids = {row["target_uniprot_id"] for row in protein_rows}
    compound_ids = set(all_master)
    form_ids = set()
    form_to_compound = {}
    form_row_count = 0
    with (CANDIDATE / "compound_form_hierarchy_v1_2.tsv").open(
        "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            form_row_count += 1
            form_ids.add(row["compound_form_id"])
            form_to_compound[row["compound_form_id"]] = row["compound_internal_id"]
    master_row_count = count_rows(CANDIDATE / "small_molecule_master_v1_2.tsv")

    errors: Counter[str] = Counter()
    evidence_ids: set[str] = set()
    evidence_default_rows = 0
    pair_counter: Counter[tuple[str, str]] = Counter()
    with gzip.open(
        CANDIDATE / "binding_evidence_master_v6_1.tsv.gz",
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            evidence_id = row["evidence_id"]
            if evidence_id in evidence_ids:
                errors["duplicate_evidence_id"] += 1
            evidence_ids.add(evidence_id)
            if row["target_uniprot_id"] not in protein_ids:
                errors["evidence_missing_target_fk"] += 1
            if row["compound_internal_id"] not in compound_ids:
                errors["evidence_missing_compound_fk"] += 1
            if row["compound_form_id"] and row["compound_form_id"] not in form_ids:
                errors["evidence_missing_form_fk"] += 1
            if row["default_release_inclusion"] == "1":
                evidence_default_rows += 1
                pair_counter[
                    (row["target_uniprot_id"], row["compound_internal_id"])
                ] += 1

    pair_rows = 0
    pair_keys: set[tuple[str, str]] = set()
    pair_evidence_sum = 0
    with gzip.open(
        CANDIDATE / "protein_compound_summary_v6_1.tsv.gz",
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["target_uniprot_id"], row["compound_internal_id"])
            if key in pair_keys:
                errors["duplicate_pair_key"] += 1
            pair_keys.add(key)
            pair_rows += 1
            pair_evidence_sum += integer(row["binding_evidence_count"])
            if pair_counter.get(key, 0) != integer(row["binding_evidence_count"]):
                errors["pair_evidence_count_mismatch"] += 1

    negative_missing_fk = 0
    with gzip.open(
        CANDIDATE / "negative_binding_evidence_v1_1.tsv.gz",
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["compound_internal_id"] not in compound_ids:
                negative_missing_fk += 1
            if row["target_uniprot_id"] not in protein_ids:
                errors["negative_missing_target_fk"] += 1
    errors["negative_missing_compound_fk"] = negative_missing_fk

    accounted = sum(reconciliation["disposition_counts"].values()) + int(
        reconciliation.get("duplicate_source_evidence_id_rows", 0)
    )
    if accounted != reconciliation["review_rows_total"]:
        errors["review_row_reconciliation_mismatch"] += 1
    if pair_evidence_sum != evidence_default_rows:
        errors["global_pair_evidence_sum_mismatch"] += 1
    if master_row_count != len(compound_ids):
        errors["duplicate_compound_primary_key"] += master_row_count - len(compound_ids)
    if form_row_count != len(form_ids):
        errors["duplicate_form_primary_key"] += form_row_count - len(form_ids)
    errors = Counter({key: value for key, value in errors.items() if value})
    return {
        "blocking_errors": dict(errors),
        "blocking_error_total": sum(errors.values()),
        "protein_primary_key_unique": len(protein_ids) == len(protein_rows),
        "compound_primary_key_unique": master_row_count == len(compound_ids),
        "form_primary_key_unique": form_row_count == len(form_ids),
        "compound_row_count": master_row_count,
        "form_row_count": form_row_count,
        "evidence_id_unique_count": len(evidence_ids),
        "default_positive_evidence_rows": evidence_default_rows,
        "pair_rows": pair_rows,
        "pair_evidence_sum": pair_evidence_sum,
        "review_rows_accounted": accounted,
        "review_rows_total": reconciliation["review_rows_total"],
    }


def write_release_docs(
    identity_counts: dict,
    reconciliation: dict,
    canonical: dict,
    master_stats: dict,
    evidence_stats: dict,
    output_stats: dict,
    validation: dict,
) -> None:
    coverage_fields = ["metric", "value"]
    coverage = {
        **{f"identity_{key}": value for key, value in identity_counts.items()},
        **{
            f"disposition_{key}": value
            for key, value in reconciliation["disposition_counts"].items()
        },
        **{f"evidence_{key}": value for key, value in evidence_stats["counts"].items()},
        **output_stats,
    }
    with atomic_text(
        CANDIDATE / "source_coverage_v6_1.tsv", coverage_fields
    ) as writer:
        for key, value in sorted(coverage.items()):
            writer.writerow({"metric": key, "value": value})

    methods = """# MemPro V6.1 methods

V6.1 is an incremental identity-resolution release over immutable V6.0.
PubChem P1/P2/P3, PDBe CCD, PDBbind, and BRENDA records were mapped by exact
full InChIKey. Connectivity-only agreement never caused an automatic merge.
RDKit cleanup, fragment-parent, and charge-parent transformations were recorded
for newly released exact forms. Stereochemically distinct full InChIKeys remain
distinct. Active/observed records with safe identities entered the positive
evidence layer; inactive records entered the negative layer; ambiguous,
inconclusive, generic, covalent-special, and unresolved records remained in
review. No quantitative values were averaged.
"""
    (CANDIDATE / "METHODS_v6_1.md").write_text(methods, encoding="utf-8")
    policy = """# MemPro V6.1 inclusion policy

- Automatic identity merge requires exact full InChIKey.
- Connectivity-only, ambiguous stereochemistry, and multiple exact candidates
  remain in review.
- Parent, salt, charge, and exact form identifiers are stored separately.
- Positive default evidence requires a safe compound identity and an
  active/observed outcome.
- Inactive results are retained as negative evidence, never converted to
  positives.
- Covalent and non-small-molecule records require dedicated review.
"""
    (CANDIDATE / "INCLUSION_POLICY_v6_1.md").write_text(
        policy, encoding="utf-8"
    )
    dictionary = """# MemPro V6.1 data dictionary

The formal release contains the membrane-protein master, protein-gene-disease
relations, small-molecule master, exact form/parent hierarchy, positive binding
evidence, binding-site instances, protein-compound summary, negative evidence,
review queue, identity conflict audit, duplicate audit, quantitative
heterogeneity audit, source coverage, validation report, manifest, and hashes.
Identifiers are stable within the release. Full provenance remains on every
evidence row.
"""
    (CANDIDATE / "DATA_DICTIONARY_v6_1.md").write_text(
        dictionary, encoding="utf-8"
    )
    readme = f"""# MemPro V6.1

Formal incremental release frozen on {RELEASE_DATE}.

- Baseline: MemPro V6.0 (immutable)
- Compound master: small-molecule-v1.2
- Review rows reconciled: {reconciliation['review_rows_total']:,}
- New canonical parent compounds: {canonical['new_parent_compounds']:,}
- New exact forms: {canonical['new_exact_forms']:,}
- Validation status: {'PASS' if validation['blocking_error_total'] == 0 else 'FAIL'}
"""
    (CANDIDATE / "README_v6_1.md").write_text(readme, encoding="utf-8")
    release_info = {
        "release_version": RELEASE_VERSION,
        "release_date": RELEASE_DATE,
        "baseline_release": "v6.0",
        "small_molecule_master_version": "v1.2",
        "status": "frozen" if validation["blocking_error_total"] == 0 else "failed",
        "identity_counts": identity_counts,
        "reconciliation": reconciliation,
        "canonicalization": {
            key: value
            for key, value in canonical.items()
            if key not in {"parent_ids", "derived_parent_structures", "records"}
        },
        "master_stats": master_stats,
        "evidence_stats": evidence_stats,
        "output_stats": output_stats,
    }
    write_json(CANDIDATE / "RELEASE_INFO_v6_1.json", release_info)


def manifest_and_hashes() -> None:
    excluded_names = {
        "RELEASE_MANIFEST_v6_1.tsv",
        "SHA256SUMS_v6_1.txt",
        "V61_VALIDATION_REPORT.json",
    }
    files = sorted(
        path
        for path in CANDIDATE.iterdir()
        if path.is_file() and path.name not in excluded_names
    )
    manifest_fields = ["file_name", "bytes", "sha256", "row_count"]
    with atomic_text(
        CANDIDATE / "RELEASE_MANIFEST_v6_1.tsv", manifest_fields
    ) as writer:
        for path in files:
            row_count = (
                count_rows(path)
                if path.suffix in {".tsv", ".gz"} and ".tsv" in path.name
                else ""
            )
            writer.writerow(
                {
                    "file_name": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "row_count": row_count,
                }
            )
    checksum_files = sorted(
        path
        for path in CANDIDATE.iterdir()
        if path.is_file() and path.name != "SHA256SUMS_v6_1.txt"
    )
    checksum_text = "\n".join(
        f"{sha256_file(path)}  {path.name}" for path in checksum_files
    )
    (CANDIDATE / "SHA256SUMS_v6_1.txt").write_text(
        checksum_text + "\n", encoding="utf-8"
    )


def main() -> None:
    require_inputs()
    if CANDIDATE.exists():
        shutil.rmtree(CANDIDATE)
    CANDIDATE.mkdir(parents=True)
    update_progress("initialize")
    master_by_id, master_by_key, form_by_key, max_compound, max_form = (
        load_baseline_identity()
    )
    con = connect_db()

    update_progress("load_identity_sources")
    identity_counts = load_identity_sources(con)
    update_progress("classify_review_queue", identity_counts=identity_counts)
    reconciliation = classify_dispositions(con)

    update_progress(
        "canonicalize_parent_form",
        reconciliation=reconciliation,
    )
    canonical = canonicalize_candidates(
        con, master_by_key, form_by_key, max_compound, max_form
    )
    update_progress(
        "build_master_and_forms",
        canonicalization={
            key: value
            for key, value in canonical.items()
            if key not in {"parent_ids", "derived_parent_structures", "records"}
        },
    )
    canonical.pop("records", None)
    gc.collect()
    all_master, master_stats = build_master_and_forms(
        con, master_by_id, canonical, set(form_by_key.values())
    )
    canonical.pop("derived_parent_structures", None)
    canonical.pop("parent_ids", None)
    gc.collect()

    update_progress("rebuild_evidence_layers", master_stats=master_stats)
    evidence_stats = build_evidence_layers(con, all_master, form_by_key)

    update_progress("build_pair_site_protein_outputs")
    output_stats = build_pair_site_protein_outputs(con, all_master)

    update_progress("validate_referential_integrity")
    validation = duplicate_and_fk_checks(
        con, all_master, reconciliation
    )
    validation.update(
        {
            "status": (
                "PASS" if validation["blocking_error_total"] == 0 else "FAIL"
            ),
            "release_version": RELEASE_VERSION,
            "release_date": RELEASE_DATE,
            "generated_utc": now(),
            "identity_counts": identity_counts,
            "reconciliation": reconciliation,
            "canonicalization": {
                key: value
                for key, value in canonical.items()
                if key
                not in {"parent_ids", "derived_parent_structures", "records"}
            },
            "master_stats": master_stats,
            "evidence_stats": evidence_stats,
            "output_stats": output_stats,
        }
    )
    write_json(CANDIDATE / "V61_VALIDATION_REPORT.json", validation)
    write_release_docs(
        identity_counts,
        reconciliation,
        canonical,
        master_stats,
        evidence_stats,
        output_stats,
        validation,
    )
    manifest_and_hashes()

    if validation["blocking_error_total"]:
        update_progress(
            "validation_failed",
            status="failed",
            validation_report=str(CANDIDATE / "V61_VALIDATION_REPORT.json"),
            blocking_errors=validation["blocking_errors"],
        )
        raise RuntimeError(
            f"V6.1 not frozen; blocking errors: {validation['blocking_errors']}"
        )

    shutil.copytree(CANDIDATE, FINAL)
    freeze = {
        "status": "frozen",
        "release_version": RELEASE_VERSION,
        "release_path": str(FINAL),
        "frozen_utc": now(),
        "validation_report_sha256": sha256_file(
            FINAL / "V61_VALIDATION_REPORT.json"
        ),
        "release_manifest_sha256": sha256_file(
            FINAL / "RELEASE_MANIFEST_v6_1.tsv"
        ),
    }
    write_json(FINAL / "FREEZE_v6_1.json", freeze)
    update_progress(
        "complete",
        status="complete",
        release_path=str(FINAL),
        validation_status="PASS",
        freeze_manifest=str(FINAL / "FREEZE_v6_1.json"),
    )
    con.close()
    print(json.dumps(freeze, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
