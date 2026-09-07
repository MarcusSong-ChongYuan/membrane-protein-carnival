#!/usr/bin/env python3
"""Build a refreshable V6.1 premerge from currently completed identity sources."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize


RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
RUN = ROOT / "runs" / "incremental_v61_20260727"
OUT = RUN / "merge_candidate" / "premerge_existing_inputs"
QA = RUN / "qa" / "V61_PREMERGE_EXISTING_INPUTS_REPORT.json"

MASTER = RELEASE / "small_molecule_master_v1_1.tsv"
REVIEW = RELEASE / "binding_evidence_review_queue_v6_0.tsv.gz"
P1_CROSSWALK = (
    RUN
    / "staging"
    / "pubchem_identity"
    / "pubchem_identity_crosswalk_p1_v61.tsv.gz"
)
P1_NEW = (
    RUN
    / "staging"
    / "pubchem_identity"
    / "pubchem_new_structure_clusters_p1_v61.tsv.gz"
)
PDBE = RUN / "staging" / "pdbe_identity" / "pdbe_ccd_identity_v61.tsv.gz"
PDBBIND = (
    RUN / "staging" / "pdbbind_identity" / "pdbbind_ligand_identity_v61.tsv.gz"
)
STRUCTURAL_NEW = (
    RUN
    / "staging"
    / "structural_identity_merge"
    / "structural_new_compound_clusters_v61.tsv.gz"
)

PUBCHEM_RE = re.compile(r"(?:PUBCHEM:|CID[: ]?)(\d+)", re.IGNORECASE)


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def provisional_id(inchikey: str) -> str:
    token = hashlib.sha1(inchikey.encode("ascii", errors="ignore")).hexdigest()[:14]
    return f"HMPD-V61-CAND-{token.upper()}"


def mol_identity(mol: Chem.Mol) -> tuple[str, str, str]:
    smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    inchi = Chem.MolToInchi(mol)
    key = Chem.InchiToInchiKey(inchi) if inchi else ""
    return smiles, inchi, key


def standardize_candidate(smiles: str, expected_key: str) -> dict[str, object]:
    result: dict[str, object] = {
        "structure_parse_status": "unresolved",
        "computed_full_inchikey": "",
        "full_key_matches_source": 0,
        "fragment_count": "",
        "fragment_parent_smiles": "",
        "fragment_parent_inchikey": "",
        "charge_parent_smiles": "",
        "charge_parent_inchikey": "",
        "parent_connectivity_key": "",
        "parent_relation": "manual_review",
        "stereocenter_count": "",
        "unassigned_stereocenter_count": "",
        "standardization_error": "",
    }
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("RDKit returned no molecule")
        cleaned = rdMolStandardize.Cleanup(mol)
        full_smiles, _full_inchi, computed_key = mol_identity(cleaned)
        fragments = len(Chem.GetMolFrags(cleaned))
        fragment_parent = rdMolStandardize.FragmentParent(cleaned)
        fragment_smiles, _fragment_inchi, fragment_key = mol_identity(fragment_parent)
        charge_parent = rdMolStandardize.ChargeParent(fragment_parent)
        charge_smiles, _charge_inchi, charge_key = mol_identity(charge_parent)
        centers = Chem.FindMolChiralCenters(
            cleaned, includeUnassigned=True, includeCIP=True
        )
        unassigned = sum(1 for _index, label in centers if label == "?")
        if fragment_key and fragment_key != computed_key and fragments > 1:
            relation = "salt_or_multicomponent_form"
        elif charge_key and charge_key != fragment_key:
            relation = "charged_or_protonation_form"
        elif centers:
            relation = "stereochemically_defined_or_partial_parent"
        else:
            relation = "neutral_parent_candidate"
        result.update(
            {
                "structure_parse_status": "ok",
                "computed_full_inchikey": computed_key,
                "full_key_matches_source": int(
                    bool(expected_key) and computed_key == expected_key
                ),
                "standardized_smiles": full_smiles,
                "fragment_count": fragments,
                "fragment_parent_smiles": fragment_smiles,
                "fragment_parent_inchikey": fragment_key,
                "charge_parent_smiles": charge_smiles,
                "charge_parent_inchikey": charge_key,
                "parent_connectivity_key": charge_key[:14] if charge_key else "",
                "parent_relation": relation,
                "stereocenter_count": len(centers),
                "unassigned_stereocenter_count": unassigned,
            }
        )
    except Exception as exc:
        result["standardization_error"] = f"{type(exc).__name__}:{exc}"
    return result


def load_master() -> tuple[dict, dict, int]:
    full: dict[str, list[tuple[str, str]]] = defaultdict(list)
    connectivity: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    rows = 0
    with MASTER.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            key = row["standard_inchikey"].strip().upper()
            conn = row["connectivity_key"].strip().upper()
            item = (row["compound_internal_id"], row["preferred_name"])
            if key:
                full[key].append(item)
            if conn:
                connectivity[conn].append((*item, key))
    return full, connectivity, rows


def add_candidate(
    candidates: dict[str, dict],
    inchikey: str,
    connectivity_key: str,
    smiles: str,
    inchi: str,
    formula: str,
    source: str,
    source_ids: list[str],
    evidence_rows: int,
) -> None:
    key = inchikey.strip().upper()
    if not key:
        return
    item = candidates.setdefault(
        key,
        {
            "inchikey": key,
            "connectivity_key": connectivity_key.strip().upper() or key[:14],
            "smiles": smiles,
            "inchi": inchi,
            "formula": formula,
            "sources": set(),
            "pubchem_cids": set(),
            "pdbe_component_ids": set(),
            "pdbbind_complex_ids": set(),
            "source_evidence_rows": 0,
        },
    )
    item["sources"].add(source)
    if not item["smiles"] and smiles:
        item["smiles"] = smiles
    if not item["inchi"] and inchi:
        item["inchi"] = inchi
    if not item["formula"] and formula:
        item["formula"] = formula
    bucket = {
        "PubChem": "pubchem_cids",
        "PDBe": "pdbe_component_ids",
        "PDBbind": "pdbbind_complex_ids",
    }[source]
    item[bucket].update(value for value in source_ids if value)
    item["source_evidence_rows"] += evidence_rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    started = now()
    required = [MASTER, REVIEW, P1_CROSSWALK, P1_NEW, PDBE, PDBBIND, STRUCTURAL_NEW]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + "; ".join(missing))

    master_full, master_connectivity, master_rows = load_master()
    candidates: dict[str, dict] = {}

    with gzip.open(P1_NEW, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            add_candidate(
                candidates,
                row["inchikey"],
                row["connectivity_key"],
                row["smiles"],
                row["inchi"],
                row["molecular_formula"],
                "PubChem",
                row["cids"].split(";"),
                int(row["cid_count"]),
            )

    with gzip.open(
        STRUCTURAL_NEW, "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if "PDBe" in row["sources"].split(";"):
                add_candidate(
                    candidates,
                    row["inchikey"],
                    row["connectivity_key"],
                    row["smiles"],
                    row["inchi"],
                    row["formula"],
                    "PDBe",
                    row["pdbe_component_ids"].split(";"),
                    int(row["evidence_row_count"]),
                )
            if "PDBbind" in row["sources"].split(";"):
                add_candidate(
                    candidates,
                    row["inchikey"],
                    row["connectivity_key"],
                    row["smiles"],
                    row["inchi"],
                    row["formula"],
                    "PDBbind",
                    row["pdbbind_complex_ids"].split(";"),
                    int(row["pdbbind_complex_count"]),
                )

    relation_counts: Counter[str] = Counter()
    parse_counts: Counter[str] = Counter()
    source_overlap_counts: Counter[str] = Counter()
    exact_master_collisions = 0
    for key, item in candidates.items():
        structure = standardize_candidate(item["smiles"], key)
        item.update(structure)
        existing = master_full.get(key, [])
        parent_key = str(item.get("charge_parent_inchikey", ""))
        parent_existing = master_full.get(parent_key, []) if parent_key else []
        parent_conn = str(item.get("parent_connectivity_key", ""))
        parent_connectivity_matches = master_connectivity.get(parent_conn, [])
        item["existing_exact_compound_ids"] = sorted(
            {value[0] for value in existing}
        )
        item["existing_parent_compound_ids"] = sorted(
            {value[0] for value in parent_existing}
        )
        item["existing_parent_connectivity_compound_ids"] = sorted(
            {value[0] for value in parent_connectivity_matches}
        )
        if existing:
            exact_master_collisions += 1
            item["candidate_action"] = "map_existing_unexpected_collision"
            item["provisional_cluster_id"] = ""
        else:
            item["candidate_action"] = "new_candidate_pending_full_refresh"
            item["provisional_cluster_id"] = provisional_id(key)
        relation_counts[str(item["parent_relation"])] += 1
        parse_counts[str(item["structure_parse_status"])] += 1
        source_overlap_counts["+".join(sorted(item["sources"]))] += 1

    candidate_path = OUT / "global_identity_clusters_current_inputs_v61.tsv.gz"
    candidate_fields = [
        "provisional_cluster_id",
        "inchikey",
        "connectivity_key",
        "smiles",
        "inchi",
        "formula",
        "sources",
        "source_count",
        "pubchem_cid_count",
        "pubchem_cids",
        "pdbe_component_count",
        "pdbe_component_ids",
        "pdbbind_complex_count",
        "pdbbind_complex_ids",
        "source_evidence_rows",
        "structure_parse_status",
        "computed_full_inchikey",
        "full_key_matches_source",
        "fragment_count",
        "fragment_parent_smiles",
        "fragment_parent_inchikey",
        "charge_parent_smiles",
        "charge_parent_inchikey",
        "parent_connectivity_key",
        "parent_relation",
        "stereocenter_count",
        "unassigned_stereocenter_count",
        "existing_exact_compound_ids",
        "existing_parent_compound_ids",
        "existing_parent_connectivity_compound_ids",
        "candidate_action",
        "standardization_error",
    ]
    with gzip.open(candidate_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=candidate_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(candidates):
            item = candidates[key]
            row = {
                **item,
                "sources": ";".join(sorted(item["sources"])),
                "source_count": len(item["sources"]),
                "pubchem_cid_count": len(item["pubchem_cids"]),
                "pubchem_cids": ";".join(
                    sorted(item["pubchem_cids"], key=lambda x: int(x))
                ),
                "pdbe_component_count": len(item["pdbe_component_ids"]),
                "pdbe_component_ids": ";".join(sorted(item["pdbe_component_ids"])),
                "pdbbind_complex_count": len(item["pdbbind_complex_ids"]),
                "pdbbind_complex_ids": ";".join(
                    sorted(item["pdbbind_complex_ids"])
                ),
                "existing_exact_compound_ids": ";".join(
                    item["existing_exact_compound_ids"]
                ),
                "existing_parent_compound_ids": ";".join(
                    item["existing_parent_compound_ids"]
                ),
                "existing_parent_connectivity_compound_ids": ";".join(
                    item["existing_parent_connectivity_compound_ids"]
                ),
            }
            writer.writerow({field: row.get(field, "") for field in candidate_fields})

    form_path = OUT / "parent_form_stereo_candidates_current_inputs_v61.tsv.gz"
    form_fields = [
        "provisional_cluster_id",
        "exact_inchikey",
        "connectivity_key",
        "fragment_parent_inchikey",
        "charge_parent_inchikey",
        "parent_connectivity_key",
        "parent_relation",
        "fragment_count",
        "stereocenter_count",
        "unassigned_stereocenter_count",
        "existing_parent_compound_ids",
        "existing_parent_connectivity_compound_ids",
        "required_release_action",
    ]
    with gzip.open(form_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=form_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(candidates):
            item = candidates[key]
            relation = item["parent_relation"]
            if item["structure_parse_status"] != "ok":
                action = "manual_structure_review"
            elif item["existing_parent_compound_ids"]:
                action = "link_form_to_existing_parent_after_audit"
            elif len(item["existing_parent_connectivity_compound_ids"]) == 1:
                action = "review_stereo_or_charge_link_to_existing_connectivity"
            elif len(item["existing_parent_connectivity_compound_ids"]) > 1:
                action = "ambiguous_existing_parent_family_review"
            elif relation in {
                "salt_or_multicomponent_form",
                "charged_or_protonation_form",
            }:
                action = "create_new_parent_and_form_if_absent_after_global_refresh"
            else:
                action = "create_new_parent_candidate_after_global_refresh"
            writer.writerow(
                {
                    "provisional_cluster_id": item["provisional_cluster_id"],
                    "exact_inchikey": key,
                    "connectivity_key": item["connectivity_key"],
                    "fragment_parent_inchikey": item["fragment_parent_inchikey"],
                    "charge_parent_inchikey": item["charge_parent_inchikey"],
                    "parent_connectivity_key": item["parent_connectivity_key"],
                    "parent_relation": relation,
                    "fragment_count": item["fragment_count"],
                    "stereocenter_count": item["stereocenter_count"],
                    "unassigned_stereocenter_count": item[
                        "unassigned_stereocenter_count"
                    ],
                    "existing_parent_compound_ids": ";".join(
                        item["existing_parent_compound_ids"]
                    ),
                    "existing_parent_connectivity_compound_ids": ";".join(
                        item["existing_parent_connectivity_compound_ids"]
                    ),
                    "required_release_action": action,
                }
            )

    lookup_rows: list[dict[str, object]] = []

    with gzip.open(P1_CROSSWALK, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["inchikey"].upper()
            candidate = candidates.get(key, {})
            lookup_rows.append(
                {
                    "source_object_type": "PUBCHEM_CID",
                    "source_object_id": row["cid"],
                    "source_database": "PubChem BioAssay",
                    "input_identity_resolution": row["identity_resolution"],
                    "release_action": row["release_action"],
                    "mapped_existing_compound_ids": row[
                        "matched_compound_internal_ids"
                    ],
                    "provisional_cluster_id": candidate.get(
                        "provisional_cluster_id", ""
                    ),
                    "exact_inchikey": key,
                    "connectivity_key": row["connectivity_key"],
                    "eligibility": (
                        "eligible_identity"
                        if row["identity_resolution"]
                        in {"exact_full_inchikey", "new_structure_candidate"}
                        else "form_or_stereo_review"
                    ),
                    "notes": (
                        "CID structure conflict"
                        if row["cid_structure_conflict"] == "1"
                        else ""
                    ),
                }
            )

    with gzip.open(PDBE, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["inchikey"].upper()
            candidate = candidates.get(key, {})
            lookup_rows.append(
                {
                    "source_object_type": "PDB_CCD",
                    "source_object_id": row["component_id"],
                    "source_database": "PDBe",
                    "input_identity_resolution": row["identity_resolution"],
                    "release_action": row["release_eligibility"],
                    "mapped_existing_compound_ids": row[
                        "matched_compound_internal_ids"
                    ],
                    "provisional_cluster_id": candidate.get(
                        "provisional_cluster_id", ""
                    ),
                    "exact_inchikey": key,
                    "connectivity_key": row["connectivity_key"],
                    "eligibility": row["release_eligibility"],
                    "notes": "",
                }
            )

    with gzip.open(PDBBIND, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["standard_inchikey"].upper()
            candidate = candidates.get(key, {})
            lookup_rows.append(
                {
                    "source_object_type": "PDBBIND_COMPLEX",
                    "source_object_id": row["complex_id"],
                    "source_database": "PDBbind",
                    "input_identity_resolution": row["identity_resolution"],
                    "release_action": row["structure_eligibility"],
                    "mapped_existing_compound_ids": row[
                        "matched_compound_internal_ids"
                    ],
                    "provisional_cluster_id": candidate.get(
                        "provisional_cluster_id", ""
                    ),
                    "exact_inchikey": key,
                    "connectivity_key": row["connectivity_key"],
                    "eligibility": row["structure_eligibility"],
                    "notes": row["standardization_error"],
                }
            )

    lookup_path = OUT / "source_identity_resolution_lookup_current_inputs_v61.tsv.gz"
    lookup_fields = [
        "source_object_type",
        "source_object_id",
        "source_database",
        "input_identity_resolution",
        "release_action",
        "mapped_existing_compound_ids",
        "provisional_cluster_id",
        "exact_inchikey",
        "connectivity_key",
        "eligibility",
        "notes",
    ]
    lookup_index: dict[tuple[str, str], dict] = {}
    with gzip.open(lookup_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=lookup_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in sorted(
            lookup_rows,
            key=lambda value: (
                str(value["source_object_type"]),
                str(value["source_object_id"]),
            ),
        ):
            writer.writerow(row)
            lookup_index[
                (str(row["source_object_type"]), str(row["source_object_id"]))
            ] = row

    review_actions: Counter[str] = Counter()
    review_sources: Counter[str] = Counter()
    review_rows = 0
    with gzip.open(REVIEW, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            review_rows += 1
            source = row["source_database"]
            review_sources[source] += 1
            lookup = None
            if source.startswith("PubChem"):
                cid = ""
                match = PUBCHEM_RE.search(row["compound_source_id"])
                if match:
                    cid = match.group(1)
                lookup = lookup_index.get(("PUBCHEM_CID", cid))
                if lookup is None:
                    review_actions["pending_P2_or_P3"] += 1
                    continue
            elif source == "PDBe":
                component = row["compound_source_id"].split(":", 1)[-1].upper()
                lookup = lookup_index.get(("PDB_CCD", component))
            elif source == "PDBbind":
                complex_id = (
                    row["source_record_id"].split(":", 1)[0].lower()
                    or row["pdb_ids"].split(";", 1)[0].lower()
                )
                lookup = lookup_index.get(("PDBBIND_COMPLEX", complex_id))
            elif source == "BRENDA":
                review_actions["pending_BRENDA"] += 1
                continue
            if lookup is None:
                review_actions["unmatched_source_object"] += 1
            elif lookup["mapped_existing_compound_ids"]:
                review_actions["map_existing_compound"] += 1
            elif lookup["provisional_cluster_id"]:
                review_actions["provisional_new_compound"] += 1
            elif "exclude" in str(lookup["eligibility"]):
                review_actions["exclude_nonrelease_identity"] += 1
            elif "review" in str(lookup["eligibility"]) or "audit" in str(
                lookup["release_action"]
            ):
                review_actions["retain_identity_review"] += 1
            else:
                review_actions["retain_unresolved"] += 1

    duplicate_provisional_ids = len(candidates) - len(
        {
            item["provisional_cluster_id"]
            for item in candidates.values()
            if item["provisional_cluster_id"]
        }
    )
    report = {
        "status": "complete_premerge_not_release",
        "started_utc": started,
        "completed_utc": now(),
        "input_scope": ["V6.0 compound master", "PubChem P1", "PDBe", "PDBbind"],
        "upstream_not_yet_included": ["PubChem P2", "PubChem P3", "BRENDA"],
        "baseline_compound_master_rows": master_rows,
        "baseline_unique_full_inchikeys": len(master_full),
        "baseline_unique_connectivity_keys": len(master_connectivity),
        "candidate_unique_full_inchikeys": len(candidates),
        "candidate_source_overlap_counts": dict(source_overlap_counts),
        "candidate_parent_relation_counts": dict(relation_counts),
        "candidate_parse_counts": dict(parse_counts),
        "candidate_exact_master_collisions": exact_master_collisions,
        "duplicate_provisional_ids": duplicate_provisional_ids,
        "source_identity_lookup_rows": len(lookup_rows),
        "review_queue_rows_scanned": review_rows,
        "review_queue_source_rows": dict(review_sources),
        "review_queue_premerge_actions": dict(review_actions),
        "outputs": {
            "global_identity_clusters": str(candidate_path),
            "parent_form_candidates": str(form_path),
            "source_identity_lookup": str(lookup_path),
        },
        "release_ready": False,
        "release_blockers": [
            "PubChem P2 identity fetch and mapping incomplete",
            "PubChem P3 identity fetch and mapping not started",
            "BRENDA identity resolution not started",
            "Final evidence rebuild, foreign-key validation, and source reconciliation pending",
        ],
    }
    write_json(QA, report)
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
