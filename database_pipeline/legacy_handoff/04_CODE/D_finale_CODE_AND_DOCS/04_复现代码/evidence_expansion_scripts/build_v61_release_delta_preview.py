#!/usr/bin/env python3
"""Build a non-release V6.1 compound/evidence delta from completed inputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import time
from collections import Counter
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
RUN = ROOT / "runs" / "incremental_v61_20260727"
PRE = RUN / "merge_candidate" / "premerge_existing_inputs"
BUILD_REPORT = RUN / "qa" / "V61_PREMERGE_EXISTING_INPUTS_REPORT.json"
OUT = RUN / "merge_candidate" / "release_delta_preview"
QA = RUN / "qa" / "V61_RELEASE_DELTA_PREVIEW_REPORT.json"
REVIEW = RELEASE / "binding_evidence_review_queue_v6_0.tsv.gz"
PUBCHEM_RE = re.compile(r"(?:PUBCHEM:|CID[: ]?)(\d+)", re.IGNORECASE)


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def form_id(inchikey: str) -> str:
    token = hashlib.sha1(("FORM:" + inchikey).encode("ascii")).hexdigest()[:14]
    return f"HMPD-V61-FORM-{token.upper()}"


def wait_for_premerge() -> dict:
    deadline = time.time() + 3 * 60 * 60
    while time.time() < deadline:
        if BUILD_REPORT.exists():
            try:
                payload = json.loads(BUILD_REPORT.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                payload = {}
            if payload.get("status") == "complete_premerge_not_release":
                return payload
        time.sleep(15)
    raise TimeoutError("Timed out waiting for the V6.1 premerge")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    premerge = wait_for_premerge()
    candidates_path = PRE / "global_identity_clusters_current_inputs_v61.tsv.gz"
    lookup_path = PRE / "source_identity_resolution_lookup_current_inputs_v61.tsv.gz"

    candidates: dict[str, dict[str, str]] = {}
    master_delta_path = OUT / "small_molecule_master_delta_preview_v61.tsv.gz"
    master_fields = [
        "provisional_compound_internal_id",
        "preferred_name_preview",
        "standard_smiles",
        "standard_inchi",
        "standard_inchikey",
        "connectivity_key",
        "molecular_formula",
        "source_databases",
        "pubchem_cids",
        "pdbe_component_ids",
        "pdbbind_complex_ids",
        "parent_relation",
        "proposed_parent_inchikey",
        "existing_parent_compound_ids",
        "identity_confidence_preview",
        "release_action",
    ]
    with gzip.open(
        candidates_path, "rt", encoding="utf-8-sig", newline=""
    ) as source, gzip.open(
        master_delta_path, "wt", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(
            destination,
            fieldnames=master_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            candidates[row["provisional_cluster_id"]] = row
            if not row["provisional_cluster_id"]:
                continue
            if row["structure_parse_status"] != "ok":
                confidence = "low_manual_review"
            elif row["source_count"] and int(row["source_count"]) > 1:
                confidence = "high_cross_source_structure"
            else:
                confidence = "medium_single_source_structure"
            if row["pubchem_cids"]:
                preferred = "PubChem CID " + row["pubchem_cids"].split(";")[0]
            elif row["pdbe_component_ids"]:
                preferred = "PDB CCD " + row["pdbe_component_ids"].split(";")[0]
            else:
                preferred = "PDBbind " + row["pdbbind_complex_ids"].split(";")[0]
            writer.writerow(
                {
                    "provisional_compound_internal_id": row[
                        "provisional_cluster_id"
                    ],
                    "preferred_name_preview": preferred,
                    "standard_smiles": row["smiles"],
                    "standard_inchi": row["inchi"],
                    "standard_inchikey": row["inchikey"],
                    "connectivity_key": row["connectivity_key"],
                    "molecular_formula": row["formula"],
                    "source_databases": row["sources"],
                    "pubchem_cids": row["pubchem_cids"],
                    "pdbe_component_ids": row["pdbe_component_ids"],
                    "pdbbind_complex_ids": row["pdbbind_complex_ids"],
                    "parent_relation": row["parent_relation"],
                    "proposed_parent_inchikey": row["charge_parent_inchikey"],
                    "existing_parent_compound_ids": row[
                        "existing_parent_compound_ids"
                    ],
                    "identity_confidence_preview": confidence,
                    "release_action": "preview_only_wait_for_P2_P3_BRENDA",
                }
            )

    form_delta_path = OUT / "compound_form_hierarchy_delta_preview_v61.tsv.gz"
    form_fields = [
        "provisional_form_id",
        "provisional_compound_internal_id",
        "exact_inchikey",
        "fragment_parent_inchikey",
        "charge_parent_inchikey",
        "parent_relation",
        "formal_release_status",
    ]
    with gzip.open(
        PRE / "parent_form_stereo_candidates_current_inputs_v61.tsv.gz",
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as source, gzip.open(
        form_delta_path, "wt", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(
            destination,
            fieldnames=form_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            if not row["provisional_cluster_id"]:
                continue
            writer.writerow(
                {
                    "provisional_form_id": form_id(row["exact_inchikey"]),
                    "provisional_compound_internal_id": row[
                        "provisional_cluster_id"
                    ],
                    "exact_inchikey": row["exact_inchikey"],
                    "fragment_parent_inchikey": row[
                        "fragment_parent_inchikey"
                    ],
                    "charge_parent_inchikey": row["charge_parent_inchikey"],
                    "parent_relation": row["parent_relation"],
                    "formal_release_status": "preview_only",
                }
            )

    lookup: dict[tuple[str, str], dict[str, str]] = {}
    with gzip.open(
        lookup_path, "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lookup[(row["source_object_type"], row["source_object_id"])] = row

    evidence_path = OUT / "binding_evidence_resolution_delta_preview_v61.tsv.gz"
    evidence_fields = [
        "source_evidence_id",
        "target_uniprot_id",
        "source_database",
        "source_record_id",
        "compound_source_id",
        "resolved_compound_internal_id_preview",
        "identity_resolution",
        "evidence_tier",
        "activity_outcome",
        "proposed_release_action",
        "proposed_default_inclusion",
        "review_reason",
    ]
    action_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    written = 0
    with gzip.open(
        REVIEW, "rt", encoding="utf-8-sig", newline=""
    ) as source, gzip.open(
        evidence_path, "wt", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(
            destination,
            fieldnames=evidence_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            source_db = row["source_database"]
            identity = None
            if source_db.startswith("PubChem"):
                match = PUBCHEM_RE.search(row["compound_source_id"])
                cid = match.group(1) if match else ""
                identity = lookup.get(("PUBCHEM_CID", cid))
            elif source_db == "PDBe":
                component = row["compound_source_id"].split(":", 1)[-1].upper()
                identity = lookup.get(("PDB_CCD", component))
            elif source_db == "PDBbind":
                complex_id = (
                    row["source_record_id"].split(":", 1)[0].lower()
                    or row["pdb_ids"].split(";", 1)[0].lower()
                )
                identity = lookup.get(("PDBBIND_COMPLEX", complex_id))
            if identity is None:
                continue

            existing_ids = [
                value
                for value in identity["mapped_existing_compound_ids"].split(";")
                if value
            ]
            preview_id = (
                existing_ids[0]
                if len(existing_ids) == 1
                else identity["provisional_cluster_id"]
            )
            outcome = row["activity_outcome"].casefold()
            eligibility = identity["eligibility"].casefold()
            notes = identity["notes"]
            if len(existing_ids) > 1 or "conflict" in notes.casefold():
                action = "identity_conflict_review"
                default = "0"
            elif "exclude" in eligibility:
                action = "exclude_nonrelease_identity"
                default = "0"
            elif not preview_id:
                action = "retain_identity_review"
                default = "0"
            elif outcome == "inactive":
                action = "negative_evidence"
                default = "0"
            elif outcome == "inconclusive":
                action = "retain_inconclusive_review"
                default = "0"
            elif outcome in {"active", "observed"}:
                if "covalent" in eligibility:
                    action = "positive_covalent_separate_layer"
                    default = "0"
                elif "review" in eligibility or "audit" in identity[
                    "release_action"
                ].casefold():
                    action = "retain_identity_review"
                    default = "0"
                else:
                    action = "positive_release_candidate"
                    default = "1"
            else:
                action = "retain_unspecified_review"
                default = "0"
            writer.writerow(
                {
                    "source_evidence_id": row["source_evidence_id"],
                    "target_uniprot_id": row["target_uniprot_id"],
                    "source_database": source_db,
                    "source_record_id": row["source_record_id"],
                    "compound_source_id": row["compound_source_id"],
                    "resolved_compound_internal_id_preview": preview_id,
                    "identity_resolution": identity[
                        "input_identity_resolution"
                    ],
                    "evidence_tier": row["evidence_tier"],
                    "activity_outcome": row["activity_outcome"],
                    "proposed_release_action": action,
                    "proposed_default_inclusion": default,
                    "review_reason": row["review_reason"],
                }
            )
            written += 1
            action_counts[action] += 1
            source_counts[source_db] += 1

    report = {
        "status": "complete_preview_not_release",
        "input_premerge_status": premerge["status"],
        "candidate_compound_rows": sum(
            1 for value in candidates.values() if value["provisional_cluster_id"]
        ),
        "evidence_resolution_preview_rows": written,
        "evidence_action_counts": dict(action_counts),
        "evidence_source_counts": dict(source_counts),
        "outputs": {
            "small_molecule_master_delta_preview": str(master_delta_path),
            "compound_form_hierarchy_delta_preview": str(form_delta_path),
            "binding_evidence_resolution_delta_preview": str(evidence_path),
        },
        "release_ready": False,
        "required_refresh_inputs": ["PubChem P2", "PubChem P3", "BRENDA"],
    }
    write_json(QA, report)
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
