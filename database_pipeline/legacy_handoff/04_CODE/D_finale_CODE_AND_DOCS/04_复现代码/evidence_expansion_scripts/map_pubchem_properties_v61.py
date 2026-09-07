#!/usr/bin/env python3
"""Map fetched PubChem structures to the frozen V6.0 compound identities."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import time
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
DEFAULT_RELEASE = ROOT / "releases" / "release_mempro_v6_0_20260727"
DEFAULT_RUN = ROOT / "runs" / "incremental_v61_20260727"


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def wait_for_fetch(progress: Path, properties: Path, max_hours: float) -> dict:
    deadline = time.time() + max_hours * 3600
    while time.time() < deadline:
        if progress.exists():
            try:
                state = json.loads(progress.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                state = {}
            if state.get("status") == "complete" and properties.exists():
                return state
            if state.get("status") == "incomplete":
                raise RuntimeError(f"Upstream fetch is incomplete: {state}")
        time.sleep(30)
    raise TimeoutError(f"Timed out waiting for {progress}")


def load_master(path: Path):
    full: dict[str, list[tuple[str, str]]] = defaultdict(list)
    connectivity: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    cids: dict[str, set[str]] = defaultdict(set)
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
            for cid in row["pubchem_cids"].replace(",", "|").split("|"):
                cid = cid.strip()
                if cid.isdigit():
                    cids[cid].add(compound_id)
    return full, connectivity, cids


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--job-tag", default="p1")
    parser.add_argument(
        "--fetch-progress-name",
        default="PUBCHEM_IDENTITY_FETCH_V61_PROGRESS.json",
        help="Progress filename under run/qa. The active P1 process uses the legacy name.",
    )
    parser.add_argument(
        "--properties-name",
        default="pubchem_properties_p1_v61.tsv.gz",
        help="Fetched properties filename under staging/pubchem_identity.",
    )
    parser.add_argument("--wait-hours", type=float, default=12)
    args = parser.parse_args()

    output_dir = args.run / "staging" / "pubchem_identity"
    qa_dir = args.run / "qa"
    fetch_progress = qa_dir / args.fetch_progress_name
    properties = output_dir / args.properties_name
    map_progress = qa_dir / f"PUBCHEM_IDENTITY_MAP_{args.job_tag.upper()}_V61_PROGRESS.json"
    output = output_dir / f"pubchem_identity_crosswalk_{args.job_tag}_v61.tsv.gz"
    cluster_output = (
        output_dir / f"pubchem_new_structure_clusters_{args.job_tag}_v61.tsv.gz"
    )
    started = now()
    write_json(
        map_progress,
        {
            "status": "waiting_for_fetch",
            "started_utc": started,
            "updated_utc": now(),
            "upstream_progress": str(fetch_progress),
        },
    )
    upstream = wait_for_fetch(fetch_progress, properties, args.wait_hours)
    write_json(
        map_progress,
        {
            "status": "running",
            "started_utc": started,
            "updated_utc": now(),
            "upstream_counts": upstream.get("counts", {}),
        },
    )

    full_map, connectivity_map, cid_map = load_master(
        args.release / "small_molecule_master_v1_1.tsv"
    )
    header = [
        "cid",
        "fetch_status",
        "molecular_formula",
        "molecular_weight",
        "smiles",
        "connectivity_smiles",
        "inchi",
        "inchikey",
        "connectivity_key",
        "formal_charge",
        "identity_resolution",
        "matched_compound_internal_ids",
        "matched_compound_names",
        "existing_cid_compound_internal_ids",
        "cid_structure_conflict",
        "release_action",
        "fetch_attempts",
        "fetch_error_message",
    ]
    counts: dict[str, int] = defaultdict(int)
    action_counts: dict[str, int] = defaultdict(int)
    conflicts = 0
    new_clusters: dict[str, dict] = {}
    processed = 0
    with gzip.open(properties, "rt", encoding="utf-8-sig", newline="") as source:
        with gzip.open(output, "wt", encoding="utf-8", newline="") as target:
            reader = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(
                target, fieldnames=header, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in reader:
                processed += 1
                cid = row["cid"]
                key = row["inchikey"].strip().upper()
                conn_key = key[:14] if key else ""
                exact = full_map.get(key, [])
                connected = connectivity_map.get(conn_key, [])
                existing_cid_ids = cid_map.get(cid, set())
                if row["status"] != "resolved" or not key:
                    resolution = "unresolved_fetch_or_structure"
                    matches = []
                    action = "retain_review"
                elif exact:
                    resolution = "exact_full_inchikey"
                    matches = exact
                    action = "map_existing_exact"
                elif connected:
                    resolution = "connectivity_only_form_or_stereo"
                    matches = [(value[0], value[1]) for value in connected]
                    action = "form_or_stereo_audit"
                else:
                    resolution = "new_structure_candidate"
                    matches = []
                    action = "create_candidate_after_global_dedup"
                    cluster = new_clusters.setdefault(
                        key,
                        {
                            "inchikey": key,
                            "connectivity_key": conn_key,
                            "smiles": row["smiles"],
                            "inchi": row["inchi"],
                            "molecular_formula": row["molecular_formula"],
                            "molecular_weight": row["molecular_weight"],
                            "cids": [],
                        },
                    )
                    cluster["cids"].append(cid)
                match_ids = {value[0] for value in matches}
                conflict = bool(
                    existing_cid_ids
                    and row["status"] == "resolved"
                    and (not match_ids or existing_cid_ids.isdisjoint(match_ids))
                )
                if conflict:
                    conflicts += 1
                    action = "identity_conflict_review"
                out = {
                    "cid": cid,
                    "fetch_status": row["status"],
                    "molecular_formula": row["molecular_formula"],
                    "molecular_weight": row["molecular_weight"],
                    "smiles": row["smiles"],
                    "connectivity_smiles": row["connectivity_smiles"],
                    "inchi": row["inchi"],
                    "inchikey": key,
                    "connectivity_key": conn_key,
                    "formal_charge": row["charge"],
                    "identity_resolution": resolution,
                    "matched_compound_internal_ids": ";".join(sorted(match_ids)),
                    "matched_compound_names": ";".join(
                        sorted({value[1] for value in matches})
                    ),
                    "existing_cid_compound_internal_ids": ";".join(
                        sorted(existing_cid_ids)
                    ),
                    "cid_structure_conflict": int(conflict),
                    "release_action": action,
                    "fetch_attempts": row["attempts"],
                    "fetch_error_message": row["error_message"],
                }
                writer.writerow(out)
                counts[resolution] += 1
                action_counts[action] += 1
                if processed % 25_000 == 0:
                    write_json(
                        map_progress,
                        {
                            "status": "running",
                            "started_utc": started,
                            "updated_utc": now(),
                            "rows_processed": processed,
                            "identity_resolution_counts": counts,
                            "release_action_counts": action_counts,
                            "cid_structure_conflicts": conflicts,
                        },
                    )

    cluster_header = [
        "inchikey",
        "connectivity_key",
        "smiles",
        "inchi",
        "molecular_formula",
        "molecular_weight",
        "cid_count",
        "cids",
    ]
    with gzip.open(cluster_output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=cluster_header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(new_clusters):
            item = new_clusters[key]
            writer.writerow(
                {
                    **item,
                    "cid_count": len(item["cids"]),
                    "cids": ";".join(sorted(item["cids"], key=int)),
                }
            )

    summary = {
        "status": "complete",
        "started_utc": started,
        "completed_utc": now(),
        "rows_processed": processed,
        "identity_resolution_counts": counts,
        "release_action_counts": action_counts,
        "cid_structure_conflicts": conflicts,
        "new_structure_cluster_count": len(new_clusters),
        "output": str(output),
        "new_cluster_output": str(cluster_output),
    }
    write_json(map_progress, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
