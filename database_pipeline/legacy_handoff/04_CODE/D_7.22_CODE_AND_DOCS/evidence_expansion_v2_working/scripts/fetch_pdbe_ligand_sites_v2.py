#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
import time
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
PROTEIN_MASTER = Path(
    r"D:\7.22\small_molecule_v1_working\releases\release_small_molecule_v1_0"
    r"\human_membrane_audit_master_v5_5.tsv"
)
RAW_DIR = ROOT / "raw" / "pdbe_202607" / "ligand_sites"
CHECKPOINT = ROOT / "raw" / "pdbe_202607" / "pdbe_ligand_sites_fetch_v2.jsonl"
INVENTORY = ROOT / "intermediate" / "pdbe_ligand_sites_inventory_v2.tsv"
REPORT = ROOT / "reports" / "PDBE_LIGAND_SITES_FETCH_REPORT.json"
LOG = ROOT / "logs" / "pdbe_ligand_sites_fetch_v2.log"
BASE = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/ligand_sites"
USER_AGENT = "MemProDB-evidence-expansion/2.0 (scientific database integration)"


def log(message: str) -> None:
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def load_targets() -> list[dict[str, str]]:
    targets = []
    with PROTEIN_MASTER.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("membrane_class_v52") not in {"A", "B", "C"}:
                continue
            pdb_ids = (row.get("pdb_ids") or "").strip()
            if not pdb_ids:
                continue
            targets.append(
                {
                    "target_uniprot_id": row["target_uniprot_id"].strip(),
                    "approved_symbol": (row.get("approved_symbol") or "").strip(),
                    "baseline_pdb_ids": pdb_ids,
                }
            )
    targets.sort(key=lambda row: row["target_uniprot_id"])
    return targets


def load_completed() -> dict[str, dict]:
    completed: dict[str, dict] = {}
    if not CHECKPOINT.exists():
        return completed
    with CHECKPOINT.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            uniprot = record.get("target_uniprot_id")
            if uniprot and record.get("status") in {"ok", "not_found", "no_data"}:
                completed[uniprot] = record
    return completed


def fetch(uniprot: str) -> tuple[str, bytes, str]:
    url = f"{BASE}/{uniprot}"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                payload = response.read()
            decoded = json.loads(payload.decode("utf-8"))
            data = decoded.get(uniprot, {})
            if not data or not data.get("data"):
                return "no_data", payload, url
            return "ok", payload, url
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return "not_found", b"", url
            if error.code in {429, 500, 502, 503, 504} and attempt < 4:
                time.sleep(5 * (attempt + 1))
                continue
            return f"http_{error.code}", b"", url
        except Exception as error:
            if attempt < 4:
                time.sleep(5 * (attempt + 1))
                continue
            return f"error_{type(error).__name__}", b"", url
    return "retry_exhausted", b"", url


def summarize_payload(uniprot: str, payload: bytes) -> dict:
    if not payload:
        return {
            "ligand_count": 0,
            "ligands_with_residues": 0,
            "residue_contact_count": 0,
            "pdb_count": 0,
        }
    decoded = json.loads(payload.decode("utf-8"))
    data = decoded.get(uniprot, {}).get("data", [])
    pdb_ids = set()
    residue_count = 0
    ligands_with_residues = 0
    for ligand in data:
        residues = ligand.get("residues") or []
        if residues:
            ligands_with_residues += 1
        residue_count += len(residues)
        additional = ligand.get("additionalData") or {}
        pdb_ids.update((pdb or "").lower() for pdb in additional.get("pdbEntries") or [])
        for residue in residues:
            pdb_ids.update(
                (pdb or "").lower() for pdb in residue.get("allPDBEntries") or []
            )
    pdb_ids.discard("")
    return {
        "ligand_count": len(data),
        "ligands_with_residues": ligands_with_residues,
        "residue_contact_count": residue_count,
        "pdb_count": len(pdb_ids),
    }


def write_outputs(targets: list[dict[str, str]], completed: dict[str, dict]) -> None:
    records = [
        completed[target["target_uniprot_id"]]
        for target in targets
        if target["target_uniprot_id"] in completed
    ]
    fields = [
        "target_uniprot_id",
        "approved_symbol",
        "status",
        "ligand_count",
        "ligands_with_residues",
        "residue_contact_count",
        "pdb_count",
        "source_url",
        "retrieved_at_utc",
    ]
    with INVENTORY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({field: record.get(field, "") for field in fields} for record in records)

    status_counts = Counter(record["status"] for record in records)
    report = {
        "status": (
            "passed"
            if len(records) == len(targets)
            and set(status_counts) <= {"ok", "not_found", "no_data"}
            else "incomplete"
        ),
        "eligible_abc_proteins_with_baseline_pdb_ids": len(targets),
        "completed_targets": len(records),
        "targets_with_ligand_sites": sum(record["ligand_count"] > 0 for record in records),
        "ligand_site_records": sum(record["ligand_count"] for record in records),
        "residue_contacts": sum(record["residue_contact_count"] for record in records),
        "fetch_status_counts": dict(status_counts),
        "raw_directory": str(RAW_DIR),
        "checkpoint_file": str(CHECKPOINT),
    }
    with REPORT.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    targets = load_targets()
    completed = load_completed()
    log(
        f"{datetime.now(timezone.utc).isoformat()} "
        f"start targets={len(targets)} completed={len(completed)}"
    )

    with CHECKPOINT.open("a", encoding="utf-8") as checkpoint:
        for index, target in enumerate(targets, 1):
            uniprot = target["target_uniprot_id"]
            if uniprot in completed:
                continue
            status, payload, url = fetch(uniprot)
            summary = summarize_payload(uniprot, payload)
            if payload:
                with gzip.open(
                    RAW_DIR / f"{uniprot}.json.gz", "wb", compresslevel=6
                ) as handle:
                    handle.write(payload)
            record = {
                **target,
                "status": status,
                **summary,
                "source_url": url,
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            }
            checkpoint.write(json.dumps(record, ensure_ascii=False) + "\n")
            checkpoint.flush()
            if status in {"ok", "not_found", "no_data"}:
                completed[uniprot] = record
            if index % 25 == 0:
                log(
                    f"{datetime.now(timezone.utc).isoformat()} "
                    f"progress={index}/{len(targets)} "
                    f"ligands={sum(r['ligand_count'] for r in completed.values())}"
                )
            time.sleep(0.4)

    write_outputs(targets, completed)
    log(
        f"{datetime.now(timezone.utc).isoformat()} "
        f"complete targets={len(completed)}/{len(targets)}"
    )


if __name__ == "__main__":
    main()
