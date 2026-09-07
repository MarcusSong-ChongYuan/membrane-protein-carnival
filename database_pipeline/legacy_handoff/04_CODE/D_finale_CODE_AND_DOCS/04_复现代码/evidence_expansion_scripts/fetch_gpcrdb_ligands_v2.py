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
PROTEINS = ROOT / "intermediate" / "protein_universe_abc_v2.tsv"
RECEPTORS = Path(
    r"D:\7.22\membrane_master_v5_working\raw\gpcrdb_receptorlist_2026-07-24.json"
)
RAW_DIR = ROOT / "raw" / "gpcrdb_202607" / "ligands"
CHECKPOINT = ROOT / "raw" / "gpcrdb_202607" / "gpcrdb_ligand_fetch_v2.jsonl"
INVENTORY = ROOT / "intermediate" / "gpcrdb_ligand_source_inventory_v2.tsv"
REPORT = ROOT / "reports" / "GPCRDB_LIGAND_FETCH_REPORT.json"
LOG = ROOT / "logs" / "gpcrdb_ligand_fetch_v2.log"
BASE = "https://gpcrdb.org/services/ligands"
USER_AGENT = "MemProDB-evidence-expansion/2.0 (scientific database integration)"


def load_completed() -> dict[str, dict]:
    completed = {}
    if not CHECKPOINT.exists():
        return completed
    with CHECKPOINT.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if record.get("entry_name") and record.get("status") in {"ok", "not_found"}:
                completed[record["entry_name"]] = record
    return completed


def fetch_json(entry_name: str) -> tuple[str, bytes, str]:
    url = f"{BASE}/{entry_name}/"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                payload = response.read()
            json.loads(payload.decode("utf-8"))
            return "ok", payload, url
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return "not_found", b"[]", url
            if error.code in {429, 500, 502, 503, 504} and attempt < 4:
                time.sleep(4 * (attempt + 1))
                continue
            return f"http_{error.code}", b"", url
        except Exception as error:
            if attempt < 4:
                time.sleep(4 * (attempt + 1))
                continue
            return f"error_{type(error).__name__}", b"", url
    return "retry_exhausted", b"", url


def log(message: str) -> None:
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    abc = set()
    with PROTEINS.open("r", encoding="utf-8-sig", newline="") as handle:
        abc.update(row["target_uniprot_id"] for row in csv.DictReader(handle, delimiter="\t"))

    with RECEPTORS.open("r", encoding="utf-8") as handle:
        receptor_rows = json.load(handle)
    targets = [
        row
        for row in receptor_rows
        if row.get("species") == "Homo sapiens" and row.get("accession") in abc
    ]
    targets.sort(key=lambda row: row["entry_name"])
    completed = load_completed()
    log(
        f"{datetime.now(timezone.utc).isoformat()} start targets={len(targets)} completed={len(completed)}"
    )

    with CHECKPOINT.open("a", encoding="utf-8") as checkpoint:
        for index, target in enumerate(targets, 1):
            entry_name = target["entry_name"]
            if entry_name in completed:
                continue
            status, payload, url = fetch_json(entry_name)
            rows = []
            if payload:
                rows = json.loads(payload.decode("utf-8"))
                with gzip.open(RAW_DIR / f"{entry_name}.json.gz", "wb", compresslevel=6) as handle:
                    handle.write(payload)

            source_counts = Counter((row.get("Source") or "(blank)").strip() for row in rows)
            ligand_type_counts = Counter(
                (row.get("Ligand type") or "(blank)").strip() for row in rows
            )
            record = {
                "entry_name": entry_name,
                "target_uniprot_id": target.get("accession", ""),
                "status": status,
                "row_count": len(rows),
                "source_counts": dict(source_counts),
                "ligand_type_counts": dict(ligand_type_counts),
                "url": url,
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            }
            checkpoint.write(json.dumps(record, ensure_ascii=False) + "\n")
            checkpoint.flush()
            completed[entry_name] = record
            if index % 25 == 0:
                log(
                    f"{datetime.now(timezone.utc).isoformat()} progress={index}/{len(targets)} rows={sum(r['row_count'] for r in completed.values())}"
                )
            time.sleep(0.35)

    records = [completed[target["entry_name"]] for target in targets if target["entry_name"] in completed]
    source_totals = Counter()
    with INVENTORY.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "entry_name",
            "target_uniprot_id",
            "status",
            "row_count",
            "source_database",
            "source_row_count",
            "retrieved_at_utc",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            if record["source_counts"]:
                for source, count in sorted(record["source_counts"].items()):
                    source_totals[source] += count
                    writer.writerow(
                        {
                            "entry_name": record["entry_name"],
                            "target_uniprot_id": record["target_uniprot_id"],
                            "status": record["status"],
                            "row_count": record["row_count"],
                            "source_database": source,
                            "source_row_count": count,
                            "retrieved_at_utc": record["retrieved_at_utc"],
                        }
                    )
            else:
                writer.writerow(
                    {
                        "entry_name": record["entry_name"],
                        "target_uniprot_id": record["target_uniprot_id"],
                        "status": record["status"],
                        "row_count": record["row_count"],
                        "source_database": "",
                        "source_row_count": 0,
                        "retrieved_at_utc": record["retrieved_at_utc"],
                    }
                )

    status_counts = Counter(record["status"] for record in records)
    report = {
        "status": "passed" if len(records) == len(targets) and set(status_counts) <= {"ok", "not_found"} else "incomplete",
        "human_gpcr_targets_in_membrane_universe": len(targets),
        "completed_targets": len(records),
        "targets_with_ligand_rows": sum(record["row_count"] > 0 for record in records),
        "total_ligand_activity_rows": sum(record["row_count"] for record in records),
        "fetch_status_counts": dict(status_counts),
        "source_row_counts": dict(source_totals),
        "raw_directory": str(RAW_DIR),
    }
    with REPORT.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    log(f"{datetime.now(timezone.utc).isoformat()} complete {json.dumps(report)}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
