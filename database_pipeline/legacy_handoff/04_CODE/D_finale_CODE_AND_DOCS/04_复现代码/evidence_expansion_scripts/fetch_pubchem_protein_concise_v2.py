#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
import shutil
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
PROTEINS = ROOT / "intermediate" / "protein_universe_abc_v2.tsv"
RAW_DIR = ROOT / "raw" / "pubchem_202607" / "protein_concise"
CHECKPOINT = (
    ROOT / "raw" / "pubchem_202607" / "pubchem_protein_concise_v2.jsonl"
)
INVENTORY = ROOT / "intermediate" / "pubchem_protein_concise_inventory_v2.tsv"
REPORT = ROOT / "reports" / "PUBCHEM_PROTEIN_CONCISE_FETCH_REPORT.json"
LOG = ROOT / "logs" / "pubchem_protein_concise_v2.log"
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/protein/accession"
USER_AGENT = "MemProDB-evidence-expansion/2.0 (scientific database integration)"
MIN_FREE_BYTES = 12 * 1024**3


def log(message: str) -> None:
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def load_targets() -> list[dict[str, str]]:
    targets = []
    with PROTEINS.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            targets.append(
                {
                    "target_uniprot_id": row["target_uniprot_id"].strip(),
                    "approved_symbol": row.get("approved_symbol", "").strip(),
                }
            )
    targets.sort(key=lambda row: row["target_uniprot_id"])
    return targets


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
            if record.get("status") in {"ok", "no_data"}:
                completed[record["target_uniprot_id"]] = record
    return completed


def fetch(uniprot: str) -> tuple[str, bytes, str, str]:
    encoded = urllib.parse.quote(uniprot, safe="")
    url = f"{BASE}/{encoded}/concise/CSV"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept-Encoding": "gzip",
            "Accept": "text/csv",
        },
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=240) as response:
                payload = response.read()
                encoding = response.headers.get("Content-Encoding", "")
            return "ok", payload, encoding, url
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return "no_data", b"", "", url
            if error.code in {408, 429, 500, 502, 503, 504} and attempt < 4:
                time.sleep(5 * (attempt + 1))
                continue
            return f"http_{error.code}", b"", "", url
        except Exception as error:
            if attempt < 4:
                time.sleep(5 * (attempt + 1))
                continue
            return f"error_{type(error).__name__}", b"", "", url
    return "retry_exhausted", b"", "", url


def save_payload(uniprot: str, payload: bytes, encoding: str) -> Path:
    destination = RAW_DIR / f"{uniprot}.csv.gz"
    temporary = RAW_DIR / f"{uniprot}.csv.gz.part"
    if encoding.lower() == "gzip" or payload[:2] == b"\x1f\x8b":
        with temporary.open("wb") as handle:
            handle.write(payload)
    else:
        with gzip.open(temporary, "wb", compresslevel=6) as handle:
            handle.write(payload)
    temporary.replace(destination)
    return destination


def summarize_csv(path: Path) -> dict:
    outcomes = Counter()
    assay_types = Counter()
    activity_types = Counter()
    aids = set()
    cids = set()
    rows = 0
    quantitative_rows = 0
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"AID", "SID", "CID", "Activity Outcome"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"unexpected_csv_header:{reader.fieldnames}")
        for row in reader:
            rows += 1
            outcomes[(row.get("Activity Outcome") or "(blank)").strip()] += 1
            assay_types[(row.get("Assay Type") or "(blank)").strip()] += 1
            activity_types[(row.get("Activity Name") or "(blank)").strip()] += 1
            aid = (row.get("AID") or "").strip()
            cid = (row.get("CID") or "").strip()
            if aid:
                aids.add(aid)
            if cid:
                cids.add(cid)
            if (row.get("Activity Value [uM]") or "").strip():
                quantitative_rows += 1
    return {
        "row_count": rows,
        "quantitative_row_count": quantitative_rows,
        "unique_aid_count": len(aids),
        "unique_cid_count": len(cids),
        "outcome_counts": dict(outcomes),
        "assay_type_counts": dict(assay_types),
        "activity_type_counts": dict(activity_types),
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
        "row_count",
        "quantitative_row_count",
        "unique_aid_count",
        "unique_cid_count",
        "compressed_bytes",
        "source_url",
        "retrieved_at_utc",
    ]
    with INVENTORY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field, "") for field in fields})

    status_counts = Counter(record["status"] for record in records)
    outcome_counts = Counter()
    assay_type_counts = Counter()
    activity_type_counts = Counter()
    for record in records:
        outcome_counts.update(record.get("outcome_counts") or {})
        assay_type_counts.update(record.get("assay_type_counts") or {})
        activity_type_counts.update(record.get("activity_type_counts") or {})
    report = {
        "status": (
            "passed"
            if len(records) == len(targets)
            and set(status_counts) <= {"ok", "no_data"}
            else "incomplete"
        ),
        "target_rows": len(targets),
        "completed_targets": len(records),
        "targets_with_concise_rows": sum(record["row_count"] > 0 for record in records),
        "total_concise_rows": sum(record["row_count"] for record in records),
        "total_quantitative_rows": sum(
            record["quantitative_row_count"] for record in records
        ),
        "compressed_bytes": sum(record["compressed_bytes"] for record in records),
        "status_counts": dict(status_counts),
        "outcome_counts": dict(outcome_counts),
        "assay_type_counts": dict(assay_type_counts),
        "activity_type_counts": dict(activity_type_counts),
        "checkpoint_file": str(CHECKPOINT),
        "raw_directory": str(RAW_DIR),
        "inventory_file": str(INVENTORY),
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
            free_bytes = shutil.disk_usage(ROOT).free
            if free_bytes < MIN_FREE_BYTES:
                log(
                    f"{datetime.now(timezone.utc).isoformat()} "
                    f"stopped_low_disk free_bytes={free_bytes}"
                )
                break
            status, payload, encoding, url = fetch(uniprot)
            summary = {
                "row_count": 0,
                "quantitative_row_count": 0,
                "unique_aid_count": 0,
                "unique_cid_count": 0,
                "outcome_counts": {},
                "assay_type_counts": {},
                "activity_type_counts": {},
            }
            compressed_bytes = 0
            if status == "ok":
                try:
                    path = save_payload(uniprot, payload, encoding)
                    summary = summarize_csv(path)
                    compressed_bytes = path.stat().st_size
                except Exception as error:
                    status = f"parse_error_{type(error).__name__}"
            record = {
                **target,
                "status": status,
                **summary,
                "compressed_bytes": compressed_bytes,
                "source_url": url,
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            }
            checkpoint.write(json.dumps(record, ensure_ascii=False) + "\n")
            checkpoint.flush()
            if status in {"ok", "no_data"}:
                completed[uniprot] = record
            if index % 25 == 0:
                log(
                    f"{datetime.now(timezone.utc).isoformat()} "
                    f"progress={index}/{len(targets)} "
                    f"ok={sum(r['status'] == 'ok' for r in completed.values())} "
                    f"rows={sum(r['row_count'] for r in completed.values())} "
                    f"compressed_bytes={sum(r['compressed_bytes'] for r in completed.values())}"
                )
            time.sleep(0.25)
    write_outputs(targets, completed)
    log(
        f"{datetime.now(timezone.utc).isoformat()} "
        f"complete completed={len(completed)}/{len(targets)}"
    )


if __name__ == "__main__":
    main()
