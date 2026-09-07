#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
PROTEINS = ROOT / "intermediate" / "protein_universe_abc_v2.tsv"
RAW_DIR = ROOT / "raw" / "pubchem_202607"
CHECKPOINT = RAW_DIR / "pubchem_target_aid_discovery_v2.jsonl"
OUTPUT = ROOT / "intermediate" / "pubchem_target_aids_v2.tsv"
REPORT = ROOT / "reports" / "PUBCHEM_AID_DISCOVERY_REPORT.json"
LOG = ROOT / "logs" / "pubchem_aid_discovery_v2.log"

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target"
USER_AGENT = "MemProDB-evidence-expansion/2.0 (scientific database integration)"


def split_numeric_ids(value: str) -> list[str]:
    ids = []
    for item in (value or "").replace("|", ";").split(";"):
        item = item.strip()
        if item.isdigit():
            ids.append(item)
    return sorted(set(ids))


def request_aids(namespace: str, identifier: str) -> tuple[str, list[int], str]:
    encoded = urllib.parse.quote(identifier, safe="")
    url = f"{BASE}/{namespace}/{encoded}/aids/TXT"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=75) as response:
                text = response.read().decode("utf-8", errors="replace")
            aids = sorted({int(line.strip()) for line in text.splitlines() if line.strip().isdigit()})
            return "ok" if aids else "empty", aids, url
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return "empty", [], url
            if error.code in {429, 500, 502, 503, 504} and attempt < 4:
                time.sleep(3 * (attempt + 1))
                continue
            return f"http_{error.code}", [], url
        except Exception as error:  # network failures are audited and resumable
            if attempt < 4:
                time.sleep(3 * (attempt + 1))
                continue
            return f"error_{type(error).__name__}", [], url
    return "retry_exhausted", [], url


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
            if record.get("target_uniprot_id"):
                completed[record["target_uniprot_id"]] = record
    return completed


def append_log(message: str) -> None:
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    completed = load_completed()
    targets = []
    with PROTEINS.open("r", encoding="utf-8-sig", newline="") as handle:
        targets.extend(csv.DictReader(handle, delimiter="\t"))

    append_log(
        f"{datetime.now(timezone.utc).isoformat()} start targets={len(targets)} completed={len(completed)}"
    )
    with CHECKPOINT.open("a", encoding="utf-8") as checkpoint:
        for index, row in enumerate(targets, 1):
            uniprot = row["target_uniprot_id"]
            if uniprot in completed:
                continue
            query_items = [("geneid", value) for value in split_numeric_ids(row.get("ncbi_gene_ids", ""))]
            if not query_items and row.get("approved_symbol"):
                query_items = [("genesymbol", row["approved_symbol"])]

            all_aids = set()
            queries = []
            statuses = []
            for namespace, identifier in query_items:
                status, aids, url = request_aids(namespace, identifier)
                all_aids.update(aids)
                queries.append(url)
                statuses.append(status)
                time.sleep(0.22)

            if not query_items:
                overall = "no_gene_query_identifier"
            elif any(status.startswith(("http_", "error_", "retry_")) for status in statuses):
                overall = "partial_or_error"
            elif all_aids:
                overall = "aids_found"
            else:
                overall = "no_aids"

            record = {
                "target_uniprot_id": uniprot,
                "approved_symbol": row.get("approved_symbol", ""),
                "ncbi_gene_ids": row.get("ncbi_gene_ids", ""),
                "query_urls": queries,
                "query_statuses": statuses,
                "aid_count": len(all_aids),
                "aids": sorted(all_aids),
                "overall_status": overall,
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            }
            checkpoint.write(json.dumps(record, ensure_ascii=False) + "\n")
            checkpoint.flush()
            completed[uniprot] = record
            if index % 100 == 0:
                append_log(
                    f"{datetime.now(timezone.utc).isoformat()} progress={index}/{len(targets)}"
                )

    records = [completed[row["target_uniprot_id"]] for row in targets if row["target_uniprot_id"] in completed]
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "target_uniprot_id",
            "approved_symbol",
            "ncbi_gene_ids",
            "aid_count",
            "aids",
            "overall_status",
            "query_statuses",
            "retrieved_at_utc",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "target_uniprot_id": record["target_uniprot_id"],
                    "approved_symbol": record["approved_symbol"],
                    "ncbi_gene_ids": record["ncbi_gene_ids"],
                    "aid_count": record["aid_count"],
                    "aids": ";".join(str(value) for value in record["aids"]),
                    "overall_status": record["overall_status"],
                    "query_statuses": ";".join(record["query_statuses"]),
                    "retrieved_at_utc": record["retrieved_at_utc"],
                }
            )

    status_counts = Counter(record["overall_status"] for record in records)
    unique_aids = {aid for record in records for aid in record["aids"]}
    report = {
        "status": "passed" if len(records) == len(targets) else "incomplete",
        "target_rows": len(targets),
        "completed_rows": len(records),
        "status_counts": dict(status_counts),
        "targets_with_aids": sum(record["aid_count"] > 0 for record in records),
        "unique_aids": len(unique_aids),
        "total_target_aid_links": sum(record["aid_count"] for record in records),
        "checkpoint_file": str(CHECKPOINT),
    }
    with REPORT.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    append_log(f"{datetime.now(timezone.utc).isoformat()} complete {json.dumps(report)}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
