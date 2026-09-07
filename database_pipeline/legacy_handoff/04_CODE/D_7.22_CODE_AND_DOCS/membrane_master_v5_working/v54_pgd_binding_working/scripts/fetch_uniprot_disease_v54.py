#!/usr/bin/env python3
"""Retrieve current UniProt disease annotations for the v5.3 protein universe.

The script is restartable. Each batch is cached as a TSV under raw/, and the
combined table is written only after every batch has succeeded.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"D:\7.22")
MASTER = (
    ROOT
    / "membrane_master_v5_working"
    / "releases"
    / "release_v5_3_sequences"
    / "human_membrane_audit_master_v5_3.tsv"
)
WORK = ROOT / "membrane_master_v5_working" / "v54_pgd_binding_working"
RAW = WORK / "raw" / "uniprot_disease_batches_75"
OUT = WORK / "raw" / "uniprot_disease_annotations_v54.tsv"
REPORT = WORK / "reports" / "UNIPROT_DISEASE_RETRIEVAL_REPORT.json"

BATCH_SIZE = 75
FIELDS = [
    "accession",
    "gene_primary",
    "cc_disease",
    "xref_mim",
    "xref_orphanet",
    "xref_disgenet",
    "xref_malacards",
    "xref_opentargets",
]
BASE_URL = "https://rest.uniprot.org/uniprotkb/search"


def read_accessions() -> list[str]:
    with MASTER.open("r", encoding="utf-8-sig", newline="") as handle:
        return [
            row["target_uniprot_id"].strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("target_uniprot_id", "").strip()
        ]


def fetch_batch(accessions: list[str], destination: Path) -> None:
    query = " OR ".join(f"accession:{accession}" for accession in accessions)
    params = urllib.parse.urlencode(
        {
            "query": f"({query})",
            "format": "tsv",
            "fields": ",".join(FIELDS),
            "size": 500,
        }
    )
    request = urllib.request.Request(
        f"{BASE_URL}?{params}",
        headers={"User-Agent": "human-membrane-database-v5.4/1.0"},
    )

    last_error: Exception | None = None
    for attempt in range(1, 7):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            if not payload.startswith(b"Entry\t"):
                raise RuntimeError("UniProt response does not contain the expected TSV header")
            destination.write_bytes(payload)
            return
        except (urllib.error.URLError, TimeoutError, RuntimeError) as exc:
            last_error = exc
            if attempt == 6:
                break
            time.sleep(min(3 * attempt, 15))
    raise RuntimeError(f"UniProt batch failed after retries: {last_error}")


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)

    accessions = read_accessions()
    batches = [
        accessions[index : index + BATCH_SIZE]
        for index in range(0, len(accessions), BATCH_SIZE)
    ]

    for index, batch in enumerate(batches, start=1):
        destination = RAW / f"batch_{index:03d}.tsv"
        if destination.exists() and destination.stat().st_size > 20:
            with destination.open("r", encoding="utf-8-sig") as handle:
                cached_rows = max(sum(1 for _ in handle) - 1, 0)
            if cached_rows == len(batch):
                continue
        fetch_batch(batch, destination)

    rows_by_accession: dict[str, dict[str, str]] = {}
    source_header: list[str] | None = None
    for batch_file in sorted(RAW.glob("batch_*.tsv")):
        with batch_file.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if source_header is None:
                source_header = reader.fieldnames
            for row in reader:
                accession = row.get("Entry", "").strip()
                if accession:
                    rows_by_accession[accession] = row

    if source_header is None:
        raise RuntimeError("No UniProt batch data found")

    output_header = [
        "target_uniprot_id",
        "gene_symbol_uniprot",
        "uniprot_disease_comment",
        "uniprot_mim_ids",
        "uniprot_orphanet_ids",
        "uniprot_disgenet_ids",
        "uniprot_malacards_ids",
        "uniprot_opentargets_ids",
        "uniprot_disease_source_url",
        "retrieval_date",
    ]
    today = datetime.now(timezone.utc).date().isoformat()
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for accession in accessions:
            source = rows_by_accession.get(accession, {})
            writer.writerow(
                {
                    "target_uniprot_id": accession,
                    "gene_symbol_uniprot": source.get("Gene Names (primary)", ""),
                    "uniprot_disease_comment": source.get(
                        "Involvement in disease", ""
                    ),
                    "uniprot_mim_ids": source.get("MIM", ""),
                    "uniprot_orphanet_ids": source.get("Orphanet", ""),
                    "uniprot_disgenet_ids": source.get("DisGeNET", ""),
                    "uniprot_malacards_ids": source.get("MalaCards", ""),
                    "uniprot_opentargets_ids": source.get("OpenTargets", ""),
                    "uniprot_disease_source_url": (
                        f"https://rest.uniprot.org/uniprotkb/{accession}.txt"
                    ),
                    "retrieval_date": today,
                }
            )

    disease_count = sum(
        bool(rows_by_accession.get(accession, {}).get("Involvement in disease", ""))
        for accession in accessions
    )
    missing = sorted(set(accessions) - set(rows_by_accession))
    report = {
        "status": "passed" if not missing else "warning",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "UniProtKB REST API",
        "source_fields": FIELDS,
        "protein_universe_rows": len(accessions),
        "returned_rows": len(rows_by_accession),
        "proteins_with_disease_comment": disease_count,
        "missing_accessions": missing,
        "batch_size": BATCH_SIZE,
        "batch_count": len(batches),
        "output": str(OUT),
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
