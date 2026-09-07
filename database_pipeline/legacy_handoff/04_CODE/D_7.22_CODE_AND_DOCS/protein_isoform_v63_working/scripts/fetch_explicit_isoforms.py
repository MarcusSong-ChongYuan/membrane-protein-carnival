#!/usr/bin/env python3
"""Freeze UniProt FASTA responses for isoform IDs explicitly named by complexes."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--complex-components", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path, required=True)
    return parser.parse_args()


def fetch(url: str, attempts: int = 5) -> tuple[int, bytes, dict[str, str]]:
    last_error = ""
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "MemPro-identity-layer/0.1"})
            with urllib.request.urlopen(request, timeout=60) as response:
                headers = {key.lower(): value for key, value in response.headers.items()}
                return response.status, response.read(), headers
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = str(exc)
            time.sleep(min(2 ** attempt, 10))
    return 0, last_error.encode("utf-8"), {}


def main() -> None:
    args = parse_args()
    opener = gzip.open if args.complex_components.suffix == ".gz" else open
    with opener(args.complex_components, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        isoform_ids = sorted({
            (row.get("component_uniprot_isoform_id") or "").strip()
            for row in reader
            if (row.get("component_uniprot_isoform_id") or "").strip()
        })
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    records = []
    fasta_chunks: list[bytes] = []
    for index, isoform_id in enumerate(isoform_ids, start=1):
        url = f"https://rest.uniprot.org/uniprotkb/{isoform_id}.fasta"
        status, body, headers = fetch(url)
        valid = status == 200 and body.startswith(b">")
        if valid:
            fasta_chunks.append(body.rstrip() + b"\n")
        records.append({
            "isoform_uniprot_accession": isoform_id,
            "url": url,
            "http_status": status,
            "valid_fasta": int(valid),
            "bytes": len(body) if valid else 0,
            "sha256": hashlib.sha256(body).hexdigest() if valid else "",
            "x_uniprot_release": headers.get("x-uniprot-release", ""),
            "x_uniprot_release_date": headers.get("x-uniprot-release-date", ""),
            "error": "" if valid else body.decode("utf-8", errors="replace")[:500],
        })
        if index % 25 == 0:
            print(f"fetched {index}/{len(isoform_ids)}", flush=True)
    args.output_fasta.write_bytes(b"".join(fasta_chunks))
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "requested_isoform_count": len(isoform_ids),
        "valid_fasta_count": sum(row["valid_fasta"] for row in records),
        "failed_count": sum(1 for row in records if not row["valid_fasta"]),
        "records": records,
    }
    args.output_manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({key: manifest[key] for key in ("requested_isoform_count", "valid_fasta_count", "failed_count")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
