#!/usr/bin/env python3
"""Freeze direct UniProt responses for external complex proteins missing in the bulk proteome."""

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


def get(url: str):
    error = ""
    for attempt in range(5):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "MemPro-identity-layer/0.2"})
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.status, response.read(), {k.lower(): v for k, v in response.headers.items()}
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return 404, b"", {}
            error = str(exc)
        except (urllib.error.URLError, TimeoutError) as exc:
            error = str(exc)
        time.sleep(min(2 ** attempt, 10))
    return 0, error.encode(), {}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--external-entities", type=Path, required=True)
    p.add_argument("--output-fasta", type=Path, required=True)
    p.add_argument("--output-manifest", type=Path, required=True)
    args = p.parse_args()
    with gzip.open(args.external_entities, "rt", encoding="utf-8", newline="") as handle:
        missing = sorted({row["canonical_uniprot_accession"] for row in csv.DictReader(handle, delimiter="\t") if row["sequence_status"] != "resolved"})
    records, chunks = [], []
    for accession in missing:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        status, body, headers = get(url)
        valid = status == 200 and body.startswith(b">")
        returned = body.split(b"|", 2)[1].decode() if valid and b"|" in body.splitlines()[0] else ""
        if valid:
            chunks.append(body.rstrip() + b"\n")
        records.append({
            "requested_accession": accession,
            "returned_accession": returned,
            "url": url,
            "http_status": status,
            "valid_fasta": int(valid),
            "bytes": len(body) if valid else 0,
            "sha256": hashlib.sha256(body).hexdigest() if valid else "",
            "x_uniprot_release": headers.get("x-uniprot-release", ""),
            "x_uniprot_release_date": headers.get("x-uniprot-release-date", ""),
        })
    args.output_fasta.write_bytes(b"".join(chunks))
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "requested_count": len(missing),
        "valid_fasta_count": sum(r["valid_fasta"] for r in records),
        "unresolved_count": sum(1 for r in records if not r["valid_fasta"]),
        "records": records,
    }
    args.output_manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("requested_count", "valid_fasta_count", "unresolved_count")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
