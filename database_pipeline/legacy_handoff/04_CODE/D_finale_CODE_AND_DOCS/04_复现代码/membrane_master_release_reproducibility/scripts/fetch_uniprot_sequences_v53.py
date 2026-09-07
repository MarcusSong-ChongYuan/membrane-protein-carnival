from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
INPUT = (
    ROOT
    / "releases"
    / "release_v5_2_evidence_tiered"
    / "human_membrane_audit_master_v5_2.tsv"
)
WORKING = ROOT / "v53_sequence_working"
RAW = WORKING / "raw" / "uniprot_sequence_batches_75"
LOOKUP = WORKING / "uniprot_sequence_lookup_v53.tsv"
REPORT = WORKING / "uniprot_sequence_retrieval_report_v53.json"
API = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,id,reviewed,sequence,length,sequence_version"
RETRIEVAL_DATE = "2026-07-24"
CHUNK_SIZE = 75
MAX_WORKERS = 4


def chunks(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def request_tsv(query: str) -> tuple[str, dict[str, str], str]:
    params = {
        "format": "tsv",
        "fields": FIELDS,
        "size": 500,
        "query": query,
    }
    url = f"{API}?{urlencode(params)}"
    request = Request(
        url,
        headers={
            "User-Agent": "HumanMembraneProteinDatabase/5.3 sequence-retrieval "
            "(research dataset; UniProt canonical sequences)"
        },
    )
    last_error: Exception | None = None
    for attempt in range(7):
        try:
            with urlopen(request, timeout=120) as response:
                headers = {
                    "uniprot_release": response.headers.get("X-UniProt-Release", ""),
                    "uniprot_release_date": response.headers.get(
                        "X-UniProt-Release-Date", ""
                    ),
                    "api_deployment_date": response.headers.get(
                        "X-API-Deployment-Date", ""
                    ),
                }
                return response.read().decode("utf-8"), headers, url
        except HTTPError as error:
            last_error = error
            if error.code not in {429, 500, 502, 503, 504} or attempt == 6:
                raise
            retry_after = error.headers.get("Retry-After")
            wait = float(retry_after) if retry_after else 0.6 * (2**attempt)
            time.sleep(min(wait, 30))
        except (URLError, TimeoutError) as error:
            last_error = error
            if attempt == 6:
                raise
            time.sleep(min(0.6 * (2**attempt), 30))
    raise RuntimeError(f"UniProt request failed: {last_error}")


def fetch_batch(batch_index: int, accessions: list[str]) -> dict[str, object]:
    RAW.mkdir(parents=True, exist_ok=True)
    raw_path = RAW / f"batch_{batch_index:03d}.tsv"
    meta_path = RAW / f"batch_{batch_index:03d}.json"
    if raw_path.exists() and meta_path.exists():
        return {
            "batch_index": batch_index,
            "accessions": accessions,
            "text": raw_path.read_text(encoding="utf-8"),
            "metadata": json.loads(meta_path.read_text(encoding="utf-8")),
            "cached": True,
        }

    query = "(" + " OR ".join(f"accession:{item}" for item in accessions) + ")"
    text, headers, url = request_tsv(query)
    raw_path.write_text(text, encoding="utf-8")
    metadata = {
        "batch_index": batch_index,
        "requested_count": len(accessions),
        "request_url": url,
        **headers,
    }
    meta_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return {
        "batch_index": batch_index,
        "accessions": accessions,
        "text": text,
        "metadata": metadata,
        "cached": False,
    }


def parse_result(text: str) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def normalize_row(
    requested_accession: str,
    result: dict[str, str] | None,
    status: str,
    metadata: dict[str, str],
) -> dict[str, str]:
    if result is None:
        return {
            "requested_accession": requested_accession,
            "resolved_accession": "",
            "entry_name": "",
            "reviewed": "",
            "canonical_sequence": "",
            "sequence_length_api": "",
            "sequence_version": "",
            "sequence_sha256": "",
            "sequence_status": status,
            "isoform_scope": "canonical",
            "sequence_source": "UniProtKB REST API",
            "sequence_source_url": "",
            "uniprot_release": metadata.get("uniprot_release", ""),
            "uniprot_release_date": metadata.get("uniprot_release_date", ""),
            "sequence_retrieval_date": RETRIEVAL_DATE,
        }
    sequence = re.sub(r"\s+", "", result.get("Sequence", "")).upper()
    resolved = result.get("Entry", "")
    return {
        "requested_accession": requested_accession,
        "resolved_accession": resolved,
        "entry_name": result.get("Entry Name", ""),
        "reviewed": result.get("Reviewed", ""),
        "canonical_sequence": sequence,
        "sequence_length_api": result.get("Length", ""),
        "sequence_version": result.get("Sequence version", ""),
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if sequence
        else "",
        "sequence_status": status,
        "isoform_scope": "canonical",
        "sequence_source": "UniProtKB REST API",
        "sequence_source_url": f"https://rest.uniprot.org/uniprotkb/{resolved}.fasta"
        if resolved
        else "",
        "uniprot_release": metadata.get("uniprot_release", ""),
        "uniprot_release_date": metadata.get("uniprot_release_date", ""),
        "sequence_retrieval_date": RETRIEVAL_DATE,
    }


def resolve_missing(accession: str) -> tuple[str, dict[str, str] | None, dict[str, str]]:
    query = f"(accession:{accession} OR sec_acc:{accession})"
    text, metadata, _ = request_tsv(query)
    results = parse_result(text)
    if not results:
        return accession, None, metadata
    return accession, results[0], metadata


def main() -> None:
    WORKING.mkdir(parents=True, exist_ok=True)
    with INPUT.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    accessions = list(dict.fromkeys(row["target_uniprot_id"] for row in rows))
    if len(accessions) != 10_997:
        raise ValueError(f"Expected 10,997 unique accessions, found {len(accessions):,}")

    batches = list(enumerate(chunks(accessions, CHUNK_SIZE), start=1))
    fetched: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(fetch_batch, index, batch): index
            for index, batch in batches
        }
        for completed, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            fetched.append(result)
            print(
                f"Fetched batch {result['batch_index']:03d}; "
                f"{completed}/{len(batches)} complete",
                flush=True,
            )

    sequence_rows: dict[str, dict[str, str]] = {}
    releases: set[str] = set()
    release_dates: set[str] = set()
    for batch in fetched:
        metadata = batch["metadata"]
        if metadata.get("uniprot_release"):
            releases.add(str(metadata["uniprot_release"]))
        if metadata.get("uniprot_release_date"):
            release_dates.add(str(metadata["uniprot_release_date"]))
        for result in parse_result(str(batch["text"])):
            accession = result.get("Entry", "")
            if accession:
                sequence_rows[accession] = normalize_row(
                    accession, result, "retrieved_primary", metadata
                )

    missing = [accession for accession in accessions if accession not in sequence_rows]
    if missing:
        print(f"Resolving {len(missing)} missing/secondary accessions", flush=True)
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {
                executor.submit(resolve_missing, accession): accession
                for accession in missing
            }
            for future in as_completed(futures):
                requested, result, metadata = future.result()
                if result is None:
                    sequence_rows[requested] = normalize_row(
                        requested, None, "not_found", metadata
                    )
                else:
                    resolved = result.get("Entry", "")
                    status = (
                        "retrieved_primary"
                        if resolved == requested
                        else "retrieved_secondary_mapped"
                    )
                    sequence_rows[requested] = normalize_row(
                        requested, result, status, metadata
                    )

    ordered = [sequence_rows[accession] for accession in accessions]
    fieldnames = list(ordered[0])
    with LOOKUP.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(ordered)

    invalid_sequences = [
        row["requested_accession"]
        for row in ordered
        if row["canonical_sequence"]
        and not re.fullmatch(r"[A-Z]+", row["canonical_sequence"])
    ]
    report = {
        "requested_accessions": len(accessions),
        "retrieved_sequences": sum(bool(row["canonical_sequence"]) for row in ordered),
        "retrieved_primary": sum(
            row["sequence_status"] == "retrieved_primary" for row in ordered
        ),
        "retrieved_secondary_mapped": sum(
            row["sequence_status"] == "retrieved_secondary_mapped" for row in ordered
        ),
        "not_found": sum(row["sequence_status"] == "not_found" for row in ordered),
        "invalid_sequence_characters": invalid_sequences,
        "uniprot_releases_seen": sorted(releases),
        "uniprot_release_dates_seen": sorted(release_dates),
        "retrieval_date": RETRIEVAL_DATE,
        "canonical_only": True,
        "batch_count": len(batches),
    }
    REPORT.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
