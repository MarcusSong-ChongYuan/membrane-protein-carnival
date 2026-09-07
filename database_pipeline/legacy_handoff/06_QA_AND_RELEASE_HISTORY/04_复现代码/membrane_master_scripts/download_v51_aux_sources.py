from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
RAW = ROOT / "v51_review_working" / "raw"
USER_AGENT = "HumanMembraneProteinMaster/5.1 (evidence review)"

SOURCES = [
    (
        "tmbed_human_2022_predictions.txt",
        "https://zenodo.org/api/records/14705941/files/"
        "human_210422_tmbed.preds/content",
        "TMbed precomputed human-proteome predictions dated 2022-04-21; "
        "used only as independent computational evidence",
    ),
    (
        "goa_human_2026-07-08.gaf.gz",
        "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz",
        "Current human Gene Ontology annotations with evidence codes",
    ),
    (
        "go_basic_2026-06-26.obo",
        "https://purl.obolibrary.org/obo/go/go-basic.obo",
        "Gene Ontology term names, namespaces, and parent relations",
    ),
    (
        "interpro_entry_list_2026-06-10.tsv",
        "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list",
        "Current InterPro accession, type, and name list",
    ),
    (
        "interpro_short_names_2026-06-10.tsv",
        "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/short_names.dat",
        "Current InterPro short names",
    ),
    (
        "interpro2go_2026-06-10.txt",
        "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go",
        "Current InterPro-to-GO mappings",
    ),
    (
        "interpro_release_notes_2026-06-10.txt",
        "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/release_notes.txt",
        "InterPro release metadata",
    ),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path, attempts: int = 6) -> dict[str, str]:
    last_error: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": USER_AGENT, "Accept": "*/*"},
        )
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                with destination.open("wb") as output:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        output.write(block)
                return {
                    "http_date": response.headers.get("Date", ""),
                    "last_modified": response.headers.get("Last-Modified", ""),
                    "etag": response.headers.get("ETag", ""),
                }
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if attempt + 1 == attempts:
                break
            time.sleep(min(30, 2 ** attempt))
    raise RuntimeError(f"Failed after {attempts} attempts: {url}") from last_error


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    records = []
    for filename, url, role in SOURCES:
        destination = RAW / filename
        headers = download(url, destination)
        record = {
            "file": filename,
            "url": url,
            "role": role,
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "bytes": destination.stat().st_size,
            "sha256": sha256(destination),
            **headers,
        }
        records.append(record)
        print(f"Downloaded {filename}: {record['bytes']} bytes")

    manifest_path = RAW / "v51_aux_source_manifest.json"
    manifest_path.write_text(
        json.dumps({"sources": records}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"manifest": str(manifest_path), "sources": len(records)}, indent=2))


if __name__ == "__main__":
    main()
