from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
V51 = ROOT / "v51_review_working"
RAW = V51 / "raw"
BASELINE = ROOT / "review" / "human_membrane_review_queue_v5.tsv"
USER_AGENT = "HumanMembraneProteinMaster/5.1 (evidence review)"

FIELDS = [
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "reviewed",
    "length",
    "sequence",
    "ft_transmem",
    "ft_intramem",
    "ft_signal",
    "ft_lipid",
    "ft_topo_dom",
    "cc_subcellular_location",
    "go_c",
    "go_f",
    "go_id",
    "keyword",
    "xref_pdb",
    "xref_alphafolddb",
    "xref_interpro",
    "xref_pfam",
    "xref_smart",
    "xref_prosite",
    "cc_function",
    "cc_subunit",
    "cc_domain",
    "cc_ptm",
    "protein_families",
    "lit_pubmed_id",
    "date_modified",
    "version",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, destination: Path, attempts: int = 6) -> dict[str, str]:
    last_error: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(
            url,
            headers={
                "User-Agent": USER_AGENT,
                "Accept": "text/tab-separated-values",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                body = response.read()
                destination.write_bytes(body)
                return {
                    "http_date": response.headers.get("Date", ""),
                    "etag": response.headers.get("ETag", ""),
                    "last_modified": response.headers.get("Last-Modified", ""),
                }
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if attempt + 1 == attempts:
                break
            time.sleep(min(30, 2 ** attempt))
    raise RuntimeError(f"Failed after {attempts} attempts: {url}") from last_error


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    with BASELINE.open("r", encoding="utf-8-sig", newline="") as handle:
        accessions = [
            row["target_uniprot_id"].strip()
            for row in csv.DictReader(handle, delimiter="\t")
        ]

    if len(accessions) != 1109 or len(set(accessions)) != 1109:
        raise ValueError("The frozen v5 R queue must contain 1,109 unique accessions")

    batch_size = 40
    batch_files: list[Path] = []
    manifest_records: list[dict[str, object]] = []
    for batch_number, start in enumerate(range(0, len(accessions), batch_size), 1):
        batch = accessions[start : start + batch_size]
        query = "(" + " OR ".join(f"accession:{item}" for item in batch) + ")"
        params = urllib.parse.urlencode(
            {
                "query": query,
                "format": "tsv",
                "fields": ",".join(FIELDS),
                "size": 500,
            }
        )
        url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
        destination = RAW / f"uniprot_r_review_batch_{batch_number:03d}.tsv"
        headers = fetch(url, destination)
        batch_files.append(destination)
        manifest_records.append(
            {
                "source": "UniProtKB REST API",
                "url": url,
                "batch": batch_number,
                "accession_count_requested": len(batch),
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                "file": destination.name,
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
                **headers,
            }
        )
        print(
            f"UniProt batch {batch_number:02d}: "
            f"{start + 1}-{min(start + batch_size, len(accessions))}"
        )

    combined = RAW / f"uniprot_r_review_{date.today().isoformat()}.tsv"
    header: str | None = None
    row_count = 0
    with combined.open("w", encoding="utf-8", newline="") as output:
        for batch_file in batch_files:
            lines = batch_file.read_text(encoding="utf-8").splitlines()
            if not lines:
                continue
            if header is None:
                header = lines[0]
                output.write(header + "\n")
            elif lines[0] != header:
                raise ValueError(f"Header mismatch in {batch_file}")
            for line in lines[1:]:
                if line.strip():
                    output.write(line + "\n")
                    row_count += 1

    if row_count != 1109:
        raise ValueError(
            f"UniProt returned {row_count} rows for 1,109 frozen R accessions"
        )

    manifest = {
        "release": "v5.1 review evidence",
        "retrieval_date": date.today().isoformat(),
        "baseline_queue": str(BASELINE),
        "baseline_queue_sha256": sha256(BASELINE),
        "requested_unique_accessions": len(accessions),
        "combined_uniprot_rows": row_count,
        "combined_file": combined.name,
        "combined_sha256": sha256(combined),
        "fields": FIELDS,
        "batches": manifest_records,
    }
    (RAW / "v51_source_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "combined_file": str(combined),
                "rows": row_count,
                "sha256": manifest["combined_sha256"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
