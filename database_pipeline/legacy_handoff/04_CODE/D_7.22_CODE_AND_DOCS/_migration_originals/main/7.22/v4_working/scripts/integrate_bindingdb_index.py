from __future__ import annotations

import csv
import hashlib
import json
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
RAW = WORK / "raw" / "bindingdb_2026_07"
INTERMEDIATE = WORK / "intermediate"
REPORTS = WORK / "reports"
URL = "https://www.bindingdb.org/rwd/bind/BindingDB_UniProt.txt"
VERSION = "2026-07"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    path = RAW / "BindingDB_UniProt.txt"
    if not path.exists():
        request = urllib.request.Request(URL, headers={"User-Agent": "mempro1-v4-build/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            path.write_bytes(response.read())
    mappings = read_tsv(path)
    membrane = read_tsv(INTERMEDIATE / "human_membrane_protein_master.tsv")
    membrane_ids = {row["target_uniprot_id"] for row in membrane}
    by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mappings:
        for accession in row["UniProt"].replace(",", ";").split(";"):
            accession = accession.strip()
            if accession:
                by_accession[accession].append(row)
    output = []
    for accession in sorted(membrane_ids & set(by_accession)):
        records = by_accession[accession]
        output.append(
            {
                "target_uniprot_id": accession,
                "bindingdb_polymer_ids": ";".join(sorted({row["polymerid"] for row in records})),
                "bindingdb_target_names": ";".join(sorted({row["BindingDB Name"] for row in records})),
                "bindingdb_mapping_count": len(records),
                "evidence_status": "present_in_bindingdb_target_index",
                "candidate_sm_level": "SM4_candidate",
                "candidate_sm_level_note": (
                    "BindingDB target presence implies quantitative binding data; "
                    "measurement-level import is required before final SM4 promotion."
                ),
                "source_database": "BindingDB",
                "source_version": VERSION,
                "source_url": URL,
            }
        )
    write_tsv(
        INTERMEDIATE / "external_bindingdb_target_index.tsv",
        output,
        list(output[0]),
    )
    report = {
        "status": "passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "BindingDB",
        "source_version": VERSION,
        "mapping_file_bytes": path.stat().st_size,
        "mapping_file_sha256": sha256(path),
        "mapping_rows": len(mappings),
        "unique_uniprot_accessions_in_mapping": len(by_accession),
        "reviewed_membrane_targets_present_in_bindingdb": len(output),
        "integration_scope": "target-presence index; no measurement rows imported",
    }
    (REPORTS / "BINDINGDB_INDEX_INTEGRATION_REPORT.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
