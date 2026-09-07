from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
RAW = WORK / "raw" / "chembl_37"
INTERMEDIATE = WORK / "intermediate"
REPORTS = WORK / "reports"
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/target.json"
VERSION = "ChEMBL 37"
BATCH_SIZE = 250


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "mempro1-v4-build/1.0"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("ChEMBL request failed")


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    membrane = read_tsv(INTERMEDIATE / "human_membrane_protein_master.tsv")
    accessions = sorted(row["target_uniprot_id"] for row in membrane)
    raw_jsonl = RAW / "human_membrane_target_index.jsonl"
    targets_by_id: dict[str, dict] = {}
    accession_to_targets: dict[str, set[str]] = defaultdict(set)

    with raw_jsonl.open("w", encoding="utf-8") as raw_handle:
        for start in range(0, len(accessions), BATCH_SIZE):
            batch = accessions[start : start + BATCH_SIZE]
            params = urllib.parse.urlencode(
                {
                    "target_components__accession__in": ",".join(batch),
                    "limit": 1000,
                }
            )
            next_url = f"{BASE_URL}?{params}"
            while next_url:
                data = fetch_json(next_url)
                for target in data.get("targets", []):
                    if target.get("organism") != "Homo sapiens":
                        continue
                    target_id = target["target_chembl_id"]
                    targets_by_id[target_id] = target
                    for component in target.get("target_components", []):
                        accession = component.get("accession", "")
                        if accession in batch:
                            accession_to_targets[accession].add(target_id)
                    raw_handle.write(json.dumps(target, ensure_ascii=False) + "\n")
                next_url = data.get("page_meta", {}).get("next") or ""
                if next_url and next_url.startswith("/"):
                    next_url = f"https://www.ebi.ac.uk{next_url}"
            print(f"chembl_batch={start // BATCH_SIZE + 1}")

    output = []
    for accession in sorted(accession_to_targets):
        target_ids = sorted(accession_to_targets[accession])
        target_types = sorted({targets_by_id[target_id].get("target_type", "") for target_id in target_ids})
        names = sorted({targets_by_id[target_id].get("pref_name", "") for target_id in target_ids})
        single_ids = [
            target_id
            for target_id in target_ids
            if targets_by_id[target_id].get("target_type") == "SINGLE PROTEIN"
        ]
        output.append(
            {
                "target_uniprot_id": accession,
                "chembl_target_ids": ";".join(target_ids),
                "chembl_single_protein_target_ids": ";".join(single_ids),
                "chembl_target_types": ";".join(target_types),
                "chembl_target_names": ";".join(names),
                "chembl_target_count": len(target_ids),
                "evidence_status": "present_in_chembl_target_index",
                "candidate_sm_level": "SM2_candidate",
                "candidate_sm_level_note": (
                    "Target-index presence is an activity-data discovery flag; "
                    "activity-level import is required before final promotion."
                ),
                "source_database": "ChEMBL",
                "source_version": VERSION,
            }
        )
    write_tsv(
        INTERMEDIATE / "external_chembl_target_index.tsv",
        output,
        list(output[0]),
    )
    digest = hashlib.sha256(raw_jsonl.read_bytes()).hexdigest()
    report = {
        "status": "passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "ChEMBL",
        "source_version": VERSION,
        "reviewed_membrane_accessions_queried": len(accessions),
        "reviewed_membrane_targets_present_in_chembl": len(output),
        "unique_chembl_targets_returned": len(targets_by_id),
        "single_protein_target_mappings": sum(
            bool(row["chembl_single_protein_target_ids"]) for row in output
        ),
        "raw_jsonl_bytes": raw_jsonl.stat().st_size,
        "raw_jsonl_sha256": digest,
        "integration_scope": "target-presence index; no new ChEMBL activity rows imported",
    }
    (REPORTS / "CHEMBL_TARGET_INDEX_INTEGRATION_REPORT.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
