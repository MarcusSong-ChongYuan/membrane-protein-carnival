#!/usr/bin/env python3
"""Query ChEMBL targets for v5.3 proteins absent from the earlier v4 universe."""

from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"D:\7.22")
V53 = (
    ROOT
    / "membrane_master_v5_working"
    / "releases"
    / "release_v5_3_sequences"
    / "human_membrane_audit_master_v5_3.tsv"
)
V4 = (
    ROOT
    / "v4_working"
    / "releases"
    / "release_v4"
    / "human_membrane_protein_master_v4.tsv"
)
WORK = ROOT / "membrane_master_v5_working" / "v54_pgd_binding_working"
RAW = WORK / "raw" / "chembl_37_v54_supplemental_targets.jsonl"
OUT = WORK / "intermediate" / "external_chembl_target_index_v54_supplement.tsv"
REPORT = WORK / "reports" / "CHEMBL_TARGET_INDEX_V54_SUPPLEMENT_REPORT.json"
API = "https://www.ebi.ac.uk/chembl/api/data/target.json"


def read_ids(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return {
            row["target_uniprot_id"].strip()
            for row in csv.DictReader(handle, delimiter="\t")
        }


def fetch(url: str) -> dict:
    request = urllib.request.Request(
        url, headers={"User-Agent": "human-membrane-database-v5.4/1.0"}
    )
    for attempt in range(1, 7):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.loads(response.read())
        except Exception:
            if attempt == 6:
                raise
            time.sleep(min(3 * attempt, 15))
    raise RuntimeError("unreachable")


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    RAW.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    accessions = sorted(read_ids(V53) - read_ids(V4))
    target_by_id: dict[str, dict] = {}
    accession_to_targets: dict[str, set[str]] = defaultdict(set)
    with RAW.open("w", encoding="utf-8") as raw_handle:
        for start in range(0, len(accessions), 250):
            batch = accessions[start : start + 250]
            params = urllib.parse.urlencode(
                {
                    "target_components__accession__in": ",".join(batch),
                    "limit": 1000,
                }
            )
            url = f"{API}?{params}"
            while url:
                payload = fetch(url)
                for target in payload.get("targets", []):
                    if target.get("organism") != "Homo sapiens":
                        continue
                    target_id = target["target_chembl_id"]
                    target_by_id[target_id] = target
                    for component in target.get("target_components", []):
                        accession = component.get("accession", "")
                        if accession in batch:
                            accession_to_targets[accession].add(target_id)
                    raw_handle.write(
                        json.dumps(target, separators=(",", ":")) + "\n"
                    )
                url = payload.get("page_meta", {}).get("next") or ""
                if url.startswith("/"):
                    url = f"https://www.ebi.ac.uk{url}"

    header = [
        "target_uniprot_id",
        "chembl_target_ids",
        "chembl_single_protein_target_ids",
        "chembl_target_types",
        "chembl_target_names",
        "chembl_target_count",
        "source_database",
        "source_version",
    ]
    rows = []
    for accession in sorted(accession_to_targets):
        ids = sorted(accession_to_targets[accession])
        rows.append(
            {
                "target_uniprot_id": accession,
                "chembl_target_ids": ";".join(ids),
                "chembl_single_protein_target_ids": ";".join(
                    target_id
                    for target_id in ids
                    if target_by_id[target_id].get("target_type")
                    == "SINGLE PROTEIN"
                ),
                "chembl_target_types": ";".join(
                    sorted(
                        {
                            target_by_id[target_id].get("target_type", "")
                            for target_id in ids
                        }
                    )
                ),
                "chembl_target_names": ";".join(
                    sorted(
                        {
                            target_by_id[target_id].get("pref_name", "")
                            for target_id in ids
                        }
                    )
                ),
                "chembl_target_count": len(ids),
                "source_database": "ChEMBL",
                "source_version": "ChEMBL 37",
            }
        )
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "status": "passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "v53_accessions_absent_from_v4_universe_queried": len(accessions),
        "accessions_with_chembl_target_hit": len(rows),
        "single_protein_target_mappings": sum(
            bool(row["chembl_single_protein_target_ids"]) for row in rows
        ),
        "unique_chembl_targets": len(target_by_id),
        "output": str(OUT),
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
