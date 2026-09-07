#!/usr/bin/env python3
"""Fetch standardized ChEMBL direct-binding activities for v5.3 proteins.

Scope is intentionally conservative: ChEMBL 37 single-protein targets already
mapped to the membrane universe, assay_type B, standard_flag 1, nM values, and
Ki/Kd endpoints. IC50/EC50 are excluded here because they are not intrinsically
direct binding measurements. The raw response objects are retained as JSONL and a
normalized TSV is produced for downstream evidence integration.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"D:\7.22")
TARGET_INDEX = (
    ROOT
    / "v4_working"
    / "releases"
    / "release_v4"
    / "external_chembl_target_index_v4.tsv"
)
WORK = ROOT / "membrane_master_v5_working" / "v54_pgd_binding_working"
RAW = WORK / "raw" / "chembl_37_direct_binding_batches"
OUT = WORK / "intermediate" / "chembl_37_binding_activities_v54.tsv"
REPORT = WORK / "reports" / "CHEMBL_BINDING_ACTIVITY_RETRIEVAL_REPORT.json"
API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
TARGET_BATCH_SIZE = 25
PAGE_LIMIT = 1000
MAX_WORKERS = 24

ONLY_FIELDS = [
    "activity_id",
    "assay_chembl_id",
    "assay_type",
    "bao_label",
    "document_chembl_id",
    "molecule_chembl_id",
    "target_chembl_id",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "data_validity_comment",
    "potential_duplicate",
    "src_id",
    "target_organism",
    "target_pref_name",
]


def read_target_map() -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = defaultdict(set)
    with TARGET_INDEX.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["target_uniprot_id"].strip()
            for target_id in row["chembl_single_protein_target_ids"].split(";"):
                target_id = target_id.strip()
                if target_id:
                    mapping[target_id].add(accession)
    return dict(mapping)


def request_json(url: str) -> dict:
    request = urllib.request.Request(
        url, headers={"User-Agent": "human-membrane-database-v5.4/1.0"}
    )
    last_error: Exception | None = None
    for attempt in range(1, 7):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                return json.loads(response.read())
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt == 6:
                break
            time.sleep(min(3 * attempt, 20))
    raise RuntimeError(f"ChEMBL request failed: {last_error}; URL={url}")


def batch_url(target_ids: list[str], offset: int) -> str:
    params = {
        "target_chembl_id__in": ",".join(target_ids),
        "assay_type": "B",
        "standard_type__in": "Ki,Kd",
        "standard_units": "nM",
        "standard_flag": "1",
        "limit": str(PAGE_LIMIT),
        "offset": str(offset),
        "only": ",".join(ONLY_FIELDS),
    }
    return f"{API}?{urllib.parse.urlencode(params)}"


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    target_map = read_target_map()
    target_ids = sorted(target_map)
    batches = [
        target_ids[index : index + TARGET_BATCH_SIZE]
        for index in range(0, len(target_ids), TARGET_BATCH_SIZE)
    ]

    def fetch_one(item: tuple[int, list[str]]) -> None:
        batch_index, target_batch = item
        done_file = RAW / f"batch_{batch_index:03d}.done.json"
        if done_file.exists():
            return
        batch_jsonl = RAW / f"batch_{batch_index:03d}.jsonl"
        offset = 0
        total = None
        written = 0
        with batch_jsonl.open("w", encoding="utf-8", newline="") as handle:
            while total is None or offset < total:
                payload = request_json(batch_url(target_batch, offset))
                activities = payload.get("activities", [])
                page_meta = payload.get("page_meta", {})
                total = int(page_meta.get("total_count", len(activities)))
                for activity in activities:
                    handle.write(json.dumps(activity, separators=(",", ":")) + "\n")
                written += len(activities)
                if not activities:
                    break
                offset += len(activities)
        done_file.write_text(
            json.dumps(
                {
                    "batch_index": batch_index,
                    "target_count": len(target_batch),
                    "activity_count": written,
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        list(executor.map(fetch_one, enumerate(batches, start=1)))

    header = [
        "chembl_activity_id",
        "target_uniprot_id",
        "chembl_target_id",
        "chembl_molecule_id",
        "chembl_assay_id",
        "chembl_document_id",
        "activity_type",
        "activity_relation",
        "activity_value_nM",
        "activity_unit",
        "pchembl_value",
        "assay_type",
        "bao_label",
        "target_assignment_status",
        "data_validity_comment",
        "potential_duplicate",
        "source_database",
        "source_version",
        "source_url",
    ]
    output_count = 0
    mapped_target_count: set[str] = set()
    with OUT.open("w", encoding="utf-8", newline="") as output_handle:
        writer = csv.DictWriter(
            output_handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for batch_jsonl in sorted(RAW.glob("batch_*.jsonl")):
            with batch_jsonl.open("r", encoding="utf-8") as input_handle:
                for line in input_handle:
                    activity = json.loads(line)
                    target_id = str(activity.get("target_chembl_id") or "")
                    accessions = sorted(target_map.get(target_id, set()))
                    assignment = (
                        "single_uniprot_mapping"
                        if len(accessions) == 1
                        else "multiple_uniprot_mappings"
                    )
                    for accession in accessions:
                        writer.writerow(
                            {
                                "chembl_activity_id": activity.get("activity_id", ""),
                                "target_uniprot_id": accession,
                                "chembl_target_id": target_id,
                                "chembl_molecule_id": activity.get(
                                    "molecule_chembl_id", ""
                                ),
                                "chembl_assay_id": activity.get(
                                    "assay_chembl_id", ""
                                ),
                                "chembl_document_id": activity.get(
                                    "document_chembl_id", ""
                                ),
                                "activity_type": activity.get("standard_type", ""),
                                "activity_relation": activity.get(
                                    "standard_relation", ""
                                ),
                                "activity_value_nM": activity.get(
                                    "standard_value", ""
                                ),
                                "activity_unit": activity.get("standard_units", ""),
                                "pchembl_value": activity.get("pchembl_value", ""),
                                "assay_type": activity.get("assay_type", ""),
                                "bao_label": activity.get("bao_label", ""),
                                "target_assignment_status": assignment,
                                "data_validity_comment": activity.get(
                                    "data_validity_comment", ""
                                ),
                                "potential_duplicate": activity.get(
                                    "potential_duplicate", ""
                                ),
                                "source_database": "ChEMBL",
                                "source_version": "ChEMBL 37",
                                "source_url": (
                                    "https://www.ebi.ac.uk/chembl/explore/activity/"
                                    f"{activity.get('activity_id', '')}"
                                ),
                            }
                        )
                        output_count += 1
                        mapped_target_count.add(accession)

    report = {
        "status": "passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "ChEMBL",
        "source_version": "ChEMBL 37",
        "target_index_rows": len(target_ids),
        "target_batches": len(batches),
        "parallel_workers": MAX_WORKERS,
        "normalized_activity_rows": output_count,
        "unique_membrane_proteins": len(mapped_target_count),
        "filters": {
            "target_type": "SINGLE PROTEIN (from existing mapped target index)",
            "assay_type": "B",
            "standard_flag": 1,
            "standard_type": ["Ki", "Kd"],
            "standard_units": "nM",
        },
        "output": str(OUT),
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
