#!/usr/bin/env python3
"""Preflight audit for the completed MemProDB incremental source staging.

All published baselines and source staging files are read-only.  The script
only writes a compact QA report into the current run's qa directory.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import pathlib
import sys
from collections import Counter


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "incremental_20260727"
QA = RUN / "qa"

PUBCHEM_PROGRESS = QA / "PUBCHEM_NORMALIZATION_V3_PROGRESS.json"
PUBCHEM_DIR = RUN / "staging" / "pubchem" / "by_target"
GPCR = ROOT / "intermediate" / "gpcrdb_default_new_evidence_v2.tsv"
PDBe = RUN / "staging" / "pdbe" / "pdbe_normalized_evidence_v3.tsv.gz"
PDBBIND = RUN / "staging" / "pdbbind" / "pdbbind_normalized_evidence_v3.tsv.gz"
BRENDA = RUN / "staging" / "brenda" / "brenda_normalized_evidence_v3.tsv.gz"


def text_reader(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def audit_table(path: pathlib.Path, source: str) -> dict[str, object]:
    counts: Counter[str] = Counter()
    mapped_compounds: set[str] = set()
    targets: set[str] = set()
    with text_reader(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["rows"] += 1
            inclusion = row.get("default_release_inclusion", "") == "1"
            outcome = (row.get("activity_outcome") or "").strip().lower()
            mapping = (row.get("compound_mapping_status") or "").strip()
            target = (row.get("target_uniprot_id") or "").strip()
            compound = (row.get("compound_internal_id") or "").strip()
            if target:
                targets.add(target)
            if compound:
                mapped_compounds.add(compound)
            counts[f"mapping_{mapping or 'blank'}"] += 1
            counts[f"outcome_{outcome or 'blank'}"] += 1
            counts["default_rows" if inclusion else "nondefault_rows"] += 1
            if inclusion and not target:
                counts["default_missing_target"] += 1
            if inclusion and not compound:
                counts["default_missing_compound"] += 1
            if inclusion and mapping not in {"mapped_core", "mapped_extended"}:
                counts["default_not_uniquely_mapped"] += 1
    return {
        "source": source,
        "path": str(path),
        "counts": dict(sorted(counts.items())),
        "unique_targets": len(targets),
        "unique_mapped_compounds": len(mapped_compounds),
    }


def audit_pubchem() -> dict[str, object]:
    counts: Counter[str] = Counter()
    mapped_compounds: set[str] = set()
    targets: set[str] = set()
    cids: set[str] = set()
    files = sorted(PUBCHEM_DIR.glob("*.tsv.gz"))
    for number, path in enumerate(files, 1):
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                counts["rows"] += 1
                bucket = (row.get("release_bucket") or "").strip()
                mapping = (row.get("compound_mapping_status") or "").strip()
                outcome = (row.get("activity_outcome") or "").strip().lower()
                tier = (row.get("evidence_tier") or "").strip()
                target = (row.get("target_uniprot_id") or "").strip()
                compound = (row.get("compound_internal_id") or "").strip()
                cid = (row.get("pubchem_cid") or "").strip()
                inclusion = row.get("default_release_inclusion", "") == "1"
                counts[f"bucket_{bucket or 'blank'}"] += 1
                counts[f"mapping_{mapping or 'blank'}"] += 1
                counts[f"outcome_{outcome or 'blank'}"] += 1
                counts[f"tier_{tier or 'blank'}"] += 1
                counts["default_rows" if inclusion else "nondefault_rows"] += 1
                if target:
                    targets.add(target)
                if compound:
                    mapped_compounds.add(compound)
                if cid:
                    cids.add(cid)
                if inclusion and not target:
                    counts["default_missing_target"] += 1
                if inclusion and not compound:
                    counts["default_missing_compound"] += 1
                if inclusion and mapping not in {"mapped_core", "mapped_extended"}:
                    counts["default_not_uniquely_mapped"] += 1
        if number % 250 == 0:
            print(
                f"pubchem files={number}/{len(files)} rows={counts['rows']}",
                flush=True,
            )
    return {
        "source": "PubChem BioAssay",
        "path": str(PUBCHEM_DIR),
        "file_count": len(files),
        "counts": dict(sorted(counts.items())),
        "unique_targets": len(targets),
        "unique_pubchem_cids": len(cids),
        "unique_mapped_compounds": len(mapped_compounds),
    }


def main() -> int:
    required = [PUBCHEM_PROGRESS, GPCR, PDBe, PDBBIND, BRENDA]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + "; ".join(missing))

    pubchem_progress = json.loads(PUBCHEM_PROGRESS.read_text(encoding="utf-8"))
    sources = [
        audit_pubchem(),
        audit_table(GPCR, "GPCRdb curated source records"),
        audit_table(PDBe, "PDBe"),
        audit_table(PDBBIND, "PDBbind"),
        audit_table(BRENDA, "BRENDA"),
    ]

    blocking: list[str] = []
    if pubchem_progress.get("status") != "passed":
        blocking.append("PubChem normalization status is not passed")
    if pubchem_progress.get("download_completed_targets") != 10556:
        blocking.append("PubChem target coverage is not 10,556")
    if pubchem_progress.get("normalized_ok_targets") != pubchem_progress.get(
        "download_ok_targets"
    ):
        blocking.append("PubChem normalization has not caught up with downloads")
    for item in sources:
        counts = item["counts"]
        for field in (
            "default_missing_target",
            "default_missing_compound",
            "default_not_uniquely_mapped",
        ):
            if counts.get(field, 0):
                blocking.append(
                    f"{item['source']}: {field}={counts[field]}"
                )

    report = {
        "status": "passed" if not blocking else "failed",
        "run_type": "incremental_supplement_not_rebuild",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "baseline_policy": "read_only_no_overwrite",
        "pubchem_progress": pubchem_progress,
        "sources": sources,
        "blocking_findings": blocking,
    }
    output = QA / "INCREMENTAL_PREFLIGHT_V4_QA.json"
    temporary = output.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not blocking else 1


if __name__ == "__main__":
    sys.exit(main())
