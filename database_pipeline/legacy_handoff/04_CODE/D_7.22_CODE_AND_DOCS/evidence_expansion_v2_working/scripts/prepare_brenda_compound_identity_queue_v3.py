#!/usr/bin/env python3
"""Prepare a prioritized BRENDA compound identity verification queue.

Exact name matches are candidates only.  They do not mutate evidence rows or
trigger a canonical merge until an identifier/structure check confirms them.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import pathlib
import re
import sys
from collections import Counter, defaultdict


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "runs" / "incremental_20260727"
INVENTORY = (
    RUN_ROOT
    / "staging"
    / "brenda"
    / "brenda_compound_name_inventory_v3.tsv.gz"
)
BASELINE_PATHS = ROOT / "config" / "baseline_paths.json"
OUT_DIR = RUN_ROOT / "staging" / "brenda"


def normalize_name(value: str) -> str:
    value = value.replace("\u00ad", "").strip().casefold()
    return re.sub(r"\s+", " ", value)


def load_name_candidates():
    paths = json.loads(BASELINE_PATHS.read_text(encoding="utf-8"))
    master = pathlib.Path(paths["compound_master_v1_0"])
    evidence = pathlib.Path(paths["binding_evidence_v4_2"])
    names: dict[str, set[str]] = defaultdict(set)
    metadata: dict[str, dict[str, str]] = {}
    with master.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            internal = row["compound_internal_id"]
            preferred = normalize_name(row.get("preferred_name") or "")
            if preferred:
                names[preferred].add(internal)
            metadata[internal] = {
                "preferred_name": row.get("preferred_name", ""),
                "pubchem_cids": row.get("pubchem_cids", ""),
                "standard_inchikey": row.get("standard_inchikey", ""),
                "compound_scope_status": row.get("compound_scope_status", ""),
                "identity_confidence": row.get("identity_confidence", ""),
            }
    with evidence.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            internal = (row.get("compound_internal_id") or "").strip()
            name = normalize_name(row.get("compound_name") or "")
            if internal and name:
                names[name].add(internal)
    return names, metadata


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    name_candidates, metadata = load_name_candidates()
    output = OUT_DIR / "brenda_compound_identity_queue_v3.tsv.gz"
    temporary = output.with_suffix(output.suffix + ".tmp")
    fields = [
        "priority_rank",
        "compound_name_normalized",
        "preferred_source_name",
        "roles",
        "brenda_fields",
        "target_count",
        "evidence_row_count",
        "generic_or_common",
        "name_candidate_count",
        "candidate_compound_internal_ids",
        "candidate_preferred_names",
        "candidate_pubchem_cids",
        "candidate_standard_inchikeys",
        "identity_queue_status",
        "required_next_check",
    ]
    rows = []
    counts: Counter[str] = Counter()
    with gzip.open(INVENTORY, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["compound_name_normalized"]
            candidates = sorted(name_candidates.get(key, set()))
            generic = row.get("generic_or_common") == "1"
            roles = set(filter(None, row.get("roles", "").split(";")))
            evidence_rows = int(row.get("evidence_row_count") or 0)
            target_count = int(row.get("target_count") or 0)
            high_value = bool(
                roles & {"inhibitor", "activator", "substrate", "cofactor"}
            )
            if generic:
                status = "excluded_generic_or_common"
                required = "none unless specifically requested"
                priority_group = 9
            elif len(candidates) == 1:
                status = "unique_exact_name_candidate_unverified"
                required = (
                    "resolve BRENDA name to PubChem/structure and compare "
                    "against candidate CID/InChIKey"
                )
                priority_group = 1 if high_value else 2
            elif len(candidates) > 1:
                status = "ambiguous_exact_name_candidates"
                required = (
                    "resolve structure then select only the matching canonical "
                    "compound/form"
                )
                priority_group = 3
            else:
                status = "no_existing_name_candidate"
                required = (
                    "resolve name to a unique PubChem CID/structure; create a "
                    "new compound only if structure identity is unique"
                )
                priority_group = 4 if high_value else 5
            priority_score = (
                priority_group,
                -min(evidence_rows, 10**9),
                -min(target_count, 10**9),
                key,
            )
            rows.append(
                {
                    "_sort": priority_score,
                    "compound_name_normalized": key,
                    "preferred_source_name": row["preferred_source_name"],
                    "roles": row["roles"],
                    "brenda_fields": row["brenda_fields"],
                    "target_count": target_count,
                    "evidence_row_count": evidence_rows,
                    "generic_or_common": int(generic),
                    "name_candidate_count": len(candidates),
                    "candidate_compound_internal_ids": ";".join(candidates),
                    "candidate_preferred_names": ";".join(
                        metadata[value]["preferred_name"] for value in candidates
                    ),
                    "candidate_pubchem_cids": ";".join(
                        metadata[value]["pubchem_cids"] for value in candidates
                    ),
                    "candidate_standard_inchikeys": ";".join(
                        metadata[value]["standard_inchikey"] for value in candidates
                    ),
                    "identity_queue_status": status,
                    "required_next_check": required,
                }
            )
            counts[status] += 1
    rows.sort(key=lambda row: row["_sort"])
    with gzip.open(
        temporary, "wt", encoding="utf-8", newline="", compresslevel=6
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for rank, row in enumerate(rows, 1):
            row["priority_rank"] = rank
            writer.writerow(row)
    temporary.replace(output)
    report = {
        "status": "passed",
        "run_type": "incremental_supplement_not_rebuild",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "compound_policy": "name matches are candidates only; no merge performed",
        "counts": dict(sorted(counts.items())),
        "queue_rows": len(rows),
        "output": str(output),
    }
    qa = RUN_ROOT / "qa" / "BRENDA_COMPOUND_IDENTITY_QUEUE_V3_QA.json"
    tmp = qa.with_suffix(".json.tmp")
    tmp.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    tmp.replace(qa)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
