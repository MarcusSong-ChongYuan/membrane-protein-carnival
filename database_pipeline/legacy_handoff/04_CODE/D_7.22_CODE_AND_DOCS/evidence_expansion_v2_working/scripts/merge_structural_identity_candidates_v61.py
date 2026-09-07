#!/usr/bin/env python3
"""Cross-source deduplication for PDBe and PDBbind identity candidates."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RUN = ROOT / "runs" / "incremental_v61_20260727"
PDBE = RUN / "staging" / "pdbe_identity" / "pdbe_ccd_identity_v61.tsv.gz"
PDBBIND = (
    RUN / "staging" / "pdbbind_identity" / "pdbbind_ligand_identity_v61.tsv.gz"
)
OUTDIR = RUN / "staging" / "structural_identity_merge"
QA = RUN / "qa" / "STRUCTURAL_IDENTITY_MERGE_V61_SUMMARY.json"


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tmp.replace(path)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    new: dict[str, dict] = {}
    forms: dict[str, dict] = {}
    exact_rows = defaultdict(int)

    with gzip.open(PDBE, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            resolution = row["identity_resolution"]
            evidence_rows = int(row["evidence_row_count"])
            if resolution == "exact_full_inchikey":
                exact_rows["PDBe"] += evidence_rows
            elif resolution == "new_structure_candidate" and row[
                "release_eligibility"
            ] == "eligible_small_molecule":
                key = row["inchikey"]
                item = new.setdefault(
                    key,
                    {
                        "inchikey": key,
                        "connectivity_key": row["connectivity_key"],
                        "smiles": row["smiles_stereo"] or row["smiles"],
                        "inchi": row["inchi"],
                        "formula": row["formula"],
                        "sources": set(),
                        "pdbe_component_ids": set(),
                        "pdbbind_complex_ids": set(),
                        "evidence_row_count": 0,
                    },
                )
                item["sources"].add("PDBe")
                item["pdbe_component_ids"].add(row["component_id"])
                item["evidence_row_count"] += evidence_rows
            elif resolution == "connectivity_only_form_or_stereo":
                key = row["connectivity_key"]
                item = forms.setdefault(
                    key,
                    {
                        "connectivity_key": key,
                        "observed_inchikeys": set(),
                        "matched_compound_internal_ids": set(),
                        "sources": set(),
                        "source_object_ids": set(),
                        "evidence_row_count": 0,
                    },
                )
                item["observed_inchikeys"].add(row["inchikey"])
                item["matched_compound_internal_ids"].update(
                    value
                    for value in row["matched_compound_internal_ids"].split(";")
                    if value
                )
                item["sources"].add("PDBe")
                item["source_object_ids"].add(f"PDBe:{row['component_id']}")
                item["evidence_row_count"] += evidence_rows

    with gzip.open(PDBBIND, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            resolution = row["identity_resolution"]
            if resolution == "exact_full_inchikey":
                exact_rows["PDBbind"] += 1
            elif resolution == "new_structure_candidate" and row[
                "structure_eligibility"
            ] in {
                "eligible_small_molecule",
                "covalent_separate_layer",
                "form_hierarchy_review",
            }:
                key = row["standard_inchikey"]
                item = new.setdefault(
                    key,
                    {
                        "inchikey": key,
                        "connectivity_key": row["connectivity_key"],
                        "smiles": row["standard_smiles"],
                        "inchi": row["standard_inchi"],
                        "formula": row["molecular_formula"],
                        "sources": set(),
                        "pdbe_component_ids": set(),
                        "pdbbind_complex_ids": set(),
                        "evidence_row_count": 0,
                    },
                )
                item["sources"].add("PDBbind")
                item["pdbbind_complex_ids"].add(row["complex_id"])
                item["evidence_row_count"] += 1
            elif resolution == "connectivity_only_form_or_stereo":
                key = row["connectivity_key"]
                item = forms.setdefault(
                    key,
                    {
                        "connectivity_key": key,
                        "observed_inchikeys": set(),
                        "matched_compound_internal_ids": set(),
                        "sources": set(),
                        "source_object_ids": set(),
                        "evidence_row_count": 0,
                    },
                )
                item["observed_inchikeys"].add(row["standard_inchikey"])
                item["matched_compound_internal_ids"].update(
                    value
                    for value in row["matched_compound_internal_ids"].split(";")
                    if value
                )
                item["sources"].add("PDBbind")
                item["source_object_ids"].add(f"PDBbind:{row['complex_id']}")
                item["evidence_row_count"] += 1

    new_path = OUTDIR / "structural_new_compound_clusters_v61.tsv.gz"
    new_header = [
        "inchikey",
        "connectivity_key",
        "smiles",
        "inchi",
        "formula",
        "sources",
        "source_count",
        "pdbe_component_count",
        "pdbe_component_ids",
        "pdbbind_complex_count",
        "pdbbind_complex_ids",
        "evidence_row_count",
        "provisional_action",
    ]
    with gzip.open(new_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=new_header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(new):
            item = new[key]
            writer.writerow(
                {
                    **item,
                    "sources": ";".join(sorted(item["sources"])),
                    "source_count": len(item["sources"]),
                    "pdbe_component_count": len(item["pdbe_component_ids"]),
                    "pdbe_component_ids": ";".join(
                        sorted(item["pdbe_component_ids"])
                    ),
                    "pdbbind_complex_count": len(item["pdbbind_complex_ids"]),
                    "pdbbind_complex_ids": ";".join(
                        sorted(item["pdbbind_complex_ids"])
                    ),
                    "provisional_action": "wait_for_pubchem_global_dedup",
                }
            )

    form_path = OUTDIR / "structural_form_stereo_audit_v61.tsv.gz"
    form_header = [
        "connectivity_key",
        "observed_inchikey_count",
        "observed_inchikeys",
        "matched_compound_count",
        "matched_compound_internal_ids",
        "sources",
        "source_object_count",
        "source_object_ids",
        "evidence_row_count",
        "required_action",
    ]
    with gzip.open(form_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=form_header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(forms):
            item = forms[key]
            writer.writerow(
                {
                    "connectivity_key": key,
                    "observed_inchikey_count": len(item["observed_inchikeys"]),
                    "observed_inchikeys": ";".join(
                        sorted(item["observed_inchikeys"])
                    ),
                    "matched_compound_count": len(
                        item["matched_compound_internal_ids"]
                    ),
                    "matched_compound_internal_ids": ";".join(
                        sorted(item["matched_compound_internal_ids"])
                    ),
                    "sources": ";".join(sorted(item["sources"])),
                    "source_object_count": len(item["source_object_ids"]),
                    "source_object_ids": ";".join(
                        sorted(item["source_object_ids"])
                    ),
                    "evidence_row_count": item["evidence_row_count"],
                    "required_action": "audit_parent_form_stereo_charge",
                }
            )

    both = sum(1 for item in new.values() if len(item["sources"]) > 1)
    summary = {
        "status": "complete",
        "completed_utc": now(),
        "exact_existing_evidence_rows": exact_rows,
        "new_structure_clusters": len(new),
        "new_clusters_seen_in_both_structural_sources": both,
        "form_or_stereo_connectivity_clusters": len(forms),
        "new_cluster_output": str(new_path),
        "form_audit_output": str(form_path),
    }
    write_json(QA, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
