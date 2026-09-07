#!/usr/bin/env python3
"""Validate and freeze the MemProDB V6 incremental release."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import pathlib
import shutil
import sys
from collections import Counter


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "incremental_20260727"
CANDIDATE = RUN / "release_candidate_v6_0"
RELEASES = ROOT / "releases"
FINAL = RELEASES / "release_mempro_v6_0_20260727"
QA_DIR = RUN / "qa"
RUN_MANIFEST = RUN / "RUN_MANIFEST.json"

META_NAMES = {
    "RELEASE_MANIFEST_v6_0.tsv",
    "SHA256SUMS_v6_0.txt",
    "RELEASE_INFO_v6_0.json",
    "V6_VALIDATION_REPORT.json",
}


def open_text(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def count_rows(path: pathlib.Path) -> int:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        return sum(1 for _ in reader)


def main() -> int:
    evidence_qa_path = QA_DIR / "INCREMENTAL_EVIDENCE_BUILD_V6_QA.json"
    master_qa_path = QA_DIR / "INCREMENTAL_MASTER_BUILD_V6_QA.json"
    for path in [evidence_qa_path, master_qa_path, RUN_MANIFEST]:
        if not path.exists():
            raise FileNotFoundError(path)
    evidence_qa = json.loads(evidence_qa_path.read_text(encoding="utf-8"))
    master_qa = json.loads(master_qa_path.read_text(encoding="utf-8"))
    run_manifest = json.loads(RUN_MANIFEST.read_text(encoding="utf-8"))
    blocking: list[str] = []
    counts: Counter[str] = Counter()
    if evidence_qa.get("status") != "passed":
        blocking.append("evidence build QA is not passed")
    if master_qa.get("status") != "passed":
        blocking.append("master build QA is not passed")

    required = [
        "human_membrane_protein_master_v6_0.tsv",
        "protein_gene_disease_relations_v6_0.tsv",
        "binding_evidence_master_v6_0.tsv.gz",
        "small_molecule_master_v1_1.tsv",
        "compound_form_hierarchy_v1_1.tsv",
        "protein_compound_summary_v6_0.tsv.gz",
        "binding_site_instances_v6_0.tsv.gz",
        "negative_binding_evidence_v1_0.tsv.gz",
        "binding_evidence_review_queue_v6_0.tsv.gz",
        "pubchem_context_conflicts_v1_0.tsv.gz",
        "quantitative_heterogeneity_audit_v1_0.tsv.gz",
        "incremental_source_coverage_v6_0.tsv",
        "README_v6_0.md",
        "METHODS_v6_0.md",
        "INCLUSION_POLICY_v6_0.md",
        "DATA_DICTIONARY_v6_0.md",
        "SOURCE_REGISTRY_v6_0.tsv",
    ]
    for name in required:
        if not (CANDIDATE / name).exists():
            blocking.append(f"missing release file: {name}")
    if blocking:
        raise RuntimeError("; ".join(blocking))

    protein_ids: set[str] = set()
    protein_path = CANDIDATE / "human_membrane_protein_master_v6_0.tsv"
    with protein_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("target_uniprot_id", "")
            counts["protein_rows"] += 1
            if not key or key in protein_ids:
                blocking.append(f"duplicate or blank protein key: {key}")
                if len(blocking) > 50:
                    break
            protein_ids.add(key)

    compound_ids: set[str] = set()
    compound_path = CANDIDATE / "small_molecule_master_v1_1.tsv"
    with compound_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("compound_internal_id", "")
            counts["compound_rows"] += 1
            if not key or key in compound_ids:
                blocking.append(f"duplicate or blank compound key: {key}")
                if len(blocking) > 50:
                    break
            compound_ids.add(key)

    form_ids: set[str] = set()
    form_path = CANDIDATE / "compound_form_hierarchy_v1_1.tsv"
    with form_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("compound_form_id", "")
            parent = row.get("compound_internal_id", "")
            counts["form_rows"] += 1
            if not key or key in form_ids:
                blocking.append(f"duplicate or blank form key: {key}")
            if parent not in compound_ids:
                counts["form_missing_parent"] += 1
            form_ids.add(key)
    if counts["form_missing_parent"]:
        blocking.append(
            f"forms with missing parent: {counts['form_missing_parent']}"
        )

    disease_ids: set[str] = set()
    disease_path = CANDIDATE / "protein_gene_disease_relations_v6_0.tsv"
    with disease_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("disease_relation_id", "")
            target = row.get("target_uniprot_id", "")
            counts["disease_rows"] += 1
            if not key or key in disease_ids:
                blocking.append(f"duplicate or blank disease relation: {key}")
            if target not in protein_ids:
                counts["disease_missing_protein"] += 1
            disease_ids.add(key)
    if counts["disease_missing_protein"]:
        blocking.append(
            "disease relations with missing protein: "
            f"{counts['disease_missing_protein']}"
        )

    evidence_ids: set[str] = set()
    binding_path = CANDIDATE / "binding_evidence_master_v6_0.tsv.gz"
    with gzip.open(
        binding_path, "rt", encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("evidence_id", "")
            target = row.get("target_uniprot_id", "")
            compound = row.get("compound_internal_id", "")
            counts["binding_rows"] += 1
            if row.get("default_release_inclusion_v4_2") == "1":
                counts["binding_default_rows"] += 1
            if not key or key in evidence_ids:
                counts["duplicate_evidence_ids"] += 1
            evidence_ids.add(key)
            if target not in protein_ids:
                counts["binding_missing_protein"] += 1
            if compound not in compound_ids:
                if (
                    not compound
                    and row.get("default_release_inclusion_v4_2") == "0"
                    and row.get("compound_mapping_status")
                    in {"unresolved", "excluded"}
                ):
                    counts[
                        "binding_allowed_nondefault_unresolved_compound"
                    ] += 1
                else:
                    counts["binding_missing_compound"] += 1
            form = row.get("compound_form_id", "")
            if form and form not in form_ids:
                counts["binding_missing_form"] += 1
    for field in [
        "duplicate_evidence_ids",
        "binding_missing_protein",
        "binding_missing_compound",
        "binding_missing_form",
    ]:
        if counts[field]:
            blocking.append(f"{field}: {counts[field]}")

    pair_keys: set[tuple[str, str]] = set()
    pair_path = CANDIDATE / "protein_compound_summary_v6_0.tsv.gz"
    with gzip.open(pair_path, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (
                row.get("target_uniprot_id", ""),
                row.get("compound_internal_id", ""),
            )
            counts["pair_rows"] += 1
            if key in pair_keys:
                counts["duplicate_pair_keys"] += 1
            pair_keys.add(key)
            if key[0] not in protein_ids:
                counts["pair_missing_protein"] += 1
            if key[1] not in compound_ids:
                counts["pair_missing_compound"] += 1
    for field in [
        "duplicate_pair_keys",
        "pair_missing_protein",
        "pair_missing_compound",
    ]:
        if counts[field]:
            blocking.append(f"{field}: {counts[field]}")

    site_ids: set[str] = set()
    site_path = CANDIDATE / "binding_site_instances_v6_0.tsv.gz"
    with gzip.open(site_path, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("binding_site_instance_id", "")
            counts["site_rows"] += 1
            if not key or key in site_ids:
                counts["duplicate_site_ids"] += 1
            site_ids.add(key)
            if row.get("target_uniprot_id", "") not in protein_ids:
                counts["site_missing_protein"] += 1
            site_compound = row.get("compound_internal_id", "")
            if site_compound not in compound_ids:
                if (
                    not site_compound
                    and row.get("compound_mapping_status")
                    in {"unresolved", "excluded"}
                    and row.get("site_merge_status_v60") == "baseline_site"
                ):
                    counts["site_allowed_unresolved_baseline_compound"] += 1
                else:
                    counts["site_missing_compound"] += 1
    for field in [
        "duplicate_site_ids",
        "site_missing_protein",
        "site_missing_compound",
    ]:
        if counts[field]:
            blocking.append(f"{field}: {counts[field]}")

    negative_path = CANDIDATE / "negative_binding_evidence_v1_0.tsv.gz"
    review_path = CANDIDATE / "binding_evidence_review_queue_v6_0.tsv.gz"
    counts["negative_rows"] = count_rows(negative_path)
    counts["review_rows"] = count_rows(review_path)
    expected = evidence_qa.get("counts", {})
    expected_binding = expected.get("baseline_binding_rows", 0) + expected.get(
        "new_binding_rows", 0
    )
    for label, actual, wanted in [
        ("binding_rows", counts["binding_rows"], expected_binding),
        ("negative_rows", counts["negative_rows"], expected.get("negative_rows")),
        ("review_rows", counts["review_rows"], expected.get("review_rows")),
        (
            "pair_rows",
            counts["pair_rows"],
            master_qa.get("counts", {}).get("protein_compound_pair_rows"),
        ),
    ]:
        if wanted is not None and actual != wanted:
            blocking.append(f"{label}: expected {wanted}, found {actual}")

    for item in run_manifest.get("baseline_validation", []):
        path = pathlib.Path(item["path"])
        actual_hash = sha256(path)
        counts["baseline_files_rehashed"] += 1
        if actual_hash != item.get("expected_sha256"):
            blocking.append(f"baseline hash changed: {path}")

    reproducibility = CANDIDATE / "reproducibility" / "scripts"
    reproducibility.mkdir(parents=True, exist_ok=True)
    for name in [
        "audit_incremental_ready_v4.py",
        "build_incremental_evidence_release_v6.py",
        "build_incremental_master_release_v6.py",
        "validate_and_freeze_release_v6.py",
    ]:
        shutil.copy2(ROOT / "scripts" / name, reproducibility / name)
    shutil.copy2(RUN_MANIFEST, CANDIDATE / "RUN_MANIFEST.json")

    report = {
        "status": "passed" if not blocking else "failed",
        "release": "MemProDB V6.0",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "counts": dict(sorted(counts.items())),
        "row_reconciliation": evidence_qa.get("row_reconciliation"),
        "baseline_integrity": (
            "all frozen input hashes unchanged" if not any(
                item.startswith("baseline hash changed") for item in blocking
            ) else "failed"
        ),
        "blocking_findings": blocking,
        "frozen": False,
        "final_release_dir": "",
    }
    qa_path = QA_DIR / "V6_VALIDATION_REPORT.json"
    qa_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (CANDIDATE / "V6_VALIDATION_REPORT.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if blocking:
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 1

    payload = sorted(
        path
        for path in CANDIDATE.rglob("*")
        if path.is_file() and path.name not in META_NAMES
    )
    manifest_rows = []
    for path in payload:
        manifest_rows.append(
            {
                "relative_path": path.relative_to(CANDIDATE).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest_path = CANDIDATE / "RELEASE_MANIFEST_v6_0.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["relative_path", "bytes", "sha256"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    checksum_path = CANDIDATE / "SHA256SUMS_v6_0.txt"
    with checksum_path.open("w", encoding="utf-8", newline="\n") as handle:
        for item in manifest_rows:
            handle.write(
                f"{item['sha256']}  {item['relative_path']}\n"
            )
    info = {
        "release": "MemProDB V6.0",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "validation_status": "passed",
        "payload_file_count": len(manifest_rows),
        "canonical_tables": {
            "protein": "human_membrane_protein_master_v6_0.tsv",
            "protein_gene_disease": (
                "protein_gene_disease_relations_v6_0.tsv"
            ),
            "binding_evidence": "binding_evidence_master_v6_0.tsv.gz",
            "small_molecule": "small_molecule_master_v1_1.tsv",
        },
    }
    (CANDIDATE / "RELEASE_INFO_v6_0.json").write_text(
        json.dumps(info, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    RELEASES.mkdir(parents=True, exist_ok=True)
    candidate_resolved = CANDIDATE.resolve()
    final_resolved = FINAL.resolve()
    if candidate_resolved.parent != RUN.resolve():
        raise RuntimeError("candidate directory escaped intended run root")
    if final_resolved.parent != RELEASES.resolve():
        raise RuntimeError("final directory escaped intended release root")
    if FINAL.exists():
        raise FileExistsError(FINAL)
    CANDIDATE.rename(FINAL)

    report.update(
        {
            "frozen": True,
            "final_release_dir": str(FINAL),
            "manifest_file_count": len(manifest_rows),
        }
    )
    qa_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (FINAL / "V6_VALIDATION_REPORT.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
