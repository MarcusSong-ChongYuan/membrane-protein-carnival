#!/usr/bin/env python3
import csv
import hashlib
import json
import re
from pathlib import Path


WORK = Path(__file__).resolve().parents[1]
RELEASE = WORK / "release"
CROSSWALKS = WORK / "crosswalks"
REVIEW = WORK / "review"
REPORTS = WORK / "reports"
RAW = WORK / "raw"
V4 = WORK.parent / "v4_working" / "releases" / "release_v4"


def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


checks = []


def check(name, condition, detail):
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


stats = json.loads((REPORTS / "V5_BUILD_STATS.json").read_text(encoding="utf-8"))
master = read_tsv(RELEASE / "human_membrane_protein_master_v5.tsv")
integral = read_tsv(RELEASE / "human_integral_membrane_view_v5.tsv")
core = read_tsv(RELEASE / "human_core_membrane_view_v5.tsv")
extended = read_tsv(RELEASE / "human_extended_membrane_view_v5.tsv")
review = read_tsv(REVIEW / "human_membrane_review_queue_v5.tsv")
evidence = read_tsv(CROSSWALKS / "source_evidence_matrix_v5.tsv")
unreviewed = read_tsv(CROSSWALKS / "unreviewed_canonical_resolution_v5.tsv")
special82 = read_tsv(REVIEW / "unreviewed_82_no_reviewed_gene_symbol_match_v5.tsv")
legacy = read_tsv(CROSSWALKS / "legacy_377_resolution_v5.tsv")
topology = read_tsv(CROSSWALKS / "single_pass_resolution_v5.tsv")
classification = read_tsv(CROSSWALKS / "classification_resolution_v5.tsv")
v4 = read_tsv(V4 / "human_membrane_protein_master_v4.tsv")

master_ids = [row["target_uniprot_id"] for row in master]
master_set = set(master_ids)
v4_set = {row["target_uniprot_id"] for row in v4}
integral_set = {row["target_uniprot_id"] for row in integral}
core_set = {row["target_uniprot_id"] for row in core}
extended_set = {row["target_uniprot_id"] for row in extended}
review_set = {row["target_uniprot_id"] for row in review}

check("master_row_count", len(master) == stats["v5_union_rows"], str(len(master)))
check("unique_canonical_accession", len(master_ids) == len(master_set), str(len(master_set)))
check("all_v4_rows_retained", v4_set <= master_set, f"missing={len(v4_set-master_set)}")
check(
    "accession_format",
    all(re.fullmatch(r"[A-Z0-9]+", accession) for accession in master_ids),
    "canonical UniProt accession characters",
)
check(
    "species_9606",
    all(row["organism_id"] == "9606" for row in master),
    "all master rows are human",
)
check(
    "required_decision_fields",
    all(
        row["release_tier_v5"]
        and row["membrane_decision_v5"]
        and row["membrane_scope_v5"]
        and row["classification_status_v5"]
        for row in master
    ),
    "tier, decision, scope and classification status are non-empty",
)
check(
    "valid_tiers",
    {row["release_tier_v5"] for row in master} <= {"A", "B", "C", "R"},
    str(sorted({row["release_tier_v5"] for row in master})),
)
check(
    "tier_A_integral_scope",
    all(row["membrane_scope_v5"] == "integral_membrane" for row in integral),
    "Tier A final scope is integral_membrane",
)
check(
    "integral_view_exact",
    integral_set
    == {row["target_uniprot_id"] for row in master if row["release_tier_v5"] == "A"},
    str(len(integral_set)),
)
check(
    "core_view_exact",
    core_set
    == {
        row["target_uniprot_id"]
        for row in master
        if row["release_tier_v5"] in {"A", "B"}
    },
    str(len(core_set)),
)
check(
    "extended_view_exact",
    extended_set
    == {
        row["target_uniprot_id"]
        for row in master
        if row["release_tier_v5"] in {"A", "B", "C"}
    },
    str(len(extended_set)),
)
check(
    "review_view_exact",
    review_set
    == {row["target_uniprot_id"] for row in master if row["release_tier_v5"] == "R"},
    str(len(review_set)),
)
check("evidence_matrix_one_per_master", len(evidence) == len(master), str(len(evidence)))
check(
    "classification_one_per_master",
    len(classification) == len(master),
    str(len(classification)),
)
check("unreviewed_candidate_count", len(unreviewed) == 9823, str(len(unreviewed)))
check(
    "unreviewed_candidate_unique",
    len({row["candidate_uniprot_id"] for row in unreviewed}) == len(unreviewed),
    "candidate accession uniqueness",
)
check("special_82_count", len(special82) == 82, str(len(special82)))
check("legacy_377_count", len(legacy) == 377, str(len(legacy)))
check(
    "legacy_identifier_resolved",
    all(row["canonical_uniprot_id_v5"] for row in legacy),
    "all legacy identifiers have a canonical result",
)
check("single_pass_680_count", len(topology) == 680, str(len(topology)))
check(
    "single_pass_allowed_outcomes",
    {
        row["resolved_topology_v5"] for row in topology
    }
    <= {
        "single_pass_type_i",
        "single_pass_type_ii",
        "single_pass_type_iii",
        "single_pass_type_i_or_iii",
        "single_pass_unresolved",
    },
    str(sorted({row["resolved_topology_v5"] for row in topology})),
)

manifest = json.loads((RAW / "source_download_manifest.json").read_text(encoding="utf-8-sig"))
hash_failures = []
for item in manifest:
    path = RAW / item["local_file"]
    if not path.exists() or sha256(path) != item["sha256"]:
        hash_failures.append(item["source_id"])
check("raw_source_hashes", not hash_failures, f"failures={hash_failures}")

passed = all(item["passed"] for item in checks)
report = {
    "validation_passed": passed,
    "check_count": len(checks),
    "passed_count": sum(item["passed"] for item in checks),
    "failed_count": sum(not item["passed"] for item in checks),
    "checks": checks,
}
(REPORTS / "V5_VALIDATION_REPORT.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(report, ensure_ascii=False, indent=2))
if not passed:
    raise SystemExit(1)
