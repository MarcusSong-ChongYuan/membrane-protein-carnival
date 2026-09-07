import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path(r"C:\github-repos\upload\normalized_tables\outputs\v3_1_completion")
FILES = {
    "binary": BASE / "drug_protein_binary_relationships_v3_1.tsv",
    "protein": BASE / "protein_database_v3_1.tsv",
    "molecule": BASE / "small_molecule_database_v3_1.tsv",
    "pocket": BASE / "pocket_instances_v3_1.tsv",
}
csv.field_size_limit(100_000_000)


def load(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def shape(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        r = csv.reader(f, delimiter="\t")
        h = next(r)
        bad = 0
        n = 0
        for row in r:
            n += 1
            bad += len(row) != len(h)
    return {"rows": n, "cols": len(h), "malformed": bad}


def main():
    binary = load(FILES["binary"])
    proteins = load(FILES["protein"])
    molecules = load(FILES["molecule"])
    pockets = load(FILES["pocket"])
    pids = {r["target_uniprot_id"] for r in proteins}
    mids = {r["drug_id"] for r in molecules}
    bpid = {r["target_uniprot_id"] for r in binary}
    bmid = {r["drug_id"] for r in binary}
    relids = [r["relationship_context_id"] for r in binary]
    pocket_counts = Counter(r["target_uniprot_id"] for r in pockets)
    pocket_mismatch = [r["target_uniprot_id"] for r in proteins if int(r["n_pockets_alphafold"] or 0) != pocket_counts.get(r["target_uniprot_id"], 0)]
    suspicious = [r for r in binary if r["record_qc_status"] != "ok"]
    orphan = [r for r in molecules if r["record_qc_status"] == "orphan_no_current_relationship"]
    report = {
        "shapes": {k: shape(v) for k, v in FILES.items()},
        "unique_keys": {
            "protein_rows_minus_unique_ids": len(proteins) - len(pids),
            "molecule_rows_minus_unique_ids": len(molecules) - len(mids),
            "pocket_rows_minus_unique_ids": len(pockets) - len({r["pocket_id"] for r in pockets}),
            "relationship_rows_minus_unique_ids": len(binary) - len(set(relids)),
        },
        "foreign_keys": {
            "binary_proteins_missing": sorted(bpid - pids),
            "binary_compounds_missing": sorted(bmid - mids),
            "pocket_proteins_missing": sorted({r["target_uniprot_id"] for r in pockets} - pids),
        },
        "pocket_count_mismatches": pocket_mismatch,
        "uniprot": {
            "returned": sum(r["record_qc_status"] != "uniprot_record_not_returned" for r in proteins),
            "family_filled": sum(bool(r["protein_family"]) for r in proteins),
            "alphafold_filled": sum(bool(r["alphafold_model_id"]) for r in proteins),
            "membrane_scope_review": sum(r["record_qc_status"] == "review_membrane_scope" for r in proteins),
        },
        "relationship_qc": {
            "flagged_rows": len(suspicious),
            "flagged_types": dict(Counter(r["activity_type"] for r in suspicious)),
            "evidence_types": dict(Counter(r["evidence_type"] for r in binary)),
            "activity_types": dict(Counter(r["activity_type"] for r in binary if r["activity_type"])),
        },
        "orphan_compounds": [{"drug_id": r["drug_id"], "drug_name": r["drug_name"], "source_row_count": r["source_row_count"], "unique_target_count": r["unique_target_count"]} for r in orphan],
        "spot_checks": {
            "P08172": next(({k: r[k] for k in ["protein_class", "protein_family", "protein_subfamily", "transmembrane_count", "transmembrane_topology_summary", "alphafold_model_id"]} for r in proteins if r["target_uniprot_id"] == "P08172"), None),
            "IC50_less_than": [{k: r[k] for k in ["target_uniprot_id", "drug_id", "activity_type", "activity_relation", "activity_value_uM", "activity_value_unit"]} for r in binary if r["activity_relation"] == "<"][:5],
        },
    }
    ok = (
        all(x["malformed"] == 0 for x in report["shapes"].values())
        and all(x == 0 for x in report["unique_keys"].values())
        and all(not x for x in report["foreign_keys"].values())
        and not pocket_mismatch
        and report["uniprot"]["returned"] == len(proteins)
    )
    report["validation_passed"] = ok
    (BASE / "V3_1_VALIDATION_REPORT.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
