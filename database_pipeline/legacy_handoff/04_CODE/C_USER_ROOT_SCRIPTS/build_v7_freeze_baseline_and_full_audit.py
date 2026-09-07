#!/usr/bin/env python3
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

BASE = Path(r"D:\finale\08_V6.4_database_freeze_and_figures_20260812")
ROOT = Path(r"D:\finale\09_V7_data_freeze_working_20260813")
for d in [ROOT/"00_baseline", ROOT/"01_rules", ROOT/"02_protein_audit", ROOT/"12_validation", ROOT/"logs"]:
    d.mkdir(parents=True, exist_ok=True)
protein = BASE/"01_release_overlay"/"human_membrane_protein_master_v6_4_candidate.tsv.gz"

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

files = list((BASE/"01_release_overlay").glob("*")) + [BASE/"RELEASE_MANIFEST_V6_4.tsv", BASE/"FREEZE_V6_4_INTERNAL.json"]
with (BASE/"RELEASE_MANIFEST_V6_4.tsv").open(encoding="utf-8-sig", newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["scope"] == "immutable_dependency":
            files.append(Path(row["relative_path"]))
files = list(dict.fromkeys(p.resolve() for p in files if p.exists()))
manifest = []
for path in files:
    manifest.append({"absolute_path": str(path), "bytes": path.stat().st_size, "sha256": sha(path), "baseline_role": "v64_overlay_or_declared_dependency", "read_only_policy": 1})
with (ROOT/"00_baseline"/"BASELINE_MANIFEST.tsv").open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=manifest[0].keys(), delimiter="\t")
    w.writeheader(); w.writerows(manifest)

def truth(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}

def num(value):
    try: return int(float(value or 0))
    except Exception: return 0

def lipid_basis(row):
    text = " ".join([row.get("membrane_evidence_basis", ""), row.get("membrane_decision_v5", ""), row.get("evidence_basis_v52", ""), row.get("subcellular_location", ""), row.get("uniprot_keywords", "")]).lower()
    return any(x in text for x in ["lipid-anchor", "lipid anchor", "gpi-anchor", "gpi anchor", "myristoyl", "palmitoyl", "prenyl"])

def peripheral_basis(row):
    text = " ".join([row.get("membrane_evidence_basis", ""), row.get("subcellular_location", ""), row.get("membrane_scope_v5", "")]).lower()
    return any(x in text for x in ["peripheral membrane", "peripheral_membrane", "membrane-associated", "membrane association"])

fields = ["target_uniprot", "approved_symbol", "protein_name", "reviewed", "organism_id", "current_class", "current_evidence", "current_release_scope", "website_default", "tm_count", "intramembrane_count", "signal_peptide_present", "lipid_anchor_basis", "peripheral_basis", "opm_integral", "pdbtm_present", "hpa_only_or_single_source", "cross_source_conflict", "isoform_scope", "docking_related", "initial_risk_priority", "initial_risk_flags", "initial_rule_status", "source_protein_table"]
out, seen, duplicates = [], set(), []
with gzip.open(protein, "rt", encoding="utf-8", newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc = row["target_uniprot_id"].strip(); flags = []
        if acc in seen: duplicates.append(acc)
        seen.add(acc)
        cls = row.get("membrane_class_v52", ""); ev = row.get("evidence_level_v52", "")
        tm = num(row.get("transmembrane_count_v5") or row.get("transmembrane_count")); intra = num(row.get("intramembrane_count"))
        lipid = lipid_basis(row); peripheral = peripheral_basis(row)
        opm = truth(row.get("opm_integral_structure_v5")); pdbtm = truth(row.get("pdbtm_present_v5"))
        sources = num(row.get("distinct_contributing_membrane_database_count_v63") or row.get("independent_membrane_source_count_v5"))
        source_text = (row.get("contributing_membrane_databases_v63", "") + row.get("independent_membrane_sources_v5", "")).lower()
        hpa_only = sources <= 1 and "hpa" in source_text
        conflict = truth(row.get("cross_source_conflict_v5"))
        if row.get("organism_id") != "9606": flags.append("NON_HUMAN_ORGANISM")
        if not re.match(r"^[A-Z0-9]{6,10}$", acc): flags.append("ACCESSION_FORMAT_REVIEW")
        integral = bool(tm or intra or opm or pdbtm)
        if cls == "A" and not integral: flags.append("A_WITHOUT_INTEGRAL_OR_INTRAMEMBRANE_EVIDENCE")
        if cls == "A" and lipid and not integral: flags.append("LIKELY_A_TO_B_LIPID_ANCHOR")
        if cls == "A" and peripheral and not integral: flags.append("LIKELY_A_TO_C_PERIPHERAL")
        if cls == "B" and not lipid: flags.append("B_WITHOUT_EXPLICIT_LIPID_ANCHOR")
        if cls == "C" and not peripheral and not lipid: flags.append("C_WITHOUT_EXPLICIT_PERIPHERAL_MODE")
        if hpa_only: flags.append("SINGLE_HPA_MEMBRANE_SUPPORT")
        if conflict: flags.append("CROSS_SOURCE_CONFLICT")
        if ev == "E3": flags.append("E3_CANDIDATE_NOT_DEFAULT")
        if row.get("reviewed", "").lower() != "reviewed": flags.append("UNREVIEWED_PROTEIN")
        if "A_WITHOUT_INTEGRAL_OR_INTRAMEMBRANE_EVIDENCE" in flags or "NON_HUMAN_ORGANISM" in flags: priority = "P0"
        elif conflict or "B_WITHOUT_EXPLICIT_LIPID_ANCHOR" in flags: priority = "P1"
        elif hpa_only or ev == "E2": priority = "P2"
        elif cls in {"B", "C"} or ev == "E3": priority = "P3"
        else: priority = "P4_CONTROL"
        out.append({"target_uniprot": acc, "approved_symbol": row.get("approved_symbol", ""), "protein_name": row.get("protein_name", ""), "reviewed": row.get("reviewed", ""), "organism_id": row.get("organism_id", ""), "current_class": cls, "current_evidence": ev, "current_release_scope": row.get("release_scope_v52", ""), "website_default": row.get("website_default_v52", ""), "tm_count": tm, "intramembrane_count": intra, "signal_peptide_present": int(bool(row.get("signal_peptide_features"))), "lipid_anchor_basis": int(lipid), "peripheral_basis": int(peripheral), "opm_integral": int(opm), "pdbtm_present": int(pdbtm), "hpa_only_or_single_source": int(hpa_only), "cross_source_conflict": int(conflict), "isoform_scope": row.get("isoform_scope_v53", ""), "docking_related": 0, "initial_risk_priority": priority, "initial_risk_flags": ";".join(flags), "initial_rule_status": "INITIAL_TRIAGE_NOT_FINAL_DECISION", "source_protein_table": str(protein)})

dest = ROOT/"02_protein_audit"/"protein_identity_audit_all_v7_initial.tsv"
with dest.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(out)
summary = {"protein_rows": len(out), "unique_accessions": len(seen), "duplicate_accessions": len(duplicates), "class_counts": dict(Counter(x["current_class"] for x in out)), "evidence_counts": dict(Counter(x["current_evidence"] for x in out)), "priority_counts": dict(Counter(x["initial_risk_priority"] for x in out)), "top_flags": dict(Counter(flag for x in out for flag in x["initial_risk_flags"].split(";") if flag)), "baseline_files": len(manifest), "blocking_errors": [] if len(out) == len(seen) == 10997 else ["PROTEIN_COUNT_OR_DUPLICATE_MISMATCH"], "status": "INITIAL_TRIAGE_COMPLETE_NOT_RELEASE"}
(ROOT/"02_protein_audit"/"FULL_PROTEIN_INITIAL_AUDIT_SUMMARY.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
(ROOT/"00_baseline"/"BASELINE_COUNTS.json").write_text(json.dumps({"protein_master_rows": len(out), "protein_master_columns": 228, "baseline_files": len(manifest)}, indent=2), encoding="utf-8")
(ROOT/"01_rules"/"V7_SCOPE_NOTICE.md").write_text("# MemPro V7 pure-data freeze\n\nV6.4 and declared dependencies are immutable. All outputs are candidates until final QA. Initial triage flags are not deletion decisions.\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
