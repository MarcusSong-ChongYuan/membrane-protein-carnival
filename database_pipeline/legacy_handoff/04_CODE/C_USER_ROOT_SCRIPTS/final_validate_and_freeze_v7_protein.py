import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(r"D:\finale\09_V7_data_freeze_working_20260813\02_protein_audit")
SRC = ROOT / "final_adjudication_v7_3" / "protein_identity_final_adjudication_all_v7_3.tsv"
OUT = ROOT / "final_validation_v7_3"
OUT.mkdir(parents=True, exist_ok=True)


def intval(value):
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def truthy(value):
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def stable_rank(accession, salt):
    return hashlib.sha256(f"{salt}|{accession}".encode()).hexdigest()


with SRC.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

blocking = []
warnings = []
class_counts = Counter()
decision_counts = Counter()
mode_counts = Counter()
entity_counts = Counter()
seen = Counter()

for row in rows:
    acc = row["target_uniprot"]
    cls = row["final_class_v7_3"]
    decision = row["final_decision_v7_3"]
    mode = row["final_mode_v7_3"]
    included = truthy(row["public_inclusion_v7_3"])
    entity = row["entity_type_v7_2"]
    seen[acc] += 1
    class_counts[cls] += 1
    decision_counts[decision] += 1
    mode_counts[mode] += 1
    entity_counts[entity] += 1

    tm = intval(row["uniprot_transmembrane_count_current"])
    intra = intval(row["uniprot_intramembrane_count_current"])
    lipid = intval(row["uniprot_lipidation_count_current"])
    loc = row["uniprot_subcellular_locations_current"].lower()
    external_integral = bool(row["external_integral_supports_v7"].strip())
    opm_integral = truthy(row["opm_integral"])
    pdbtm = truthy(row["pdbtm_present"])
    has_anchor_text = "lipid-anchor" in loc or "gpi-anchor" in loc
    has_integral_mechanism = tm > 0 or intra > 0 or external_integral or opm_integral or pdbtm
    has_anchor_mechanism = lipid > 0 or has_anchor_text

    if included != (cls in {"A", "B", "C"}):
        blocking.append((acc, "INCLUSION_CLASS_MISMATCH", cls, decision))
    if cls == "A" and not has_integral_mechanism:
        blocking.append((acc, "A_WITHOUT_INTEGRAL_MECHANISM", cls, decision))
    if cls == "B" and not has_anchor_mechanism:
        blocking.append((acc, "B_WITHOUT_ANCHOR_MECHANISM", cls, decision))
    if cls == "C" and mode in {"", "location_without_binding_mode", "soluble_secreted_or_lumenal"}:
        blocking.append((acc, "C_WITHOUT_PERIPHERAL_MODE", cls, decision))
    if entity == "immune_gene_segment" and included:
        blocking.append((acc, "IMMUNE_GENE_SEGMENT_IN_PUBLIC_PROTEIN_LAYER", cls, decision))
    if cls == "EXCLUDED" and decision.startswith("CONFIRMED_"):
        blocking.append((acc, "EXCLUDED_WITH_CONFIRMED_MEMBRANE_DECISION", cls, decision))
    if cls in {"A", "B", "C"} and decision.startswith("EXCLUDED"):
        blocking.append((acc, "INCLUDED_WITH_EXCLUSION_DECISION", cls, decision))
    if row["final_validation_status_v7_3"] != "T09_APPLIED_NO_UNRESOLVED_STATUS":
        blocking.append((acc, "NONFINAL_VALIDATION_STATUS", cls, decision))
    if cls == "EXCLUDED" and (tm > 0 or intra > 0):
        warnings.append((acc, "EXCLUDED_DESPITE_CURRENT_UNIPROT_TM_OR_INTRAMEMBRANE", cls, decision))
    if cls == "EXCLUDED" and has_anchor_mechanism:
        warnings.append((acc, "EXCLUDED_DESPITE_CURRENT_ANCHOR_ANNOTATION", cls, decision))

duplicates = sorted(acc for acc, count in seen.items() if count != 1)
for acc in duplicates:
    blocking.append((acc, "NONUNIQUE_UNIPROT_PRIMARY_KEY", "", str(seen[acc])))

# Deterministic stratified pressure-test sample. This is an audit artifact, not a
# new classifier: it includes every decision stratum and oversamples rare modes.
strata = defaultdict(list)
for row in rows:
    key = (row["final_class_v7_3"], row["final_decision_v7_3"], row["final_mode_v7_3"])
    strata[key].append(row)

sample = []
for key, group in sorted(strata.items()):
    group = sorted(group, key=lambda r: stable_rank(r["target_uniprot"], "V7.3-final-pressure-test"))
    n = min(len(group), 12 if len(group) >= 12 else len(group))
    for row in group[:n]:
        item = dict(row)
        item["validation_stratum"] = " | ".join(key)
        item["sample_reason"] = "deterministic_decision_mode_stratified_sample"
        sample.append(item)

sample_fields = ["validation_stratum", "sample_reason"] + list(rows[0].keys())
with (OUT / "V7_3_FINAL_STRATIFIED_PRESSURE_TEST.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=sample_fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(sample)

issue_fields = ["target_uniprot", "issue", "final_class", "final_decision"]
for name, issues in [("BLOCKING", blocking), ("WARNINGS", warnings)]:
    with (OUT / f"V7_3_FINAL_VALIDATION_{name}.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(issue_fields)
        writer.writerows(issues)

# Freeze protein-layer outputs only if blocking validation passes.
freeze_dir = ROOT / "V7_PROTEIN_FREEZE_CANDIDATE"
if not blocking:
    freeze_dir.mkdir(parents=True, exist_ok=True)
    public_rows = [r for r in rows if truthy(r["public_inclusion_v7_3"])]
    excluded_rows = [r for r in rows if not truthy(r["public_inclusion_v7_3"])]
    fields = list(rows[0].keys())
    for filename, subset in [
        ("human_membrane_protein_master_V7.tsv", public_rows),
        ("protein_excluded_or_other_entity_audit_V7.tsv", excluded_rows),
        ("protein_identity_adjudication_all_V7.tsv", rows),
    ]:
        path = freeze_dir / filename
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(subset)

summary = {
    "input_records": len(rows),
    "unique_uniprot_accessions": len(seen),
    "class_counts": dict(class_counts),
    "public_ABC_records": sum(class_counts[c] for c in ("A", "B", "C")),
    "decision_counts": dict(decision_counts),
    "mode_counts": dict(mode_counts),
    "entity_counts": dict(entity_counts),
    "blocking_error_count": len(blocking),
    "warning_count": len(warnings),
    "stratified_sample_records": len(sample),
    "strata_count": len(strata),
    "freeze_candidate_created": not blocking,
    "status": "PASS_PROTEIN_LAYER_FREEZE_CANDIDATE" if not blocking else "FAIL_BLOCKING_ERRORS",
}
(OUT / "V7_3_FINAL_VALIDATION_SUMMARY.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(summary, ensure_ascii=False, indent=2))
