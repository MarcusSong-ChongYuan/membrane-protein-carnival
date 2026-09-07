import csv
import json
from pathlib import Path

root = Path(r"D:\finale\09_V7_data_freeze_working_20260813\02_protein_audit")
src = root / "final_adjudication_v7_3" / "protein_identity_final_adjudication_all_v7_3.tsv"
outdir = root / "final_adjudication_v7_3_1"
outdir.mkdir(parents=True, exist_ok=True)
dst = outdir / "protein_identity_final_adjudication_all_v7_3_1.tsv"

with src.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
fields = list(rows[0].keys())
changed = []
for row in rows:
    if row["target_uniprot"] != "P01880":
        continue
    before = {k: row.get(k, "") for k in ["entity_type_v7_2", "form_specificity_v7_2", "final_class_v7_3", "final_mode_v7_3", "final_decision_v7_3", "public_inclusion_v7_3"]}
    row["entity_type_v7_2"] = "canonical_protein"
    row["form_specificity_v7_2"] = "membrane_isoform_A;secreted_isoform_nonmembrane"
    row["final_class_v7_3"] = "A"
    row["final_mode_v7_3"] = "immunoglobulin_heavy_chain_membrane_isoform"
    row["final_decision_v7_3"] = "CONFIRMED_A_FORM_SPECIFIC"
    row["final_evidence_v7_3"] = "Current UniProt records an explicit transmembrane helix (407-427) and both secreted and single-pass type I membrane forms; classification A applies to the membrane form"
    row["public_inclusion_v7_3"] = "1"
    row["final_validation_status_v7_3"] = "CORRECTED_BOUNDARY_CASE_NO_UNRESOLVED_STATUS"
    after = {k: row.get(k, "") for k in before}
    changed.append({"target_uniprot": row["target_uniprot"], "before": before, "after": after})
if len(changed) != 1:
    raise RuntimeError(f"Expected exactly one P01880 row, changed={len(changed)}")
with dst.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
(outdir / "V7_3_1_BOUNDARY_CORRECTION.json").write_text(json.dumps(changed, ensure_ascii=False, indent=2), encoding="utf-8")
print(dst)
