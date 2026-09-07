from pathlib import Path

source = Path(r"C:\Users\Administrator\final_validate_and_freeze_v7_protein.py").read_text(encoding="utf-8")
source = source.replace(
    'SRC = ROOT / "final_adjudication_v7_3" / "protein_identity_final_adjudication_all_v7_3.tsv"',
    'SRC = ROOT / "final_adjudication_v7_3_1" / "protein_identity_final_adjudication_all_v7_3_1.tsv"',
)
source = source.replace(
    'OUT = ROOT / "final_validation_v7_3"',
    'OUT = ROOT / "final_validation_v7_3_1"',
)
source = source.replace(
    'if row["final_validation_status_v7_3"] != "T09_APPLIED_NO_UNRESOLVED_STATUS":',
    'if row["final_validation_status_v7_3"] not in {"T09_APPLIED_NO_UNRESOLVED_STATUS", "CORRECTED_BOUNDARY_CASE_NO_UNRESOLVED_STATUS"}:',
)
exec(compile(source, "final_validate_and_freeze_v7_protein_v7_3_1.py", "exec"))
