import csv
import gzip
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

root = Path(r"D:\finale\09_V7_data_freeze_working_20260813\02_protein_audit")
freeze = root / "V7_PROTEIN_FREEZE_CANDIDATE"
validation = root / "final_validation_v7_3_1" / "V7_3_FINAL_VALIDATION_SUMMARY.json"
correction = root / "final_adjudication_v7_3_1" / "V7_3_1_BOUNDARY_CORRECTION.json"

for tsv in sorted(freeze.glob("*.tsv")):
    gz = tsv.with_suffix(tsv.suffix + ".gz")
    with tsv.open("rb") as src, gzip.open(gz, "wb", compresslevel=9) as dst:
        shutil.copyfileobj(src, dst)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

manifest_rows = []
for path in sorted(freeze.iterdir()):
    if path.is_file() and path.name not in {"MANIFEST.tsv", "V7_PROTEIN_FREEZE_REPORT.json"}:
        manifest_rows.append({
            "filename": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
with (freeze / "MANIFEST.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["filename", "bytes", "sha256"], delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest_rows)

summary = json.loads(validation.read_text(encoding="utf-8"))
report = {
    "release_name": "MemPro V7 protein-layer freeze candidate",
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "source_release": "V6.4.1 official public candidate plus current UniProt and cross-source T09 adjudication",
    "classification_policy": {
        "A": "integral membrane: transmembrane/intramembrane or strong structure-supported integral evidence",
        "B": "direct membrane insertion by covalent lipid/GPI anchor",
        "C": "experimentally/structurally supported peripheral or complex-mediated membrane association",
        "EXCLUDED": "no validated membrane-binding mechanism, non-membrane, or incomplete immune gene-segment entity",
    },
    "validation": summary,
    "boundary_correction": json.loads(correction.read_text(encoding="utf-8")),
    "files": manifest_rows,
    "release_scope": "protein layer only; relation tables must be rebuilt against the 7,800-accession whitelist before complete V7 freeze",
}
(freeze / "V7_PROTEIN_FREEZE_REPORT.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps({"files": len(manifest_rows), "status": summary["status"], "public_ABC_records": summary["public_ABC_records"]}, ensure_ascii=False))
