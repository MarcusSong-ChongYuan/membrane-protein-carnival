import csv
import hashlib
import json
from pathlib import Path

root = Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813")
manifest = root / "MANIFEST.tsv"
checked = missing = mismatched = 0
details = []
with manifest.open("r", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        path = root / row["relative_path"]
        checked += 1
        if not path.is_file():
            missing += 1
            details.append({"file": row["relative_path"], "status": "MISSING"})
            continue
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
        actual = digest.hexdigest()
        if actual != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            mismatched += 1
            details.append({"file": row["relative_path"], "status": "MISMATCH"})

report = {
    "status": "PASS" if missing == 0 and mismatched == 0 else "FAIL",
    "manifest_entries_checked": checked,
    "missing": missing,
    "mismatched": mismatched,
    "details": details,
}
(root / "MANIFEST_VERIFICATION.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(report, ensure_ascii=False, indent=2))
