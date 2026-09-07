from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

BASE = Path(r"D:\finale\16_MemPro_V7.0.2_final_20260814")
OUT = Path(r"D:\finale\17_MemPro_V7.1_candidate_20260814")


def opener(path: Path):
    return gzip.open(path, "rt", encoding="utf-8-sig", newline="") if path.suffix == ".gz" else path.open("r", encoding="utf-8-sig", newline="")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inspect(path: Path) -> dict:
    rel = path.relative_to(BASE).as_posix()
    rec = {"relative_path": rel, "bytes": path.stat().st_size, "sha256": sha256(path)}
    if path.name.endswith((".tsv", ".tsv.gz")):
        with opener(path) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, [])
            rows = sum(1 for _ in reader)
        rec.update({"format": "tsv", "row_count": rows, "column_count": len(header), "columns": header})
    elif path.suffix == ".json":
        rec.update({"format": "json", "row_count": None, "column_count": None, "columns": []})
    else:
        rec.update({"format": path.suffix.lstrip("."), "row_count": None, "column_count": None, "columns": []})
    return rec


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "00_baseline").mkdir(exist_ok=True)
    records = [inspect(p) for p in sorted(BASE.rglob("*")) if p.is_file()]
    with (OUT / "00_baseline" / "V702_BASELINE_INVENTORY.json").open("w", encoding="utf-8") as f:
        json.dump({"baseline": str(BASE), "file_count": len(records), "files": records}, f, ensure_ascii=False, indent=2)
    with (OUT / "00_baseline" / "V702_BASELINE_INVENTORY.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["relative_path", "bytes", "sha256", "format", "row_count", "column_count", "columns"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for rec in records:
            row = dict(rec)
            row["columns"] = "|".join(row["columns"])
            w.writerow(row)
    print(json.dumps({"file_count": len(records), "out": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
