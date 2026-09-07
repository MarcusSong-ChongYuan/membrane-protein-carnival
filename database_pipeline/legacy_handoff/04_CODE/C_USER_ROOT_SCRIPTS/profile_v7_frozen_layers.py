import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813")
FILES = {
    "sites": ROOT / "02_frozen_review" / "binding_site_coordinate_frozen_v7.tsv.gz",
    "isoforms": ROOT / "02_frozen_review" / "isoform_frozen_review_v7.tsv.gz",
    "processed": ROOT / "02_frozen_review" / "processed_form_frozen_review_v7.tsv.gz",
}

def op(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="")

def summarize(path, name):
    with op(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = []
        nonempty = Counter()
        values = defaultdict(Counter)
        likely_status = [f for f in fields if any(k in f.lower() for k in ("status", "reason", "action", "resolution", "mapping", "disposition", "review"))]
        total = 0
        for row in reader:
            total += 1
            if len(rows) < 8:
                rows.append(row)
            for f, v in row.items():
                if v:
                    nonempty[f] += 1
            for f in likely_status:
                values[f][row.get(f, "") or "<EMPTY>"] += 1
        return {
            "name": name,
            "path": str(path),
            "rows": total,
            "fields": fields,
            "status_distributions": {f: dict(c.most_common(30)) for f, c in values.items()},
            "nonempty_counts": dict(nonempty),
            "examples": rows,
        }

out = {name: summarize(path, name) for name, path in FILES.items()}
dest = Path(r"D:\finale\12_V7_final_completion_20260813\07_qa\FROZEN_LAYER_PROFILE.json")
dest.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
for name, data in out.items():
    print("\n###", name, data["rows"])
    print("FIELDS:", data["fields"])
    for field, dist in data["status_distributions"].items():
        print(field, dist)
    print("EXAMPLES:")
    for row in data["examples"][:3]:
        print({k: v for k, v in row.items() if v})
print("\nREPORT", dest)
