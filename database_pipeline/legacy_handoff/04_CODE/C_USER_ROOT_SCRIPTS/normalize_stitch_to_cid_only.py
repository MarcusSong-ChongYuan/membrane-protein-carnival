import csv
import hashlib
import json
import re
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")
BINARY = BASE / "drug_protein_binary_relationships.tsv"
PROTEIN = BASE / "protein_database.tsv"
MANIFEST = BASE / "MANIFEST.json"
REPORT = BASE / "stitch_cid_only_normalization_report.json"

csv.field_size_limit(100_000_000)

CHEMBL_CID_SCORE = re.compile(r"\bCHEMBL\d+\(CID(\d+)\)(:\d+(?:\.\d+)?)")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def update_manifest(paths):
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    files = manifest.get("files", [])
    by_name = {Path(item["path"]).name: item for item in files}
    for path in paths:
        item = by_name.get(path.name)
        if item is None:
            item = {"path": f"upload/normalized_tables/{path.name}"}
            files.append(item)
        item["sha256"] = sha256(path)
        item["bytes"] = path.stat().st_size
    manifest["files"] = files
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def split_preserve_separator(value):
    if "|" in value:
        return "|", [part.strip() for part in value.split("|") if part.strip()]
    return ";", [part.strip() for part in value.split(";") if part.strip()]


def normalize_value(value):
    if not value:
        return value, 0, 0
    separator, parts = split_preserve_separator(value)
    output = []
    seen = set()
    replacements = 0
    duplicates_removed = 0
    for part in parts:
        new_part, n = CHEMBL_CID_SCORE.subn(r"CID\1\2", part)
        replacements += n
        if new_part in seen:
            duplicates_removed += 1
            continue
        seen.add(new_part)
        output.append(new_part)
    return separator.join(output), replacements, duplicates_removed


def normalize_file(path, columns):
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)

    changed_rows = 0
    replacements = 0
    duplicates_removed = 0
    for row in rows:
        row_changed = False
        for column in columns:
            before = row.get(column, "")
            after, n_replacements, n_duplicates = normalize_value(before)
            if after != before:
                row[column] = after
                row_changed = True
            replacements += n_replacements
            duplicates_removed += n_duplicates
        if row_changed:
            changed_rows += 1

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    return {
        "path": str(path),
        "rows": len(rows),
        "changed_rows": changed_rows,
        "replacements": replacements,
        "duplicates_removed": duplicates_removed,
    }


def main():
    results = [
        normalize_file(BINARY, ["stitch_matched_compounds"]),
        normalize_file(PROTEIN, ["stitch_compounds_cleaned"]),
    ]
    REPORT.write_text(json.dumps({"files": results}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    update_manifest([BINARY, PROTEIN, REPORT])
    print(json.dumps({"files": results, "report": str(REPORT)}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
