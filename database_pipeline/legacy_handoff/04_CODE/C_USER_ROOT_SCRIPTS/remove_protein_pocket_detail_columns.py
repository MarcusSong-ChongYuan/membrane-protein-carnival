import csv
import hashlib
import json
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")
PROTEIN = BASE / "protein_database.tsv"
MANIFEST = BASE / "MANIFEST.json"

REMOVE_COLUMNS = {
    "pocket_ids_alphafold",
    "pocket_centers_alphafold",
    "pocket_scores_alphafold",
    "pocket_residues_all_alphafold",
}

csv.field_size_limit(100_000_000)


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


def main():
    with PROTEIN.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        original_fields = list(reader.fieldnames or [])
        rows = list(reader)

    missing = sorted(REMOVE_COLUMNS - set(original_fields))
    if missing:
        raise SystemExit(f"Columns not found in {PROTEIN}: {missing}")

    output_fields = [field for field in original_fields if field not in REMOVE_COLUMNS]
    with PROTEIN.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=output_fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in output_fields})

    update_manifest([PROTEIN])
    print(json.dumps({
        "protein_rows": len(rows),
        "columns_before": len(original_fields),
        "columns_after": len(output_fields),
        "removed_columns": sorted(REMOVE_COLUMNS),
        "kept_summary_column": "n_pockets_alphafold",
        "output": str(PROTEIN),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
