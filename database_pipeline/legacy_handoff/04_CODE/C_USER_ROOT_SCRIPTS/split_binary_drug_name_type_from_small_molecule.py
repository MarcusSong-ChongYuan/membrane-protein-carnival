from __future__ import annotations

import csv
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")


def main() -> None:
    small_path = BASE / "small_molecule_database.tsv"
    role_by_drug: dict[str, tuple[str, str]] = {}
    with small_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            role_by_drug[row["drug_id"]] = (
                row["compound_name_category"],
                row["compound_biological_role"],
            )

    binary_path = BASE / "drug_protein_binary_relationships.tsv"
    temp = binary_path.with_suffix(".tmp")
    with binary_path.open("r", encoding="utf-8", newline="") as src, temp.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        old_fields = reader.fieldnames or []
        fields: list[str] = []
        for field in old_fields:
            if field == "drug_name_type":
                fields.extend(["compound_name_category", "compound_biological_role"])
            else:
                fields.append(field)
        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        missing = 0
        for row in reader:
            category, role = role_by_drug.get(row["drug_id"], ("unknown_name_type", "bioactive_research_ligand"))
            if row["drug_id"] not in role_by_drug:
                missing += 1
            out: dict[str, str] = {}
            for field in old_fields:
                if field == "drug_name_type":
                    out["compound_name_category"] = category
                    out["compound_biological_role"] = role
                else:
                    out[field] = row.get(field, "")
            writer.writerow(out)
    temp.replace(binary_path)
    print({"missing_drug_ids": missing, "new_columns": len(fields)})


if __name__ == "__main__":
    main()
