from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")


def main() -> None:
    rels: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "drugs": set(),
            "sources": set(),
            "has_matched": False,
        }
    )
    with (BASE / "drug_protein_binary_relationships.tsv").open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            uid = row["target_uniprot_id"]
            item = rels[uid]
            item["drugs"].add(row["drug_id"])
            item["sources"].update(x.strip() for x in row["source_database"].split(";") if x.strip())
            if any(row.get(col) == "1" for col in [
                "has_pdb_biolip_matched",
                "has_pdb_scpdb_matched",
                "has_pdbbind_matched",
                "has_stitch_matched",
                "has_exp_binding_site",
            ]):
                item["has_matched"] = True

    path = BASE / "protein_database.tsv"
    temp = path.with_suffix(".tmp")
    with path.open("r", encoding="utf-8", newline="") as src, temp.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fields = reader.fieldnames or []
        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            uid = row["target_uniprot_id"]
            item = rels.get(uid)
            if item:
                drugs = sorted(item["drugs"])
                row["unique_drug_count"] = str(len(drugs))
                row["unique_drug_cid_sample"] = ";".join(drugs[:50])
                row["source_databases"] = ";".join(sorted(item["sources"]))
                row["has_matched_evidence"] = "1" if item["has_matched"] else "0"
            else:
                row["unique_drug_count"] = "0"
                row["unique_drug_cid_sample"] = ""
                row["source_databases"] = ""
                row["has_matched_evidence"] = "0"
            writer.writerow(row)
    temp.replace(path)


if __name__ == "__main__":
    main()
