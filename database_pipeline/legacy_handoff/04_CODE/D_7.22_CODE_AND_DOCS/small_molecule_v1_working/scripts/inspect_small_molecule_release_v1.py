from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path


DB = Path(r"D:\7.22\small_molecule_v1_working\intermediate\small_molecule_v1.sqlite")
EVIDENCE = Path(
    r"D:\7.22\small_molecule_v1_working\releases\release_small_molecule_v1_0"
    r"\binding_evidence_master_v4_2.tsv"
)


def main() -> None:
    database = sqlite3.connect(DB)
    result = {
        "scope_reasons": database.execute(
            """
            SELECT final_scope_status, final_scope_reason, COUNT(*)
            FROM source_curated
            GROUP BY final_scope_status, final_scope_reason
            ORDER BY COUNT(*) DESC
            LIMIT 30
            """
        ).fetchall(),
        "excluded_molecule_types": database.execute(
            """
            SELECT chembl_molecule_type, COUNT(*)
            FROM source_curated
            WHERE final_scope_status='excluded'
            GROUP BY chembl_molecule_type
            ORDER BY COUNT(*) DESC
            """
        ).fetchall(),
        "unresolved_source_types": database.execute(
            """
            SELECT source_compound_id_type, COUNT(*)
            FROM source_curated
            WHERE final_scope_status='unresolved'
            GROUP BY source_compound_id_type
            ORDER BY COUNT(*) DESC
            """
        ).fetchall(),
    }
    missing = []
    with EVIDENCE.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["compound_mapping_status"] == "not_in_v4_1_compound_index":
                missing.append(
                    {
                        "compound_id": row["compound_id"],
                        "compound_id_type": row["compound_id_type"],
                        "compound_name": row["compound_name"],
                        "source_database": row["source_database"],
                        "evidence_id": row["evidence_id"],
                    }
                )
    result["binding_rows_not_in_v4_1_index"] = missing
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
