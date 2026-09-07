#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sqlite3
from collections import Counter
from pathlib import Path

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
BASE = Path(r"D:\7.22\small_molecule_v1_working\releases\release_small_molecule_v1_0")
PROTEINS = ROOT / "intermediate" / "protein_universe_abc_v2.tsv"
XREFS = BASE / "compound_source_xref_v1_0.tsv"
PAIRS = BASE / "protein_compound_summary_v4_2.tsv"
DB = ROOT / "intermediate" / "integration_index_v2.sqlite"


def split_ids(value: str) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.replace("|", ";").split(";") if item.strip()]


def add_alias(cursor: sqlite3.Cursor, alias: str, uniprot: str, alias_type: str) -> None:
    if alias:
        cursor.execute(
            "INSERT OR IGNORE INTO protein_alias(alias, target_uniprot_id, alias_type) VALUES (?, ?, ?)",
            (alias, uniprot, alias_type),
        )


def main() -> None:
    for path in (PROTEINS, XREFS, PAIRS):
        if not path.exists():
            raise FileNotFoundError(path)
    if DB.exists():
        DB.unlink()

    connection = sqlite3.connect(DB)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=MEMORY")
    cursor = connection.cursor()
    cursor.executescript(
        """
        CREATE TABLE protein (
          target_uniprot_id TEXT PRIMARY KEY,
          approved_symbol TEXT,
          protein_name TEXT,
          primary_gene_id TEXT,
          membrane_class TEXT,
          evidence_level TEXT,
          website_default INTEGER,
          functional_primary_class TEXT,
          functional_subclass TEXT,
          membrane_topology TEXT
        );
        CREATE TABLE protein_alias (
          alias TEXT NOT NULL,
          target_uniprot_id TEXT NOT NULL,
          alias_type TEXT NOT NULL,
          PRIMARY KEY(alias, target_uniprot_id, alias_type)
        );
        CREATE INDEX protein_alias_lookup ON protein_alias(alias);

        CREATE TABLE compound_xref (
          id_type TEXT NOT NULL,
          id_value TEXT NOT NULL,
          compound_internal_id TEXT,
          compound_form_id TEXT,
          mapping_status TEXT,
          scope_status TEXT,
          source_record_id TEXT,
          preferred_source_name TEXT,
          PRIMARY KEY(id_type, id_value, compound_internal_id, compound_form_id, source_record_id)
        );
        CREATE INDEX compound_xref_lookup ON compound_xref(id_type, id_value);

        CREATE TABLE existing_pair (
          target_uniprot_id TEXT NOT NULL,
          compound_internal_id TEXT NOT NULL,
          best_evidence_tier TEXT,
          evidence_count INTEGER,
          independent_sources TEXT,
          PRIMARY KEY(target_uniprot_id, compound_internal_id)
        );
        CREATE INDEX existing_pair_compound ON existing_pair(compound_internal_id);
        """
    )

    metrics = Counter()
    with PROTEINS.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            uniprot = row["target_uniprot_id"].strip()
            cursor.execute(
                """
                INSERT INTO protein VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    uniprot,
                    row.get("approved_symbol", ""),
                    row.get("protein_name", ""),
                    row.get("primary_gene_id", ""),
                    row.get("membrane_class_v52", ""),
                    row.get("evidence_level_v52", ""),
                    int(row.get("website_default_v52") or 0),
                    row.get("functional_primary_class_v5", ""),
                    row.get("functional_subclass_v5", ""),
                    row.get("membrane_topology_v5", ""),
                ),
            )
            add_alias(cursor, uniprot, uniprot, "UniProt")
            add_alias(cursor, row.get("approved_symbol", "").strip(), uniprot, "gene_symbol")
            add_alias(cursor, row.get("primary_gene_id", "").strip(), uniprot, "primary_gene_id")
            for gene_id in split_ids(row.get("ncbi_gene_ids", "")):
                add_alias(cursor, gene_id, uniprot, "NCBI_Gene")
            for ensembl_id in split_ids(row.get("ensembl_gene_ids", "")):
                add_alias(cursor, ensembl_id.split(".")[0], uniprot, "Ensembl")
            metrics["protein_rows"] += 1
    connection.commit()

    xref_batch = []
    with XREFS.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            base_values = (
                row.get("compound_internal_id", ""),
                row.get("compound_form_id", ""),
                row.get("mapping_status", ""),
                row.get("final_scope_status", ""),
                row.get("source_record_id", ""),
                row.get("source_compound_name", ""),
            )
            native_type = row.get("source_compound_id_type", "").strip()
            native_id = row.get("source_compound_id", "").strip()
            if native_type and native_id:
                xref_batch.append((native_type, native_id, *base_values))
            for cid in split_ids(row.get("pubchem_cids", "")):
                xref_batch.append(("PubChem CID", cid, *base_values))
            for chembl_id in split_ids(row.get("chembl_ids", "")):
                xref_batch.append(("ChEMBL ID", chembl_id, *base_values))
            if len(xref_batch) >= 50000:
                cursor.executemany(
                    "INSERT OR IGNORE INTO compound_xref VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    xref_batch,
                )
                metrics["compound_xref_attempted"] += len(xref_batch)
                xref_batch.clear()
    if xref_batch:
        cursor.executemany(
            "INSERT OR IGNORE INTO compound_xref VALUES (?, ?, ?, ?, ?, ?, ?, ?)", xref_batch
        )
        metrics["compound_xref_attempted"] += len(xref_batch)
    connection.commit()

    pair_batch = []
    with PAIRS.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pair_batch.append(
                (
                    row["target_uniprot_id"],
                    row["compound_internal_id"],
                    row.get("best_binding_evidence_level", ""),
                    int(row.get("binding_evidence_count") or 0),
                    row.get("independent_sources", ""),
                )
            )
            if len(pair_batch) >= 50000:
                cursor.executemany("INSERT INTO existing_pair VALUES (?, ?, ?, ?, ?)", pair_batch)
                metrics["existing_pair_rows"] += len(pair_batch)
                pair_batch.clear()
    if pair_batch:
        cursor.executemany("INSERT INTO existing_pair VALUES (?, ?, ?, ?, ?)", pair_batch)
        metrics["existing_pair_rows"] += len(pair_batch)
    connection.commit()

    metrics["protein_alias_rows"] = cursor.execute(
        "SELECT COUNT(*) FROM protein_alias"
    ).fetchone()[0]
    metrics["compound_xref_rows"] = cursor.execute(
        "SELECT COUNT(*) FROM compound_xref"
    ).fetchone()[0]
    metrics["pubchem_cid_xref_rows"] = cursor.execute(
        "SELECT COUNT(*) FROM compound_xref WHERE id_type='PubChem CID'"
    ).fetchone()[0]
    metrics["chembl_id_xref_rows"] = cursor.execute(
        "SELECT COUNT(*) FROM compound_xref WHERE id_type='ChEMBL ID'"
    ).fetchone()[0]
    metrics["database_size_bytes"] = DB.stat().st_size

    integrity = cursor.execute("PRAGMA integrity_check").fetchone()[0]
    report = {"status": "passed" if integrity == "ok" else "failed", "integrity": integrity}
    report.update(metrics)
    with (ROOT / "reports" / "INTEGRATION_INDEX_V2_REPORT.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(report, handle, indent=2)
    connection.close()
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
