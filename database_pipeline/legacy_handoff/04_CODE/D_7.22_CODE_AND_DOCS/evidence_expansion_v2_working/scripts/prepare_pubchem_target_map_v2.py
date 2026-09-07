#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
PROTEINS = ROOT / "intermediate" / "protein_universe_abc_v2.tsv"
IDENTIFIERS = ROOT / "raw" / "pubchem_202607" / "protein-identifiers.tsv.gz"
OUTPUT = ROOT / "intermediate" / "pubchem_protein_accession_map_v2.tsv"
REPORT = ROOT / "reports" / "PUBCHEM_TARGET_MAPPING_REPORT.json"


def main() -> None:
    proteins = {}
    with PROTEINS.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            proteins[row["target_uniprot_id"]] = row

    accession_map: dict[str, set[str]] = defaultdict(set)
    with gzip.open(IDENTIFIERS, "rt", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            protein_accession, external_id, id_type = parts
            if id_type == "UniProt ID" and external_id in proteins:
                accession_map[external_id].add(protein_accession)

    counts = defaultdict(int)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "target_uniprot_id",
            "approved_symbol",
            "ncbi_gene_ids",
            "pubchem_protein_accessions",
            "pubchem_protein_accession_count",
            "pubchem_mapping_status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for uniprot, row in proteins.items():
            accessions = sorted(accession_map.get(uniprot, set()))
            if len(accessions) == 1:
                status = "single_pubchem_protein_accession"
            elif accessions:
                status = "multiple_pubchem_protein_accessions"
            else:
                status = "no_pubchem_protein_accession"
            counts[status] += 1
            writer.writerow(
                {
                    "target_uniprot_id": uniprot,
                    "approved_symbol": row.get("approved_symbol", ""),
                    "ncbi_gene_ids": row.get("ncbi_gene_ids", ""),
                    "pubchem_protein_accessions": ";".join(accessions),
                    "pubchem_protein_accession_count": len(accessions),
                    "pubchem_mapping_status": status,
                }
            )

    report = {
        "status": "passed",
        "protein_universe_rows": len(proteins),
        "mapping_counts": dict(counts),
        "source_file": str(IDENTIFIERS),
        "output_file": str(OUTPUT),
    }
    with REPORT.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
