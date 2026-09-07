from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path

from standardize_compound_structures_v1 import OUTPUT_FIELDS, standardize_row


ROOT = Path(r"D:\7.22")
WORK = ROOT / "small_molecule_v1_working"
CHUNKS = WORK / "intermediate" / "standardized_chunks"
RAW = WORK / "raw" / "pubchem"
REPORT = WORK / "reports" / "BINDING_COMPOUND_INDEX_GAP_REPAIR_REPORT.json"
OUTPUT = CHUNKS / "standardized_chunk_0061_supplement.tsv"

CIDS = ["17107496", "131801155", "135450187", "139033340", "148568163"]
PROPERTIES = ",".join(
    [
        "Title",
        "ConnectivitySMILES",
        "SMILES",
        "InChI",
        "InChIKey",
        "MolecularFormula",
        "MolecularWeight",
    ]
)


def main() -> None:
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{','.join(CIDS)}/property/{PROPERTIES}/JSON"
    )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MemProDB compound curation/1.0"},
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        properties = json.load(response)["PropertyTable"]["Properties"]
    pubchem = {str(row["CID"]): row for row in properties}

    records = []
    for offset, cid in enumerate(CIDS, start=302_880):
        property_row = pubchem[cid]
        source_row = {
            "compound_id": cid,
            "compound_id_type": "PubChem CID",
            "compound_name": property_row.get("Title", ""),
            "canonical_smiles": property_row.get("SMILES")
            or property_row.get("ConnectivitySMILES")
            or "",
            "inchikey": property_row.get("InChIKey", ""),
            "pubchem_cids": cid,
            "chembl_ids": "",
            "small_molecule_scope_status": "accepted_core",
            "small_molecule_scope_class": "research_ligand",
            "is_core_small_molecule": "1",
            "source_databases": "IUPHAR/BPS Guide to PHARMACOLOGY",
            "source_versions": "GtoPdb 2026.2",
            "record_qc_status": "repaired_missing_v4_1_index_record",
        }
        records.append(standardize_row(source_row, offset, {}, pubchem))

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(records)

    RAW.mkdir(parents=True, exist_ok=True)
    raw_out = RAW / "pubchem_binding_index_gap_repair_v1.json"
    raw_out.write_text(
        json.dumps(properties, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    report = {
        "missing_distinct_pubchem_cids": len(CIDS),
        "retrieved_records": len(records),
        "repaired_evidence_rows": 8,
        "pubchem_cids": CIDS,
        "source": "PubChem PUG REST",
        "source_url": "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest",
        "retrieval_date": "2026-07-27",
        "output_chunk": str(OUTPUT),
    }
    REPORT.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
