from __future__ import annotations

import csv
import http.client
import json
import socket
import time
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(r"D:\7.22")
INDEX = (
    ROOT
    / "membrane_master_v5_working"
    / "releases"
    / "release_v5_4_pgd_binding"
    / "small_molecule_index_v4_1.tsv"
)
WORK = ROOT / "small_molecule_v1_working"
RAW = WORK / "raw" / "pubchem"
REPORTS = WORK / "reports"
RAW.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

OUT_JSONL = RAW / "pubchem_missing_properties_v1.jsonl"
FAILED_TSV = RAW / "pubchem_missing_failed_v1.tsv"
REPORT = REPORTS / "PUBCHEM_MISSING_RETRIEVAL_REPORT.json"

PROPERTIES = ",".join(
    [
        "Title",
        "ConnectivitySMILES",
        "SMILES",
        "InChI",
        "InChIKey",
        "MolecularFormula",
        "MolecularWeight",
        "ExactMass",
        "MonoisotopicMass",
        "XLogP",
        "TPSA",
        "HBondDonorCount",
        "HBondAcceptorCount",
        "RotatableBondCount",
        "Complexity",
        "Charge",
    ]
)


def load_missing_cids() -> list[str]:
    missing = set()
    with INDEX.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("compound_id_type") != "PubChem CID":
                continue
            if row.get("canonical_smiles", "").strip():
                continue
            cid = row.get("compound_id", "").strip().replace("PUBCHEM:", "")
            if cid.isdigit():
                missing.add(cid)
    return sorted(missing, key=int)


def request_batch(cids: list[str]) -> list[dict]:
    joined = ",".join(cids)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{joined}/property/{PROPERTIES}/JSON"
    )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MemProDB compound curation/1.0"},
    )
    last_error = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                payload = json.load(response)
            time.sleep(0.26)
            return payload.get("PropertyTable", {}).get("Properties", [])
        except (
            urllib.error.HTTPError,
            urllib.error.URLError,
            http.client.RemoteDisconnected,
            socket.timeout,
            TimeoutError,
            OSError,
        ) as exc:
            last_error = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in {400, 404}:
                break
            time.sleep(2**attempt)
    raise RuntimeError(str(last_error))


def fetch_recursive(cids: list[str], output: dict[str, dict], failed: dict[str, str]) -> None:
    if not cids:
        return
    try:
        rows = request_batch(cids)
        returned = set()
        for row in rows:
            cid = str(row.get("CID", ""))
            if cid:
                output[cid] = row
                returned.add(cid)
        for cid in set(cids) - returned:
            failed[cid] = "not_returned"
    except RuntimeError as exc:
        if len(cids) == 1:
            failed[cids[0]] = str(exc)
            return
        midpoint = len(cids) // 2
        fetch_recursive(cids[:midpoint], output, failed)
        fetch_recursive(cids[midpoint:], output, failed)


def save_output(output: dict[str, dict]) -> None:
    with OUT_JSONL.open("w", encoding="utf-8", newline="\n") as handle:
        for cid in sorted(output, key=int):
            handle.write(json.dumps(output[cid], ensure_ascii=False) + "\n")


def main() -> None:
    missing = load_missing_cids()
    output: dict[str, dict] = {}
    failed: dict[str, str] = {}

    if OUT_JSONL.exists():
        with OUT_JSONL.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    row = json.loads(line)
                    output[str(row["CID"])] = row

    remaining = [cid for cid in missing if cid not in output]
    for start in range(0, len(remaining), 100):
        batch = remaining[start : start + 100]
        fetch_recursive(batch, output, failed)
        if start % 500 == 0:
            save_output(output)
        if start % 1000 == 0:
            print(
                json.dumps(
                    {
                        "processed": min(start + len(batch), len(remaining)),
                        "remaining_initial": len(remaining),
                        "retrieved_total": len(output),
                        "failed_total": len(failed),
                    }
                ),
                flush=True,
            )

    save_output(output)

    with FAILED_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pubchem_cid", "failure_reason"])
        for cid in sorted(failed, key=int):
            writer.writerow([cid, failed[cid]])

    report = {
        "requested_missing_pubchem_cids": len(missing),
        "retrieved_records": len(set(missing) & set(output)),
        "failed_records": len(failed),
        "output_jsonl": str(OUT_JSONL),
        "failed_tsv": str(FAILED_TSV),
        "source": "PubChem PUG REST",
        "source_url": "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest",
        "retrieval_date": "2026-07-26",
    }
    REPORT.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
