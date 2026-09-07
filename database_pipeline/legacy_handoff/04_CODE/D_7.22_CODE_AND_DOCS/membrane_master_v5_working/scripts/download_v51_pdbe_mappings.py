from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
V51 = ROOT / "v51_review_working"
RAW = V51 / "raw" / "pdbe_uniprot_mappings"
NORMALIZED = V51 / "normalized"
QUEUE = ROOT / "review" / "human_membrane_review_queue_v5.tsv"
PDBTM_XML = ROOT / "raw" / "unitmp_pdbtm_all_2026-07-24.xml"
USER_AGENT = "HumanMembraneProteinMaster/5.1 (chain-level membrane review)"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_one(pdb_id: str, attempts: int = 6) -> dict[str, object]:
    url = f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{pdb_id}"
    destination = RAW / f"{pdb_id}.json"
    if destination.exists() and destination.stat().st_size > 0:
        try:
            json.loads(destination.read_text(encoding="utf-8"))
            return {
                "pdb_id": pdb_id,
                "status": "cached",
                "url": url,
                "http_date": "",
                "last_modified": "",
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
            }
        except (json.JSONDecodeError, UnicodeDecodeError):
            pass
    last_error: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                body = response.read()
                destination.write_bytes(body)
                return {
                    "pdb_id": pdb_id,
                    "status": "downloaded",
                    "url": url,
                    "http_date": response.headers.get("Date", ""),
                    "last_modified": response.headers.get("Last-Modified", ""),
                    "bytes": destination.stat().st_size,
                    "sha256": sha256(destination),
                }
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                destination.write_text("{}\n", encoding="utf-8")
                return {
                    "pdb_id": pdb_id,
                    "status": "not_found",
                    "url": url,
                    "http_date": "",
                    "last_modified": "",
                    "bytes": destination.stat().st_size,
                    "sha256": sha256(destination),
                }
            last_error = exc
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
        if attempt + 1 < attempts:
            time.sleep(min(20, 2 ** attempt))
    raise RuntimeError(f"Failed PDBe mapping request for {pdb_id}") from last_error


def parse_pdbtm_chains(relevant_pdb_ids: set[str]) -> dict[str, dict[str, dict[str, str]]]:
    result: dict[str, dict[str, dict[str, str]]] = {}
    for _, element in ET.iterparse(PDBTM_XML, events=("end",)):
        if element.tag.split("}")[-1] != "pdbtm":
            continue
        pdb_id = element.attrib.get("ID", "").strip().lower()
        if pdb_id in relevant_pdb_ids:
            chains: dict[str, dict[str, str]] = {}
            for child in list(element):
                if child.tag.split("}")[-1] != "CHAIN":
                    continue
                chain_id = child.attrib.get("CHAINID", "").strip()
                chains[chain_id] = {
                    "num_tm": child.attrib.get("NUM_TM", "").strip(),
                    "chain_type": child.attrib.get("TYPE", "").strip(),
                }
            result[pdb_id] = chains
        element.clear()
    return result


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    NORMALIZED.mkdir(parents=True, exist_ok=True)

    queue_rows = []
    relevant_pdb_ids: set[str] = set()
    with QUEUE.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            queue_rows.append(row)
            for field in ("opm_pdb_ids_v5", "pdbtm_pdb_ids_v5"):
                relevant_pdb_ids.update(
                    value.strip().lower()
                    for value in row.get(field, "").split(";")
                    if value.strip()
                )

    records = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(fetch_one, pdb_id): pdb_id
            for pdb_id in sorted(relevant_pdb_ids)
        }
        for completed, future in enumerate(as_completed(futures), 1):
            record = future.result()
            records.append(record)
            if completed % 25 == 0 or completed == len(futures):
                print(f"PDBe mappings: {completed}/{len(futures)}")

    pdbtm_chains = parse_pdbtm_chains(relevant_pdb_ids)
    output_rows = []
    for row in queue_rows:
        accession = row["target_uniprot_id"].strip().upper()
        pdb_ids = sorted(
            {
                value.strip().lower()
                for field in ("opm_pdb_ids_v5", "pdbtm_pdb_ids_v5")
                for value in row.get(field, "").split(";")
                if value.strip()
            }
        )
        for pdb_id in pdb_ids:
            mapping_path = RAW / f"{pdb_id}.json"
            payload = json.loads(mapping_path.read_text(encoding="utf-8"))
            uniprot_payload = (
                payload.get(pdb_id, {})
                .get("UniProt", {})
                .get(accession, {})
            )
            mappings = uniprot_payload.get("mappings", [])
            target_chains = sorted(
                {
                    str(mapping.get("chain_id", "")).strip()
                    for mapping in mappings
                    if str(mapping.get("chain_id", "")).strip()
                }
            )
            topology = pdbtm_chains.get(pdb_id, {})
            chain_details = []
            target_tm_counts = []
            target_chain_types = []
            for chain_id in target_chains:
                info = topology.get(chain_id, {})
                num_tm = info.get("num_tm", "")
                chain_type = info.get("chain_type", "")
                chain_details.append(f"{chain_id}:{num_tm or '?'}:{chain_type or '?'}")
                if num_tm:
                    target_tm_counts.append(num_tm)
                if chain_type:
                    target_chain_types.append(chain_type)

            target_has_tm = any(
                value.isdigit() and int(value) > 0 for value in target_tm_counts
            )
            target_all_non_tm = bool(target_chains) and bool(target_tm_counts) and all(
                value.isdigit() and int(value) == 0 for value in target_tm_counts
            )
            output_rows.append(
                {
                    "target_uniprot_id": accession,
                    "pdb_id": pdb_id,
                    "pdbe_target_mapped": "1" if target_chains else "0",
                    "pdbe_target_chain_ids": ";".join(target_chains),
                    "pdbtm_target_chain_details": ";".join(chain_details),
                    "pdbtm_target_tm_counts": ";".join(target_tm_counts),
                    "pdbtm_target_chain_types": ";".join(target_chain_types),
                    "pdbtm_target_has_tm": "1" if target_has_tm else "0",
                    "pdbtm_target_all_non_tm": "1" if target_all_non_tm else "0",
                    "pdbe_mapping_url": (
                        f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{pdb_id}"
                    ),
                    "source_versions": (
                        "PDBe/SIFTS current retrieval; "
                        "UniTmp PDBTM current 2026-07-24"
                    ),
                }
            )

    columns = list(output_rows[0])
    output_path = NORMALIZED / "pdbtm_target_chain_topology_v51.tsv"
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    manifest = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "queue_sha256": sha256(QUEUE),
        "pdbtm_xml_sha256": sha256(PDBTM_XML),
        "requested_pdb_ids": len(relevant_pdb_ids),
        "mapping_records": sorted(records, key=lambda item: item["pdb_id"]),
        "normalized_file": str(output_path),
        "normalized_rows": len(output_rows),
        "normalized_sha256": sha256(output_path),
    }
    manifest_path = V51 / "raw" / "pdbe_mapping_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "pdb_ids": len(relevant_pdb_ids),
                "normalized_rows": len(output_rows),
                "target_chain_tm_hits": sum(
                    row["pdbtm_target_has_tm"] == "1" for row in output_rows
                ),
                "target_chain_non_tm_hits": sum(
                    row["pdbtm_target_all_non_tm"] == "1" for row in output_rows
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
