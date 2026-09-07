from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
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
RAW = WORK / "raw" / "chembl37"
REPORTS = WORK / "reports"
RAW.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

OUT_JSONL = RAW / "chembl37_selected_metadata.jsonl"
FAILED_TSV = RAW / "chembl37_metadata_failed.tsv"
REPORT = REPORTS / "CHEMBL37_METADATA_RETRIEVAL_REPORT.json"
ONLY_FIELDS = ",".join(
    [
        "molecule_chembl_id",
        "pref_name",
        "molecule_type",
        "max_phase",
        "first_approval",
        "chirality",
        "inorganic_flag",
        "natural_product",
        "therapeutic_flag",
        "prodrug",
        "dosed_ingredient",
        "oral",
        "parenteral",
        "topical",
        "black_box_warning",
        "chemical_probe",
        "helm_notation",
        "molecule_hierarchy",
        "molecule_properties",
        "atc_classifications",
    ]
)


def extract_chembl_ids(value: str) -> set[str]:
    ids = set()
    for token in (value or "").replace("|", ";").replace(",", ";").split(";"):
        token = token.strip().upper()
        if token.startswith("CHEMBL") and token[6:].isdigit():
            ids.add(token)
    return ids


def load_wanted() -> list[str]:
    wanted = set()
    with INDEX.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            wanted.update(extract_chembl_ids(row.get("compound_id", "")))
            wanted.update(extract_chembl_ids(row.get("chembl_ids", "")))
    if OUT_JSONL.exists():
        with OUT_JSONL.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                parent_id = str(json.loads(line).get("parent_chembl_id") or "").upper()
                if parent_id.startswith("CHEMBL") and parent_id[6:].isdigit():
                    wanted.add(parent_id)
    return sorted(wanted, key=lambda item: int(item[6:]))


def compact_molecule(row: dict) -> dict:
    structures = row.get("molecule_structures") or {}
    properties = row.get("molecule_properties") or {}
    hierarchy = row.get("molecule_hierarchy") or {}
    synonyms = []
    for item in row.get("molecule_synonyms") or []:
        synonym = str(item.get("molecule_synonym") or "").strip()
        if synonym and synonym.casefold() not in {x.casefold() for x in synonyms}:
            synonyms.append(synonym)
        if len(synonyms) >= 50:
            break
    cross_references = []
    for item in row.get("cross_references") or []:
        source = str(item.get("xref_src") or "").strip()
        value = str(item.get("xref_id") or "").strip()
        if source and value:
            cross_references.append(f"{source}:{value}")
    return {
        "molecule_chembl_id": row.get("molecule_chembl_id"),
        "pref_name": row.get("pref_name"),
        "molecule_type": row.get("molecule_type"),
        "max_phase": row.get("max_phase"),
        "first_approval": row.get("first_approval"),
        "chirality": row.get("chirality"),
        "inorganic_flag": row.get("inorganic_flag"),
        "natural_product": row.get("natural_product"),
        "therapeutic_flag": row.get("therapeutic_flag"),
        "prodrug": row.get("prodrug"),
        "dosed_ingredient": row.get("dosed_ingredient"),
        "oral": row.get("oral"),
        "parenteral": row.get("parenteral"),
        "topical": row.get("topical"),
        "black_box_warning": row.get("black_box_warning"),
        "chemical_probe": row.get("chemical_probe"),
        "helm_notation": row.get("helm_notation"),
        "parent_chembl_id": hierarchy.get("parent_chembl_id"),
        "active_chembl_id": hierarchy.get("active_chembl_id"),
        "canonical_smiles": structures.get("canonical_smiles"),
        "standard_inchi": structures.get("standard_inchi"),
        "standard_inchi_key": structures.get("standard_inchi_key"),
        "full_molformula": properties.get("full_molformula"),
        "full_mwt": properties.get("full_mwt"),
        "mw_freebase": properties.get("mw_freebase"),
        "alogp": properties.get("alogp"),
        "psa": properties.get("psa"),
        "hba": properties.get("hba"),
        "hbd": properties.get("hbd"),
        "rtb": properties.get("rtb"),
        "ro5_violations": properties.get("num_ro5_violations"),
        "qed_weighted": properties.get("qed_weighted"),
        "np_likeness_score": properties.get("np_likeness_score"),
        "atc_classifications": "|".join(row.get("atc_classifications") or []),
        "synonyms": "|".join(synonyms),
        "cross_references": "|".join(sorted(set(cross_references))),
    }


def fetch_batch(ids: list[str]) -> tuple[list[dict], str | None]:
    encoded = urllib.parse.quote(";".join(ids), safe=";")
    url = (
        f"https://www.ebi.ac.uk/chembl/api/data/molecule/set/{encoded}.json"
        f"?only={urllib.parse.quote(ONLY_FIELDS, safe=',')}"
    )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MemProDB compound curation/1.0"},
    )
    last_error = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.load(response)
            return [compact_molecule(row) for row in payload.get("molecules", [])], None
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
            last_error = str(exc)
            time.sleep(min(20, 2**attempt))
    return [], last_error or "unknown_error"


def save_jsonl(records: dict[str, dict]) -> None:
    with OUT_JSONL.open("w", encoding="utf-8", newline="\n") as handle:
        for chembl_id in sorted(records, key=lambda item: int(item[6:])):
            handle.write(json.dumps(records[chembl_id], ensure_ascii=False) + "\n")


def main() -> None:
    wanted = load_wanted()
    records: dict[str, dict] = {}
    if OUT_JSONL.exists():
        with OUT_JSONL.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    row = json.loads(line)
                    records[row["molecule_chembl_id"]] = row

    remaining = [item for item in wanted if item not in records]
    chunks = [remaining[index : index + 100] for index in range(0, len(remaining), 100)]
    failures: dict[str, str] = {}
    completed = 0

    with ThreadPoolExecutor(max_workers=12) as executor:
        future_to_ids = {
            executor.submit(fetch_batch, chunk): chunk
            for chunk in chunks
        }
        for future in as_completed(future_to_ids):
            ids = future_to_ids[future]
            rows, error = future.result()
            returned = set()
            for row in rows:
                chembl_id = row.get("molecule_chembl_id")
                if chembl_id:
                    records[chembl_id] = row
                    returned.add(chembl_id)
            if error:
                for chembl_id in ids:
                    failures[chembl_id] = error
            else:
                for chembl_id in set(ids) - returned:
                    failures[chembl_id] = "not_returned"
            completed += 1
            if completed % 25 == 0:
                save_jsonl(records)
                print(
                    json.dumps(
                        {
                            "completed_chunks": completed,
                            "total_chunks": len(chunks),
                            "retrieved_total": len(records),
                            "failure_total": len(failures),
                        }
                    ),
                    flush=True,
                )

    save_jsonl(records)
    with FAILED_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chembl_id", "failure_reason"])
        for chembl_id in sorted(failures, key=lambda item: int(item[6:])):
            if chembl_id not in records:
                writer.writerow([chembl_id, failures[chembl_id]])

    report = {
        "wanted_chembl_ids": len(wanted),
        "retrieved_metadata_records": len(set(wanted) & set(records)),
        "unretrieved_records": len(set(wanted) - set(records)),
        "unretrieved_examples": sorted(set(wanted) - set(records))[:100],
        "source": "ChEMBL 37 REST API",
        "source_url": "https://www.ebi.ac.uk/chembl/api/data/docs",
        "retrieval_date": "2026-07-26",
        "output_jsonl": str(OUT_JSONL),
        "failed_tsv": str(FAILED_TSV),
    }
    REPORT.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
