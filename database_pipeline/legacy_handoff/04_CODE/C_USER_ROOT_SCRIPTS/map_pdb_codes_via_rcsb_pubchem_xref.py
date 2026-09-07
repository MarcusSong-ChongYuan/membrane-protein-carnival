import csv
import json
import re
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = Path(r"C:\github-repos\upload\normalized_tables")
REVIEW = BASE / "protein_structure_ligand_context_review.tsv"
CACHE = BASE / "xlsx_build" / "pubchem_ligand_lookup_cache.json"
RCSB_XREF_CACHE = BASE / "xlsx_build" / "rcsb_pubchem_xref_cache.json"

csv.field_size_limit(100_000_000)


def fetch_rcsb_pubchem(code):
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                data = json.load(response)
            related = data.get("rcsb_chem_comp_related") or []
            for item in related:
                if item.get("resource_name") == "PubChem":
                    accession = str(item.get("resource_accession_code", "")).strip()
                    if accession.isdigit():
                        return accession
            return ""
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return ""
        except Exception:
            pass
        time.sleep(0.3 * (attempt + 1))
    return ""


def main():
    rows = []
    with REVIEW.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    codes = sorted({
        row["pdb_ligand_code"]
        for row in rows
        if row.get("cid_mapping_status") != "mapped" and row.get("pdb_ligand_code")
    })

    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    xref_cache = json.loads(RCSB_XREF_CACHE.read_text(encoding="utf-8")) if RCSB_XREF_CACHE.exists() else {}

    pending = [code for code in codes if code not in xref_cache]
    print(json.dumps({"unique_unmapped_codes": len(codes), "pending_rcsb": len(pending)}, indent=2))

    done = 0
    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(fetch_rcsb_pubchem, code): code for code in pending}
        for future in as_completed(futures):
            code = futures[future]
            cid = future.result()
            xref_cache[code] = cid
            if cid:
                cache[f"rcsb_pubchem_xref:{code}"] = cid
            done += 1
            if done % 100 == 0:
                RCSB_XREF_CACHE.write_text(json.dumps(xref_cache, indent=2, ensure_ascii=False), encoding="utf-8")
                CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
                print(json.dumps({"done": done, "mapped_xref": sum(1 for v in xref_cache.values() if v)}, indent=2))

    RCSB_XREF_CACHE.write_text(json.dumps(xref_cache, indent=2, ensure_ascii=False), encoding="utf-8")
    CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({
        "done": done,
        "mapped_xref_total": sum(1 for v in xref_cache.values() if v),
    }, indent=2))


if __name__ == "__main__":
    main()
