import csv
import json
import re
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = Path(r"C:\github-repos\upload\normalized_tables")
PROTEIN = BASE / "protein_database.tsv"
CACHE = BASE / "xlsx_build" / "pubchem_ligand_lookup_cache.json"

csv.field_size_limit(100_000_000)


def lookup_chembl(chembl_id):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{chembl_id}/cids/TXT"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                text = response.read().decode("utf-8", errors="ignore")
            for line in text.splitlines():
                line = line.strip()
                if line.isdigit():
                    return line
            return ""
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return ""
        except Exception:
            pass
        time.sleep(0.3 * (attempt + 1))
    return ""


def main():
    chembl_ids = set()
    with PROTEIN.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for col in ["pdb_biolip_other", "pdb_scpdb_other", "pdbbind_other"]:
                chembl_ids.update(m.upper() for m in re.findall(r"CHEMBL\d+", row.get(col, ""), re.I))

    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    pending = [chembl for chembl in sorted(chembl_ids) if f"pubchem_name_chembl:{chembl}" not in cache]
    print(json.dumps({"unique_chembl": len(chembl_ids), "pending": len(pending)}, indent=2))

    done = 0
    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(lookup_chembl, chembl): chembl for chembl in pending}
        for future in as_completed(futures):
            chembl = futures[future]
            cid = future.result()
            cache[f"pubchem_name_chembl:{chembl}"] = cid
            done += 1
            if done % 100 == 0:
                CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
                print(json.dumps({"done": done, "mapped": sum(1 for v in cache.values() if v)}, indent=2))

    CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"done": done, "mapped_total": sum(1 for v in cache.values() if v)}, indent=2))


if __name__ == "__main__":
    main()
