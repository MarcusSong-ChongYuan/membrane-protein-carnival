import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = Path(r"C:\github-repos\upload\normalized_tables")
PROTEIN = BASE / "protein_database.tsv"
CACHE = BASE / "xlsx_build" / "pubchem_ligand_lookup_cache.json"

csv.field_size_limit(100_000_000)

SKIP = {"", "HOH", "WAT", "DOD"}


def lookup_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name, safe='')}/cids/TXT"
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
    codes = set()
    with PROTEIN.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for col in ["pdb_biolip_other", "pdb_scpdb_other", "pdbbind_other"]:
                for label in re.findall(r"([^;@]+)@", row.get(col, "")):
                    matches = re.findall(r"\(([A-Za-z0-9]{1,8})\)", label)
                    if matches:
                        code = matches[-1].upper()
                        if code not in SKIP:
                            codes.add(code)

    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    pending = [code for code in sorted(codes) if f"pubchem_name_pdb_ligand_code:{code}" not in cache]
    print(json.dumps({"unique_codes": len(codes), "pending": len(pending)}, indent=2))

    done = 0
    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(lookup_name, code): code for code in pending}
        for future in as_completed(futures):
            code = futures[future]
            cache[f"pubchem_name_pdb_ligand_code:{code}"] = future.result()
            done += 1
            if done % 100 == 0:
                CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
                print(json.dumps({"done": done, "mapped_total": sum(1 for v in cache.values() if v)}, indent=2))

    CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"done": done, "mapped_total": sum(1 for v in cache.values() if v)}, indent=2))


if __name__ == "__main__":
    main()
