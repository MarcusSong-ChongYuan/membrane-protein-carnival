import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = Path(r"C:\github-repos\upload\normalized_tables")
PROTEIN = BASE / "protein_database.tsv"
CACHE = BASE / "xlsx_build" / "pubchem_ligand_lookup_cache.json"
RCSB_CACHE = BASE / "xlsx_build" / "rcsb_chemcomp_cache.json"

csv.field_size_limit(100_000_000)

FUNCTIONAL = {"ATP", "ADP", "AMP", "ANP", "ACP", "GDP", "GTP", "GNP", "GMP", "NAD", "NAP", "FAD", "FMN", "HEM", "HEC", "PLP", "SAM", "SAH", "COA", "NCT", "ADN", "ACH", "CMP", "NMN"}


def fetch_json(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                return json.load(response)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
        except Exception:
            pass
        time.sleep(0.4 * (attempt + 1))
    return None


def pubchem_name(name):
    if not name:
        return ""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name, safe='')}/cids/TXT"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                text = response.read().decode("utf-8", errors="ignore")
            for line in text.splitlines():
                if line.strip().isdigit():
                    return line.strip()
            return ""
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return ""
        except Exception:
            pass
        time.sleep(0.4 * (attempt + 1))
    return ""


def map_code(code, rcsb_cache):
    if code in rcsb_cache:
        name = rcsb_cache[code].get("name", "")
    else:
        data = fetch_json(f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}")
        name = ""
        if data:
            name = data.get("chem_comp", {}).get("name", "") or ""
            synonyms = data.get("rcsb_chem_comp_synonyms") or []
            rcsb_cache[code] = {
                "name": name,
                "synonyms": [s.get("name", "") for s in synonyms if isinstance(s, dict)],
                "inchikey": data.get("rcsb_chem_comp_descriptor", {}).get("InChIKey", ""),
            }
        else:
            rcsb_cache[code] = {"name": "", "synonyms": [], "inchikey": ""}
    candidates = [rcsb_cache[code].get("name", "")] + rcsb_cache[code].get("synonyms", [])
    for candidate in candidates:
        cid = pubchem_name(candidate)
        if cid:
            return cid
    return ""


def main():
    counts = Counter()
    with PROTEIN.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for col in ["pdb_biolip_other", "pdb_scpdb_other", "pdbbind_other"]:
                for label in re.findall(r"([^;@]+)@", row.get(col, "")):
                    matches = re.findall(r"\(([A-Za-z0-9]{1,8})\)", label)
                    if matches:
                        counts[matches[-1].upper()] += 1

    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    rcsb_cache = json.loads(RCSB_CACHE.read_text(encoding="utf-8")) if RCSB_CACHE.exists() else {}
    codes = sorted(code for code, count in counts.items() if count >= 3 or code in FUNCTIONAL)
    pending = [code for code in codes if not cache.get(f"pubchem_name_pdb_ligand_code:{code}")]
    print(json.dumps({"candidate_codes": len(codes), "pending": len(pending)}, indent=2))

    done = 0
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(map_code, code, rcsb_cache): code for code in pending}
        for future in as_completed(futures):
            code = futures[future]
            cache[f"pubchem_name_pdb_ligand_code:{code}"] = future.result()
            done += 1
            if done % 50 == 0:
                CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
                RCSB_CACHE.write_text(json.dumps(rcsb_cache, indent=2, ensure_ascii=False), encoding="utf-8")
                print(json.dumps({"done": done, "mapped_total": sum(1 for v in cache.values() if v)}, indent=2))

    CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
    RCSB_CACHE.write_text(json.dumps(rcsb_cache, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"done": done, "mapped_total": sum(1 for v in cache.values() if v)}, indent=2))


if __name__ == "__main__":
    main()
