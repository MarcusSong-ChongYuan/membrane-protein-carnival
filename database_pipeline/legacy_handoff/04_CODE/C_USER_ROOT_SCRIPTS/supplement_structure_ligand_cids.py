import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


BASE = Path(r"C:\github-repos\upload\normalized_tables")
REVIEW = BASE / "protein_structure_ligand_context_review.tsv"
BUILD = BASE / "xlsx_build"
PUBCHEM_CACHE = BUILD / "pubchem_ligand_lookup_cache.json"
RCSB_CHEMCOMP_CACHE = BUILD / "rcsb_chemcomp_cache.json"
REPORT = BASE / "structure_ligand_cid_supplement_report.json"


def load_json(path):
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {}


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


def url_json(url, timeout=25):
    req = urllib.request.Request(url, headers={"User-Agent": "codex-ligand-cid-curation/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8", errors="ignore"))


def url_text(url, timeout=35):
    req = urllib.request.Request(url, headers={"User-Agent": "codex-ligand-cid-curation/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="ignore").strip()


def pubchem_get_cid(namespace, value):
    if not value:
        return ""
    encoded = urllib.parse.quote(value, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{namespace}/{encoded}/cids/TXT"
    for attempt in range(2):
        try:
            text = url_text(url, timeout=12)
            for line in text.splitlines():
                line = line.strip()
                if line.isdigit():
                    return line
            return ""
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return ""
            if exc.code in {400, 405}:
                return ""
        except Exception:
            pass
        time.sleep(0.8 * (attempt + 1))
    return ""


def pubchem_get_cid_by_smiles(smiles):
    if not smiles:
        return ""
    data = urllib.parse.urlencode({"smiles": smiles}).encode("utf-8")
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT"
    for attempt in range(2):
        try:
            req = urllib.request.Request(
                url,
                data=data,
                headers={
                    "User-Agent": "codex-ligand-cid-curation/1.0",
                    "Content-Type": "application/x-www-form-urlencoded",
                },
            )
            with urllib.request.urlopen(req, timeout=12) as response:
                text = response.read().decode("utf-8", errors="ignore").strip()
            for line in text.splitlines():
                line = line.strip()
                if line.isdigit():
                    return line
            return ""
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return ""
            if exc.code in {400, 405}:
                return ""
        except Exception:
            pass
        time.sleep(0.8 * (attempt + 1))
    return ""


def unichem_get_pubchem_cid_by_inchikey(inchikey):
    if not inchikey:
        return ""
    url = f"https://www.ebi.ac.uk/unichem/rest/inchikey/{urllib.parse.quote(inchikey)}"
    for attempt in range(2):
        try:
            data = url_json(url, timeout=18)
            for item in data:
                if str(item.get("src_id", "")) == "22":
                    accession = str(item.get("src_compound_id", "")).strip()
                    if accession.isdigit():
                        return accession
            return ""
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return ""
        except Exception:
            pass
        time.sleep(0.8 * (attempt + 1))
    return ""


def get_chemcomp(code, cache):
    if code in cache:
        return cache[code]
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{urllib.parse.quote(code)}"
    try:
        data = url_json(url)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            cache[code] = {}
            return {}
        raise
    except Exception:
        return {}

    comp = data.get("chem_comp", {}) or {}
    desc = data.get("rcsb_chem_comp_descriptor", {}) or {}
    synonyms = []
    for item in data.get("rcsb_chem_comp_synonyms", []) or []:
        value = item.get("name") or item.get("comp_id") or ""
        if value and value not in synonyms:
            synonyms.append(value)
    name = comp.get("name", "")
    if name and name not in synonyms:
        synonyms.insert(0, name)

    related_pubchem = []
    for item in data.get("rcsb_chem_comp_related", []) or []:
        if str(item.get("resource_name", "")).lower() == "pubchem":
            accession = str(item.get("resource_accession_code", "")).strip()
            if accession.isdigit() and accession not in related_pubchem:
                related_pubchem.append(accession)

    cache[code] = {
        "name": name,
        "synonyms": synonyms,
        "inchikey": desc.get("InChIKey", ""),
        "inchi": desc.get("InChI", ""),
        "smiles": desc.get("SMILES_stereo") or desc.get("SMILES", ""),
        "pubchem_xrefs": related_pubchem,
    }
    return cache[code]


def main():
    rows = []
    with REVIEW.open(newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if (
                row.get("cid_mapping_status") == "unmapped"
                and row.get("ligand_category") != "empty_or_unparsed"
                and (row.get("pdb_ligand_code") or row.get("chembl_id") or row.get("ligand_name_raw"))
            ):
                rows.append(row)

    codes = sorted({row["pdb_ligand_code"] for row in rows if row.get("pdb_ligand_code")})
    chembl_ids = sorted({row["chembl_id"] for row in rows if row.get("chembl_id")})
    pubchem_cache = load_json(PUBCHEM_CACHE)
    chemcomp_cache = load_json(RCSB_CHEMCOMP_CACHE)

    # First fetch/cache RCSB details. This is cheap and gives us InChIKeys for
    # the second-pass cross-reference lookup.
    for index, code in enumerate(codes, 1):
        get_chemcomp(code, chemcomp_cache)
        if index % 25 == 0:
            save_json(RCSB_CHEMCOMP_CACHE, chemcomp_cache)
            print(f"cached RCSB chemcomp {index}/{len(codes)}")
    save_json(RCSB_CHEMCOMP_CACHE, chemcomp_cache)

    code_results = {}
    method_counts = {}

    def resolve_one_code(code):
        detail = chemcomp_cache.get(code, {})
        cid = ""
        method = ""
        if detail.get("pubchem_xrefs"):
            return code, detail["pubchem_xrefs"][0], "rcsb_pubchem_xref"
        if detail.get("inchikey"):
            cache_key = f"unichem_inchikey_pubchem:{detail['inchikey']}"
            cached = pubchem_cache.get(cache_key)
            if cached is None:
                cached = unichem_get_pubchem_cid_by_inchikey(detail["inchikey"])
                pubchem_cache[cache_key] = cached
            if cached:
                return code, cached, "unichem_inchikey_pubchem"
        return code, cid, method

    # UniChem is less congested than PubChem PUG for InChIKey lookups here.
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(resolve_one_code, code): code for code in codes}
        for index, future in enumerate(as_completed(futures), 1):
            code, cid, method = future.result()
            if cid:
                pubchem_cache[f"{method}:{code}"] = cid
                pubchem_cache[f"rcsb_pubchem_xref:{code}"] = cid
                method_counts[method] = method_counts.get(method, 0) + 1
            detail = chemcomp_cache.get(code, {})
            code_results[code] = {
                "cid": cid,
                "method": method,
                "inchikey": detail.get("inchikey", ""),
                "name": detail.get("name", ""),
                "synonym_count": len(detail.get("synonyms", [])),
            }
            if index % 10 == 0:
                save_json(PUBCHEM_CACHE, pubchem_cache)
                print(f"resolved UniChem {index}/{len(codes)} codes; mapped so far {sum(1 for x in code_results.values() if x['cid'])}")

    save_json(PUBCHEM_CACHE, pubchem_cache)

    for index, code in enumerate(codes, 1):
        if code_results.get(code, {}).get("cid"):
            continue
        detail = get_chemcomp(code, chemcomp_cache)
        cid = ""
        method = ""
        if detail.get("pubchem_xrefs"):
            cid = detail["pubchem_xrefs"][0]
            method = "rcsb_pubchem_xref"
        if not cid and detail.get("inchikey"):
            cache_key = f"pubchem_inchikey_pdb_ligand_code:{code}"
            cid = pubchem_cache.get(cache_key, "")
            if cache_key not in pubchem_cache:
                cid = pubchem_get_cid("inchikey", detail["inchikey"])
                pubchem_cache[cache_key] = cid
            if cid:
                method = "pubchem_inchikey_pdb_ligand_code"
        if not cid and detail.get("smiles"):
            cache_key = f"pubchem_smiles_pdb_ligand_code:{code}"
            cid = pubchem_cache.get(cache_key, "")
            if cache_key not in pubchem_cache:
                cid = pubchem_get_cid_by_smiles(detail["smiles"])
                pubchem_cache[cache_key] = cid
            if cid:
                method = "pubchem_smiles_pdb_ligand_code"
        if not cid:
            for name in detail.get("synonyms", [])[:5]:
                cache_key = f"pubchem_name_rcsb_chemcomp:{code}:{name}"
                cid = pubchem_cache.get(cache_key, "")
                if cache_key not in pubchem_cache:
                    cid = pubchem_get_cid("name", name)
                    pubchem_cache[cache_key] = cid
                if cid:
                    method = "pubchem_name_rcsb_chemcomp"
                    break

        if cid:
            # Make the existing review-builder pick this up without broad refactoring.
            pubchem_cache[f"{method}:{code}"] = cid
            pubchem_cache[f"rcsb_pubchem_xref:{code}"] = cid
            method_counts[method] = method_counts.get(method, 0) + 1
        code_results[code] = {
            "cid": cid,
            "method": method,
            "inchikey": detail.get("inchikey", ""),
            "name": detail.get("name", ""),
            "synonym_count": len(detail.get("synonyms", [])),
        }
        if index % 5 == 0:
            save_json(PUBCHEM_CACHE, pubchem_cache)
            save_json(RCSB_CHEMCOMP_CACHE, chemcomp_cache)
            print(f"processed {index}/{len(codes)} codes; mapped so far {sum(1 for x in code_results.values() if x['cid'])}")
            time.sleep(0.2)

    chembl_results = {}
    for chembl_id in chembl_ids:
        cache_key = f"pubchem_name_chembl:{chembl_id}"
        cid = pubchem_cache.get(cache_key, "")
        if not cid:
            cid = pubchem_get_cid("name", chembl_id)
            pubchem_cache[cache_key] = cid
        chembl_results[chembl_id] = cid

    save_json(PUBCHEM_CACHE, pubchem_cache)
    save_json(RCSB_CHEMCOMP_CACHE, chemcomp_cache)

    mapped_codes = {code: item for code, item in code_results.items() if item["cid"]}
    report = {
        "input_unmapped_non_empty_records": len(rows),
        "unique_pdb_ligand_codes": len(codes),
        "unique_chembl_ids": len(chembl_ids),
        "mapped_unique_pdb_ligand_codes": len(mapped_codes),
        "mapping_method_counts_by_code": method_counts,
        "chembl_results": chembl_results,
        "unmapped_unique_pdb_ligand_codes": {
            code: item for code, item in code_results.items() if not item["cid"]
        },
        "mapped_unique_pdb_ligand_codes_sample": dict(list(mapped_codes.items())[:50]),
    }
    save_json(REPORT, report)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
