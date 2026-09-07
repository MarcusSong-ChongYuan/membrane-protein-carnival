"""
Query PubChem substance synonyms for all 963 unmapped SIDs.
Extract ChEMBL IDs, DrugBank IDs, etc. for cross-referencing.
"""
import csv, gzip, json, subprocess, time, threading, re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from collections import Counter

OUT_DIR = Path(r"D:\finale\negative_upgrade")
UNMAPPED = r"D:\finale\01_正式数据_V6.2\negative_binding_evidence_unmapped_review_v1_2.tsv.gz"

SYNONYMS_FILE = OUT_DIR / "sid_synonyms.jsonl"
SYNONYMS_ERRORS = OUT_DIR / "sid_synonyms_errors.txt"

MAX_WORKERS = 10
os = __import__('os')
lock = threading.Lock()


def extract_sids():
    sids = set()
    with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sid = row.get('pubchem_sid', '').strip()
            if sid:
                sids.add(sid)
    return sorted(sids, key=lambda x: int(x))


def query_synonyms(sid):
    """Get synonyms for a single SID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/synonyms/JSON"

    for attempt in range(3):
        try:
            result = subprocess.run(
                ["curl", "-s", "--connect-timeout", "10", "--max-time", "30",
                 "-H", "Accept: application/json", url],
                capture_output=True, text=True, timeout=35
            )

            if result.returncode != 0:
                if attempt < 2: time.sleep(1)
                continue

            stdout = result.stdout
            if not stdout or len(stdout) < 5:
                if attempt < 2: time.sleep(2)
                continue

            if "Too Many Requests" in stdout or "Service Unavailable" in stdout:
                if attempt < 2: time.sleep(3)
                continue

            data = json.loads(stdout)
            info_list = data.get('InformationList', {}).get('Information', [])
            if info_list:
                synonyms = info_list[0].get('Synonym', [])
                return synonyms, True

            return [], True  # No synonyms, not an error

        except json.JSONDecodeError:
            if attempt < 2: time.sleep(1)
            continue
        except Exception:
            if attempt < 2: time.sleep(1)
            continue

    return [], False


def save_synonyms(sid, synonyms):
    with lock:
        with open(SYNONYMS_FILE, 'a', encoding='utf-8') as f:
            f.write(json.dumps({"sid": sid, "synonyms": synonyms}, ensure_ascii=False) + "\n")


def save_error(sid, msg):
    with lock:
        with open(SYNONYMS_ERRORS, 'a') as f:
            f.write(f"{sid}\t{msg}\n")


def main():
    print("=" * 60)
    print("Step 5c: SID Synonym Extraction")
    print("=" * 60)

    sids = extract_sids()
    print(f"  Unique SIDs: {len(sids):,}")

    # Clear outputs
    SYNONYMS_FILE.unlink(missing_ok=True)
    SYNONYMS_ERRORS.unlink(missing_ok=True)

    success = 0
    failed = 0
    with_synonyms = 0
    t0 = time.time()

    def process_sid(sid):
        synonyms, ok = query_synonyms(sid)
        if ok:
            save_synonyms(sid, synonyms)
        else:
            save_error(sid, "query_failed")
        return ok, len(synonyms)

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_sid, s): s for s in sids}
        for future in as_completed(futures):
            ok, syn_count = future.result()
            if ok:
                success += 1
                if syn_count > 0:
                    with_synonyms += 1
            else:
                failed += 1
            total = success + failed
            if total % 100 == 0 or total == len(sids):
                elapsed = time.time() - t0
                rate = total / elapsed if elapsed > 0 else 0
                eta = (len(sids) - total) / rate if rate > 0 else 0
                print(f"  {total}/{len(sids)} | ok={success} fail={failed} "
                      f"syn={with_synonyms} | {rate:.0f}/s ETA {eta:.0f}s")

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.0f}s")
    print(f"  Queried: {success:,}, Failed: {failed:,}")
    print(f"  SIDs with synonyms: {with_synonyms:,}")

    # Analyze synonym types
    print("\n" + "=" * 60)
    print("Synonym Analysis")
    print("=" * 60)

    # Patterns to extract
    patterns = {
        'CHEMBL': re.compile(r'^CHEMBL\d+$'),
        'CHEBI': re.compile(r'^CHEBI:\d+$'),
        'DRUGBANK': re.compile(r'^DB\d{5}$'),
        'UNII': re.compile(r'^[A-Z0-9]{10}$'),
        'CAS': re.compile(r'^\d{2,7}-\d{2}-\d$'),
        'ZINC': re.compile(r'^ZINC\d+$'),
        'DRUGCENTRAL': re.compile(r'^[A-Z]{2,4}\d{4,8}$'),  # DrugCentral IDs
        'InChIKey': re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'),
        'InChI': re.compile(r'^InChI=1S?/'),
        'SMILES': re.compile(r'^[A-Za-z0-9@+\-\[\]\(\)\\\/\=#\.%]+$'),  # Fuzzy
    }

    extracted = {k: {} for k in patterns}  # {type: {sid: [values]}}
    total_found = Counter()
    sid_extracted_count = 0

    if SYNONYMS_FILE.exists():
        with open(SYNONYMS_FILE, encoding='utf-8') as f:
            for line in f:
                entry = json.loads(line)
                sid = entry['sid']
                syns = entry['synonyms']
                extracted_any = False
                for syn in syns:
                    for ptype, pat in patterns.items():
                        if pat.match(syn.strip()):
                            if sid not in extracted[ptype]:
                                extracted[ptype][sid] = []
                            extracted[ptype][sid].append(syn.strip())
                            total_found[ptype] += 1
                            extracted_any = True
                if extracted_any:
                    sid_extracted_count += 1

    print(f"\n  SIDs with extracted identifiers: {sid_extracted_count}/{len(sids)}")
    for ptype in patterns:
        count = len(extracted[ptype])
        if count > 0:
            print(f"  {ptype}: {count} SIDs, {total_found[ptype]} total values")
            # Show examples
            for sid, vals in list(extracted[ptype].items())[:3]:
                print(f"    SID {sid}: {vals[:3]}")

    # Save extracted data for Step 5d
    extracted_file = OUT_DIR / "sid_extracted_ids.json"
    # Convert to serializable format
    serializable = {}
    for ptype, sid_vals in extracted.items():
        if sid_vals:
            serializable[ptype] = {sid: vals for sid, vals in sid_vals.items()}

    with open(extracted_file, 'w', encoding='utf-8') as f:
        json.dump(serializable, f, ensure_ascii=False, indent=2)
    print(f"\n  Saved extracted IDs to: {extracted_file}")

    return serializable


if __name__ == "__main__":
    main()
