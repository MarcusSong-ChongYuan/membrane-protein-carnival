"""
Step 5: Resolve PubChem SID → CID for the remaining 2,641 records.
Uses curl subprocess (urllib broken on this machine).
"""
import csv, gzip, json, os, subprocess, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

OUT_DIR = Path(r"D:\finale\negative_upgrade")
UNMAPPED = r"D:\finale\01_正式数据_V6.2\negative_binding_evidence_unmapped_review_v1_2.tsv.gz"

SID_FILE = OUT_DIR / "sids_to_query.txt"
SID_RESULTS = OUT_DIR / "sid_to_cid.jsonl"
SID_ERRORS = OUT_DIR / "sid_errors.txt"

MAX_WORKERS = 10
BATCH_SIZE = 100  # Try batch SID→CID
os.makedirs(OUT_DIR, exist_ok=True)

lock = threading.Lock()


def extract_sids():
    """Extract unique SIDs from unmapped records."""
    sids = set()
    sid_data = {}  # sid -> list of record data for backfill

    with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sid = row.get('pubchem_sid', '').strip()
            if sid:
                sids.add(sid)
                if sid not in sid_data:
                    sid_data[sid] = []
                sid_data[sid].append({
                    'source_evidence_id': row['source_evidence_id'],
                    'compound_source_id': row.get('compound_source_id', ''),
                    'target_uniprot_id': row.get('target_uniprot_id', ''),
                })

    with open(SID_FILE, 'w') as f:
        for s in sorted(sids, key=lambda x: int(x)):
            f.write(s + '\n')

    return list(sids), sid_data


def query_sid_batch(sids_batch):
    """Query PubChem SID→CID for a batch of SIDs."""
    sids_str = ",".join(sids_batch)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sids_str}/cids/JSON"

    for attempt in range(3):
        try:
            result = subprocess.run(
                ["curl", "-s", "--connect-timeout", "15", "--max-time", "60",
                 "-H", "Accept: application/json", url],
                capture_output=True, text=True, timeout=70
            )

            if result.returncode != 0:
                if attempt < 2:
                    time.sleep(2 ** attempt * 2)
                continue

            stdout = result.stdout
            if not stdout or len(stdout) < 5:
                if attempt < 2:
                    time.sleep(3)
                continue

            if "Too Many Requests" in stdout or "Service Unavailable" in stdout:
                if attempt < 2:
                    time.sleep(5)
                continue

            data = json.loads(stdout)

            # Response format: {"IdentifierList": {"CID": [123, 456]}}
            id_list = data.get("IdentifierList", {})
            cids = id_list.get("CID", [])

            # Build mapping: CID at index i corresponds to SID at index i
            results = {}
            for i, sid in enumerate(sids_batch):
                if i < len(cids) and cids[i] != 0:  # 0 means no mapping
                    results[sid] = str(cids[i])

            return results, True

        except json.JSONDecodeError:
            if attempt < 2:
                time.sleep(2)
            continue
        except subprocess.TimeoutExpired:
            if attempt < 2:
                time.sleep(3)
            continue
        except Exception:
            if attempt < 2:
                time.sleep(2)
            continue

    return {}, False


def save_result(sid, cid):
    with lock:
        with open(SID_RESULTS, 'a') as f:
            f.write(json.dumps({"sid": sid, "cid": cid}) + "\n")


def save_error(sid, msg):
    with lock:
        with open(SID_ERRORS, 'a') as f:
            f.write(f"{sid}\t{msg}\n")


def main():
    print("=" * 60)
    print("Step 5: SID → CID Resolution")
    print("=" * 60)

    # Extract
    print("\nExtracting SIDs from remaining unmapped records...")
    sids_list, sid_data = extract_sids()
    print(f"  Unique SIDs: {len(sids_list):,}")
    print(f"  Total records: {sum(len(v) for v in sid_data.values()):,}")

    # Query in batches
    batches = [sids_list[i:i + BATCH_SIZE] for i in range(0, len(sids_list), BATCH_SIZE)]
    print(f"\nQuerying {len(batches):,} batches (size={BATCH_SIZE}, workers={MAX_WORKERS})...")

    # Clear output files
    SID_RESULTS.unlink(missing_ok=True)
    SID_ERRORS.unlink(missing_ok=True)

    success = 0
    failed = 0
    t0 = time.time()

    def process_batch(batch):
        results, ok = query_sid_batch(batch)
        local_success = 0
        local_failed = 0
        for sid in batch:
            if sid in results:
                save_result(sid, results[sid])
                local_success += 1
            else:
                save_error(sid, "no_cid_found")
                local_failed += 1
        return local_success, local_failed

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_batch, b): i for i, b in enumerate(batches)}
        for future in as_completed(futures):
            s, f = future.result()
            success += s
            failed += f
            total_done = success + failed
            if total_done % 100 == 0 or total_done == len(sids_list):
                elapsed = time.time() - t0
                print(f"  {total_done}/{len(sids_list)} | ok={success} fail={failed} | {elapsed:.0f}s")

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.0f}s: {success:,} mapped, {failed:,} failed")
    print(f"  Coverage: {success/len(sids_list)*100:.1f}%")

    # Load results
    sid_to_cid = {}
    if SID_RESULTS.exists():
        with open(SID_RESULTS) as f:
            for line in f:
                r = json.loads(line)
                sid_to_cid[r['sid']] = r['cid']

    # Calculate how many records can now be backfilled
    records_savable = sum(len(sid_data[s]) for s in sid_to_cid)
    records_unsavable = sum(len(sid_data[s]) for s in sids_list if s not in sid_to_cid)
    print(f"\n  Records now mappable: {records_savable:,}")
    print(f"  Records still unmappable: {records_unsavable:,}")

    print("\nDone. Run 05b_backfill_sid.py to backfill these records.")
    return sid_to_cid


if __name__ == "__main__":
    main()
