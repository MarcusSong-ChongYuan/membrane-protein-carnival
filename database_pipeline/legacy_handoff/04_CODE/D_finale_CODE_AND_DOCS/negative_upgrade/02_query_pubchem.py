"""
Step 2: Batch-query PubChem PUG REST API for compound properties (v3).
- Uses curl via subprocess (urllib broken on this machine)
- BATCH_SIZE=200, MAX_WORKERS=10
- Auto-detects rate-limiting and backs off
- Retries failed CIDs in smaller batches
"""
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

OUT_DIR = Path(r"D:\finale\negative_upgrade")
CIDS_FILE = OUT_DIR / "cids_to_query.txt"
CHECKPOINT = OUT_DIR / "cids_queried.txt"
OUTPUT = OUT_DIR / "pubchem_properties.jsonl"
ERROR_LOG = OUT_DIR / "pubchem_errors.txt"
RETRY_FILE = OUT_DIR / "cids_to_retry.txt"

PROPERTIES = "IsomericSMILES,InChI,InChIKey,MolecularWeight,MolecularFormula,IUPACName,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,Charge,ExactMass"

BATCH_SIZE = 200
RETRY_BATCH_SIZE = 20
MAX_WORKERS = 10
RATE_LIMIT_BACKOFF = 5.0

os.makedirs(OUT_DIR, exist_ok=True)

import threading
lock = threading.Lock()
rate_limit_until = [0.0]  # mutable for cross-thread sharing


def load_cids():
    cids = []
    with open(CIDS_FILE) as f:
        for line in f:
            line = line.strip()
            if line.isdigit():
                cids.append(line)
    return cids


def load_checkpoint():
    if not CHECKPOINT.exists():
        return set()
    done = set()
    with open(CHECKPOINT) as f:
        for line in f:
            done.add(line.strip())
    return done


def save_result(cid, props):
    with lock:
        with open(OUTPUT, "a", encoding="utf-8") as f:
            f.write(json.dumps({"cid": cid, "properties": props}, ensure_ascii=False) + "\n")
        with open(CHECKPOINT, "a") as cf:
            cf.write(cid + "\n")


def save_error(cid, msg):
    with lock:
        with open(ERROR_LOG, "a") as f:
            f.write(f"{cid}\t{msg}\n")
        with open(CHECKPOINT, "a") as cf:
            cf.write(cid + "\n")


def add_to_retry(cid):
    with lock:
        with open(RETRY_FILE, "a") as f:
            f.write(cid + "\n")


def wait_if_rate_limited():
    """If we've been rate-limited, wait until the backoff expires."""
    while time.time() < rate_limit_until[0]:
        time.sleep(0.5)


def query_batch(cids_batch):
    """Query PubChem via curl. Returns (results_dict, success, rate_limited)."""
    if not cids_batch:
        return {}, True, False

    cids_str = ",".join(cids_batch)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids_str}/property/{PROPERTIES}/JSON"

    for attempt in range(3):
        try:
            result = subprocess.run(
                ["curl", "-s", "--connect-timeout", "15", "--max-time", "60",
                 "-H", "Accept: application/json", url],
                capture_output=True, text=True, timeout=70
            )

            if result.returncode != 0:
                if attempt < 2:
                    time.sleep(2 ** attempt * 3)
                continue

            # Check for rate limiting in response
            stdout = result.stdout
            if not stdout or len(stdout) < 10:
                # Empty response = likely rate limited
                rate_limit_until[0] = time.time() + RATE_LIMIT_BACKOFF * (attempt + 1)
                if attempt < 2:
                    time.sleep(RATE_LIMIT_BACKOFF * (attempt + 1))
                continue

            if "PUGREST.NotFound" in stdout or "Too Many Requests" in stdout or "Service Unavailable" in stdout:
                rate_limit_until[0] = time.time() + RATE_LIMIT_BACKOFF * (attempt + 1)
                if attempt < 2:
                    time.sleep(RATE_LIMIT_BACKOFF * (attempt + 1))
                continue

            data = json.loads(stdout)
            props_list = data.get("PropertyTable", {}).get("Properties", [])
            results = {}
            for entry in props_list:
                cid = str(entry.get("CID", ""))
                # PubChem returns 'SMILES' key when IsomericSMILES is requested,
                # 'ConnectivitySMILES' when CanonicalSMILES is requested
                smiles = entry.get("SMILES", entry.get("ConnectivitySMILES", "")) or ""
                results[cid] = {
                    "SMILES": smiles,
                    "InChI": entry.get("InChI", "") or "",
                    "InChIKey": entry.get("InChIKey", "") or "",
                    "MolecularWeight": entry.get("MolecularWeight", ""),
                    "MolecularFormula": entry.get("MolecularFormula", "") or "",
                    "IUPACName": entry.get("IUPACName", "") or "",
                    "XLogP": entry.get("XLogP", ""),
                    "TPSA": entry.get("TPSA", ""),
                    "HBondDonorCount": entry.get("HBondDonorCount", ""),
                    "HBondAcceptorCount": entry.get("HBondAcceptorCount", ""),
                    "RotatableBondCount": entry.get("RotatableBondCount", ""),
                    "HeavyAtomCount": entry.get("HeavyAtomCount", ""),
                    "Charge": entry.get("Charge", ""),
                    "ExactMass": entry.get("ExactMass", ""),
                }
            return results, True, False

        except json.JSONDecodeError:
            if attempt < 2:
                time.sleep(2 ** attempt * 3)
            continue
        except subprocess.TimeoutExpired:
            if attempt < 2:
                time.sleep(2 ** attempt * 3)
            continue
        except Exception:
            if attempt < 2:
                time.sleep(2 ** attempt * 3)
            continue

    return {}, False, False


def process_batches(cids, done_set, label="initial", batch_size=BATCH_SIZE):
    """Process CIDs in batches."""
    remaining = [c for c in cids if c not in done_set]
    if not remaining:
        return

    batches = [remaining[i:i + batch_size] for i in range(0, len(remaining), batch_size)]
    print(f"\n  [{label}] {len(remaining):,} CIDs in {len(batches):,} batches (size={batch_size})")

    completed = 0
    success = 0
    failed = 0
    t0 = time.time()
    last_report = 0

    def process_one_batch(batch):
        wait_if_rate_limited()
        results, ok, rate_limited = query_batch(batch)
        local_success = 0
        local_failed = 0
        for cid in batch:
            if cid in results:
                save_result(cid, results[cid])
                local_success += 1
            else:
                save_error(cid, "no_data_returned")
                add_to_retry(cid)  # queue for retry
                local_failed += 1
        return local_success, local_failed

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_one_batch, b): i for i, b in enumerate(batches)}
        for future in as_completed(futures):
            try:
                s, f = future.result()
                success += s
                failed += f
                completed += (s + f)
                now = time.time()
                if now - last_report > 20:
                    elapsed = now - t0
                    rate = completed / elapsed if elapsed > 0 else 0
                    eta = (len(remaining) - completed) / rate if rate > 0 else 0
                    print(f"    {completed:,}/{len(remaining):,} ({100*completed/len(remaining):.0f}%) "
                          f"ok={success} fail={failed} | {rate:.0f}/s | ETA {eta/60:.0f}m")
                    last_report = now
            except Exception as e:
                print(f"  [ERROR] {e}")

    elapsed = time.time() - t0
    print(f"  [{label}] Done: {success:,} ok, {failed:,} failed in {elapsed:.0f}s ({completed/elapsed:.0f}/s)")


def main():
    print("=" * 60)
    print("Step 2: Query PubChem (v3 - curl, optimized)")
    print("=" * 60)

    cids = load_cids()
    print(f"  CIDs: {len(cids):,}")

    # Clean up bad previous runs
    done_set = load_checkpoint()
    result_count = sum(1 for _ in open(OUTPUT, encoding="utf-8")) if OUTPUT.exists() else 0
    error_count = sum(1 for _ in open(ERROR_LOG)) if ERROR_LOG.exists() else 0
    print(f"  Existing: {result_count:,} results, {error_count:,} errors, {len(done_set):,} checkpointed")

    if result_count == 0:
        print("  Starting fresh...")
        CHECKPOINT.unlink(missing_ok=True)
        if ERROR_LOG.exists():
            ERROR_LOG.unlink(missing_ok=True)
        if RETRY_FILE.exists():
            RETRY_FILE.unlink(missing_ok=True)
        done_set = set()

    # Phase 1: initial batch query
    process_batches(cids, done_set, label="Phase 1", batch_size=BATCH_SIZE)

    # Phase 2: retry failed CIDs with smaller batches
    if RETRY_FILE.exists():
        retry_cids = []
        with open(RETRY_FILE) as f:
            for line in f:
                c = line.strip()
                if c.isdigit():
                    retry_cids.append(c)

        # Reset: remove retry CIDs from checkpoint so they get re-queried
        if retry_cids:
            print(f"\n  Retrying {len(retry_cids):,} failed CIDs with batch_size={RETRY_BATCH_SIZE}...")
            # Clear the retry file for fresh retry tracking
            RETRY_FILE.unlink()
            # Remove from done_set so they get re-processed
            retry_set = set(retry_cids)
            old_done = load_checkpoint()
            # Write new checkpoint without retry CIDs
            new_done = old_done - retry_set
            CHECKPOINT.unlink()
            for c in new_done:
                with open(CHECKPOINT, "a") as cf:
                    cf.write(c + "\n")
            done_set = new_done

            process_batches(retry_cids, done_set, label="Phase 2 retry", batch_size=RETRY_BATCH_SIZE)

    # Final stats
    final_results = sum(1 for _ in open(OUTPUT, encoding="utf-8")) if OUTPUT.exists() else 0
    final_errors = sum(1 for _ in open(ERROR_LOG)) if ERROR_LOG.exists() else 0
    print(f"\n  FINAL: {final_results:,} with properties, {final_errors:,} without")
    print(f"  Coverage: {final_results/len(cids)*100:.1f}%")
    print("\nDone.")


if __name__ == "__main__":
    main()
