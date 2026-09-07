import csv, gzip, hashlib, json, time, urllib.error, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites")
IDS = [x.strip().lower() for x in (ROOT / "pdb_ids_for_frozen_rescue.txt").read_text().splitlines() if x.strip()]
CP = ROOT / "sifts_residue_mapping_v7.checkpoint.tsv"
STATUS = ROOT / "sifts_download_status_v7.tsv"
SUMMARY = ROOT / "SIFTS_FROZEN_RESCUE_DOWNLOAD_SUMMARY.json"

def lname(tag): return tag.split("}")[-1]

def fetch(pid):
    url = f"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{pid}.xml.gz"
    last = ""
    for attempt, delay in enumerate((0, 1, 5), 1):
        if delay: time.sleep(delay)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MemPro-V7-SIFTS/1.1"})
            with urllib.request.urlopen(req, timeout=60) as response: blob = response.read()
            raw = gzip.decompress(blob) if blob[:2] == b"\x1f\x8b" else blob
            sha = hashlib.sha256(raw).hexdigest(); root = ET.fromstring(raw); rows = []
            for entity in root.iter():
                if lname(entity.tag) != "entity": continue
                asym = entity.get("entityId", "")
                for residue in entity.iter():
                    if lname(residue.tag) != "residue": continue
                    pdb = {}; uni = {}; details = []
                    for child in residue.iter():
                        name = lname(child.tag)
                        if name == "crossRefDb":
                            src = child.get("dbSource", "").lower()
                            if src in {"pdb", "pdbe"} and not pdb: pdb = dict(child.attrib)
                            elif src == "uniprot" and not uni: uni = dict(child.attrib)
                        elif name == "residueDetail": details.append((child.text or "") + " " + child.get("property", ""))
                    if not uni: continue
                    chain = pdb.get("dbChainId") or residue.get("dbChainId", "")
                    rn = pdb.get("dbResNum") or residue.get("dbResNum", ""); icode = ""
                    if rn and rn[-1:].isalpha(): icode = rn[-1]; rn = rn[:-1]
                    detail = " ".join(details).lower(); observed = "unobserved" if "not_observed" in detail or "unobserved" in detail else "observed"
                    rows.append([pid.upper(),chain,asym,uni.get("dbAccessionId",""),uni.get("dbResNum",""),rn,icode,pdb.get("dbResName") or residue.get("dbResName",""),observed,"mapped" if uni.get("dbResNum") else "no_uniprot_number",sha])
            return pid, "ok" if rows else "no_uniprot_mapping", rows, sha, attempt, ""
        except urllib.error.HTTPError as exc:
            last = f"HTTP {exc.code}"
            if exc.code == 404: return pid, "not_found", [], "", attempt, last
        except Exception as exc: last = repr(exc)[:300]
    return pid, "failed", [], "", 3, last

done = set()
with STATUS.open(encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"): done.add(row["pdb_id"].lower())
pending = [pid for pid in IDS if pid not in done]
print(json.dumps({"requested":len(IDS),"already_done":len(done & set(IDS)),"pending":len(pending)}), flush=True)

batch_size = 200
for start in range(0, len(pending), batch_size):
    batch = pending[start:start+batch_size]; results = []
    with ThreadPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(fetch, pid) for pid in batch]
        for future in as_completed(futures): results.append(future.result())
    with CP.open("a", encoding="utf-8", newline="") as mapping_out, STATUS.open("a", encoding="utf-8", newline="") as status_out:
        wm = csv.writer(mapping_out, delimiter="\t"); ws = csv.writer(status_out, delimiter="\t")
        for pid, status, rows, sha, attempts, message in sorted(results):
            wm.writerows(rows); ws.writerow([pid.upper(),status,len(rows),sha,attempts,message])
    print(f"checkpoint {min(start+len(batch),len(pending))}/{len(pending)}", flush=True)

status_rows = list(csv.DictReader(STATUS.open(encoding="utf-8"), delimiter="\t")); counts = {}
for row in status_rows: counts[row["status"]] = counts.get(row["status"], 0) + 1
complete = {row["pdb_id"].lower() for row in status_rows}
report = {"requested_pdb":len(IDS),"completed_requested":len(set(IDS)&complete),"pending":len(set(IDS)-complete),"status_counts_all_cache":counts,"source":"EBI SIFTS XML","checkpoint_batch_size":batch_size}
SUMMARY.write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2), flush=True)
