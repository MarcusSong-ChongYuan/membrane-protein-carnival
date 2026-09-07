"""
Step 1: Extract unique CIDs from unmapped_review that are NOT in small_molecule_master.
Outputs:
  - cids_to_query.txt: one CID per line (444,467 CIDs)
  - cid_rows_map.json: {cid: [row_indices]} for backfill
  - existing_cids.txt: CIDs already in master (29, for reference)
"""
import csv
import gzip
import json
import os
import sys

UNMAPPED = r"D:\finale\01_正式数据_V6.2\negative_binding_evidence_unmapped_review_v1_2.tsv.gz"
MASTER = r"D:\finale\01_正式数据_V6.2\small_molecule_master_v1_3.tsv"
OUT_DIR = r"D:\finale\negative_upgrade"
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("Step 1a: Load all CIDs from small_molecule_master")
print("=" * 60)

master_cids = set()
with open(MASTER, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for i, row in enumerate(reader):
        cids_str = row.get("pubchem_cids", "")
        if cids_str:
            for cid in cids_str.split(";"):
                cid = cid.strip()
                if cid.isdigit():
                    master_cids.add(cid)
        if (i + 1) % 500000 == 0:
            print(f"  {i+1:,} rows scanned, {len(master_cids):,} unique CIDs found")

print(f"  Total: {len(master_cids):,} unique CIDs in master")

print("\n" + "=" * 60)
print("Step 1b: Extract CIDs from unmapped_review, collect row mappings")
print("=" * 60)

cid_rows = {}  # cid -> list of row dicts (for backfill)
cid_empty = 0
total = 0

with gzip.open(UNMAPPED, "rt", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames
    for row in reader:
        total += 1
        cid = row.get("pubchem_cid", "").strip()
        if not cid:
            cid_empty += 1
            continue
        if cid not in cid_rows:
            cid_rows[cid] = []
        cid_rows[cid].append(row)
        if total % 500000 == 0:
            print(f"  {total:,} records scanned, {len(cid_rows):,} unique CIDs")

print(f"  Total records: {total:,}")
print(f"  Records with empty CID: {cid_empty:,}")
print(f"  Unique CIDs: {len(cid_rows):,}")

# Split into existing (already in master) and new (need to query)
existing = {c for c in cid_rows if c in master_cids}
new = {c for c in cid_rows if c not in master_cids}

print(f"\n  CIDs already in master: {len(existing):,}")
print(f"  CIDs NOT in master (need query): {len(new):,}")

# Save new CIDs list
new_cids_file = os.path.join(OUT_DIR, "cids_to_query.txt")
with open(new_cids_file, "w") as f:
    for cid in sorted(new, key=lambda x: int(x)):
        f.write(cid + "\n")
print(f"\n  Saved {len(new):,} CIDs to: {new_cids_file}")

# Save existing CIDs (for reference)
existing_file = os.path.join(OUT_DIR, "existing_cids.txt")
with open(existing_file, "w") as f:
    for cid in sorted(existing, key=lambda x: int(x)):
        f.write(cid + "\n")
print(f"  Saved {len(existing):,} existing CIDs to: {existing_file}")

# Save CID->rows mapping (only for NEW CIDs) — we'll backfill after query
# Save as JSONL for memory efficiency
map_file = os.path.join(OUT_DIR, "cid_rows_map.jsonl")
# Also save record count per CID
count_file = os.path.join(OUT_DIR, "cid_record_count.tsv")

with open(map_file, "w", encoding="utf-8") as fmap, open(count_file, "w", encoding="utf-8") as fcnt:
    fcnt.write("cid\trecord_count\n")
    for cid in sorted(new, key=lambda x: int(x)):
        rows = cid_rows[cid]
        fcnt.write(f"{cid}\t{len(rows)}\n")
        # Store minimal info for backfill: source_evidence_id -> compound_source_id, compound_name
        for row in rows:
            rec = {
                "source_evidence_id": row["source_evidence_id"],
                "compound_source_id": row.get("compound_source_id", ""),
                "compound_name": row.get("compound_name", ""),
                "target_uniprot_id": row.get("target_uniprot_id", ""),
                "approved_symbol": row.get("approved_symbol", ""),
            }
            fmap.write(json.dumps({"cid": cid, "record": rec}, ensure_ascii=False) + "\n")

print(f"  Saved CID→rows mapping to: {map_file}")
print(f"  Saved CID record counts to: {count_file}")

# Stats
total_records_to_upgrade = sum(len(cid_rows[c]) for c in new)
print(f"\n  Total records that can be upgraded: {total_records_to_upgrade:,}")
print(f"  Total records that remain unmapped (no CID): {cid_empty:,}")

print("\nDone.")
