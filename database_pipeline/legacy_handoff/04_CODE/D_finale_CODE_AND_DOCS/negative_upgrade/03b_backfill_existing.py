"""Step 3b: Backfill the ~207 records from 29 CIDs already existing in master."""
import csv, gzip, json, os, shutil, subprocess

OUT_DIR = r'D:\finale\negative_upgrade'
MASTER = r'D:\finale\01_正式数据_V6.2\small_molecule_master_v1_3.tsv'
UNMAPPED = r'D:\finale\01_正式数据_V6.2\negative_binding_evidence_unmapped_review_v1_2.tsv.gz'
NEG_EVIDENCE = r'D:\finale\01_正式数据_V6.2\negative_binding_evidence_v1_2.tsv.gz'
UPGRADED_FILE = os.path.join(OUT_DIR, 'upgraded_negative_evidence.tsv.gz')
NEW_LOOKUP_FILE = os.path.join(OUT_DIR, 'cid_to_internal_ids.json')

# Step 1: Build CID->internal_id from master (for existing CIDs)
print("Building master CID lookup...")
master_lookup = {}
with open(MASTER, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cids_str = row.get('pubchem_cids', '')
        internal = row['compound_internal_id']
        if cids_str:
            for cid in cids_str.split(';'):
                cid = cid.strip()
                if cid and cid.isdigit():
                    master_lookup[cid] = internal
print(f"  {len(master_lookup):,} CIDs from master")

# Step 2: Load newly generated lookup
with open(NEW_LOOKUP_FILE) as f:
    new_lookup = json.load(f)
print(f"  {len(new_lookup):,} CIDs from new query")

# Step 3: Find records for existing CIDs
print("\nScanning unmapped_review for existing-CID records...")
leftover_backfill = []
empty_cid_records = []

with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = reader.fieldnames
    for row in reader:
        cid = row.get('pubchem_cid', '').strip()
        if not cid:
            empty_cid_records.append(row)
        elif cid not in new_lookup and cid in master_lookup:
            row['compound_internal_id'] = master_lookup[cid]
            row['compound_form_id'] = ''
            row['compound_mapping_status'] = 'mapped_core_existing'
            row['release_action_v62'] = 'mapped_negative_release_existing'
            row['unmapped_reason_v62'] = ''
            row['release_version_v62'] = 'v6.2'
            leftover_backfill.append(row)

print(f"  Existing-CID records: {len(leftover_backfill):,}")
print(f"  Empty-CID records: {len(empty_cid_records):,}")

# Step 4: Append to upgraded evidence
print("\nUpdating upgraded evidence...")
existing_upgraded = []
with gzip.open(UPGRADED_FILE, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        existing_upgraded.append(row)
print(f"  Previously upgraded: {len(existing_upgraded):,}")

all_upgraded = existing_upgraded + leftover_backfill
shutil.copy(UPGRADED_FILE, os.path.join(OUT_DIR, 'backup', 'upgraded_negative_evidence_v2.tsv.gz'))
with gzip.open(UPGRADED_FILE, 'wt', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(all_upgraded)
print(f"  Total upgraded now: {len(all_upgraded):,}")

# Step 5: Append to main negative evidence (streaming)
print("\nAppending to main negative evidence (streaming)...")
tmp_neg = os.path.join(OUT_DIR, 'neg_evidence_merged.tsv')

with gzip.open(NEG_EVIDENCE, 'rt', encoding='utf-8') as fin:
    reader = csv.DictReader(fin, delimiter='\t')
    neg_fieldnames = reader.fieldnames
    with open(tmp_neg, 'w', encoding='utf-8', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=neg_fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        cnt = 0
        for row in reader:
            writer.writerow(row)
            cnt += 1
            if cnt % 2000000 == 0:
                print(f"  Copied {cnt:,}...")
        print(f"  Copied {cnt:,} existing negative records")
        for row in leftover_backfill:
            clean_row = {k: row.get(k, '') for k in neg_fieldnames}
            writer.writerow(clean_row)
        print(f"  Appended {len(leftover_backfill):,} existing-CID records")

print("Compressing...")
subprocess.run(['gzip', '-f', tmp_neg], check=True)
shutil.move(tmp_neg + '.gz', NEG_EVIDENCE)
print("Done.")

# Verify
total_neg = 0
with gzip.open(NEG_EVIDENCE, 'rt', encoding='utf-8') as f:
    for _ in f:
        total_neg += 1
print(f"\nFinal negative evidence lines: {total_neg:,}")

# Summary
print(f"\n{'='*60}")
print(f"FINAL SUMMARY")
print(f"{'='*60}")
print(f"Total upgraded: {len(all_upgraded):,}")
print(f"  - Newly queried CIDs (444,466): {len(existing_upgraded):,}")
print(f"  - Existing master CIDs (29):    {len(leftover_backfill):,}")
print(f"Remaining unmapped (empty CID):   {len(empty_cid_records):,}")
print(f"Total negative evidence:          {total_neg:,}")
