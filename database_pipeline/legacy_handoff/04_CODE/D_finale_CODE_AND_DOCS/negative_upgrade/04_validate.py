"""Final validation of the negative evidence upgrade."""
import csv, gzip, json, os

OUT_DIR = r'D:\finale\negative_upgrade'
MASTER_DIR = r'D:\finale\01_正式数据_V6.2'

print("=" * 60)
print("FINAL VALIDATION REPORT")
print("=" * 60)

# 1. Master table counts
print("\n--- 1. Master Table Row Counts ---")
for name, path in [
    ("small_molecule_master", os.path.join(MASTER_DIR, "small_molecule_master_v1_3.tsv")),
    ("compound_form_hierarchy", os.path.join(MASTER_DIR, "compound_form_hierarchy_v1_3.tsv")),
]:
    with open(path, encoding='utf-8') as f:
        data_rows = sum(1 for _ in f) - 1
    print(f"  {name}: {data_rows:,} data rows")

# 2. Negative evidence counts
print("\n--- 2. Negative Evidence Counts ---")
neg_file = os.path.join(MASTER_DIR, "negative_binding_evidence_v1_2.tsv.gz")
unmapped_file = os.path.join(MASTER_DIR, "negative_binding_evidence_unmapped_review_v1_2.tsv.gz")

with gzip.open(neg_file, 'rt', encoding='utf-8') as f:
    neg_total = sum(1 for _ in f) - 1
print(f"  Formal negative evidence: {neg_total:,}")

with gzip.open(unmapped_file, 'rt', encoding='utf-8') as f:
    unmapped_total = sum(1 for _ in f) - 1
print(f"  Unmapped review: {unmapped_total:,}")

# 3. FK integrity: compound_internal_id in negative evidence -> master
print("\n--- 3. FK Integrity Check ---")
# Load all compound_internal_ids from master
master_ids = set()
with open(os.path.join(MASTER_DIR, "small_molecule_master_v1_3.tsv"), encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        master_ids.add(row['compound_internal_id'])
print(f"  Master compound IDs: {len(master_ids):,}")

# Check negative evidence
missing_ids = set()
neg_with_id = 0
neg_without_id = 0
with gzip.open(neg_file, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cid_val = row.get('compound_internal_id', '').strip()
        if cid_val:
            neg_with_id += 1
            if cid_val not in master_ids:
                missing_ids.add(cid_val)
        else:
            neg_without_id += 1

print(f"  Negative records with compound_internal_id: {neg_with_id:,}")
print(f"  Negative records without compound_internal_id: {neg_without_id:,}")
if missing_ids:
    print(f"  ** FK VIOLATIONS: {len(missing_ids)} IDs not in master **")
    for mid in list(missing_ids)[:5]:
        print(f"    {mid}")
else:
    print(f"  ** FK integrity: PASS (0 violations) **")

# 4. Check unmapped review
print("\n--- 4. Unmapped Review Check ---")
empty_cid = 0
has_cid_but_no_id = 0
with gzip.open(unmapped_file, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cid = row.get('pubchem_cid', '').strip()
        internal = row.get('compound_internal_id', '').strip()
        if not cid:
            empty_cid += 1
        elif cid and not internal:
            has_cid_but_no_id += 1

print(f"  Empty CID: {empty_cid:,}")
print(f"  Has CID but no internal ID: {has_cid_but_no_id:,}")
if has_cid_but_no_id == 0:
    print(f"  ** All remaining unmapped are empty-CID only (PASS) **")

# 5. Compound ID ranges
print("\n--- 5. New Compound ID Ranges ---")
# Check last new compound
with open(os.path.join(MASTER_DIR, "small_molecule_master_v1_3.tsv"), encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        pass  # iterate to last
    last_cpd = row['compound_internal_id']
    last_row_src = row.get('source_databases', '')
print(f"  Last compound_internal_id: {last_cpd}")

# Count newly added compounds
new_cpds = sum(1 for cid in master_ids if cid.startswith('HMPD-CMPD-') and int(cid.split('-')[-1]) >= 2016065)
print(f"  Newly added compounds (>=2016065): {new_cpds:,}")

# 6. Upgrade stats
print("\n--- 6. Upgrade Stats ---")
# Read upgrade report
with open(os.path.join(OUT_DIR, 'upgrade_report.json')) as f:
    report = json.load(f)
print(f"  New compounds created: {report['new_compounds']:,}")
print(f"  New forms created: {report['new_forms']:,}")
print(f"  Records upgraded: {report['records_upgraded']:,}")
print(f"  CIDs still unmapped: {report['cids_still_unmapped']}")

# Count all upgraded (new + existing)
upgraded_file = os.path.join(OUT_DIR, 'upgraded_negative_evidence.tsv.gz')
with gzip.open(upgraded_file, 'rt', encoding='utf-8') as f:
    upgraded_total = sum(1 for _ in f) - 1
print(f"  Total upgraded (incl existing-CID backfill): {upgraded_total:,}")

print("\n" + "=" * 60)
print("VALIDATION COMPLETE")
print("=" * 60)
