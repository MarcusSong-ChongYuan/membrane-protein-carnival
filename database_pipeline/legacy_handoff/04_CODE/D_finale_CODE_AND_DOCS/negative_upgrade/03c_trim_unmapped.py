"""Trim unmapped_review to only remaining 2,641 empty-CID records."""
import csv, gzip, os, shutil

UNMAPPED = r'D:\finale\01_正式数据_V6.2\negative_binding_evidence_unmapped_review_v1_2.tsv.gz'
OUT_DIR = r'D:\finale\negative_upgrade'

remaining = []
with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = reader.fieldnames
    for row in reader:
        cid = row.get('pubchem_cid', '').strip()
        if not cid:
            remaining.append(row)

print(f"Remaining unmapped (empty CID): {len(remaining):,}")

# Backup full version and update
shutil.copy(UNMAPPED, os.path.join(OUT_DIR, 'backup', 'unmapped_review_v1_2_full.tsv.gz'))
with gzip.open(UNMAPPED, 'wt', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(remaining)

# Verify
with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
    verify = sum(1 for _ in f)
print(f"Verify: {verify:,} lines (expected 2,642)")
print("Done.")
