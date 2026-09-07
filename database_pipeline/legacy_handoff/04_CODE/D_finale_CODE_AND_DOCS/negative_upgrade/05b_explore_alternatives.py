"""
Explore alternative approaches for the 2,641 unmappable records.
"""
import csv, gzip, subprocess, json, time
from collections import Counter

UNMAPPED = r"D:\finale\01_正式数据_V6.2\negative_binding_evidence_unmapped_review_v1_2.tsv.gz"

# 1. Analyze source_record_id patterns
print("=" * 60)
print("1. Source Record ID Analysis")
print("=" * 60)
src_records = set()
with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        src_id = row.get('source_record_id', '').strip()
        src_records.add(src_id)
print(f"  Unique source_record_ids: {len(src_records)}")
print(f"  Samples:")
for s in list(src_records)[:10]:
    print(f"    '{s}'")

# 2. Check AID patterns and try AID→CID for a few representative AIDs
print("\n" + "=" * 60)
print("2. AID Analysis")
print("=" * 60)
aids = Counter()
aid_sids = {}
with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        aid = row.get('pubchem_aid', '').strip()
        sid = row.get('pubchem_sid', '').strip()
        if aid:
            aids[aid] += 1
            if aid not in aid_sids:
                aid_sids[aid] = set()
            aid_sids[aid].add(sid)

print(f"  Unique AIDs: {len(aids)}")
top_aids = aids.most_common(5)
for aid, cnt in top_aids:
    print(f"  AID {aid}: {cnt} records, {len(aid_sids[aid])} unique SIDs")

# 3. Try AID→SID→CID mapping for top AIDs
print("\n" + "=" * 60)
print("3. AID → SID → CID Test (top AID)")
print("=" * 60)

test_aid = top_aids[0][0]
print(f"  Testing AID {test_aid}...")

# Get all SIDs for this AID, then try SID→CID
r = subprocess.run(["curl", "-s", "--connect-timeout", "15", "--max-time", "30",
    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{test_aid}/sids/JSON"],
    capture_output=True, text=True)
try:
    d = json.loads(r.stdout)
    sids_in_assay = d.get('InformationList', {}).get('Information', [])
    assay_sids = [str(s['SID']) for s in sids_in_assay]
    print(f"  Total SIDs in assay: {len(assay_sids)}")

    # How many of our target SIDs are in this assay?
    target_in_assay = [s for s in aid_sids[test_aid] if s in assay_sids]
    print(f"  Our target SIDs found in assay SID list: {len(target_in_assay)}/{len(aid_sids[test_aid])}")

    # Now try getting CIDs for the SIDs that are in the assay
    if target_in_assay:
        test_sids = target_in_assay[:5]
        for sid in test_sids:
            r2 = subprocess.run(["curl", "-s", "--connect-timeout", "10", "--max-time", "30",
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/JSON"],
                capture_output=True, text=True)
            try:
                d2 = json.loads(r2.stdout)
                if 'IdentifierList' in d2:
                    cids = d2['IdentifierList'].get('CID', [])
                    print(f"    SID {sid} -> CIDs: {cids}")
                elif 'Fault' in d2:
                    pass  # No CID found
            except:
                pass
except Exception as e:
    print(f"  Error querying assay: {e}")
    print(f"  Raw: {r.stdout[:300]}")

# 4. Alternative: try to query PubChem by the exact SID list but in smaller batches
print("\n" + "=" * 60)
print("4. Try individual SID queries (small sample)")
print("=" * 60)

all_sids = []
with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        sid = row.get('pubchem_sid', '').strip()
        if sid:
            all_sids.append(sid)

# Try a random sample of 20 SIDs individually
import random
random.seed(42)
sample = random.sample(all_sids, min(20, len(all_sids)))

found = 0
for sid in sample:
    r = subprocess.run(["curl", "-s", "--connect-timeout", "10", "--max-time", "30",
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/JSON"],
        capture_output=True, text=True)
    try:
        d = json.loads(r.stdout)
        if 'IdentifierList' in d:
            cids = d['IdentifierList'].get('CID', [])
            if cids and cids[0] != 0:
                print(f"  SID {sid} -> CID {cids[0]}")
                found += 1
    except:
        pass

if found == 0:
    print(f"  0/{len(sample)} randomly sampled SIDs mapped to CIDs")
    print(f"  These SIDs genuinely have no CID in PubChem.")

# 5. Try to get substance synonyms — maybe they contain InChIKeys we can use
print("\n" + "=" * 60)
print("5. Substance Synonyms Check")
print("=" * 60)
test_sid = all_sids[0]
r = subprocess.run(["curl", "-s", "--connect-timeout", "10", "--max-time", "30",
    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{test_sid}/synonyms/JSON"],
    capture_output=True, text=True)
try:
    d = json.loads(r.stdout)
    info = d.get('InformationList', {}).get('Information', [])
    if info:
        syns = info[0].get('Synonym', [])
        print(f"  SID {test_sid} synonyms ({len(syns)}):")
        for s in syns[:10]:
            print(f"    {s}")
    else:
        print(f"  No synonyms. Response: {r.stdout[:200]}")
except Exception as e:
    print(f"  Error: {e}")
    print(f"  Raw: {r.stdout[:200]}")

print("\nDone.")
