"""
Resolve extracted identifiers (ChEMBL, CAS) to PubChem CIDs.
Then backfill any newly mappable records.
"""
import csv, gzip, json, subprocess, time, os, shutil
from collections import defaultdict

OUT_DIR = r"D:\finale\negative_upgrade"
MASTER_DIR = r"D:\finale\01_正式数据_V6.2"
UNMAPPED = os.path.join(MASTER_DIR, "negative_binding_evidence_unmapped_review_v1_2.tsv.gz")
EXTRACTED_FILE = os.path.join(OUT_DIR, "sid_extracted_ids.json")
SYNONYMS_FILE = os.path.join(OUT_DIR, "sid_synonyms.jsonl")
NEG_EVIDENCE = os.path.join(MASTER_DIR, "negative_binding_evidence_v1_2.tsv.gz")

# Load extracted identifiers
with open(EXTRACTED_FILE, encoding='utf-8') as f:
    extracted = json.load(f)

# Load SID synonyms raw
sid_synonyms = {}
if os.path.exists(SYNONYMS_FILE):
    with open(SYNONYMS_FILE, encoding='utf-8') as f:
        for line in f:
            entry = json.loads(line)
            sid_synonyms[entry['sid']] = entry['synonyms']

print("=" * 60)
print("Step 5d: Resolve Identifiers to CIDs")
print("=" * 60)

# ====== Approach 1: ChEMBL ID → CID ======
print("\n--- 1. ChEMBL ID → CID ---")
chembl_sids = extracted.get('CHEMBL', {})
print(f"  ChEMBL IDs: {len(chembl_sids)}")
for sid, chembl_ids in chembl_sids.items():
    print(f"  SID {sid}: {chembl_ids}")

# Resolve ChEMBL → CID via PubChem
chembl_to_cid = {}
for sid, chembl_ids in chembl_sids.items():
    for chembl_id in chembl_ids:
        # Query PubChem: name/CHEMBLxxx/cids
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chembl_id}/cids/JSON"
        r = subprocess.run(["curl", "-s", "--connect-timeout", "10", "--max-time", "30",
            "-H", "Accept: application/json", url],
            capture_output=True, text=True, timeout=35)
        try:
            d = json.loads(r.stdout)
            cids = d.get('IdentifierList', {}).get('CID', [])
            if cids:
                chembl_to_cid[sid] = str(cids[0])
                print(f"    {chembl_id} -> CID {cids[0]}")
        except:
            pass
        time.sleep(0.2)  # Be gentle

# ====== Approach 2: CAS Number → CID ======
print("\n--- 2. CAS Number → CID ---")
cas_sids = extracted.get('CAS', {})
print(f"  CAS numbers: {len(cas_sids)}")
for sid, cas_nums in cas_sids.items():
    for cas in cas_nums[:2]:
        print(f"  SID {sid}: CAS {cas}")

cas_to_cid = {}
for sid, cas_nums in cas_sids.items():
    for cas in cas_nums:
        # URL-encode the CAS number
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}/cids/JSON"
        r = subprocess.run(["curl", "-s", "--connect-timeout", "10", "--max-time", "30",
            "-H", "Accept: application/json", url],
            capture_output=True, text=True, timeout=35)
        try:
            d = json.loads(r.stdout)
            cids = d.get('IdentifierList', {}).get('CID', [])
            if cids:
                cas_to_cid[sid] = str(cids[0])
                print(f"    CAS {cas} -> CID {cids[0]}")
                break  # First successful mapping
        except:
            pass
        time.sleep(0.2)

# ====== Approach 3: Check if synonyms contain actual InChIKeys ======
print("\n--- 3. InChIKey/SMILES in Synonyms ---")
import re
inchikey_pat = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
smiles_like_pat = re.compile(r'^[A-Za-z0-9@+\-\[\]()\\/=#.%]{10,}$')  # Longer alphanumeric

synth_to_cid = {}
for sid, syns in sid_synonyms.items():
    if sid in chembl_to_cid or sid in cas_to_cid:
        continue  # Already mapped

    for syn in syns:
        if inchikey_pat.match(syn):
            # Look up by InChIKey
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{syn}/cids/JSON"
            r = subprocess.run(["curl", "-s", "--connect-timeout", "10", "--max-time", "30",
                "-H", "Accept: application/json", url],
                capture_output=True, text=True, timeout=35)
            try:
                d = json.loads(r.stdout)
                cids = d.get('IdentifierList', {}).get('CID', [])
                if cids:
                    synth_to_cid[sid] = str(cids[0])
                    print(f"  SID {sid}: InChIKey {syn} -> CID {cids[0]}")
                    break
            except:
                pass
            time.sleep(0.2)
            break  # Only try first InChIKey

# ====== Approach 4: Try some synonyms as compound names ======
print("\n--- 4. Try synonyms as PubChem compound names ---")
# Focus on synonyms that look like names (not IDs)
name_pat = re.compile(r'^[A-Za-z][A-Za-z0-9\s\-_,;:()\[\]+\']{3,}$')
chemical_like = re.compile(r'^\d', re.IGNORECASE)  # Number-like, probably not a name

name_to_cid = {}
for sid, syns in sid_synonyms.items():
    if sid in chembl_to_cid or sid in cas_to_cid or sid in synth_to_cid:
        continue

    for syn in syns[:5]:  # Try first 5 synonyms
        # Skip IDs
        if re.match(r'^(NCGC|NSC|SMR|NCI|AIDS|SCHEMBL|AKOS|MFCD|AK Scientific|MolPort|ZINC|BDBM|CS|HY|LS|EU|FT|AK|AO|DTXSID|CID|SID|AID)', syn):
            continue
        if syn.startswith('UNII-') or syn.startswith('CAS-'):
            continue
        if inchikey_pat.match(syn):
            continue

        # If it looks like a chemical name (contains spaces, letters)
        if ' ' in syn and len(syn) > 10 and len(syn) < 200:
            # URL-encode the name for curl
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{syn}/cids/JSON"
            r = subprocess.run(["curl", "-s", "--connect-timeout", "10", "--max-time", "30",
                "--data-urlencode", f"name={syn}",
                "-G", "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/cids/JSON"],
                capture_output=True, text=True, timeout=35)
            try:
                d = json.loads(r.stdout)
                cids = d.get('IdentifierList', {}).get('CID', [])
                if cids:
                    name_to_cid[sid] = str(cids[0])
                    print(f"  SID {sid}: '{syn[:80]}' -> CID {cids[0]}")
                    break
            except:
                pass
            time.sleep(0.3)
            break  # One attempt per SID for now

# ====== Combine all mappings ======
print("\n" + "=" * 60)
print("COMBINED RESULTS")
print("=" * 60)

all_mappings = {}
all_mappings.update(chembl_to_cid)
all_mappings.update(cas_to_cid)
all_mappings.update(synth_to_cid)
all_mappings.update(name_to_cid)

print(f"  ChEMBL → CID: {len(chembl_to_cid)}")
print(f"  CAS → CID:    {len(cas_to_cid)}")
print(f"  InChIKey → CID: {len(synth_to_cid)}")
print(f"  Name → CID:   {len(name_to_cid)}")
print(f"  Total new mappings: {len(all_mappings)}")

# Save mappings
with open(os.path.join(OUT_DIR, 'sid_to_cid_via_synonyms.json'), 'w') as f:
    json.dump(all_mappings, f, ensure_ascii=False, indent=2)

# ====== Backfill with new mappings ======
if all_mappings:
    print(f"\n--- Backfilling with {len(all_mappings)} new SID→CID mappings ---")

    # Load CID→internal_id from master
    print("  Loading master CID lookup...")
    master_lookup = {}
    with open(os.path.join(MASTER_DIR, "small_molecule_master_v1_3.tsv"), encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cids_str = row.get('pubchem_cids', '')
            internal = row['compound_internal_id']
            if cids_str:
                for cid in cids_str.split(';'):
                    cid = cid.strip()
                    if cid and cid.isdigit():
                        master_lookup[cid] = internal
    print(f"  Loaded {len(master_lookup):,} CIDs from master")

    # Backfill unmapped records
    newly_mapped = []
    still_unmapped = []
    with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        for row in reader:
            sid = row.get('pubchem_sid', '').strip()
            if sid in all_mappings:
                cid = all_mappings[sid]
                row['pubchem_cid'] = cid
                row['compound_source_id'] = f"PubChem CID:{cid}"
                if cid in master_lookup:
                    row['compound_internal_id'] = master_lookup[cid]
                    row['compound_mapping_status'] = 'mapped_via_synonym'
                else:
                    row['compound_mapping_status'] = 'cid_found_not_in_master'
                row['release_action_v62'] = 'synonym_resolved'
                row['unmapped_reason_v62'] = ''
                newly_mapped.append(row)
            else:
                still_unmapped.append(row)

    print(f"  Newly mapped: {len(newly_mapped):,}")
    print(f"  Still unmapped: {len(still_unmapped):,}")

    # How many have compound_internal_id?
    with_id = sum(1 for r in newly_mapped if r.get('compound_internal_id', '').strip())
    without_id = len(newly_mapped) - with_id
    print(f"    With internal ID (in master): {with_id}")
    print(f"    Without internal ID (need master entry): {without_id}")

    # Update unmapped_review
    shutil.copy(UNMAPPED, os.path.join(OUT_DIR, 'backup', 'unmapped_review_v1_2_pre_synonym.tsv.gz'))
    with gzip.open(UNMAPPED, 'wt', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(still_unmapped)

    with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
        remaining_count = sum(1 for _ in f) - 1
    print(f"  Updated unmapped_review: {remaining_count:,} records remaining")

    # Append newly mapped (with internal ID) to negative evidence
    with_id_records = [r for r in newly_mapped if r.get('compound_internal_id', '').strip()]
    if with_id_records:
        print(f"\n  Appending {len(with_id_records):,} records to negative evidence...")
        tmp_neg = os.path.join(OUT_DIR, 'neg_evidence_final.tsv')
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
                print(f"    Copied {cnt:,} existing records")
                for row in with_id_records:
                    clean_row = {k: row.get(k, '') for k in neg_fieldnames}
                    writer.writerow(clean_row)
                print(f"    Appended {len(with_id_records):,} synonym-resolved records")

        subprocess.run(['gzip', '-f', tmp_neg], check=True)
        shutil.move(tmp_neg + '.gz', NEG_EVIDENCE)

        with gzip.open(NEG_EVIDENCE, 'rt', encoding='utf-8') as f:
            final_count = sum(1 for _ in f) - 1
        print(f"    Final negative evidence: {final_count:,} records")

    # Final summary
    print(f"\n{'='*60}")
    print(f"FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"  SID→CID via synonyms: {len(all_mappings):,}")
    print(f"  Records newly mapped: {len(newly_mapped):,}")
    print(f"  Records added to negative evidence: {len(with_id_records):,}")
    print(f"  Records still unmapped: {len(still_unmapped):,}")
else:
    print("\n  No new mappings found. All 2,641 records remain unmapped.")

print("\nDone.")
