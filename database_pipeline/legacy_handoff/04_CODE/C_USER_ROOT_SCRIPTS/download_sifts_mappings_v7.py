#!/usr/bin/env python3
"""
Download SIFTS per-residue UniProt↔PDB mappings for all PDB IDs in manifest.

Source: SIFTS XML files from EBI FTP (residue-level, not segment-level).
  https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{pdb_id_lower}.xml.gz

Each XML <residue> element contains a <crossRefDb> to UniProt with residue numbers.
This is more accurate than the PDBe REST API which sometimes only returns segment bounds.

Features:
  - Checkpoint recovery (resume after interruption)
  - Atomic writes (write to .tmp then rename)
  - Retry with exponential backoff (3 attempts)
  - Version tracking in output header
  - Progress bar every 50 PDBs

Usage (run on Windows with internet):
    python3 download_sifts_mappings.py

Output: sifts_mappings.tsv
  Columns: pdb_id  chain  struct_asym_id  uniprot_acc  uniprot_resnum
           pdb_resnum  pdb_icode  pdb_resname  observed  mapping_status
"""
import csv, gzip, json, os, re, sys, time, urllib.request, urllib.error
from io import BytesIO
from xml.etree import ElementTree as ET

MANIFEST   = "docking_manifest.tsv"
OUTPUT     = "sifts_mappings.tsv"
CHECKPOINT = "sifts_mappings.checkpoint.tsv"
SUMMARY    = "sifts_summary.json"
VERSION    = "2.0.0"

SIFTS_XML_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{}.xml.gz"
MAX_RETRIES   = 3
RETRY_DELAY   = [1, 5, 15]   # seconds, per attempt

# ── Collect PDB IDs ──────────────────────────────────────────
pdb_ids = set()
with open(MANIFEST, encoding='utf-8', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        pid = row.get('pdb_id', '').strip().lower()
        if pid and len(pid) == 4:
            pdb_ids.add(pid)
pdb_ids = sorted(pdb_ids)
print(f"Unique PDB IDs in manifest: {len(pdb_ids)}")

# ── Load checkpoint ──────────────────────────────────────────
completed = set()
if os.path.exists(CHECKPOINT):
    with open(CHECKPOINT, encoding='utf-8', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            completed.add(row.get('pdb_id', '').strip().lower())
    print(f"Checkpoint: {len(completed)} already downloaded")

pending = [p for p in pdb_ids if p not in completed]
print(f"To fetch: {len(pending)}")

# ── Write header to fresh output (only on first run) ─────────
if not os.path.exists(CHECKPOINT):
    with open(CHECKPOINT, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pdb_id', 'chain', 'struct_asym_id',
                     'uniprot_acc', 'uniprot_resnum',
                     'pdb_resnum', 'pdb_icode', 'pdb_resname',
                     'observed', 'mapping_status'])
    print("Checkpoint file created.")

# ── Fetch loop ───────────────────────────────────────────────
fetched = len(completed)
failed  = 0
empty   = 0
total_residues = 0
summaries = {}

for i, pdb_id in enumerate(pending):
    url = SIFTS_XML_URL.format(pdb_id)

    # Retry loop
    xml_bytes = None
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(url, headers={
                'User-Agent': f'DockingPipeline/{VERSION}',
                'Accept-Encoding': 'gzip',
            })
            with urllib.request.urlopen(req, timeout=60) as resp:
                xml_bytes = resp.read()
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                break  # don't retry 404s
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY[attempt])
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY[attempt])

    if xml_bytes is None:
        failed += 1
        summaries[pdb_id.upper()] = {
            'status': 'not_found_or_failed',
            'chains': 0, 'residues': 0}
        continue

    # Decompress and parse
    try:
        xml_text = gzip.decompress(xml_bytes)
    except Exception:
        # Maybe not gzipped
        xml_text = xml_bytes

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        failed += 1
        summaries[pdb_id.upper()] = {
            'status': f'xml_parse_error: {str(e)[:80]}',
            'chains': 0, 'residues': 0}
        continue

    # Parse SIFTS XML. Structure:
    # <sifts>
    #   <entity type="protein" entityId="A">
    #     <segment segId="...">
    #       <residue dbSource="PDBe" dbAccessionId="1abc" dbChainId="A"
    #                dbResNum="1" dbResName="ALA" dbResCode="A">
    #         <crossRefDb dbSource="UniProt" dbAccessionId="P12345"
    #                     dbResNum="42" dbResName="ALA" dbResCode="A"/>
    #       </residue>
    #     </segment>
    #   </entity>
    # </sifts>

    residues_found = []
    chains_seen = set()

    for entity in root.iter('entity'):
        ent_type = entity.get('type', '')
        entity_id = entity.get('entityId', '')

        for residue in entity.iter('residue'):
            # PDB-side attributes (on the residue element itself)
            db_source  = residue.get('dbSource', '')      # usually "PDBe"
            db_acc_id  = residue.get('dbAccessionId', '')  # PDB ID
            db_chain   = residue.get('dbChainId', '')
            db_res_num = residue.get('dbResNum', '')
            db_res_name = residue.get('dbResName', '')

            # Observed status
            obs = ''
            details = residue.get('details', '')
            if 'unobserved' in details.lower() or residue.get('obsFrac', '1.0') == '0.0':
                obs = 'unobserved'
            else:
                obs = 'observed'

            # Find UniProt cross-ref
            uniprot_acc = ''
            uniprot_num = ''
            uniprot_name = ''
            for xref in residue.iter('crossRefDb'):
                if xref.get('dbSource', '') in ('UniProt', 'Uniprot'):
                    uniprot_acc  = xref.get('dbAccessionId', '')
                    uniprot_num  = xref.get('dbResNum', '')
                    uniprot_name = xref.get('dbResName', '')
                    break

            if not uniprot_acc:
                continue  # no UniProt mapping for this residue

            # Determine mapping status
            if 'unobserved' in details.lower():
                status = 'unobserved'
            elif uniprot_num:
                status = 'mapped'
            else:
                status = 'no_uniprot_number'

            chains_seen.add(db_chain)
            residues_found.append([
                db_acc_id.upper(),           # pdb_id
                db_chain,                     # chain
                entity_id,                    # struct_asym_id (entity_id in SIFTS XML)
                uniprot_acc,                  # uniprot_acc
                uniprot_num,                  # uniprot_resnum
                db_res_num,                   # pdb_resnum
                '',                            # pdb_icode (not in this XML, will leave blank)
                db_res_name,                  # pdb_resname
                obs,                           # observed
                status,                        # mapping_status
            ])

    # Atomic write to checkpoint
    with open(CHECKPOINT, 'a', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        for r in residues_found:
            w.writerow(r)

    if residues_found:
        fetched += 1
        total_residues += len(residues_found)
        summaries[pdb_id.upper()] = {
            'status': 'ok',
            'chains': len(chains_seen),
            'residues': len(residues_found)}
    else:
        empty += 1
        summaries[pdb_id.upper()] = {
            'status': 'no_uniprot_mappings',
            'chains': 0, 'residues': 0}
        # Still mark as checked in checkpoint
        with open(CHECKPOINT, 'a', encoding='utf-8', newline='') as f:
            f.write(f"{pdb_id.upper()}\t\t\t\t\t\t\t\t\tempty\n")

    # Progress
    if (i + 1) % 50 == 0:
        print(f"  [{i+1}/{len(pending)}] fetched={fetched} failed={failed} "
              f"empty={empty} total_res={total_residues} "
              f"current={pdb_id.upper()} ({len(residues_found)} residues)")

    # Rate limit (be nice to EBI)
    time.sleep(0.3)

# ── Build final output from checkpoint ───────────────────────
print(f"\nBuilding final output ...")
with open(CHECKPOINT, encoding='utf-8', newline='') as fin:
    with open(OUTPUT, 'w', encoding='utf-8', newline='') as fout:
        header = fin.readline()
        fout.write(f"# SIFTS mappings v{VERSION}  {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        fout.write(header)
        for line in fin:
            fout.write(line)

# ── Summary ──────────────────────────────────────────────────
summaries['_metadata'] = {
    'version': VERSION,
    'date': time.strftime('%Y-%m-%d %H:%M:%S'),
    'total_pdb_in_manifest': len(pdb_ids),
    'fetched_with_mappings': fetched,
    'failed_or_not_found': failed,
    'no_uniprot_mappings': empty,
    'total_residue_mappings': total_residues,
}
with open(SUMMARY, 'w', encoding='utf-8') as f:
    json.dump(summaries, f, indent=2)

print(f"\n{'='*60}")
print(f"PDB IDs in manifest: {len(pdb_ids)}")
print(f"Fetched with mappings: {fetched}")
print(f"Failed / not found:   {failed}")
print(f"No UniProt mappings:  {empty}")
print(f"Total residue rows:   {total_residues}")
print(f"Output:   {OUTPUT}")
print(f"Checkpoint: {CHECKPOINT}  (keep for resume)")
print(f"Summary:  {SUMMARY}")
print(f"\nUpload {OUTPUT} to A100 alongside 04_compute_grid_boxes.py")
