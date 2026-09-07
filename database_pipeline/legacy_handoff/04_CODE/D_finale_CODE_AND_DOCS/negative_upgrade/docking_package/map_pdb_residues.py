#!/usr/bin/env python3
"""
Map UniProt residue numbers to PDB residue numbers using DBREF records.

PDB files contain DBREF lines like:
DBREF  7R8D B 3252 3252  UNP  P29274  ADA2A_HUMAN   3205   3205

This means:
  PDB chain B, PDB seq positions 3252-3252
  ← UniProt P29274 positions 3205-3205
  Offset = 3205 - 3252 = -47  (UniProt number = PDB number - 47)

Output: pdb_uniprot_offset_map.tsv — lookup table for every PDB/chain/UniProt residue

Usage: python3 map_pdb_residues.py
"""
import csv
import os
import re
import glob

MANIFEST = "docking_manifest.tsv"
PDB_DIR = "pdb_structures"
OUTPUT_MANIFEST = "docking_manifest_with_grid.tsv"  # output from 04
OFFSET_MAP = "pdb_uniprot_offset_map.tsv"

# DBREF line parser
# DBREF  PDB_ID CHAIN_ID  PDB_SEQ_BEGIN PDB_SEQ_BEGIN_INS  PDB_SEQ_END PDB_SEQ_END_INS
#        DATABASE  DB_ACCESSION  DB_ID_CODE  DB_SEQ_BEGIN DB_SEQ_BEGIN_INS  DB_SEQ_END DB_SEQ_END_INS
# Cols:  1-6      8     12-15  15-16          18-21   22      27     33-41      44     48-52    53-54

DBREF_RE = re.compile(
    r'^DBREF\s+'
    r'\S+\s+'           # PDB ID
    r'(?P<chain>\S)'    # chain
    r'\s*'
    r'(?P<pdb_begin>-?\d+)(?P<pdb_begin_ins>[A-Za-z]?)\s*'
    r'(?P<pdb_end>-?\d+)(?P<pdb_end_ins>[A-Za-z]?)\s*'
    r'\S+\s+'            # database
    r'(?P<db_acc>\S+)\s+' # accession
    r'\S+\s+'             # db_id_code
    r'(?P<db_begin>-?\d+)(?P<db_begin_ins>[A-Za-z]?)\s*'
    r'(?P<db_end>-?\d+)(?P<db_end_ins>[A-Za-z]?)'
)

# 1-letter → 3-letter AA code
AA1TO3 = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
}

# ── Step 1: Extract DBREF offsets from all PDB files ──────────

print("=== Step 1: Extracting DBREF offsets ===")
offsets = {}  # (pdb_id, chain) → {uniprot_num: pdb_num, ..., '_offset': delta}

pdb_files = sorted(glob.glob(os.path.join(PDB_DIR, '*.pdb')))
print(f"PDB files: {len(pdb_files)}")

db_count = 0
for pdb_path in pdb_files:
    pdb_id = os.path.basename(pdb_path).replace('.pdb', '')
    with open(pdb_path, encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.startswith('DBREF'):
                continue
            m = DBREF_RE.match(line)
            if not m:
                continue

            chain = m.group('chain').strip().upper()
            try:
                pdb_begin = int(m.group('pdb_begin'))
                pdb_end = int(m.group('pdb_end'))
                db_begin = int(m.group('db_begin'))
                db_end = int(m.group('db_end'))
            except (ValueError, TypeError):
                continue

            # Compute per-residue mapping (linear interpolation)
            key = (pdb_id, chain)
            if key not in offsets:
                offsets[key] = {}

            pdb_range = pdb_end - pdb_begin
            db_range = db_end - db_begin
            if pdb_range == 0:
                # Single residue
                offsets[key][db_begin] = pdb_begin
            else:
                # Check for insertion codes that would break linear mapping
                pdb_begin_ins = (m.group('pdb_begin_ins') or '').strip()
                pdb_end_ins = (m.group('pdb_end_ins') or '').strip()
                db_begin_ins = (m.group('db_begin_ins') or '').strip()
                db_end_ins = (m.group('db_end_ins') or '').strip()

                if pdb_begin_ins or pdb_end_ins or db_begin_ins or db_end_ins:
                    # Has insertion codes — can't simple linear map
                    offsets[key]['_has_icode'] = True

                # Store offset for linear mapping
                if db_range == pdb_range:
                    offsets[key]['_offset'] = db_begin - pdb_begin
                elif db_range != 0:
                    # Proportional mapping
                    offsets[key]['_offset_type'] = 'proportional'
                    offsets[key]['_db_begin'] = db_begin
                    offsets[key]['_db_end'] = db_end
                    offsets[key]['_pdb_begin'] = pdb_begin
                    offsets[key]['_pdb_end'] = pdb_end

            db_count += 1

print(f"DBREF records: {db_count}")
chains_mapped = len(offsets)
print(f"Chains with DBREF: {chains_mapped}")

# ── Step 2: Apply offsets to manifest residues ─────────────────

print("\n=== Step 2: Applying offsets to residue numbers ===")

tasks = []
with open(MANIFEST, encoding='utf-8', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        tasks.append(row)

# Add new column for PDB-mapped residues
# Reuse existing grid output or create new
fieldnames = list(tasks[0].keys())
for col in ['pdb_mapped_residues', 'residue_remap_status', 'remap_details']:
    if col not in fieldnames:
        fieldnames.append(col)

remapped = 0
no_dbref = 0
already_ok = 0
still_failed = 0

MANIFEST_REMAPPED = "docking_manifest_remapped.tsv"
with open(MANIFEST_REMAPPED, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                            extrasaction='ignore')
    writer.writeheader()

    for task in tasks:
        pdb_id = task.get('pdb_id', '').strip()
        chain = task.get('chain', '').strip()
        residue_desc = (
            task.get('binding_residues', '').strip() or
            task.get('site_residues', '').strip() or
            task.get('evidence_residues', '').strip()
        )

        task['pdb_mapped_residues'] = ''
        task['residue_remap_status'] = 'unchanged'
        task['remap_details'] = ''

        if not residue_desc or not pdb_id:
            writer.writerow(task)
            continue

        key = (pdb_id, chain)
        dbref = offsets.get(key)

        if not dbref:
            # Try chainless key
            chainless_key = (pdb_id, '')
            dbref = offsets.get(chainless_key)
            if not dbref:
                # Try any chain match
                for k, v in offsets.items():
                    if k[0] == pdb_id:
                        dbref = v
                        break

        if not dbref:
            no_dbref += 1
            task['residue_remap_status'] = 'no_dbref'
            task['remap_details'] = f'No DBREF for {pdb_id} chain {chain}'
            writer.writerow(task)
            continue

        # Parse residues
        desc = residue_desc
        # Extract after PDBbind '@'
        if '@' in desc:
            desc = desc.split('@')[-1]
        if '=' in desc:
            _, desc = desc.split('=', 1)
        elif ':' in desc:
            parts = desc.split(':', 1)
            if len(parts) > 1:
                desc = parts[-1]

        parts = [p.strip() for p in desc.replace(';', ' ').split() if p.strip()]

        mapped_parts = []
        offsets_used = []
        for part in parts:
            # Try to extract residue number
            # 3-letter: GLN71, GLN71A
            # 1-letter: Q71, Q71A
            m_aa = re.match(r'([A-Za-z]{1,3})(-?\d+)([A-Za-z]?)', part)
            if not m_aa:
                mapped_parts.append(part)
                continue

            aa_code = m_aa.group(1).upper()
            resnum = int(m_aa.group(2))
            icode = m_aa.group(3)

            # Convert 1-letter to 3-letter if needed
            if len(aa_code) == 1:
                aa3 = AA1TO3.get(aa_code, aa_code + 'AA')
            else:
                aa3 = aa_code

            # Try DBREF mapping
            if resnum in dbref:
                pdb_num = dbref[resnum]
                offsets_used.append((resnum, pdb_num))
                if icode:
                    mapped_parts.append(f'{aa3}{pdb_num}{icode}')
                else:
                    mapped_parts.append(f'{aa3}{pdb_num}')
            elif '_offset' in dbref:
                offset = dbref['_offset']
                pdb_num = resnum - offset
                offsets_used.append((resnum, pdb_num))
                if icode:
                    mapped_parts.append(f'{aa3}{pdb_num}{icode}')
                else:
                    mapped_parts.append(f'{aa3}{pdb_num}')
            elif '_offset_type' in dbref and dbref['_offset_type'] == 'proportional':
                # Proportional mapping
                db_range = dbref['_db_end'] - dbref['_db_begin']
                pdb_range = dbref['_pdb_end'] - dbref['_pdb_begin']
                if db_range != 0:
                    fraction = (resnum - dbref['_db_begin']) / db_range
                    pdb_num = int(dbref['_pdb_begin'] + fraction * pdb_range)
                    offsets_used.append((resnum, pdb_num))
                    if icode:
                        mapped_parts.append(f'{aa3}{pdb_num}{icode}')
                    else:
                        mapped_parts.append(f'{aa3}{pdb_num}')
                else:
                    mapped_parts.append(part)
            else:
                mapped_parts.append(part)

        if offsets_used:
            remapped += 1
            task['pdb_mapped_residues'] = ';'.join(mapped_parts)
            task['residue_remap_status'] = 'remapped'
            task['remap_details'] = '; '.join(
                f'UniProt{n}→PDB{p}' for n, p in offsets_used[:10]
            )
            if len(offsets_used) > 10:
                task['remap_details'] += f'... ({len(offsets_used)} total)'
        else:
            task['pdb_mapped_residues'] = residue_desc
            task['residue_remap_status'] = 'unchanged'
            task['remap_details'] = 'No offsets applied'

        writer.writerow(task)

print(f"\n=== Summary ===")
print(f"Total tasks: {len(tasks)}")
print(f"Remapped via DBREF: {remapped}")
print(f"No DBREF available: {no_dbref}")
print(f"Output: {MANIFEST_REMAPPED}")

# ── Step 3: Save offset map for audit ──
with open(OFFSET_MAP, 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['pdb_id', 'chain', 'offset_type', 'offset_value', 'details'])
    for (pdb_id, chain), d in sorted(offsets.items()):
        if '_offset' in d:
            writer.writerow([pdb_id, chain, 'linear', d['_offset'], ''])
        elif '_offset_type' in d:
            writer.writerow([pdb_id, chain, d['_offset_type'], '',
                             f"UniProt{d['_db_begin']}-{d['_db_end']}→PDB{d['_pdb_begin']}-{d['_pdb_end']}"])
        elif '_has_icode' in d:
            writer.writerow([pdb_id, chain, 'per_residue_map', '',
                             f"{len([k for k in d if not k.startswith('_')])} residues mapped"])
        else:
            writer.writerow([pdb_id, chain, 'per_residue_map', '',
                             f"{len([k for k in d if not k.startswith('_')])} residues mapped"])

print(f"Offset map: {OFFSET_MAP}")
