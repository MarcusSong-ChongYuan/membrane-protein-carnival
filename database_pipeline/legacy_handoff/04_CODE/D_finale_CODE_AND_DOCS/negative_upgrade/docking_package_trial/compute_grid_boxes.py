#!/usr/bin/env python3
"""
Compute Vina grid boxes from binding site residues.
Run AFTER downloading PDB structures (01_download_pdbs.sh).

Supports three residue formats from the MemPro database:
  1. "PHE144;PHE151;THR147"           -> chain from manifest chain column
  2. "A:ALA767;A:ALA771;B:THR200"     -> chain embedded per residue
  3. "CID907|(LP5)@4g8a:B=S415 F440" -> chain after @PDBID: prefix
"""
import csv, os, re

MANIFEST = 'docking_manifest.tsv'
PDB_DIR = 'pdb_structures'
MARGIN = 8.0


def parse_residues(desc, default_chain='A'):
    """Parse residue description into (chain, resname, resnum) tuples."""
    results = []

    # Format 3: contains @PDBID:CHAIN prefix
    if '@' in desc:
        after_at = desc.split('@')[-1]
        after_at = re.sub(r'\([^)]*\)', '', after_at).strip()
        chain = default_chain
        if ':' in after_at:
            part = after_at.split(':')[-1]
            if '=' in part:
                eq_parts = part.split('=')
                if len(eq_parts[0].strip()) == 1 and eq_parts[0].strip().isalpha():
                    chain = eq_parts[0].strip().upper()
                    part = '='.join(eq_parts[1:])
            tokens = part.replace(';', ' ').split()
            for tok in tokens:
                m = re.match(r'([A-Za-z]{3})(\d+)', tok)
                if m:
                    results.append((chain, m.group(1).upper(), int(m.group(2))))
        return results

    # Format 2: "A:ALA767;A:ALA771;..." chain prefix per residue
    if ':' in desc:
        parts = desc.split(';')
        has_chain_prefix = all(
            re.match(r'[A-Za-z]:[A-Za-z]{3}\d+', p.strip())
            for p in parts if p.strip()
        )
        if has_chain_prefix:
            for part in parts:
                part = part.strip()
                m = re.match(r'([A-Za-z]):([A-Za-z]{3})(\d+)', part)
                if m:
                    results.append((m.group(1).upper(), m.group(2).upper(), int(m.group(3))))
            return results

    # Format 1: "PHE144;PHE151;THR147"
    for part in desc.replace(' ', ';').split(';'):
        part = part.strip()
        if ':' in part:
            part = part.split(':')[-1]
        m = re.match(r'([A-Za-z]{3})(\d+)', part)
        if m:
            results.append((default_chain, m.group(1).upper(), int(m.group(2))))

    return results


# Main
print(f'Loading manifest: {MANIFEST}')
tasks = []
with open(MANIFEST, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = reader.fieldnames
    for row in reader:
        tasks.append(row)
print(f'  {len(tasks)} tasks')

updated = 0
no_residues = 0
no_pdb_file = 0
no_atom_match = 0

for task in tasks:
    pdb_id = task.get('pdb_id', '').strip()
    residues_desc = task.get('binding_residues', '').strip()
    chain = task.get('chain', 'A').strip()

    if not residues_desc:
        no_residues += 1
        continue

    residues = parse_residues(residues_desc, chain)
    if not residues:
        no_residues += 1
        continue

    pdb_path = os.path.join(PDB_DIR, f'{pdb_id}.pdb')
    if not os.path.exists(pdb_path):
        no_pdb_file += 1
        continue

    coords = []
    matched = set()
    with open(pdb_path) as pf:
        for line in pf:
            if not line.startswith('ATOM') and not line.startswith('HETATM'):
                continue
            if line[12:16].strip() != 'CA':
                continue
            pdb_chain = line[21].strip()
            pdb_resname = line[17:20].strip()
            try:
                pdb_resnum = int(line[22:26])
            except ValueError:
                continue

            for rc, rn, ri in residues:
                if (rc, rn, ri) in matched:
                    continue
                if pdb_chain == rc and pdb_resname == rn and pdb_resnum == ri:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
                        matched.add((rc, rn, ri))
                    except ValueError:
                        pass

    if not coords:
        no_atom_match += 1
        continue

    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    span_x = max(c[0] for c in coords) - min(c[0] for c in coords)
    span_y = max(c[1] for c in coords) - min(c[1] for c in coords)
    span_z = max(c[2] for c in coords) - min(c[2] for c in coords)

    sx = max(span_x + 2 * MARGIN, 18.0)
    sy = max(span_y + 2 * MARGIN, 18.0)
    sz = max(span_z + 2 * MARGIN, 18.0)

    task['grid_center_x'] = f'{cx:.3f}'
    task['grid_center_y'] = f'{cy:.3f}'
    task['grid_center_z'] = f'{cz:.3f}'
    task['grid_size_x'] = f'{sx:.1f}'
    task['grid_size_y'] = f'{sy:.1f}'
    task['grid_size_z'] = f'{sz:.1f}'
    updated += 1

    if updated % 500 == 0:
        print(f'  Computed {updated} grid boxes...')

with open(MANIFEST, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(tasks)

print(f'\nGrid box results:')
print(f'  Computed:             {updated}')
print(f'  No residues parsable: {no_residues}')
print(f'  No PDB file:          {no_pdb_file}')
print(f'  No CA atom match:     {no_atom_match}')
total = updated + no_residues + no_pdb_file + no_atom_match
if total > 0:
    print(f'  Success rate:         {updated}/{total} ({updated*100/total:.1f}%)')

print('\nUpdating Vina config files...')
conf_dir = 'vina_configs'
conf_updated = 0
for task in tasks:
    tid = task['task_id']
    cx = task.get('grid_center_x', '').strip()
    if not cx:
        continue
    conf_path = os.path.join(conf_dir, f'{tid}.conf')
    if not os.path.exists(conf_path):
        continue
    with open(conf_path) as cf:
        lines = cf.readlines()
    with open(conf_path, 'w') as cf:
        for line in lines:
            if line.startswith('center_x ='):
                cf.write(f"center_x = {task['grid_center_x']}\n")
            elif line.startswith('center_y ='):
                cf.write(f"center_y = {task['grid_center_y']}\n")
            elif line.startswith('center_z ='):
                cf.write(f"center_z = {task['grid_center_z']}\n")
            elif line.startswith('size_x ='):
                cf.write(f"size_x = {task['grid_size_x']}\n")
            elif line.startswith('size_y ='):
                cf.write(f"size_y = {task['grid_size_y']}\n")
            elif line.startswith('size_z ='):
                cf.write(f"size_z = {task['grid_size_z']}\n")
            else:
                cf.write(line)
    conf_updated += 1

print(f'  Config files updated: {conf_updated}')
print('Done.')
