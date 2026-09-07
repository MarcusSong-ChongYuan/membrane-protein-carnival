#!/usr/bin/env python3
"""
Audit UNK/UNX atoms removed from receptor PDBQT files.

For each removed atom, records:
  - PDB ID, chain, residue name, atom name, coordinates
  - Distance to grid box center (if available in manifest)
  - Removal reason

Output: receptor_cleanup_logs/unk_audit.tsv
"""
import csv
import os
import math
import re
import glob

RECEPTOR_DIR = "receptors_pdbqt"
MANIFEST_WITH_GRID = "docking_manifest_with_grid.tsv"
AUDIT_LOG = "receptor_cleanup_logs/unk_audit.tsv"
POCKET_DISTANCE_THRESHOLD = 15.0  # Angstrom — flag if within this distance

# Load grid box centers from manifest
print("Loading grid box centers...")
grid_centers = {}
with open(MANIFEST_WITH_GRID, encoding='utf-8', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        pdb_id = row.get('pdb_id', '').strip()
        try:
            cx = float(row.get('grid_center_x', 0) or 0)
            cy = float(row.get('grid_center_y', 0) or 0)
            cz = float(row.get('grid_center_z', 0) or 0)
        except ValueError:
            cx = cy = cz = 0.0
        task_id = row.get('task_id', '').strip()
        chain = row.get('chain', '').strip()
        # Key by (pdb_id, chain) for lookup
        key = (pdb_id, chain)
        if key not in grid_centers:
            grid_centers[key] = []
        grid_centers[key].append((task_id, cx, cy, cz))

# Parse PDBQT atom line for coordinates
ATOM_RE = re.compile(
    r'^(ATOM  |HETATM)'
    r'.{6}'             # columns 7-11: serial
    r'(?P<atom_name>.{4})'   # 13-16: atom name
    r'(?P<res_name>.{3})'    # 18-20: residue name
    r'.{1}'             # 22: chain
    r'(?P<chain>.)'
    r'(?P<res_num>.{4})'     # 23-26: residue number
    r'.{4}'             # 27-30
    r'(?P<x>.{8})'      # 31-38: x
    r'(?P<y>.{8})'      # 39-46: y
    r'(?P<z>.{8})'      # 47-54: z
)

def parse_atom_coords(line):
    """Extract coordinates from a PDBQT ATOM line."""
    if len(line) < 54:
        return None
    try:
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        # Extract chain and residue info
        chain = line[21:22].strip() if len(line) >= 22 else ''
        res_name = line[17:20].strip() if len(line) >= 20 else ''
        res_num = line[22:26].strip() if len(line) >= 26 else ''
        atom_name = line[12:16].strip() if len(line) >= 16 else ''
        return {
            'x': x, 'y': y, 'z': z,
            'chain': chain, 'res_name': res_name,
            'res_num': res_num, 'atom_name': atom_name,
        }
    except (ValueError, IndexError):
        return None

def distance(ax, ay, az, bx, by, bz):
    return math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)

# ── Main ──
os.makedirs('receptor_cleanup_logs', exist_ok=True)

print("Auditing UNK/UNX atoms...")
audit_rows = []
total_removed = 0
receptors_affected = 0

for pdbqt_file in sorted(glob.glob(os.path.join(RECEPTOR_DIR, '*.pdbqt'))):
    pdb_id = os.path.basename(pdbqt_file).replace('.pdbqt', '')
    removed_atoms = []

    with open(pdbqt_file, encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if (line.startswith('ATOM  ') or line.startswith('HETATM')):
                atom_type = line[77:79].strip() if len(line) >= 79 else ''
                if not atom_type:
                    coords = parse_atom_coords(line)
                    if coords:
                        removed_atoms.append(coords)

    if not removed_atoms:
        continue

    receptors_affected += 1
    total_removed += len(removed_atoms)

    # Compute distances to grid box centers for this PDB
    pdb_grids = []
    for key, centers in grid_centers.items():
        g_pdb_id, g_chain = key
        if g_pdb_id == pdb_id:
            pdb_grids.extend(centers)

    for atom in removed_atoms:
        min_dist = None
        nearest_task = ''
        for task_id, cx, cy, cz in pdb_grids:
            d = distance(atom['x'], atom['y'], atom['z'], cx, cy, cz)
            if min_dist is None or d < min_dist:
                min_dist = d
                nearest_task = task_id

        audit_rows.append({
            'pdb_id': pdb_id,
            'chain': atom['chain'],
            'res_name': atom['res_name'],
            'res_num': atom['res_num'],
            'atom_name': atom['atom_name'],
            'x': f"{atom['x']:.3f}",
            'y': f"{atom['y']:.3f}",
            'z': f"{atom['z']:.3f}",
            'distance_to_pocket': f"{min_dist:.1f}" if min_dist is not None else '',
            'nearest_task': nearest_task,
            'in_pocket_zone': 'YES' if (min_dist is not None and min_dist < POCKET_DISTANCE_THRESHOLD) else '',
            'removal_reason': 'blank_autodock_type',
        })

# Write audit log
fieldnames = [
    'pdb_id', 'chain', 'res_name', 'res_num', 'atom_name',
    'x', 'y', 'z', 'distance_to_pocket', 'nearest_task',
    'in_pocket_zone', 'removal_reason'
]

with open(AUDIT_LOG, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(audit_rows)

# Summary
pocket_atoms = [r for r in audit_rows if r['in_pocket_zone'] == 'YES']
print(f"\n=== UNK/UNX Audit Summary ===")
print(f"Receptors affected: {receptors_affected}")
print(f"Total atoms removed: {total_removed}")
print(f"Atoms within {POCKET_DISTANCE_THRESHOLD}A of pocket: {len(pocket_atoms)}")
if pocket_atoms:
    print(f"⚠️  These {len(pocket_atoms)} atoms need manual review:")
    for r in pocket_atoms[:20]:
        print(f"  {r['pdb_id']} chain {r['chain']} {r['res_name']}{r['res_num']} "
              f"({r['atom_name']}) at {r['distance_to_pocket']}A from {r['nearest_task']}")
    if len(pocket_atoms) > 20:
        print(f"  ... and {len(pocket_atoms) - 20} more")

print(f"\nFull audit: {AUDIT_LOG}")
