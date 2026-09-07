#!/usr/bin/env python3
import csv, math, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID = ROOT / 'docking_manifest_with_grid.tsv'
RAW = ROOT / 'receptors_pdbqt'
CLEAN = ROOT / 'receptors_pdbqt_clean'
META = ROOT / 'metadata'
AUDIT = ROOT / 'receptor_cleanup_logs' / 'unsupported_atoms_v3.tsv'
MAX_AXIS = 30.0
SUPPORTED_AUTODOCK_TYPES = {'C','A','N','NA','OA','SA','S','P','F','Cl','Br','I','HD','Mg','Mn','Zn','Ca','Fe'}

def atom(line):
    if not line.startswith(('ATOM  ', 'HETATM')):
        return None
    typ = line[77:79].strip() if len(line) >= 79 else ''
    if typ in SUPPORTED_AUTODOCK_TYPES:
        return None
    try:
        return dict(atom_name=line[12:16].strip(), res_name=line[17:20].strip(),
                    chain=line[21:22].strip(), res_num=line[22:26].strip(),
                    x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
                    atom_type=typ or 'BLANK', raw=line.rstrip('\r\n'))
    except ValueError:
        return dict(atom_name='', res_name='', chain='', res_num='', x=None, y=None, z=None,
                    atom_type=typ or 'BLANK', raw=line.rstrip('\r\n'))

def inside(a, row):
    if a['x'] is None:
        return True
    try:
        return all(abs(a[k] - float(row[c])) <= float(row[s]) / 2.0
                   for k, c, s in [('x','grid_center_x','grid_size_x'),
                                   ('y','grid_center_y','grid_size_y'),
                                   ('z','grid_center_z','grid_size_z')])
    except (ValueError, TypeError):
        return True

with GRID.open(encoding='utf-8', newline='') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
    fields = list(rows[0]) if rows else []

by_pdb = {}
for r in rows:
    by_pdb.setdefault((r.get('pdb_id') or '').lower(), []).append(r)

unsupported = {}
for pdb, tasks in by_pdb.items():
    path = RAW / f'{pdb}.pdbqt'
    atoms = []
    if path.exists():
        with path.open(encoding='utf-8', errors='ignore') as f:
            atoms = [a for line in f if (a := atom(line))]
    unsupported[pdb] = atoms

audit_rows = []
qc = {}
for pdb, tasks in by_pdb.items():
    atoms = unsupported[pdb]
    for r in tasks:
        tid = r.get('task_id','')
        grid_ok = (r.get('grid_confidence','').upper() in {'HIGH','MEDIUM'} and
                   all((r.get(k) or '').strip() for k in
                       ['grid_center_x','grid_center_y','grid_center_z','grid_size_x','grid_size_y','grid_size_z']))
        if not grid_ok:
            qc[tid] = ('BLOCKED_GRID_QC', len(atoms), 0)
            continue
        sizes = [float(r[f'grid_size_{a}']) for a in 'xyz']
        if any(v > MAX_AXIS for v in sizes):
            qc[tid] = ('REVIEW_OVERSIZED_BOX', len(atoms), 0)
            continue
        critical = [a for a in atoms if inside(a, r)]
        qc[tid] = ('PREP_FAILED_UNSUPPORTED_CRITICAL_ATOM' if critical else 'PASS',
                   len(atoms), len(critical))
        for a in atoms:
            cx,cy,cz = (float(r['grid_center_x']),float(r['grid_center_y']),float(r['grid_center_z']))
            d = '' if a['x'] is None else f"{math.dist((a['x'],a['y'],a['z']),(cx,cy,cz)):.3f}"
            audit_rows.append(dict(task_id=tid,pdb_id=pdb,chain=a['chain'],res_name=a['res_name'],
                res_num=a['res_num'],atom_name=a['atom_name'],atom_type=a['atom_type'],x=a['x'],y=a['y'],z=a['z'],
                distance_to_box_center=d,inside_grid='YES' if a in critical else 'NO',
                action='BLOCK_TASK' if a in critical else 'REMOVE_FROM_TASK_SAFE_RECEPTOR_COPY'))

CLEAN.mkdir(exist_ok=True)
for old in CLEAN.glob('*.pdbqt'):
    old.unlink()
for pdb, tasks in by_pdb.items():
    src = RAW / f'{pdb}.pdbqt'; dst = CLEAN / f'{pdb}.pdbqt'
    if not src.exists():
        continue
    safe_exists = any(qc.get(r.get('task_id',''),('','',0))[0] == 'PASS' for r in tasks)
    if not safe_exists:
        continue
    with src.open(encoding='utf-8', errors='ignore') as fi, dst.open('w', encoding='utf-8', newline='\n') as fo:
        for line in fi:
            if atom(line) is None:
                fo.write(line.rstrip('\r\n') + '\n')

META.mkdir(exist_ok=True); AUDIT.parent.mkdir(exist_ok=True)
with (META/'task_prep_qc.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['task_id','prep_qc_status','unsupported_atom_count_total','unsupported_atom_count_in_grid'])
    for r in rows:
        w.writerow([r.get('task_id',''),*qc.get(r.get('task_id',''),('MISSING_RECEPTOR',0,0))])
with AUDIT.open('w',encoding='utf-8',newline='') as f:
    names=['task_id','pdb_id','chain','res_name','res_num','atom_name','atom_type','x','y','z','distance_to_box_center','inside_grid','action']
    w=csv.DictWriter(f,fieldnames=names,delimiter='\t');w.writeheader();w.writerows(audit_rows)
from collections import Counter
print(dict(Counter(v[0] for v in qc.values())))
print(f'audit_rows={len(audit_rows)} clean_receptors={len(list(CLEAN.glob("*.pdbqt")))}')
