import json, os, math
from collections import defaultdict, Counter

drug_codes = {
    'GEFITINIB': 'IRE', 'IMATINIB': 'STI', 'ERLOTINIB': 'ERL',
    'DASATINIB': 'DAS', 'LAPATINIB': 'LAP', 'SORAFENIB': 'SOR',
    'SUNITINIB': 'SUN', 'CRIZOTINIB': 'CRZ', 'VEMURAFENIB': 'VEM',
    'DABRAFENIB': 'DAB', 'TRAMETINIB': 'TRA', 'IBRUTINIB': 'IBR',
    'IDELALISIB': 'IDE', 'VENETOCLAX': 'VEN', 'OSIMERTINIB': 'OSI',
    'OLAPARIB': 'OLA', 'PALBOCICLIB': 'PAL', 'RIBOCICLIB': 'RIB',
    'METHOTREXATE': 'MTX', 'IBUPROFEN': 'IBP', 'NAPROXEN': 'NPX',
    'ASPIRIN': 'ASA', 'CAFFEINE': 'CFF', 'THEOPHYLLINE': 'TEP',
    'WARFARIN': 'WAR', 'DIAZEPAM': 'DZP', 'FLUOXETINE': 'FLX',
    'ISOBUTYLMETHYLXANTHINE': 'IBM', 'IBMX': 'IBM',
    'PP2': 'PP2', 'STAUROSPORINE': 'STU',
}

EXCLUDE = {'HOH','WAT','DOD','SO4','PO4','GOL','EDO','PEG','DMS','CIT','ACT','CL','NA','MG','ZN','CA','K'}

# Load docking results
with open('/home/csong/docking/docking_package_trial/batch_dock_log_full.jsonl') as f:
    ok_set = {r['task']: r for r in [json.loads(l) for l in f if l.strip()] if r['status'] == 'ok'}

# Load manifest
with open('/home/csong/docking/docking_package_trial/docking_manifest_full_with_grid_v4.tsv') as f:
    header = f.readline().strip().split('\t')
    tid_i = header.index('task_id')
    name_i = header.index('preferred_name')
    pdb_i = header.index('pdb_id')
    cx_i = header.index('grid_center_x')
    cy_i = header.index('grid_center_y')
    cz_i = header.index('grid_center_z')
    manifest = {}
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) <= max(tid_i, name_i, pdb_i): continue
        try:
            manifest[cols[tid_i]] = {
                'name': cols[name_i].strip().upper(),
                'pdb': cols[pdb_i].lower(),
                'gc': (float(cols[cx_i]), float(cols[cy_i]), float(cols[cz_i])),
            }
        except: pass

# Load mapping source
with open('/home/csong/docking/docking_package_trial/grid_audit_report.tsv') as f:
    h = f.readline().strip().split('\t')
    ti, si = h.index('task_id'), h.index('mapping_source')
    audit = {}
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) > max(ti, si):
            audit[cols[ti]] = cols[si]

def centroid(atoms):
    n = len(atoms)
    return (sum(a[0] for a in atoms)/n, sum(a[1] for a in atoms)/n, sum(a[2] for a in atoms)/n)

def dist3d(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

def get_pdb_ligands(pdb_id, target_codes):
    path = '/home/csong/docking/docking_package_trial/pdb_structures/%s.pdb' % pdb_id
    if not os.path.exists(path): return []
    groups = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('HETATM'):
                resn = line[17:20].strip()
                if resn in EXCLUDE: continue
                chain = line[21:22]
                try:
                    groups[(resn, chain)].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except: pass
    result = []
    for (resn, chain), atoms in groups.items():
        if len(atoms) < 5: continue
        is_target = (resn in target_codes)
        for name, code in drug_codes.items():
            if code == resn:
                is_target = True
                break
        if is_target:
            result.append((resn, chain, atoms, centroid(atoms)))
    return result

def extract_best_pose(task_id):
    path = '/home/csong/docking/docking_package_trial/results_full/%s_out.pdbqt' % task_id
    if not os.path.exists(path): return None
    atoms = []
    in_m1 = False
    with open(path) as f:
        for line in f:
            if line.startswith('MODEL 1'): in_m1 = True; continue
            if in_m1 and line.startswith('MODEL'): break
            if in_m1 and line.startswith('ENDMDL'): break
            if in_m1 and (line.startswith('ATOM') or line.startswith('HETATM')):
                try: atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except: pass
    return atoms if atoms else None

redock = []
for tid in ok_set:
    if audit.get(tid) != 'direct': continue
    if tid not in manifest: continue
    m = manifest[tid]

    target_codes = set()
    if m['name'] in drug_codes:
        target_codes.add(drug_codes[m['name']])
    target_codes.add(m['name'])

    ref_ligands = get_pdb_ligands(m['pdb'], target_codes)
    if not ref_ligands: continue

    docked = extract_best_pose(tid)
    if not docked: continue
    dock_c = centroid(docked)

    best_dist = float('inf')
    best_ref = None
    for resn, chain, atoms, ref_c in ref_ligands:
        d = dist3d(dock_c, ref_c)
        if d < best_dist:
            best_dist = d
            best_ref = (resn, chain, atoms, ref_c)

    redock.append({
        'tid': tid, 'pdb': m['pdb'], 'name': m['name'],
        'ref_resn': best_ref[0], 'ref_chain': best_ref[1],
        'ref_n': len(best_ref[2]), 'dock_n': len(docked),
        'cent_dist': best_dist, 'energy': ok_set[tid]['best_energy'],
        'ref_c': best_ref[3], 'dock_c': dock_c,
        'grid_c': m['gc'],
    })

print('=' * 60)
print('TRUE REDOCKING (chain-aware): same compound, best chain match')
print('=' * 60)
print('Found: %d tasks' % len(redock))

if not redock:
    print('No matches. Checking a few PDBs for any ligand content...')
    for tid in list(ok_set)[:5]:
        if tid in manifest and audit.get(tid) == 'direct':
            pdb = manifest[tid]['pdb']
            path = '/home/csong/docking/docking_package_trial/pdb_structures/%s.pdb' % pdb
            if os.path.exists(path):
                het = set()
                with open(path) as f:
                    for line in f:
                        if line.startswith('HETATM'):
                            r = line[17:20].strip()
                            if r not in ('HOH','WAT','DOD',''): het.add(r)
                print('  %s (%s): HETATM residues: %s' % (tid, pdb, sorted(het)[:20]))
    exit()

dists = [r['cent_dist'] for r in redock]
sd = sorted(dists)
energies = [r['energy'] for r in redock]

print()
print('--- Centroid Distance (after chain matching) ---')
print('  N:      %d' % len(dists))
print('  Mean:   %.1f A' % (sum(dists)/len(dists)))
print('  Median: %.1f A' % sd[len(sd)//2])
print('  Min:    %.1f A' % min(dists))
print('  Max:    %.1f A' % max(dists))
print()

bins = [(0,2), (2,5), (5,10), (10,20), (20,50), (50,999)]
for lo, hi in bins:
    c = sum(1 for d in dists if lo <= d < hi)
    if c:
        bar = '#' * max(1, c * 40 // max(len(dists), 1))
        print('  %2d-%2d A: %3d (%5.1f%%) %s' % (lo, hi, c, 100.0*c/len(dists), bar))

print()
print('--- Energy ---')
print('  Mean: %.1f  Median: %.1f  Min: %.1f  Max: %.1f' % (
    sum(energies)/len(energies), sorted(energies)[len(energies)//2], min(energies), max(energies)))

print()
print('--- Top 10 (smallest centroid dist) ---')
for r in sorted(redock, key=lambda x: x['cent_dist'])[:10]:
    print('  %s  %s %s/%s  cent=%.1fA  E=%.1f' % (
        r['tid'], r['pdb'], r['ref_resn'], r['ref_chain'], r['cent_dist'], r['energy']))

print()
print('--- Worst 5 ---')
for r in sorted(redock, key=lambda x: -x['cent_dist'])[:5]:
    print('  %s  %s %s/%s  cent=%.1fA  E=%.1f  ref_n=%d dock_n=%d' % (
        r['tid'], r['pdb'], r['ref_resn'], r['ref_chain'], r['cent_dist'], r['energy'],
        r['ref_n'], r['dock_n']))

print()
print('--- Grid Center vs Reference Centroid ---')
gc_dists = [dist3d(r['grid_c'], r['ref_c']) for r in redock]
print('  Mean: %.1f A  Median: %.1f A  Min: %.1f A  Max: %.1f A' % (
    sum(gc_dists)/len(gc_dists), sorted(gc_dists)[len(gc_dists)//2], min(gc_dists), max(gc_dists)))

chain_counts = Counter(r['ref_chain'] for r in redock)
print()
print('--- Chains ---')
for ch, c in chain_counts.most_common():
    print('  Chain %s: %d' % (ch, c))

pdb_counts = Counter(r['pdb'] for r in redock)
print()
print('--- Top PDBs ---')
for pdb, c in pdb_counts.most_common(10):
    print('  %s: %d' % (pdb, c))

name_set = set(r['name'] for r in redock)
print()
print('--- Unique compounds: %d ---' % len(name_set))
for n in sorted(name_set):
    c = sum(1 for r in redock if r['name'] == n)
    if c >= 3: print('  %s: %d' % (n, c))
