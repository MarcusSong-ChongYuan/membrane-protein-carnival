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

COFACTORS = {'ADP','ATP','GTP','GDP','CTP','UTP','TTP','AMP','GMP','CMP','UMP',
             'NAD','NADH','NADP','NADPH','FAD','FADH','FMN','FMNH',
             'HEM','HEME','HEA','HEB','HEC','FES','F3S','F4S',
             'SAM','SAH','COA','PLP','TPP','THF','BH4',
             'OLE','OLA','OLC','STE','PAM','LDA','MYR','PLM'}
EXCLUDE = {'HOH','WAT','DOD','SO4','PO4','GOL','EDO','PEG','DMS','CIT','ACT',
           'CL','NA','MG','ZN','CA','K','MN','FE','CU','NI','CO','CD','BR','I',
           'ACY','GLC','MAN','GAL','NAG','BMA','FUC','SIA','MNA'}

with open('/home/csong/docking/docking_package_trial/batch_dock_log_full.jsonl') as f:
    ok_set = {r['task']: r for r in [json.loads(l) for l in f if l.strip()] if r['status'] == 'ok'}

with open('/home/csong/docking/docking_package_trial/docking_manifest_full_with_grid_v4.tsv') as f:
    header = f.readline().strip().split('\t')
    idx = {h: i for i, h in enumerate(header)}
    manifest = {}
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < len(header): continue
        try:
            manifest[cols[idx['task_id']]] = {
                'name': cols[idx['preferred_name']].strip().upper(),
                'pdb': cols[idx['pdb_id']].lower(),
                'gc': (float(cols[idx['grid_center_x']]),
                       float(cols[idx['grid_center_y']]),
                       float(cols[idx['grid_center_z']])),
                'gs': (float(cols[idx['grid_size_x']]),
                       float(cols[idx['grid_size_y']]),
                       float(cols[idx['grid_size_z']])),
            }
        except: pass

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

def kabsch_rmsd(ref, mobile):
    rc = centroid(ref)
    mc = centroid(mobile)
    rc_ref = [(a[0]-rc[0], a[1]-rc[1], a[2]-rc[2]) for a in ref]
    rc_mob = [(a[0]-mc[0], a[1]-mc[1], a[2]-mc[2]) for a in mobile]
    # Build x matrix
    x = [[0.0]*3 for _ in range(3)]
    for k in range(len(ref)):
        for i in range(3):
            for j in range(3):
                x[i][j] += rc_mob[k][i] * rc_ref[k][j]
    # Build F matrix
    F = [[0.0]*4 for _ in range(4)]
    F[0][0] = x[0][0] + x[1][1] + x[2][2]
    F[0][1] = x[1][2] - x[2][1]
    F[0][2] = x[2][0] - x[0][2]
    F[0][3] = x[0][1] - x[1][0]
    F[1][0] = x[1][2] - x[2][1]
    F[1][1] = x[0][0] - x[1][1] - x[2][2]
    F[1][2] = x[0][1] + x[1][0]
    F[1][3] = x[0][2] + x[2][0]
    F[2][0] = x[2][0] - x[0][2]
    F[2][1] = x[0][1] + x[1][0]
    F[2][2] = -x[0][0] + x[1][1] - x[2][2]
    F[2][3] = x[1][2] + x[2][1]
    F[3][0] = x[0][1] - x[1][0]
    F[3][1] = x[0][2] + x[2][0]
    F[3][2] = x[1][2] + x[2][1]
    F[3][3] = -x[0][0] - x[1][1] + x[2][2]
    q = [1.0, 0.0, 0.0, 0.0]
    for _ in range(20):
        nq = [0.0]*4
        for i in range(4):
            for j in range(4):
                nq[i] += F[i][j] * q[j]
        norm = math.sqrt(sum(v*v for v in nq))
        if norm < 1e-10: break
        q = [v/norm for v in nq]
    q0, q1, q2, q3 = q
    R = [
        [q0*q0+q1*q1-q2*q2-q3*q3, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3, 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0*q0-q1*q1-q2*q2+q3*q3],
    ]
    dists = []
    for a in rc_mob:
        rx = R[0][0]*a[0] + R[0][1]*a[1] + R[0][2]*a[2]
        ry = R[1][0]*a[0] + R[1][1]*a[1] + R[1][2]*a[2]
        rz = R[2][0]*a[0] + R[2][1]*a[1] + R[2][2]*a[2]
        min_d = min((rx-r[0])**2 + (ry-r[1])**2 + (rz-r[2])**2 for r in rc_ref)
        dists.append(min_d)
    return math.sqrt(sum(dists)/len(dists))

def get_pdb_ligands_by_residue(pdb_id, target_codes):
    path = '/home/csong/docking/docking_package_trial/pdb_structures/%s.pdb' % pdb_id
    if not os.path.exists(path): return []
    groups = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('HETATM'):
                resn = line[17:20].strip()
                if resn in EXCLUDE or resn in COFACTORS: continue
                chain = line[21:22]
                resnum = line[22:26].strip()
                try:
                    groups[(resn, chain, resnum)].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except: pass
    result = []
    for (resn, chain, resnum), atoms in groups.items():
        if len(atoms) < 8: continue
        is_target = False
        if resn in target_codes:
            is_target = True
        for name, code in drug_codes.items():
            if code == resn and name in target_codes:
                is_target = True
                break
        if is_target:
            result.append((resn, chain, resnum, atoms, centroid(atoms)))
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

# Collect drug-only redocking data
results = []
for tid in ok_set:
    if audit.get(tid) != 'direct': continue
    if tid not in manifest: continue
    m = manifest[tid]
    cname = m['name']
    target_codes = set()
    if cname in drug_codes:
        target_codes.add(drug_codes[cname])
    target_codes.add(cname)

    ref_ligands = get_pdb_ligands_by_residue(m['pdb'], target_codes)
    if not ref_ligands: continue

    docked = extract_best_pose(tid)
    if not docked: continue
    dock_c = centroid(docked)

    best_ref_dist = float('inf')
    best_ref = None
    for resn, chain, resnum, atoms, ref_c in ref_ligands:
        d = dist3d(m['gc'], ref_c)
        if d < best_ref_dist:
            best_ref_dist = d
            best_ref = (resn, chain, resnum, atoms, ref_c)

    if not best_ref: continue
    cent_dist = dist3d(dock_c, best_ref[4])
    n = min(len(docked), len(best_ref[3]))
    rmsd = kabsch_rmsd(best_ref[3][:n], docked[:n]) if n >= 8 else None

    results.append({
        'tid': tid, 'pdb': m['pdb'], 'name': cname,
        'ref_resn': best_ref[0], 'ref_chain': best_ref[1],
        'ref_n': len(best_ref[3]), 'dock_n': len(docked),
        'cent_dist': cent_dist, 'rmsd': rmsd,
        'grid2ref': best_ref_dist,
        'energy': ok_set[tid]['best_energy'],
        'grid_vol': m['gs'][0] * m['gs'][1] * m['gs'][2],
    })

print('=' * 70)
print('        DOCKING ACCURACY REPORT')
print('  Redocking same compound into its own crystal pocket')
print('  Evidence level: G1 (direct binding residue mapping)')
print('=' * 70)
print()
print('Total redocking cases found: %d' % len(results))

rmsds = [r['rmsd'] for r in results if r['rmsd'] is not None]
cds = [r['cent_dist'] for r in results]
energies = [r['energy'] for r in results]
g2r = [r['grid2ref'] for r in results]

print()
print('--- RMSD Distribution (Kabsch-aligned) ---')
print('  N=%d  Mean=%.1fA  Median=%.1fA  Std=%.1fA' % (
    len(rmsds), sum(rmsds)/len(rmsds), sorted(rmsds)[len(rmsds)//2],
    math.sqrt(sum((x-sum(rmsds)/len(rmsds))**2 for x in rmsds)/(len(rmsds)-1))))

thresholds = [(0, 1.0, 'Excellent (<1A)'), (1.0, 2.0, 'Good (1-2A)'),
              (2.0, 3.0, 'Moderate (2-3A)'), (3.0, 5.0, 'Poor (3-5A)'),
              (5.0, 99, 'Failed (>5A)')]
for lo, hi, label in thresholds:
    c = sum(1 for d in rmsds if lo <= d < hi)
    bar = '#' * max(1, c * 40 // len(rmsds))
    print('  %-22s: %3d (%5.1f%%) %s' % (label, c, 100.0*c/len(rmsds), bar))

print()
print('--- By Compound (min 2 cases) ---')
comp_data = defaultdict(list)
for r in results:
    if r['rmsd'] is not None:
        comp_data[r['name']].append(r)

for name in sorted(comp_data.keys(), key=lambda n: -len(comp_data[n])):
    entries = comp_data[name]
    if len(entries) < 2: continue
    e_rmsds = [e['rmsd'] for e in entries]
    e_energies = [e['energy'] for e in entries]
    print('  %-30s N=%2d  RMSD: mean=%.1f med=%.1f range=%.1f-%.1f  Energy: mean=%.1f' % (
        name[:30], len(entries), sum(e_rmsds)/len(e_rmsds),
        sorted(e_rmsds)[len(e_rmsds)//2], min(e_rmsds), max(e_rmsds),
        sum(e_energies)/len(e_energies)))

print()
print('--- Best Results (RMSD < 1.5A) ---')
for r in sorted([r for r in results if r['rmsd'] is not None and r['rmsd'] < 1.5], key=lambda x: x['rmsd']):
    print('  %s  %s  %-25s  RMSD=%.1fA  E=%.1f  grid2ref=%.1fA' % (
        r['tid'], r['pdb'], r['name'][:25], r['rmsd'], r['energy'], r['grid2ref']))

print()
print('--- Energy vs RMSD Correlation ---')
if len(rmsds) > 2:
    # Pearson correlation
    rmsd_vals = [r['rmsd'] for r in results if r['rmsd'] is not None]
    e_vals = [r['energy'] for r in results if r['rmsd'] is not None]
    n = len(rmsd_vals)
    mean_r = sum(rmsd_vals)/n
    mean_e = sum(e_vals)/n
    cov = sum((rmsd_vals[i]-mean_r)*(e_vals[i]-mean_e) for i in range(n))
    std_r = math.sqrt(sum((x-mean_r)**2 for x in rmsd_vals))
    std_e = math.sqrt(sum((x-mean_e)**2 for x in e_vals))
    r = cov/(std_r*std_e) if std_r*std_e > 0 else 0
    print('  Pearson r = %.3f  (positive = lower energy = better RMSD)' % r)
    if r > 0.3:
        print('  Moderate correlation: energy is a reasonable quality indicator')
    elif r > 0.1:
        print('  Weak correlation: energy alone is NOT reliable for quality')
    else:
        print('  No meaningful correlation: energy does not predict pose quality')

    # By energy bins
    print()
    print('  RMSD by energy bin:')
    for lo, hi in [(-999, -8), (-8, -7), (-7, -6), (-6, -5), (-5, 999)]:
        in_bin = [(e_vals[i], rmsd_vals[i]) for i in range(n) if lo <= e_vals[i] < hi]
        if in_bin:
            br = [x[1] for x in in_bin]
            print('    E %+4.0f to %+4.0f: N=%2d  RMSD mean=%.1f  med=%.1f' % (
                lo, hi, len(in_bin), sum(br)/len(br), sorted(br)[len(br)//2]))

print()
print('--- Context: How many direct-mapped tasks have co-crystal ligands? ---')
# Count all direct-mapped tasks
direct_tasks = [tid for tid in ok_set if audit.get(tid) == 'direct']
print('  Total direct-mapped + OK tasks: %d' % len(direct_tasks))
# Count how many have ANY ligand in PDB
pdb_lig_cache = {}
have_lig = 0
for tid in direct_tasks[:500]:
    if tid not in manifest: continue
    pdb = manifest[tid]['pdb']
    if pdb not in pdb_lig_cache:
        path = '/home/csong/docking/docking_package_trial/pdb_structures/%s.pdb' % pdb
        if os.path.exists(path):
            ligs = set()
            with open(path) as f:
                for line in f:
                    if line.startswith('HETATM'):
                        r = line[17:20].strip()
                        if r not in EXCLUDE and r not in ('HOH','WAT','DOD',''):
                            ligs.add(r)
            pdb_lig_cache[pdb] = len(ligs) > 0
        else:
            pdb_lig_cache[pdb] = False
    if pdb_lig_cache[pdb]:
        have_lig += 1
print('  Sampled 500 direct tasks: %d (%.1f%%) have co-crystal ligand in PDB' % (
    have_lig, 100.0*have_lig/min(500, len(direct_tasks))))
print()
print('  Note: Only %d cases had the SAME compound docked as the co-crystal' % len(results))
print('  Most tasks dock a DIFFERENT compound than the PDB co-crystal ligand')
print('  True redocking is thus limited to ~%d validated cases' % len(results))
