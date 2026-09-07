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

# Cofactors/nucleotides/solvents to exclude from drug analysis
COFACTORS = {'ADP','ATP','GTP','GDP','CTP','UTP','TTP','AMP','GMP','CMP','UMP',
             'NAD','NAD+','NADH','NADP','NADPH','FAD','FADH','FMN','FMNH',
             'HEM','HEME','HEA','HEB','HEC','FES','F3S','F4S',
             'SAM','SAH','COA','COA-SH','PLP','TPP','THF','BH4',
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
    """Compute RMSD after optimal rigid-body alignment (Kabsch algorithm)"""
    # Compute centroids
    rc = centroid(ref)
    mc = centroid(mobile)
    # Center both
    rc_ref = [(a[0]-rc[0], a[1]-rc[1], a[2]-rc[2]) for a in ref]
    rc_mob = [(a[0]-mc[0], a[1]-mc[1], a[2]-mc[2]) for a in mobile]
    # Compute covariance matrix
    C = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(len(ref)):
                C[i][j] += rc_ref[k][i] * rc_mob[k][j]
    # SVD via simple approach for 3x3
    # Compute C^T * C eigenvalues (squared singular values)
    # For 3x3, use direct formula
    # Actually, use the Kabsch formula: R = V * U^T where C = U * S * V^T
    # For simplicity, compute rotation using quaternion method
    # Build the F matrix (quaternion-based)
    # Sum over atom pairs
    x = [[0.0]*3 for _ in range(3)]
    for k in range(len(ref)):
        for i in range(3):
            for j in range(3):
                x[i][j] += rc_mob[k][i] * rc_ref[k][j]
    # Build the 4x4 F matrix for quaternion method
    # This is the standard approach
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
    # Find eigenvector of largest eigenvalue (power iteration)
    q = [1.0, 0.0, 0.0, 0.0]
    for _ in range(20):
        nq = [0.0]*4
        for i in range(4):
            for j in range(4):
                nq[i] += F[i][j] * q[j]
        norm = math.sqrt(sum(v*v for v in nq))
        if norm < 1e-10: break
        q = [v/norm for v in nq]
    # Quaternion to rotation matrix
    q0, q1, q2, q3 = q
    R = [
        [q0*q0+q1*q1-q2*q2-q3*q3, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3, 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0*q0-q1*q1-q2*q2+q3*q3],
    ]
    # Apply rotation to mobile atoms
    dists = []
    for a in rc_mob:
        rx = R[0][0]*a[0] + R[0][1]*a[1] + R[0][2]*a[2]
        ry = R[1][0]*a[0] + R[1][1]*a[1] + R[1][2]*a[2]
        rz = R[2][0]*a[0] + R[2][1]*a[1] + R[2][2]*a[2]
        # Find nearest reference atom
        min_d = min((rx-r[0])**2 + (ry-r[1])**2 + (rz-r[2])**2 for r in rc_ref)
        dists.append(min_d)
    return math.sqrt(sum(dists)/len(dists))

def get_pdb_ligands_by_residue(pdb_id, target_codes):
    """Return list of (resn, chain, resnum, atoms, centroid)"""
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
        if len(atoms) < 8: continue  # real ligands have >8 heavy atoms
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

results_all = []
results_drug = []
no_match_pdbs = Counter()
no_pose = 0

for tid in ok_set:
    if audit.get(tid) != 'direct': continue
    if tid not in manifest: continue
    m = manifest[tid]

    # Only check drug-like compounds
    cname = m['name']
    target_codes = set()
    if cname in drug_codes:
        target_codes.add(drug_codes[cname])
    target_codes.add(cname)

    # Find ALL individual ligand copies in PDB
    ref_ligands = get_pdb_ligands_by_residue(m['pdb'], target_codes)
    if not ref_ligands:
        if cname in drug_codes:
            no_match_pdbs[m['pdb']] += 1
        continue

    docked = extract_best_pose(tid)
    if not docked:
        no_pose += 1
        continue
    dock_c = centroid(docked)

    # Match: find the reference copy CLOSEST TO GRID CENTER
    # (not closest to docked pose, to avoid cherry-picking)
    best_ref_dist = float('inf')
    best_ref = None
    for resn, chain, resnum, atoms, ref_c in ref_ligands:
        d = dist3d(m['gc'], ref_c)
        if d < best_ref_dist:
            best_ref_dist = d
            best_ref = (resn, chain, resnum, atoms, ref_c)

    if not best_ref: continue

    cent_dist = dist3d(dock_c, best_ref[4])

    # Compute Kabsch RMSD for same-compound cases
    # (docked pose atoms vs reference atoms)
    rmsd = None
    if len(docked) >= 8 and len(best_ref[3]) >= 8:
        try:
            # Use all atoms from both (may have different counts)
            # Take the smaller count for proper RMSD
            n = min(len(docked), len(best_ref[3]))
            if n >= 8:
                rmsd = kabsch_rmsd(best_ref[3][:n], docked[:n])
        except:
            pass

    entry = {
        'tid': tid, 'pdb': m['pdb'], 'name': cname,
        'ref_resn': best_ref[0], 'ref_chain': best_ref[1],
        'ref_n': len(best_ref[3]), 'dock_n': len(docked),
        'cent_dist': cent_dist, 'rmsd': rmsd,
        'grid2ref': best_ref_dist,
        'energy': ok_set[tid]['best_energy'],
    }
    results_all.append(entry)

    # Drug-only subset: exclude cofactors, nucleotides, lipids
    if cname not in COFACTORS and cname not in {'ADP','ATP','GTP','GDP','CTP','OLA','OLE','OLC'} \
       and not cname.startswith('OLEIC') and not cname.startswith('CHOLESTEROL') \
       and not cname.startswith('SODIUM') and not cname.startswith('ZINC') \
       and not cname.startswith('PUBCHEM') and best_ref[0] not in COFACTORS:
        results_drug.append(entry)

def print_stats(label, data):
    print()
    print('=' * 65)
    print('  %s (N=%d)' % (label, len(data)))
    print('=' * 65)
    if not data:
        print('  (no data)')
        return

    cd = [r['cent_dist'] for r in data]
    scd = sorted(cd)
    energies = [r['energy'] for r in data]
    g2r = [r['grid2ref'] for r in data]

    print('  Centroid Distance:  Mean=%.1f  Median=%.1f  Min=%.1f  Max=%.1f' % (
        sum(cd)/len(cd), scd[len(scd)//2], min(cd), max(cd)))

    bins = [(0,3), (3,5), (5,10), (10,20), (20,999)]
    for lo, hi in bins:
        c = sum(1 for d in cd if lo <= d < hi)
        if c:
            bar = '#' * max(1, c * 40 // len(cd))
            print('    %2d-%2d A: %3d (%5.1f%%) %s' % (lo, hi, c, 100.0*c/len(cd), bar))

    rmsds = [r['rmsd'] for r in data if r['rmsd'] is not None]
    if rmsds:
        sr = sorted(rmsds)
        print('  Kabsch RMSD:        Mean=%.1f  Median=%.1f  Min=%.1f  Max=%.1f' % (
            sum(rmsds)/len(rmsds), sr[len(sr)//2], min(rmsds), max(rmsds)))
        for lo, hi in [(0,2), (2,5), (5,10), (10,20)]:
            c = sum(1 for d in rmsds if lo <= d < hi)
            if c:
                print('    %2d-%2d A: %3d (%5.1f%%)' % (lo, hi, c, 100.0*c/len(rmsds)))

    print('  Grid2Ref:           Mean=%.1f  Median=%.1f  Min=%.1f  Max=%.1f' % (
        sum(g2r)/len(g2r), sorted(g2r)[len(g2r)//2], min(g2r), max(g2r)))
    print('  Energy:             Mean=%.1f  Median=%.1f  Min=%.1f  Max=%.1f' % (
        sum(energies)/len(energies), sorted(energies)[len(energies)//2], min(energies), max(energies)))

    # Top 10
    print()
    print('  --- Best 10 (smallest centroid dist) ---')
    for r in sorted(data, key=lambda x: x['cent_dist'])[:10]:
        rstr = '  N/A' if r['rmsd'] is None else '%5.1fA' % r['rmsd']
        print('  %s  %s  %s  cent=%.1fA  rmsd=%s  E=%.1f  grid2ref=%.1f' % (
            r['tid'], r['pdb'], r['name'][:25], r['cent_dist'], rstr, r['energy'], r['grid2ref']))

    # Per-compound breakdown
    print()
    print('  --- Per-compound breakdown ---')
    comp_stats = defaultdict(list)
    for r in data:
        comp_stats[r['name']].append(r['cent_dist'])
    for name in sorted(comp_stats.keys(), key=lambda n: -len(comp_stats[n])):
        ds = comp_stats[name]
        if len(ds) >= 3:
            print('  %-30s  N=%3d  Mean=%.1f  Med=%.1f  Min=%.1f  Max=%.1f' % (
                name[:30], len(ds), sum(ds)/len(ds), sorted(ds)[len(ds)//2], min(ds), max(ds)))

print_stats('ALL REDOCKING (residue-level grouping)', results_all)
print_stats('DRUG-ONLY REDOCKING (no cofactors/lipids/nucleotides)', results_drug)

if no_match_pdbs:
    print()
    print('--- PDBs where drug code not found in structure (%d tasks, %d unique PDBs) ---' % (
        sum(no_match_pdbs.values()), len(no_match_pdbs)))
    for pdb, c in no_match_pdbs.most_common(15):
        print('  %s: %d tasks' % (pdb, c))
