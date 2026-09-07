import json, os, math, re
from collections import defaultdict, Counter
import multiprocessing as mp

CUTOFF = 4.0  # Angstroms for contact
RECEPTOR_DIR = '/home/csong/docking/docking_package_trial/receptors_pdbqt_clean'
RESULTS_DIR = '/home/csong/docking/docking_package_trial/results_full'
MANIFEST_V4 = '/home/csong/docking/docking_package_trial/docking_manifest_full_with_grid_v4.tsv'
DOCK_LOG = '/home/csong/docking/docking_package_trial/batch_dock_log_full.jsonl'
AUDIT_FILE = '/home/csong/docking/docking_package_trial/grid_audit_report.tsv'

# ---- Load data ----
print('Loading docking log...')
with open(DOCK_LOG) as f:
    ok_set = {}
    for line in f:
        if not line.strip(): continue
        r = json.loads(line)
        if r['status'] == 'ok':
            ok_set[r['task']] = r

print('Loading manifest V4...')
tasks = []
with open(MANIFEST_V4) as f:
    header = f.readline().strip().split('\t')
    idx = {h: i for i, h in enumerate(header)}
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < len(header): continue
        tid = cols[idx['task_id']]
        if tid not in ok_set: continue

        br_str = cols[idx['binding_residues']]
        if not br_str or not br_str.strip(): continue

        # Parse experimental residues: "A:ILE16;A:VAL17;..."
        exp_residues = set()
        for token in br_str.split(';'):
            token = token.strip()
            if not token: continue
            parts = token.split(':')
            if len(parts) >= 2:
                chain = parts[0].strip()
                # Extract residue number (digits after residue name)
                resn_num = parts[1].strip()
                # Separate residue name from number: "ILE16" -> ("ILE", "16")
                m = re.match(r'([A-Z]{1,3})(\d+[A-Za-z]?)', resn_num)
                if m:
                    resname = m.group(1)
                    resnum = m.group(2)
                    exp_residues.add((chain, resnum))

        if not exp_residues: continue

        tasks.append({
            'tid': tid,
            'pdb': cols[idx['pdb_id']].lower(),
            'chain': cols[idx['chain']].strip(),
            'gc': (float(cols[idx['grid_center_x']]),
                   float(cols[idx['grid_center_y']]),
                   float(cols[idx['grid_center_z']])),
            'gs': (float(cols[idx['grid_size_x']]),
                   float(cols[idx['grid_size_y']]),
                   float(cols[idx['grid_size_z']])),
            'exp_residues': exp_residues,
            'energy': ok_set[tid]['best_energy'],
            'evidence': cols[idx['evidence_tier']],
        })

print('OK tasks with binding residues: %d' % len(tasks))

# Load mapping source
with open(AUDIT_FILE) as f:
    h = f.readline().strip().split('\t')
    ti, si = h.index('task_id'), h.index('mapping_source')
    audit = {}
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) > max(ti, si):
            audit[cols[ti]] = cols[si]
for t in tasks:
    t['mapping'] = audit.get(t['tid'], 'unknown')

# ---- Group by PDB+chain for receptor reuse ----
groups = defaultdict(list)
for t in tasks:
    groups[(t['pdb'], t['chain'])].append(t)
print('Unique (PDB, chain) groups: %d' % len(groups))

# ---- Process one (PDB, chain) group ----
def process_group(args):
    (pdb, chain), task_list = args
    receptor_path = os.path.join(RECEPTOR_DIR, '%s.pdbqt' % pdb)
    if not os.path.exists(receptor_path):
        return [{'tid': t['tid'], 'error': 'no_receptor'} for t in task_list]

    # Load receptor atoms for this chain
    rec_atoms = []  # (x, y, z, chain, resname, resnum)
    with open(receptor_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atom_chain = line[21:22]
            if atom_chain != chain: continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resname = line[17:20].strip()
                resnum = line[22:26].strip()
                rec_atoms.append((x, y, z, chain, resname, resnum))
            except ValueError:
                continue

    if not rec_atoms:
        return [{'tid': t['tid'], 'error': 'no_chain_atoms'} for t in task_list]

    results = []
    for t in task_list:
        result_path = os.path.join(RESULTS_DIR, '%s_out.pdbqt' % t['tid'])
        if not os.path.exists(result_path):
            results.append({'tid': t['tid'], 'error': 'no_result'})
            continue

        # Extract MODEL 1 atoms
        lig_atoms = []
        in_m1 = False
        with open(result_path) as f:
            for line in f:
                if line.startswith('MODEL 1'):
                    in_m1 = True
                    continue
                if in_m1 and line.startswith('MODEL'):
                    break
                if in_m1 and line.startswith('ENDMDL'):
                    break
                if in_m1 and (line.startswith('ATOM') or line.startswith('HETATM')):
                    try:
                        lig_atoms.append((
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ))
                    except ValueError:
                        pass

        if not lig_atoms:
            results.append({'tid': t['tid'], 'error': 'no_ligand_atoms'})
            continue

        # Filter receptor atoms: only those within grid box + CUTOFF margin
        gc_x, gc_y, gc_z = t['gc']
        gs_x, gs_y, gs_z = t['gs']
        half_x, half_y, half_z = gs_x/2 + CUTOFF, gs_y/2 + CUTOFF, gs_z/2 + CUTOFF
        nearby_rec = []
        for rx, ry, rz, rch, rname, rnum in rec_atoms:
            dx, dy, dz = rx - gc_x, ry - gc_y, rz - gc_z
            if abs(dx) <= half_x and abs(dy) <= half_y and abs(dz) <= half_z:
                nearby_rec.append((rx, ry, rz, rch, rname, rnum))

        if not nearby_rec:
            results.append({'tid': t['tid'], 'error': 'no_nearby_receptor'})
            continue

        # Compute contacts: receptor atoms within CUTOFF of any ligand atom
        cut2 = CUTOFF * CUTOFF
        pred_residues = set()
        for rx, ry, rz, rch, rname, rnum in nearby_rec:
            for lx, ly, lz in lig_atoms:
                d2 = (rx-lx)**2 + (ry-ly)**2 + (rz-lz)**2
                if d2 <= cut2:
                    pred_residues.add((rch, rnum))
                    break

        exp = t['exp_residues']
        intersection = exp & pred_residues
        recall = len(intersection) / len(exp) if exp else 0
        precision = len(intersection) / len(pred_residues) if pred_residues else 0
        f1 = 2 * recall * precision / (recall + precision) if (recall + precision) > 0 else 0

        results.append({
            'tid': t['tid'],
            'pdb': pdb,
            'chain': chain,
            'mapping': t['mapping'],
            'evidence': t['evidence'],
            'energy': t['energy'],
            'exp_n': len(exp),
            'pred_n': len(pred_residues),
            'intersect_n': len(intersection),
            'recall': recall,
            'precision': precision,
            'f1': f1,
        })

    return results

# ---- Run in parallel ----
print('Processing %d groups with %d workers...' % (len(groups), mp.cpu_count()))
group_items = list(groups.items())

# Process in batches to show progress
batch_size = 50
all_results = []
for i in range(0, len(group_items), batch_size):
    batch = group_items[i:i+batch_size]
    with mp.Pool(32) as pool:
        batch_results = pool.map(process_group, batch)
    for rlist in batch_results:
        all_results.extend(rlist)
    print('  Progress: %d/%d groups, %d results' % (
        min(i+batch_size, len(group_items)), len(group_items), len(all_results)))

# ---- Analysis ----
print()
print('=' * 70)
print('  BINDING RESIDUE RECALL ANALYSIS')
print('  %d tasks processed' % len(all_results))
print('=' * 70)

errors = [r for r in all_results if 'error' in r]
valid = [r for r in all_results if 'error' not in r]
print('Errors: %d, Valid: %d' % (len(errors), len(valid)))

if errors:
    err_types = Counter(r['error'] for r in errors)
    print('Error types: %s' % dict(err_types))

if not valid:
    print('No valid results!')
    exit()

# Overall stats
recalls = [r['recall'] for r in valid]
precisions = [r['precision'] for r in valid]
f1s = [r['f1'] for r in valid]

sr = sorted(recalls)
sp = sorted(precisions)
sf = sorted(f1s)

print()
print('--- Overall (N=%d) ---' % len(valid))
print('  Recall:   Mean=%.3f  Median=%.3f  Q1=%.3f  Q3=%.3f' % (
    sum(recalls)/len(recalls), sr[len(sr)//2], sr[len(sr)//4], sr[3*len(sr)//4]))
print('  Precision:Mean=%.3f  Median=%.3f' % (
    sum(precisions)/len(precisions), sp[len(sp)//2]))
print('  F1:       Mean=%.3f  Median=%.3f' % (
    sum(f1s)/len(f1s), sf[len(sf)//2]))

# Recall distribution
print()
print('--- Recall Distribution ---')
bins = [(0, 0.1), (0.1, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 0.9), (0.9, 1.01)]
for lo, hi in bins:
    c = sum(1 for r in recalls if lo <= r < hi)
    bar = '#' * max(1, c * 40 // max(len(recalls), 1))
    print('  %.0f-%.0f%%: %5d (%5.1f%%) %s' % (lo*100, hi*100, c, 100.0*c/len(recalls), bar))

# By mapping source
print()
print('--- By Mapping Source ---')
for ms in ['direct', 'dbref_segment', 'homologous_chain']:
    subset = [r for r in valid if r['mapping'] == ms]
    if not subset: continue
    rs = [r['recall'] for r in subset]
    ps = [r['precision'] for r in subset]
    print('  %-20s N=%5d  Recall: mean=%.3f med=%.3f  Precision: mean=%.3f med=%.3f' % (
        ms, len(subset), sum(rs)/len(rs), sorted(rs)[len(rs)//2],
        sum(ps)/len(ps), sorted(ps)[len(ps)//2]))

# By evidence tier
print()
print('--- By Evidence Tier ---')
for et in ['BE1', 'BE2']:
    subset = [r for r in valid if r['evidence'] == et]
    if not subset: continue
    rs = [r['recall'] for r in subset]
    ps = [r['precision'] for r in subset]
    print('  %-5s N=%5d  Recall: mean=%.3f med=%.3f  Precision: mean=%.3f med=%.3f' % (
        et, len(subset), sum(rs)/len(rs), sorted(rs)[len(rs)//2],
        sum(ps)/len(ps), sorted(ps)[len(ps)//2]))

# Recall vs Energy
print()
print('--- Recall by Energy Bin ---')
for lo, hi in [(-999, -8), (-8, -7), (-7, -6), (-6, -999)]:
    bin_r = [r for r in valid if lo <= r['energy'] < hi]
    if bin_r:
        rs = [r['recall'] for r in bin_r]
        print('  E %+4.0f to %+4.0f: N=%5d  Recall: mean=%.3f med=%.3f' % (
            lo, hi, len(bin_r), sum(rs)/len(rs), sorted(rs)[len(rs)//2]))

# Top and worst
print()
print('--- Best 20 (highest F1) ---')
for r in sorted(valid, key=lambda x: -x['f1'])[:20]:
    print('  %s  %s/%s  R=%.2f P=%.2f F1=%.2f  exp=%d pred=%d int=%d  map=%s' % (
        r['tid'], r['pdb'], r['chain'],
        r['recall'], r['precision'], r['f1'],
        r['exp_n'], r['pred_n'], r['intersect_n'], r['mapping']))

# Save detailed results
with open('/home/csong/docking/docking_package_trial/residue_recall_results.jsonl', 'w') as f:
    for r in valid:
        f.write(json.dumps(r) + '\n')
print()
print('Full results saved to residue_recall_results.jsonl')
