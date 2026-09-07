import json, os, math, re
from collections import defaultdict, Counter
import multiprocessing as mp

CUTOFF = 4.0
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

def parse_binding_residues(br_str, default_chain):
    """
    Parse binding residues from various formats:
    - 'LEU46;VAL54;ALA67' (no chain prefix, use default_chain)
    - 'A:ILE16;A:VAL17;B:LEU20' (with chain prefix, format chain:3LETTER_RESNUM)
    - 'CID91801202|CHEMBL4297590(NS2)@6tp3:A=S103 P123 Q126' (complex compound@PDB:chain=RES format)
    Returns set of (chain, resnum) tuples
    """
    residues = set()
    br_str = br_str.strip()
    if not br_str: return residues

    # Split by semicolon first
    tokens = br_str.split(';')

    # Check if there are no semicolons but there are spaces (space-separated format)
    if len(tokens) == 1 and ' ' in br_str and ':' not in br_str:
        tokens = br_str.split()

    for token in tokens:
        token = token.strip()
        if not token: continue

        # Strip compound info prefix: "CID91801202|CHEMBL4297590(NS2)@6tp3:A="
        if '@' in token:
            idx = token.rfind('=')
            if idx >= 0:
                token = token[idx+1:]
            else:
                idx = token.rfind(':')
                if idx >= 0:
                    token = token[idx+1:]

        # Handle format: "chain:RESNAME_RESNUM" or just "RESNAME_RESNUM"
        chain = default_chain

        if ':' in token:
            parts = token.split(':')
            chain = parts[0].strip()
            rest = parts[1].strip() if len(parts) > 1 else ''
        else:
            rest = token

        # Parse RESNAME_RESNUM: "ILE16", "16A", "VAL17", "123"
        # Can be: 3-letter + number + optional letter
        m = re.match(r'([A-Za-z]{1,3})?(\d+[A-Za-z]?)', rest)
        if m:
            resnum = m.group(2)
            if resnum:
                residues.add((chain, resnum))
                continue

        # Try just numbers
        m2 = re.match(r'(\d+[A-Za-z]?)', rest)
        if m2:
            residues.add((chain, m2.group(1)))

    return residues


tasks = []
parse_samples = []
with open(MANIFEST_V4) as f:
    header = f.readline().strip().split('\t')
    idx = {h: i for i, h in enumerate(header)}
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < len(header): continue
        tid = cols[idx['task_id']]
        if tid not in ok_set: continue

        br_str = cols[idx['binding_residues']]
        chain = cols[idx['chain']].strip()

        exp_residues = parse_binding_residues(br_str, chain)

        if not exp_residues:
            if len(parse_samples) < 5:
                parse_samples.append((tid, br_str[:200]))
            continue

        tasks.append({
            'tid': tid,
            'pdb': cols[idx['pdb_id']].lower(),
            'chain': chain,
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
if parse_samples:
    print('Parse samples (failed):')
    for tid, bs in parse_samples:
        print('  %s -> %s...' % (tid, bs[:150]))

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
    rec_atoms = []
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
                resnum = line[22:26].strip()
                rec_atoms.append((x, y, z, resnum))
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

        # Filter receptor atoms by bounding box
        gc_x, gc_y, gc_z = t['gc']
        gs_x, gs_y, gs_z = t['gs']
        half_x, half_y, half_z = gs_x/2 + CUTOFF, gs_y/2 + CUTOFF, gs_z/2 + CUTOFF
        nearby_rec = []
        for rx, ry, rz, rnum in rec_atoms:
            dx, dy, dz = rx - gc_x, ry - gc_y, rz - gc_z
            if abs(dx) <= half_x and abs(dy) <= half_y and abs(dz) <= half_z:
                nearby_rec.append((rx, ry, rz, rnum))

        if not nearby_rec:
            results.append({'tid': t['tid'], 'error': 'no_nearby_receptor'})
            continue

        # Compute contacts
        cut2 = CUTOFF * CUTOFF
        pred_residues = set()
        for rx, ry, rz, rnum in nearby_rec:
            for lx, ly, lz in lig_atoms:
                d2 = (rx-lx)**2 + (ry-ly)**2 + (rz-lz)**2
                if d2 <= cut2:
                    pred_residues.add((chain, rnum))
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
print('Processing %d groups with 32 workers...' % len(groups))
group_items = list(groups.items())

batch_size = 50
all_results = []
for i in range(0, len(group_items), batch_size):
    batch = group_items[i:i+batch_size]
    with mp.Pool(32) as pool:
        batch_results = pool.map(process_group, batch)
    for rlist in batch_results:
        all_results.extend(rlist)
    print('  Progress: %d/%d groups' % (
        min(i+batch_size, len(group_items)), len(group_items)))

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
    print('Error types:')
    for k, v in err_types.most_common():
        print('  %s: %d' % (k, v))

if not valid:
    print('No valid results!')
    # Try to debug a few
    for t in tasks[:3]:
        print('Task: %s pdb=%s chain=%s gc=%s gs=%s' % (t['tid'], t['pdb'], t['chain'], t['gc'], t['gs']))
        rp = os.path.join(RECEPTOR_DIR, '%s.pdbqt' % t['pdb'])
        print('  Receptor exists:', os.path.exists(rp))
        xp = os.path.join(RESULTS_DIR, '%s_out.pdbqt' % t['tid'])
        print('  Result exists:', os.path.exists(xp))
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
print('  Recall:   Mean=%.3f  Median=%.3f  Q1=%.3f  Q3=%.3f  Min=%.3f  Max=%.3f' % (
    sum(recalls)/len(recalls), sr[len(sr)//2], sr[len(sr)//4], sr[3*len(sr)//4],
    min(recalls), max(recalls)))
print('  Precision:Mean=%.3f  Median=%.3f' % (
    sum(precisions)/len(precisions), sp[len(sp)//2]))
print('  F1:       Mean=%.3f  Median=%.3f' % (
    sum(f1s)/len(f1s), sf[len(sf)//2]))

print()
print('--- Recall Distribution ---')
bins = [(0, 0.10), (0.10, 0.25), (0.25, 0.50), (0.50, 0.75), (0.75, 0.90), (0.90, 1.01)]
for lo, hi in bins:
    c = sum(1 for r in recalls if lo <= r < hi)
    bar = '#' * max(1, c * 40 // max(len(recalls), 1))
    print('  %3.0f-%3.0f%%: %6d (%5.1f%%) %s' % (lo*100, hi*100, c, 100.0*c/len(recalls), bar))

print()
print('--- By Mapping Source ---')
for ms in ['direct', 'dbref_segment', 'homologous_chain']:
    subset = [r for r in valid if r['mapping'] == ms]
    if not subset: continue
    rs = [r['recall'] for r in subset]
    ps = [r['precision'] for r in subset]
    fs = [r['f1'] for r in subset]
    print('  %-20s N=%6d  Recall: mean=%.3f med=%.3f  Prec: mean=%.3f med=%.3f  F1: mean=%.3f' % (
        ms, len(subset), sum(rs)/len(rs), sorted(rs)[len(rs)//2],
        sum(ps)/len(ps), sorted(ps)[len(ps)//2], sum(fs)/len(fs)))

print()
print('--- By Evidence Tier ---')
for et in ['BE1', 'BE2']:
    subset = [r for r in valid if r['evidence'] == et]
    if not subset: continue
    rs = [r['recall'] for r in subset]
    print('  %-5s N=%6d  Recall: mean=%.3f med=%.3f' % (
        et, len(subset), sum(rs)/len(rs), sorted(rs)[len(rs)//2]))

print()
print('--- Recall vs Energy ---')
for lo, hi in [(-999, -8), (-8, -7), (-7, -6), (-6, 999)]:
    bin_r = [r for r in valid if lo <= r['energy'] < hi]
    if bin_r:
        rs = [r['recall'] for r in bin_r]
        print('  E %+4.0f to %+4.0f: N=%6d  Recall: mean=%.3f med=%.3f' % (
            lo, hi, len(bin_r), sum(rs)/len(rs), sorted(rs)[len(rs)//2]))

print()
print('--- Best 30 (highest F1) ---')
for r in sorted(valid, key=lambda x: -x['f1'])[:30]:
    print('  %s  %s/%s  R=%.3f P=%.3f F1=%.3f  exp=%d pred=%d int=%d  map=%s E=%.1f' % (
        r['tid'], r['pdb'], r['chain'],
        r['recall'], r['precision'], r['f1'],
        r['exp_n'], r['pred_n'], r['intersect_n'], r['mapping'], r['energy']))

print()
print('--- Worst 10 (lowest recall, non-zero pred) ---')
worst = sorted([r for r in valid if r['pred_n'] > 0], key=lambda x: x['recall'])[:10]
for r in worst:
    print('  %s  %s/%s  R=%.3f P=%.3f  exp=%d pred=%d int=%d  map=%s' % (
        r['tid'], r['pdb'], r['chain'],
        r['recall'], r['precision'],
        r['exp_n'], r['pred_n'], r['intersect_n'], r['mapping']))

# Save
with open('/home/csong/docking/docking_package_trial/residue_recall_v2.jsonl', 'w') as f:
    for r in valid:
        f.write(json.dumps(r) + '\n')
print()
print('Saved to residue_recall_v2.jsonl')
