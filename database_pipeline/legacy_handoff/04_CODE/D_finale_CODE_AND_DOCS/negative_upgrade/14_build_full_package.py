"""
Build COMPLETE docking package:
- ALL compounds (not just bio-status)
- Filter: BE1+BE2 + E1+E2
- Priority scoring
- Rankings from SMILES availability + PDB
"""
import csv, gzip, os, math, shutil
from collections import Counter, defaultdict

EVIDENCE = r'D:\finale\01_正式数据_V6.2\binding_evidence_master_v6_2.tsv.gz'
PROTEIN = r'D:\finale\01_正式数据_V6.2\human_membrane_protein_master_v6_2.tsv'
MASTER = r'D:\finale\01_正式数据_V6.2\small_molecule_master_v1_3.tsv'
FORM = r'D:\finale\01_正式数据_V6.2\compound_form_hierarchy_v1_3.tsv'
DOCK_DIR = r'D:\finale\negative_upgrade\docking_package'
os.makedirs(DOCK_DIR, exist_ok=True)

# ============================================================
# 1. Load protein data
# ============================================================
print('Loading proteins...')
prot_info = {}  # uniprot -> dict
with open(PROTEIN, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        up = row['target_uniprot_id']
        prot_info[up] = {
            'e_tier': row.get('evidence_level_v52', ''),
            'symbol': row.get('approved_symbol', ''),
            'name': row.get('preferred_name', ''),
            'scope': row.get('membrane_scope', ''),
            'has_be1': row.get('has_BE1_structural_binding', '0') == '1',
            'has_be2': row.get('has_BE2_direct_binding', '0') == '1',
        }
print(f'  {len(prot_info):,} membrane proteins')

# ============================================================
# 2. Load compound SMILES (from both master and form tables)
# ============================================================
print('Loading compound SMILES...')
cpd_smiles = {}
cpd_name = {}
cpd_flags = {}

with open(MASTER, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        cid = row['compound_internal_id']
        smi = row.get('standard_smiles', '')
        cpd_name[cid] = row.get('preferred_name', '')[:80]
        if smi:
            cpd_smiles[cid] = smi
        flags = []
        for flag in ['is_approved_drug', 'is_clinical_candidate', 'is_endogenous_ligand',
                      'is_natural_product', 'is_chemical_probe']:
            if row.get(flag, '0') == '1':
                flags.append(flag)
        cpd_flags[cid] = ';'.join(flags)

# Supplement from form table for missing SMILES
missing = [cid for cid in cpd_name if cid not in cpd_smiles or not cpd_smiles[cid]]
if missing:
    with open(FORM, encoding='utf-8') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            cid = row.get('compound_internal_id', '')
            if cid in missing and cid not in cpd_smiles:
                smi = row.get('exact_smiles', '') or row.get('canonical_smiles', '')
                if smi:
                    cpd_smiles[cid] = smi

print(f'  Compounds with SMILES: {len(cpd_smiles):,} / {len(cpd_name):,}')

# ============================================================
# 3. Single pass: build pairs, apply filter, collect scoring
# ============================================================
print('Scanning evidence, building pairs with filter BE1+BE2+E1+E2...')

pair_data = {}  # (cpd, up) -> best row data
prot_pdb = {}  # up -> best PDB
row_count = 0
tier_order = {'BE1': 3, 'BE2': 2, 'BE3': 1}

with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        row_count += 1
        if row_count % 1000000 == 0:
            print(f'  {row_count:,} rows...')

        cpd = row.get('compound_internal_id', '')
        up = row.get('target_uniprot_id', '')
        if up not in prot_info:
            continue
        if cpd not in cpd_smiles:
            continue

        # Get E and BE tiers
        e_tier = prot_info[up]['e_tier']
        be = row.get('evidence_tier', '')

        # Apply filter early
        if be not in ('BE1', 'BE2'):
            continue
        if e_tier not in ('E1', 'E2'):
            continue

        key = (cpd, up)
        be_val = tier_order.get(be, 0)

        if key not in pair_data or be_val > tier_order.get(pair_data[key].get('evidence_tier', ''), 0):
            pair_data[key] = row
        elif be_val == tier_order.get(pair_data[key].get('evidence_tier', ''), 0):
            # Same tier: prefer lower nM
            try:
                cur_nm = float(pair_data[key].get('standard_value_nM', '') or '')
                new_nm = float(row.get('standard_value_nM', '') or '')
                if new_nm > 0 and (cur_nm <= 0 or new_nm < cur_nm):
                    pair_data[key] = row
            except ValueError:
                pass

        # Collect protein-level best PDB
        pdb_str = row.get('pdb_ids', '').strip()
        if pdb_str and up not in prot_pdb:
            for p in pdb_str.split(';'):
                p = p.strip().lower()
                if len(p) == 4 and p[0].isdigit():
                    prot_pdb[up] = p
                    if row.get('binding_site_residues', '').strip():
                        break  # prefer residue-bearing PDB

print(f'  {row_count:,} rows scanned')
print(f'  Filtered pairs: {len(pair_data):,}')

# ============================================================
# 4. Build manifest with scoring
# ============================================================
print('Computing priority scores and building manifest...')

def compute_score(ev_row, prot, cpd_id):
    score = 0.0
    be = ev_row.get('evidence_tier', '')
    if be == 'BE1': score += 35
    elif be == 'BE2': score += 20

    try:
        snm = float(ev_row.get('standard_value_nM', '') or '')
        if snm > 0:
            p_act = 9 - math.log10(snm)
            score += max(0, min(20, p_act * 2))
    except ValueError:
        pass

    flags = cpd_flags.get(cpd_id, '')
    for flag, pts in [('is_approved_drug', 15), ('is_clinical_candidate', 12),
                       ('is_chemical_probe', 10), ('is_endogenous_ligand', 8),
                       ('is_natural_product', 5)]:
        if flag in flags:
            score += pts
            break

    ev_type = ev_row.get('evidence_type', '')
    if 'structure' in ev_type.lower(): score += 5
    elif 'curated' in ev_type.lower(): score += 4
    elif 'quantitative' in ev_type.lower(): score += 3

    # PDB and residue bonus added later when assigning PDB
    return round(score, 1)

tasks = []
stats = Counter()

for (cpd, up), ev in pair_data.items():
    pinfo = prot_info[up]

    # Determine PDB
    direct_pdb_str = ev.get('pdb_ids', '').strip()
    has_res = bool(ev.get('binding_site_residues', '').strip())
    pdb_id = ''
    pdb_source = 'none'

    if direct_pdb_str:
        for p in direct_pdb_str.split(';'):
            p = p.strip().lower()
            if len(p) == 4 and p[0].isdigit():
                pdb_id = p
                pdb_source = 'direct'
                break

    if not pdb_id and up in prot_pdb:
        pdb_id = prot_pdb[up]
        pdb_source = 'fallback'

    score = compute_score(ev, pinfo, cpd)
    if pdb_source == 'direct':
        score += 15
    elif pdb_source == 'fallback':
        score += 8
    if has_res:
        score += 10

    tasks.append({
        'compound_internal_id': cpd,
        'preferred_name': cpd_name.get(cpd, '')[:80],
        'smiles': cpd_smiles[cpd],
        'target_uniprot': up,
        'target_symbol': pinfo['symbol'],
        'target_name': pinfo['name'],
        'pdb_id': pdb_id,
        'chain': 'A',
        'pdb_source': pdb_source,
        'has_binding_residues': str(has_res),
        'evidence_residues': ev.get('binding_site_residues', ''),
        'priority_score': str(round(score, 1)),
        'best_evidence_tier': ev.get('evidence_tier', ''),
        'best_activity_nM': ev.get('standard_value_nM', ''),
        'activity_relation': ev.get('activity_relation', ''),
        'activity_outcome': ev.get('activity_outcome_v60', ''),
        'evidence_type': ev.get('evidence_type', ''),
        'evidence_directness': ev.get('evidence_directness', ''),
        'pubmed_ids': ev.get('pubmed_ids', ''),
        'doi': ev.get('doi', ''),
        'source_database': ev.get('source_database', ''),
        'protein_e_tier': pinfo['e_tier'],
        'bio_flags': cpd_flags.get(cpd, ''),
        'grid_center_x': '', 'grid_center_y': '', 'grid_center_z': '',
        'grid_size_x': '20', 'grid_size_y': '20', 'grid_size_z': '20',
    })
    stats['total'] += 1
    stats[f'BE={ev.get("evidence_tier","")}'] += 1
    stats[f'E={pinfo["e_tier"]}'] += 1
    stats[f'PDB={pdb_source}'] += 1

# Sort by priority score descending
tasks.sort(key=lambda t: float(t['priority_score']), reverse=True)

# Re-assign task IDs
for i, t in enumerate(tasks):
    t['task_id'] = f'DOCK-{i+1:06d}'

print(f'  Total tasks: {len(tasks):,}')
print(f'  BE1: {stats.get("BE=BE1",0):,}  BE2: {stats.get("BE=BE2",0):,}')
print(f'  E1: {stats.get("E=E1",0):,}  E2: {stats.get("E=E2",0):,}')
print(f'  PDB direct: {stats.get("PDB=direct",0):,}  fallback: {stats.get("PDB=fallback",0):,}  none: {stats.get("PDB=none",0):,}')
print(f'  Has residues: {sum(1 for t in tasks if t["has_binding_residues"]=="True"):,}')

# ============================================================
# 5. Write manifest
# ============================================================
print('\nWriting manifest...')
cols = ['task_id', 'priority_score', 'best_evidence_tier', 'best_activity_nM',
        'activity_relation', 'activity_outcome', 'evidence_type', 'evidence_directness',
        'pubmed_ids', 'doi', 'source_database', 'protein_e_tier',
        'compound_internal_id', 'preferred_name', 'smiles', 'bio_flags',
        'target_uniprot', 'target_symbol', 'target_name',
        'pdb_id', 'chain', 'pdb_source', 'has_binding_residues', 'evidence_residues',
        'grid_center_x', 'grid_center_y', 'grid_center_z',
        'grid_size_x', 'grid_size_y', 'grid_size_z']

manifest_path = os.path.join(DOCK_DIR, 'docking_manifest.tsv')
with open(manifest_path, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(tasks)

# ============================================================
# 6. Score distribution
# ============================================================
print('\nScore distribution:')
score_bins = Counter()
for t in tasks:
    s = float(t['priority_score'])
    if s >= 80: score_bins['80-100'] += 1
    elif s >= 60: score_bins['60-80'] += 1
    elif s >= 40: score_bins['40-60'] += 1
    elif s >= 20: score_bins['20-40'] += 1
    else: score_bins['0-20'] += 1

for b in ['80-100', '60-80', '40-60', '20-40', '0-20']:
    cnt = score_bins.get(b, 0)
    bar = '#' * (cnt // 5000)
    print(f'  {b}: {cnt:>7,}  {bar}')

# ============================================================
# 7. Ligands + PDB list
# ============================================================
print('\nWriting ligands.smi...')
unique_cpds = {}
for t in tasks:
    cid = t['compound_internal_id']
    if cid not in unique_cpds:
        unique_cpds[cid] = t

smi_file = os.path.join(DOCK_DIR, 'ligands.smi')
with open(smi_file, 'w', encoding='utf-8') as f:
    for cid in sorted(unique_cpds.keys()):
        t = unique_cpds[cid]
        smi = t['smiles'].replace('\t', ' ').replace('\n', '')
        name = t['preferred_name'].replace('\t', ' ').replace('\n', '')[:80]
        f.write(f"{smi}\t{cid}\t{name}\n")

print(f'  {len(unique_cpds):,} unique compounds')

print('Writing pdb_download_list.txt...')
pdb_counts = Counter()
for t in tasks:
    pdb = t['pdb_id'].strip()
    if pdb:
        pdb_counts[pdb] += 1

pdb_file = os.path.join(DOCK_DIR, 'pdb_download_list.txt')
with open(pdb_file, 'w') as f:
    for pdb in sorted(pdb_counts.keys()):
        f.write(f"{pdb}\t{pdb_counts[pdb]} tasks\n")
print(f'  {len(pdb_counts):,} unique PDBs')

# ============================================================
# 8. Vina configs
# ============================================================
print('\nGenerating Vina configs...')
CONF_DIR = os.path.join(DOCK_DIR, 'vina_configs')
if os.path.exists(CONF_DIR):
    shutil.rmtree(CONF_DIR)
os.makedirs(CONF_DIR)

n_configs = 0
n_skip = 0
for t in tasks:
    pdb = t['pdb_id'].strip()
    if not pdb:
        n_skip += 1
        continue
    if not t['smiles'].strip():
        n_skip += 1
        continue

    tid = t['task_id']
    cid = t['compound_internal_id']
    name = t['preferred_name'][:50]
    sym = t['target_symbol']
    score = t['priority_score']

    conf_path = os.path.join(CONF_DIR, f'{tid}.conf')
    with open(conf_path, 'w', encoding='utf-8') as f:
        f.write(f"# {tid}: {name} -> {sym} ({pdb})  score={score}\n")
        f.write(f"receptor = ../receptors_pdbqt/{pdb}.pdbqt\n")
        f.write(f"ligand = ../ligands_pdbqt/{cid}.pdbqt\n")
        f.write(f"out = ../docking_results/{tid}_out.pdbqt\n")
        f.write(f"log = ../docking_results/{tid}_log.txt\n\n")
        f.write(f"center_x = 0.0\ncenter_y = 0.0\ncenter_z = 0.0\n\n")
        f.write(f"size_x = 20\nsize_y = 20\nsize_z = 20\n\n")
        f.write(f"num_modes = 9\n")
        f.write(f"exhaustiveness = 32\n")
        f.write(f"energy_range = 3\n")
    n_configs += 1

print(f'  {n_configs:,} .conf files ({n_skip:,} skipped)')

# ============================================================
# 9. Shell scripts
# ============================================================
print('\nWriting shell scripts...')

with open(os.path.join(DOCK_DIR, '01_download_pdbs.sh'), 'w') as f:
    f.write('#!/bin/bash\n')
    f.write(f'# Download {len(pdb_counts):,} PDBs from PDBe\n')
    f.write('set -e\nPDB_DIR="pdb_structures"\nmkdir -p "$PDB_DIR"\n')
    f.write('COUNT=0\nTOTAL=$(wc -l < pdb_download_list.txt)\n')
    f.write('while IFS=$\'\\t\' read -r pdb rest; do\n')
    f.write('    if [ ! -f "$PDB_DIR/${pdb}.pdb" ]; then\n')
    f.write('        echo "[$((++COUNT))/$TOTAL] $pdb"\n')
    f.write('        curl -sf --connect-timeout 10 --max-time 60 -o "$PDB_DIR/${pdb}.pdb" \\\n')
    f.write('            "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${pdb}.ent" || \\\n')
    f.write('        curl -sf --connect-timeout 10 --max-time 60 -o "$PDB_DIR/${pdb}.pdb" \\\n')
    f.write('            "https://files.rcsb.org/download/${pdb}.pdb"\n')
    f.write('    else\n        ((COUNT++))\n    fi\n')
    f.write('done < pdb_download_list.txt\n')
    f.write('echo "Downloaded $(ls $PDB_DIR/*.pdb 2>/dev/null | wc -l) PDB files"\n')

with open(os.path.join(DOCK_DIR, '06_submit_slurm_array.sh'), 'w') as f:
    f.write('#!/bin/bash\n')
    f.write(f'#SBATCH --job-name=vina_dock\n')
    f.write(f'#SBATCH --gres=gpu:A100:1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=32G\n')
    f.write(f'#SBATCH --time=7-00:00:00\n#SBATCH --array=1-{n_configs}%1000\n')
    f.write(f'#SBATCH --output=logs/vina_%A_%a.out\n#SBATCH --error=logs/vina_%A_%a.err\n\n')
    f.write(f'set -e\ncd "{DOCK_DIR}"\nmkdir -p logs docking_results\n\n')
    f.write(f'CONF_FILE=$(printf "vina_configs/DOCK-%06d.conf" $SLURM_ARRAY_TASK_ID)\n')
    f.write(f'if [ ! -f "$CONF_FILE" ]; then\n    echo "SKIP: $CONF_FILE not found"\n    exit 0\nfi\n')
    f.write(f'echo "=== $(date) Task $SLURM_ARRAY_TASK_ID ==="\nhead -1 "$CONF_FILE"\n\n')
    f.write(f'vina_gpu --config "$CONF_FILE" --gpu_id 0\n')
    f.write(f'echo "Done: $(date)"\n')

# ============================================================
# 10. Quick validation
# ============================================================
print(f'\n{"="*60}')
print(f'PACKAGE READY')
print(f'{"="*60}')
print(f'Filter: BE1+BE2 + E1+E2 (ALL compounds)')
print(f'')
print(f'  Total tasks:     {len(tasks):>10,}')
print(f'  With PDB:        {n_configs:>10,} ({n_configs/len(tasks)*100:.1f}%)')
print(f'  No PDB:          {n_skip:>10,}')
print(f'  Unique compounds:{len(unique_cpds):>10,}')
print(f'  Unique PDBs:     {len(pdb_counts):>10,}')
print(f'  Unique targets:  {len(set(t["target_uniprot"] for t in tasks)):>10,}')
print(f'')

# Top 10 by priority
print('  Top 10:')
for t in tasks[:10]:
    print(f'    {t["task_id"]} score={t["priority_score"]:>5}  '
          f'{t["target_symbol"]:10s} x {t["preferred_name"][:40]:40s}  '
          f'BE={t["best_evidence_tier"]}  nM={t["best_activity_nM"]}  PDB={t["pdb_source"]}')

# Top targets
tgt = Counter(t['target_symbol'] for t in tasks)
print(f'\n  Top 15 targets:')
for sym, cnt in tgt.most_common(15):
    avg_s = sum(float(t['priority_score']) for t in tasks if t['target_symbol'] == sym) / cnt
    print(f'    {sym:20s}: {cnt:5d} tasks  avg_score={avg_s:.1f}')

# File sizes
print(f'\n  Files:')
for fname in sorted(os.listdir(DOCK_DIR)):
    fpath = os.path.join(DOCK_DIR, fname)
    if os.path.isfile(fpath):
        sz = os.path.getsize(fpath)
        print(f'    {fname:40s} {sz/1024:8.1f} KB')
    elif os.path.isdir(fpath) and fname not in ('__pycache__',):
        n = len(os.listdir(fpath))
        print(f'    {fname:40s} {n:8d} items')

print('\nDone.')
