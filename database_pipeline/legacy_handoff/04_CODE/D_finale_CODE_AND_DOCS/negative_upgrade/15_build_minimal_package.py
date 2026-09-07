"""
Build minimal docking package for A100 trial run.
Only pairs with binding site residue data = highest quality, all grid boxes computable.
"""
import csv, os, shutil
from collections import Counter

DETAIL = r'D:\finale\negative_upgrade\biological_status_binding_sites_detail.tsv'
DOCK_DIR = r'D:\finale\negative_upgrade\docking_package_trial'

if os.path.exists(DOCK_DIR):
    shutil.rmtree(DOCK_DIR)
os.makedirs(DOCK_DIR)

print('='*60)
print('Loading source data...')
print('='*60)

with open(DETAIL, encoding='utf-8') as f:
    all_rows = list(csv.DictReader(f, delimiter='\t'))
print(f'  Total rows: {len(all_rows):,}')

# ============================================================
# Filter: must have BOTH PDB and residues
# ============================================================
print('\nFiltering: PDB + residues required...')
valid_rows = []
for r in all_rows:
    pdb_str = r.get('pdb_ids_site', '').strip()
    residues = r.get('residue_or_site_description', '').strip()
    smiles = r.get('standard_smiles', '').strip()
    cpd = r.get('compound_internal_id', '').strip()
    if not pdb_str or not residues or not smiles or not cpd:
        continue
    pdbs = [p.strip().lower() for p in pdb_str.split(',')
            if len(p.strip()) == 4 and p.strip()[0].isdigit()]
    if not pdbs:
        continue
    valid_rows.append(r)

print(f'  Valid rows: {len(valid_rows):,}')
print(f'  Excluded: {len(all_rows) - len(valid_rows):,}')

# ============================================================
# Build unique docking tasks
# ============================================================
print('\nBuilding docking tasks (compound x PDB x chain)...')
tasks = {}
site_type_counts = Counter()
be_counts = Counter()

for r in valid_rows:
    cpd = r['compound_internal_id'].strip()
    smi = r['standard_smiles'].strip()
    name = (r.get('preferred_name', '') or '')[:80]
    sym = r.get('approved_symbol', '').strip()
    uniprot = r.get('target_uniprot_id', '').strip()
    site_type = r.get('site_type', '').strip()
    residues_desc = r.get('residue_or_site_description', '').strip()
    evidence_tier = r.get('evidence_tier', '').strip()
    bio_flags = r.get('biological_status_flags', '').strip()
    pdb_str = r.get('pdb_ids_site', '').strip().lower()

    pdbs = [p.strip() for p in pdb_str.split(',')
            if len(p.strip()) == 4 and p.strip()[0].isdigit()]

    chains_str = r.get('pdb_chain_ids_v60', '').strip()
    if chains_str:
        chains = [c.strip() for c in chains_str.split(';')]
    else:
        chains = []

    site_type_counts[site_type] += 1
    be_counts[evidence_tier] += 1

    for i, pdb in enumerate(pdbs):
        if i < len(chains) and chains[i]:
            chain = chains[i]
        else:
            chain = 'A'
        key = (cpd, pdb, chain)
        if key not in tasks:
            tasks[key] = {
                'compound_internal_id': cpd,
                'preferred_name': name,
                'smiles': smi,
                'pdb_id': pdb,
                'chain': chain,
                'target_symbol': sym,
                'target_uniprot': uniprot,
                'site_type': site_type,
                'binding_residues': residues_desc,
                'evidence_tier': evidence_tier,
                'bio_flags': bio_flags,
            }

print(f'  Unique tasks: {len(tasks):,}')
print(f'  Unique compounds: {len(set(t["compound_internal_id"] for t in tasks.values())):,}')
print(f'  Unique PDBs: {len(set(t["pdb_id"] for t in tasks.values())):,}')
print(f'  Unique targets: {len(set(t["target_uniprot"] for t in tasks.values())):,}')
print(f'\n  Site types:')
for st, cnt in site_type_counts.most_common():
    print(f'    {st}: {cnt:,}')
print(f'\n  Evidence tiers:')
for be, cnt in be_counts.most_common():
    print(f'    {be}: {cnt:,}')

# Convert to sorted list
task_list = sorted(tasks.values(), key=lambda t: (t['target_symbol'], t['pdb_id'], t['compound_internal_id']))
for i, t in enumerate(task_list):
    t['task_id'] = f'DOCK-{i+1:05d}'

n_tasks = len(task_list)
n_cpds = len(set(t['compound_internal_id'] for t in task_list))
n_pdbs = len(set(t['pdb_id'] for t in task_list))
n_targets = len(set(t['target_uniprot'] for t in task_list))

# ============================================================
# Write manifest
# ============================================================
print('\nWriting docking_manifest.tsv...')
manifest_cols = [
    'task_id', 'compound_internal_id', 'preferred_name', 'smiles',
    'pdb_id', 'chain', 'target_symbol', 'target_uniprot',
    'site_type', 'evidence_tier', 'bio_flags', 'binding_residues',
    'grid_center_x', 'grid_center_y', 'grid_center_z',
    'grid_size_x', 'grid_size_y', 'grid_size_z',
]

manifest_path = os.path.join(DOCK_DIR, 'docking_manifest.tsv')
with open(manifest_path, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=manifest_cols, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    for t in task_list:
        row = dict(t)
        for col in ['grid_center_x','grid_center_y','grid_center_z']:
            row[col] = ''
        for col in ['grid_size_x','grid_size_y','grid_size_z']:
            row[col] = '20'
        writer.writerow(row)

# ============================================================
# Write ligands.smi (unique compounds)
# ============================================================
print('Writing ligands.smi...')
cpd_map = {}
for t in task_list:
    cid = t['compound_internal_id']
    if cid not in cpd_map:
        cpd_map[cid] = t

smi_path = os.path.join(DOCK_DIR, 'ligands.smi')
with open(smi_path, 'w', encoding='utf-8') as f:
    for cid in sorted(cpd_map):
        t = cpd_map[cid]
        f.write(f"{t['smiles']}\t{cid}\t{t['preferred_name']}\n")

# ============================================================
# Write PDB download list
# ============================================================
print('Writing pdb_download_list.txt...')
pdb_task_counts = Counter(t['pdb_id'] for t in task_list)
pdb_list_path = os.path.join(DOCK_DIR, 'pdb_download_list.txt')
with open(pdb_list_path, 'w') as f:
    for pdb_id in sorted(pdb_task_counts):
        f.write(f"{pdb_id}\n")

# ============================================================
# Generate Vina configs
# ============================================================
print(f'Generating {n_tasks:,} Vina config files...')
conf_dir = os.path.join(DOCK_DIR, 'vina_configs')
os.makedirs(conf_dir)

for t in task_list:
    tid = t['task_id']
    pdb = t['pdb_id']
    cid = t['compound_internal_id']
    name = t['preferred_name'][:50]
    sym = t['target_symbol']

    conf_path = os.path.join(conf_dir, f'{tid}.conf')
    with open(conf_path, 'w', encoding='utf-8') as f:
        f.write(f"# {tid}: {name} -> {sym} ({pdb})\n")
        f.write(f"receptor = ../receptors_pdbqt/{pdb}.pdbqt\n")
        f.write(f"ligand = ../ligands_pdbqt/{cid}.pdbqt\n")
        f.write(f"out = ../results/{tid}_out.pdbqt\n")
        f.write(f"log = ../results/{tid}_log.txt\n\n")
        f.write("center_x = 0.0\ncenter_y = 0.0\ncenter_z = 0.0\n\n")
        f.write("size_x = 20.0\nsize_y = 20.0\nsize_z = 20.0\n\n")
        f.write("num_modes = 9\n")
        f.write("exhaustiveness = 32\n")
        f.write("energy_range = 3\n")

# ============================================================
# Final summary
# ============================================================
print(f'\n{"="*60}')
print(f'PACKAGE BUILT: {DOCK_DIR}')
print(f'{"="*60}')
print(f'  Tasks:       {n_tasks:,}')
print(f'  Compounds:   {n_cpds:,}')
print(f'  PDBs:        {n_pdbs:,}')
print(f'  Targets:     {n_targets:,}')
print(f'')

total_size_kb = 0
for root, dirs, files in os.walk(DOCK_DIR):
    for fname in files:
        fpath = os.path.join(root, fname)
        total_size_kb += os.path.getsize(fpath) / 1024
print(f'  Total package size: {total_size_kb/1024:.1f} MB')

est_gb = n_pdbs*2.5/1024 + n_pdbs*4/1024 + n_cpds*0.01/1024 + n_tasks*0.05/1024
print(f'  Estimated disk needed on A100: ~{est_gb:.0f} GB')
print(f'    PDB raw:     ~{n_pdbs*2.5/1024:.0f} GB')
print(f'    Receptors:   ~{n_pdbs*4/1024:.0f} GB')
print(f'    Ligands:     ~{n_cpds*0.01:.0f} MB')
print(f'    Results:     ~{n_tasks*0.05/1024:.0f} GB')
print(f'\nDone.')
