"""
Build final docking package:
1. Filter manifest: BE1+BE2 + E1+E2
2. Regenerate Vina configs
3. Update ligand list, PDB list, shell scripts
"""
import csv, os, glob
from collections import Counter

DOCK_DIR = r'D:\finale\negative_upgrade\docking_package'
MANIFEST_RANKED = os.path.join(DOCK_DIR, 'docking_manifest_ranked.tsv')

print('Loading ranked manifest...')
all_tasks = []
with open(MANIFEST_RANKED, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    cols = reader.fieldnames
    for row in reader:
        all_tasks.append(row)
print(f'  Total: {len(all_tasks):,}')

# ============================================================
# FILTER: BE1+BE2 + E1+E2
# ============================================================
print('\nFiltering: BE1+BE2 + E1+E2...')
filtered = []
excluded = []
for t in all_tasks:
    be = t.get('best_evidence_tier', '')
    e = t.get('protein_e_tier', '')
    if be in ('BE1', 'BE2') and e in ('E1', 'E2'):
        filtered.append(t)
    else:
        excluded.append(t)

print(f'  Included: {len(filtered):,}')
print(f'  Excluded: {len(excluded):,}')

# Excluded breakdown
print(f'\n  Excluded by:')
by_be = Counter(f'BE={t.get("best_evidence_tier","?")}' for t in excluded)
by_e = Counter(f'E={t.get("protein_e_tier","?")}' for t in excluded)
for k, c in by_be.most_common():
    print(f'    {k}: {c:,}')
for k, c in by_e.most_common():
    print(f'    {k}: {c:,}')

# ============================================================
# Re-assign task IDs
# ============================================================
print(f'\nRe-assigning task IDs...')
for i, t in enumerate(filtered):
    t['task_id'] = f'DOCK-{i+1:06d}'

n_total = len(filtered)

# ============================================================
# Write filtered manifest
# ============================================================
manifest_out = os.path.join(DOCK_DIR, 'docking_manifest.tsv')
print(f'\nWriting final manifest: {manifest_out}')
with open(manifest_out, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(filtered)

# ============================================================
# Score distribution of filtered set
# ============================================================
score_bins = Counter()
for t in filtered:
    s = float(t.get('priority_score', 0))
    if s >= 80: score_bins['80-100'] += 1
    elif s >= 60: score_bins['60-80'] += 1
    elif s >= 40: score_bins['40-60'] += 1
    elif s >= 20: score_bins['20-40'] += 1
    else: score_bins['0-20'] += 1

print(f'\nFiltered set score distribution:')
for b in ['80-100', '60-80', '40-60', '20-40', '0-20']:
    cnt = score_bins.get(b, 0)
    bar = '#' * (cnt // 2000)
    print(f'  {b}: {cnt:>7,}  {bar}')

# ============================================================
# Regenerate ligands.smi
# ============================================================
print(f'\nWriting ligands.smi...')
unique_cpds = {}
for t in filtered:
    cid = t['compound_internal_id']
    if cid not in unique_cpds:
        unique_cpds[cid] = t

smi_file = os.path.join(DOCK_DIR, 'ligands.smi')
with open(smi_file, 'w', encoding='utf-8') as f:
    for cid in sorted(unique_cpds.keys()):
        t = unique_cpds[cid]
        f.write(f"{t['smiles']}\t{cid}\t{t['preferred_name']}\n")
print(f'  {len(unique_cpds):,} unique compounds')

# ============================================================
# Regenerate PDB download list
# ============================================================
print(f'Writing pdb_download_list.txt...')
pdb_cpds = Counter()
pdb_res = Counter()
for t in filtered:
    pdb = t.get('pdb_id', '').strip()
    if pdb:
        pdb_cpds[pdb] += 1
        if t.get('has_binding_residues', '') == 'True':
            pdb_res[pdb] += 1

pdb_list_file = os.path.join(DOCK_DIR, 'pdb_download_list.txt')
with open(pdb_list_file, 'w') as f:
    for pdb in sorted(pdb_cpds.keys()):
        has_r = '*' if pdb in pdb_res else ' '
        f.write(f"{pdb}\t{pdb_cpds[pdb]} tasks\t{has_r}\n")
print(f'  {len(pdb_cpds):,} unique PDBs ({"*"}=has residues)')

# ============================================================
# Regenerate Vina configs
# ============================================================
print(f'\nRegenerating Vina configs...')
CONF_DIR = os.path.join(DOCK_DIR, 'vina_configs')
# Clean old
import shutil
if os.path.exists(CONF_DIR):
    shutil.rmtree(CONF_DIR)
os.makedirs(CONF_DIR)

n_configs = 0
n_skip_no_pdb = 0
for t in filtered:
    pdb = t.get('pdb_id', '').strip()
    if not pdb:
        n_skip_no_pdb += 1
        continue

    task_id = t['task_id']
    cpd_id = t['compound_internal_id']

    cx = t.get('grid_center_x', '').strip() or '0.0'
    cy = t.get('grid_center_y', '').strip() or '0.0'
    cz = t.get('grid_center_z', '').strip() or '0.0'
    sx = t.get('grid_size_x', '').strip() or '20'
    sy = t.get('grid_size_y', '').strip() or '20'
    sz = t.get('grid_size_z', '').strip() or '20'

    name = t.get('preferred_name', '')[:50]
    symbol = t.get('target_symbol', '')

    conf_path = os.path.join(CONF_DIR, f'{task_id}.conf')
    with open(conf_path, 'w') as f:
        f.write(f"# {task_id}: {name} -> {symbol} ({pdb})  score={t['priority_score']}\n")
        f.write(f"receptor = ../receptors_pdbqt/{pdb}.pdbqt\n")
        f.write(f"ligand = ../ligands_pdbqt/{cpd_id}.pdbqt\n")
        f.write(f"out = ../docking_results/{task_id}_out.pdbqt\n")
        f.write(f"log = ../docking_results/{task_id}_log.txt\n\n")
        f.write(f"center_x = {cx}\n")
        f.write(f"center_y = {cy}\n")
        f.write(f"center_z = {cz}\n\n")
        f.write(f"size_x = {sx}\n")
        f.write(f"size_y = {sy}\n")
        f.write(f"size_z = {sz}\n\n")
        f.write("num_modes = 9\n")
        f.write("exhaustiveness = 32\n")
        f.write("energy_range = 3\n")
    n_configs += 1

print(f'  Generated {n_configs:,} .conf files')
print(f'  Skipped (no PDB): {n_skip_no_pdb:,}')

# ============================================================
# Update shell scripts
# ============================================================
print(f'\nUpdating shell scripts...')

# 01_download_pdbs.sh
with open(os.path.join(DOCK_DIR, '01_download_pdbs.sh'), 'w') as f:
    f.write('#!/bin/bash\n')
    f.write(f'# Download {len(pdb_cpds):,} PDB structures from PDBe/RCSB\n')
    f.write('set -e\n\n')
    f.write('PDB_DIR="pdb_structures"\n')
    f.write('mkdir -p "$PDB_DIR"\n\n')
    f.write('COUNT=0\n')
    f.write('TOTAL=$(wc -l < pdb_download_list.txt)\n')
    f.write('echo "Downloading $TOTAL PDBs..."\n')
    f.write('while IFS=$\'\\t\' read -r pdb rest; do\n')
    f.write('    if [ ! -f "$PDB_DIR/${pdb}.pdb" ]; then\n')
    f.write('        echo "[$((++COUNT))/$TOTAL] $pdb"\n')
    f.write('        curl -sf --connect-timeout 10 --max-time 60 \\\n')
    f.write('            -o "$PDB_DIR/${pdb}.pdb" \\\n')
    f.write('            "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${pdb}.ent" || \\\n')
    f.write('        curl -sf --connect-timeout 10 --max-time 60 \\\n')
    f.write('            -o "$PDB_DIR/${pdb}.pdb" \\\n')
    f.write('            "https://files.rcsb.org/download/${pdb}.pdb"\n')
    f.write('    else\n')
    f.write('        ((COUNT++))\n')
    f.write('    fi\n')
    f.write('done < pdb_download_list.txt\n\n')
    f.write('echo "Downloaded $(ls $PDB_DIR/*.pdb 2>/dev/null | wc -l) PDB files"\n')

# 02_prepare_proteins.sh
with open(os.path.join(DOCK_DIR, '02_prepare_proteins.sh'), 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# Prepare protein PDBQT files using ADFR suite\n')
    f.write('set -e\n\n')
    f.write('PDB_DIR="pdb_structures"\n')
    f.write('RECEPTOR_DIR="receptors_pdbqt"\n')
    f.write('mkdir -p "$RECEPTOR_DIR"\n\n')
    f.write('for pdb_file in "$PDB_DIR"/*.pdb; do\n')
    f.write('    pdb_id=$(basename "$pdb_file" .pdb)\n')
    f.write('    if [ ! -f "$RECEPTOR_DIR/${pdb_id}.pdbqt" ]; then\n')
    f.write('        echo "Preparing $pdb_id..."\n')
    f.write('        prepare_receptor -r "$pdb_file" -o "$RECEPTOR_DIR/${pdb_id}.pdbqt" \\\n')
    f.write('            -A hydrogens -U nphs\n')
    f.write('    fi\n')
    f.write('done\n')
    f.write('echo "Prepared $(ls $RECEPTOR_DIR/*.pdbqt 2>/dev/null | wc -l) receptors"\n')

# 03_prepare_ligands.sh
with open(os.path.join(DOCK_DIR, '03_prepare_ligands.sh'), 'w') as f:
    f.write('#!/bin/bash\n')
    f.write(f'# Prepare {len(unique_cpds):,} ligand PDBQT files\n')
    f.write('set -e\n\n')
    f.write('LIGAND_DIR="ligands_pdbqt"\n')
    f.write('mkdir -p "$LIGAND_DIR"\n\n')
    f.write('echo "Converting SMILES to PDBQT..."\n')
    f.write('COUNT=0\n')
    f.write('while IFS=$\'\\t\' read -r smiles cpd_id name; do\n')
    f.write('    out="$LIGAND_DIR/${cpd_id}.pdbqt"\n')
    f.write('    if [ ! -f "$out" ]; then\n')
    f.write('        echo "$smiles" | obabel -ismi -opdbqt -O "$out" --gen3d --conformers --nconf 1 2>/dev/null\n')
    f.write('        ((COUNT++))\n')
    f.write('        if [ $((COUNT % 500)) -eq 0 ]; then\n')
    f.write('            echo "  $COUNT ligands prepared..."\n')
    f.write('        fi\n')
    f.write('    fi\n')
    f.write('done < ligands.smi\n')
    f.write('echo "Prepared $COUNT ligands"\n')

# 06_submit_slurm_array.sh
with open(os.path.join(DOCK_DIR, '06_submit_slurm_array.sh'), 'w') as f:
    f.write('#!/bin/bash\n')
    f.write(f'#SBATCH --job-name=vina_dock\n')
    f.write(f'#SBATCH --gres=gpu:A100:1\n')
    f.write(f'#SBATCH --cpus-per-task=8\n')
    f.write(f'#SBATCH --mem=32G\n')
    f.write(f'#SBATCH --time=7-00:00:00\n')
    f.write(f'#SBATCH --array=1-{n_configs}%1000\n')
    f.write(f'#SBATCH --output=logs/vina_%A_%a.out\n')
    f.write(f'#SBATCH --error=logs/vina_%A_%a.err\n\n')
    f.write(f'set -e\n')
    f.write(f'cd "{DOCK_DIR}"\n')
    f.write(f'mkdir -p logs docking_results\n\n')
    f.write(f'CONF_FILE=$(printf "vina_configs/DOCK-%06d.conf" $SLURM_ARRAY_TASK_ID)\n\n')
    f.write(f'if [ ! -f "$CONF_FILE" ]; then\n')
    f.write(f'    echo "SKIP: $CONF_FILE not found (no PDB for this task)"\n')
    f.write(f'    exit 0\n')
    f.write(f'fi\n\n')
    f.write(f'echo "=== $(date) Task: $SLURM_ARRAY_TASK_ID ==="\n')
    f.write(f'head -1 "$CONF_FILE"\n\n')
    f.write(f'vina_gpu --config "$CONF_FILE" --gpu_id 0\n')
    f.write(f'echo "Done: $(date)"\n')

# ============================================================
# Final summary
# ============================================================
print(f'\n{"="*60}')
print(f'FINAL PACKAGE SUMMARY')
print(f'{"="*60}')
print(f'Filter: BE1+BE2 + E1+E2')
print(f'')
print(f'  Docking tasks:        {n_total:,}')
print(f'  Vina configs:         {n_configs:,}')
print(f'  No PDB (skipped):     {n_skip_no_pdb:,}')
print(f'  Unique compounds:     {len(unique_cpds):,}')
print(f'  Unique PDBs:          {len(pdb_cpds):,}')
print(f'  Unique targets:       {len(set(t["target_uniprot"] for t in filtered)):,}')
print(f'')

# Top targets
tgt_counts = Counter(t['target_symbol'] for t in filtered)
print(f'  Top 15 targets:')
for sym, cnt in tgt_counts.most_common(15):
    avg_score = sum(float(t['priority_score']) for t in filtered if t['target_symbol'] == sym) / cnt
    print(f'    {sym:20s}: {cnt:5d} tasks  avg_score={avg_score:.1f}')

# Top compounds by priority
cpd_scores = {}
for t in filtered:
    cid = t['compound_internal_id']
    s = float(t['priority_score'])
    if cid not in cpd_scores or s > cpd_scores[cid][0]:
        cpd_scores[cid] = (s, t['preferred_name'], t['bio_flags'])

top_cpds = sorted(cpd_scores.items(), key=lambda x: -x[1][0])[:10]
print(f'\n  Top 10 compounds by priority:')
for cid, (s, name, flags) in top_cpds:
    print(f'    {name[:50]:50s}  score={s:.1f}  {flags[:50]}')

# PDB coverage
with_pdb = sum(1 for t in filtered if t.get('pdb_id', '').strip())
with_res = sum(1 for t in filtered if t.get('has_binding_residues', '') == 'True')
print(f'\n  PDB coverage: {with_pdb}/{n_total} ({with_pdb/n_total*100:.1f}%)')
print(f'  Residue coverage: {with_res}/{n_total} ({with_res/n_total*100:.1f}%)')

# File listing
print(f'\n  Files:')
for f in sorted(os.listdir(DOCK_DIR)):
    fpath = os.path.join(DOCK_DIR, f)
    if os.path.isfile(fpath):
        size = os.path.getsize(fpath)
        print(f'    {f:40s} {size/1024:8.1f} KB')
    elif os.path.isdir(fpath) and f != '__pycache__':
        n = len(os.listdir(fpath))
        print(f'    {f:40s} {n:8d} items')

print('\nDone.')
