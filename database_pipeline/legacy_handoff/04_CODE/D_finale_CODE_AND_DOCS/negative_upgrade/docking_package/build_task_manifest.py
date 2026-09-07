#!/usr/bin/env python3
"""
Build a task-level mapping table explaining how unique (protein, compound)
pairs expand into DOCK-XXXXX tasks.

Output: docking_task_map.tsv with columns:
  docking_task_id  - unique task ID (DOCK-XXXXXX)
  source_pair_id   - unique (protein_id, compound_id) pair key
  protein_id       - UniProt accession
  compound_id      - internal compound ID
  pdb_id           - PDB structure used
  chain_id         - chain identifier
  binding_site_id  - derived from evidence source
  config_file      - corresponding vina config path
  grid_method      - how grid box was determined
  grid_qc_status   - QC status of grid box
"""
import csv
import os
import hashlib
from collections import defaultdict

MANIFEST = "docking_manifest.tsv"
OUTPUT = "docking_task_map.tsv"
CONFIG_DIR = "vina_configs"

# Read manifest
tasks = []
with open(MANIFEST, encoding='utf-8', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        tasks.append(row)

print(f"Manifest rows: {len(tasks)}")

# Build pair → tasks mapping
pair_tasks = defaultdict(list)
for t in tasks:
    protein_id = t.get('target_uniprot', '').strip()
    compound_id = t.get('compound_internal_id', '').strip()
    if protein_id and compound_id:
        pair_key = f"{protein_id}__{compound_id}"
        pair_tasks[pair_key].append(t)

unique_pairs = len(pair_tasks)
total_tasks = len(tasks)
print(f"Unique (protein, compound) pairs: {unique_pairs}")
print(f"Total tasks: {total_tasks}")
print(f"Avg tasks per pair: {total_tasks / unique_pairs:.1f}")

# Build expansion summary
expansion_examples = []
for pair_key, pair_task_list in sorted(pair_tasks.items(), key=lambda x: -len(x[1]))[:5]:
    protein_id, compound_id = pair_key.split('__')
    pdbs = sorted(set(t.get('pdb_id', '') for t in pair_task_list))
    chains = sorted(set(t.get('chain', '') for t in pair_task_list))
    expansion_examples.append({
        'pair': pair_key,
        'protein': protein_id,
        'compound': compound_id,
        'task_count': len(pair_task_list),
        'pdbs': pdbs,
        'chains': chains,
    })

print("\n=== Top 5 most-expanded pairs ===")
for ex in expansion_examples:
    print(f"  {ex['pair']}: {ex['task_count']} tasks")
    print(f"    PDBs: {ex['pdbs']}")
    print(f"    Chains: {ex['chains']}")

# Determine grid method and QC status for each task
# (will be filled by 04_compute_grid_boxes.py)
fieldnames = [
    'docking_task_id', 'source_pair_id', 'protein_id', 'compound_id',
    'pdb_id', 'chain_id', 'binding_site_id',
    'config_file', 'grid_method', 'mapped_residue_count',
    'requested_residue_count', 'grid_qc_status'
]

with open(OUTPUT, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()

    for t in tasks:
        task_id = t.get('task_id', '').strip()
        protein_id = t.get('target_uniprot', '').strip()
        compound_id = t.get('compound_internal_id', '').strip()
        pdb_id = t.get('pdb_id', '').strip()
        chain = t.get('chain', '').strip()

        pair_key = f"{protein_id}__{compound_id}"

        # Determine binding site ID
        site_residues = t.get('site_residues', '').strip()
        evidence_residues = t.get('evidence_residues', '').strip()
        if site_residues:
            site_hash = hashlib.md5(site_residues.encode()).hexdigest()[:8]
            site_id = f"site_{site_hash}"
        elif evidence_residues:
            site_hash = hashlib.md5(evidence_residues.encode()).hexdigest()[:8]
            site_id = f"evi_{site_hash}"
        else:
            site_id = "unknown"

        # Config file
        config_file = f"vina_configs/{task_id}.conf"
        config_exists = os.path.exists(os.path.join(CONFIG_DIR, f"{task_id}.conf"))

        # Count residues
        residues_used = site_residues or evidence_residues
        residue_count = 0
        if residues_used:
            residue_count = len([r for r in residues_used.split(';') if r.strip()])

        writer.writerow({
            'docking_task_id': task_id,
            'source_pair_id': pair_key,
            'protein_id': protein_id,
            'compound_id': compound_id,
            'pdb_id': pdb_id,
            'chain_id': chain,
            'binding_site_id': site_id,
            'config_file': config_file if config_exists else f"{config_file} [MISSING]",
            'grid_method': '',  # filled by 04
            'mapped_residue_count': '',  # filled by 04
            'requested_residue_count': residue_count,
            'grid_qc_status': '',  # filled by 04
        })

# Verify: task count == config count
config_count = len([f for f in os.listdir(CONFIG_DIR) if f.endswith('.conf')])
print(f"\n=== Verification ===")
print(f"Tasks in map: {total_tasks}")
print(f"Config files: {config_count}")
if total_tasks == config_count:
    print("✅ MATCH: task count == config count")
else:
    missing = total_tasks - config_count
    print(f"⚠️  MISMATCH: {abs(missing)} {'missing' if missing > 0 else 'extra'} configs")
    # List some missing tasks
    if missing > 0:
        config_ids = set(f.split('.')[0] for f in os.listdir(CONFIG_DIR) if f.endswith('.conf'))
        task_ids = set(t['task_id'].strip() for t in tasks)
        missing_ids = task_ids - config_ids
        print(f"  Missing config examples: {sorted(list(missing_ids))[:10]}")

print(f"\nOutput: {OUTPUT}")
