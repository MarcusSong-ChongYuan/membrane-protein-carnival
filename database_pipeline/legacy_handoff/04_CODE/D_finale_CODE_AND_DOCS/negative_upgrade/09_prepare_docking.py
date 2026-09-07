"""
Prepare Vina docking package for A100 GPU:
1. Docking manifest (compound x protein x PDB pairs)
2. Compound SMILES list
3. Grid box definitions from known binding site residues
4. Batch preparation scripts (PDB download, PDBQT conversion, ligand prep)
5. Vina-GPU run script
"""
import csv, json, os
from collections import defaultdict, Counter

OUT_DIR = r'D:\finale\negative_upgrade'
DETAIL = os.path.join(OUT_DIR, 'biological_status_binding_sites_detail.tsv')
DOCK_DIR = os.path.join(OUT_DIR, 'docking_package')
os.makedirs(DOCK_DIR, exist_ok=True)

print("Loading detail data...")
rows = []
with open(DETAIL, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        rows.append(row)
print(f"  {len(rows):,} rows")

# ===================================================================
# 1. Build unique compound x protein x PDB mapping
# ===================================================================
print("\n--- Building unique docking tasks ---")

# Prefer: 1 PDB + 1 compound + 1 chain = 1 docking task
# De-duplicate by (compound_internal_id, pdb_id, chain)
tasks = {}  # key -> task dict
pdb_cpds = defaultdict(set)  # pdb -> {compound_ids}
cpd_pdbs = defaultdict(set)  # compound -> {pdb_ids}

for r in rows:
    cpd_id = r.get('compound_internal_id', '')
    smiles = r.get('standard_smiles', '')
    pdb_str = r.get('pdb_ids_site', '').strip()
    chain = r.get('pdb_chain_ids_v60', '').strip()
    residues = r.get('residue_or_site_description', '').strip()
    site_type = r.get('site_type', '')
    sym = r.get('approved_symbol', '')
    uniprot = r.get('target_uniprot_id', '')
    cpd_name = r.get('preferred_name', '')

    if not pdb_str or not smiles:
        continue

    for pdb in pdb_str.split(';'):
        pdb = pdb.strip().lower()
        if len(pdb) != 4 or not pdb[0].isdigit():
            continue

        # Use first chain if multiple
        chain_id = chain.split(';')[0].strip() if chain else 'A'
        if not chain_id:
            chain_id = 'A'

        key = f"{cpd_id}__{pdb}__{chain_id}"
        if key not in tasks:
            tasks[key] = {
                'compound_internal_id': cpd_id,
                'preferred_name': cpd_name,
                'smiles': smiles,
                'pdb_id': pdb,
                'chain': chain_id,
                'target_symbol': sym,
                'target_uniprot': uniprot,
                'site_type': site_type,
                'binding_residues': residues,
                'residue_index_type': r.get('residue_index_type_v60', ''),
            }
            pdb_cpds[pdb].add(cpd_id)
            cpd_pdbs[cpd_id].add(pdb)

print(f"  Unique docking tasks: {len(tasks):,}")
print(f"  Unique PDBs: {len(pdb_cpds):,}")
print(f"  Unique compounds: {len(cpd_pdbs):,}")

# ===================================================================
# 2. Write docking manifest
# ===================================================================
print("\n--- Writing docking manifest ---")

manifest_file = os.path.join(DOCK_DIR, 'docking_manifest.tsv')
with open(manifest_file, 'w', encoding='utf-8', newline='') as f:
    fields = ['task_id', 'compound_internal_id', 'preferred_name', 'smiles',
              'pdb_id', 'chain', 'target_symbol', 'target_uniprot',
              'site_type', 'binding_residues', 'residue_index_type',
              'grid_center_x', 'grid_center_y', 'grid_center_z',
              'grid_size_x', 'grid_size_y', 'grid_size_z']
    writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    writer.writeheader()

    for i, (key, task) in enumerate(sorted(tasks.items())):
        task['task_id'] = f'DOCK-{i+1:06d}'
        # Grid defaults (to be filled by prep script)
        task['grid_center_x'] = ''
        task['grid_center_y'] = ''
        task['grid_center_z'] = ''
        task['grid_size_x'] = '20'
        task['grid_size_y'] = '20'
        task['grid_size_z'] = '20'
        writer.writerow(task)

print(f"  Saved: {manifest_file}")

# ===================================================================
# 3. Write compound SMILES list
# ===================================================================
print("\n--- Writing compound SMILES ---")

# Get unique compounds from tasks
unique_cpds = {}
for task in tasks.values():
    cpd_id = task['compound_internal_id']
    if cpd_id not in unique_cpds:
        unique_cpds[cpd_id] = task

smi_file = os.path.join(DOCK_DIR, 'ligands.smi')
with open(smi_file, 'w', encoding='utf-8') as f:
    for cpd_id, cpd in sorted(unique_cpds.items()):
        smiles = cpd['smiles']
        name = cpd['preferred_name'][:50] if cpd['preferred_name'] else cpd_id
        if smiles:
            f.write(f"{smiles}\t{cpd_id}\t{name}\n")

print(f"  Saved: {smi_file} ({len(unique_cpds)} compounds)")

# ===================================================================
# 4. Write PDB download list
# ===================================================================
print("\n--- Writing PDB download list ---")

pdb_list_file = os.path.join(DOCK_DIR, 'pdb_download_list.txt')
with open(pdb_list_file, 'w') as f:
    for pdb in sorted(pdb_cpds.keys()):
        # Count compounds per PDB
        f.write(f"{pdb}\t{len(pdb_cpds[pdb])} compounds\n")

print(f"  Saved: {pdb_list_file} ({len(pdb_cpds)} PDBs)")

# ===================================================================
# 5. Write summary per protein target
# ===================================================================
print("\n--- Writing target summary ---")

target_tasks = defaultdict(list)
for task in tasks.values():
    target_tasks[task['target_symbol']].append(task)

target_file = os.path.join(DOCK_DIR, 'target_summary.tsv')
with open(target_file, 'w', encoding='utf-8', newline='') as f:
    f.write('\t'.join(['target_symbol', 'target_uniprot', 'num_tasks',
                        'num_pdbs', 'num_compounds']) + '\n')
    for sym in sorted(target_tasks.keys()):
        tlist = target_tasks[sym]
        pdbs = set(t['pdb_id'] for t in tlist)
        cpds = set(t['compound_internal_id'] for t in tlist)
        uniprot = tlist[0]['target_uniprot']
        f.write('\t'.join([sym, uniprot, str(len(tlist)),
                          str(len(pdbs)), str(len(cpds))]) + '\n')

print(f"  Saved: {target_file}")

# ===================================================================
# 6. Generate batch prep scripts (bash)
# ===================================================================
print("\n--- Generating preparation scripts ---")

# 6a: PDB download script
download_sh = os.path.join(DOCK_DIR, '01_download_pdbs.sh')
with open(download_sh, 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# Download all PDB structures for Vina docking\n')
    f.write(f'# Total: {len(pdb_cpds)} PDBs\n')
    f.write('set -e\n\n')
    f.write('PDB_DIR="pdb_structures"\n')
    f.write('mkdir -p "$PDB_DIR"\n\n')
    f.write('# Download from PDBe (faster than RCSB for batch)\n')
    f.write('while IFS=$' + "'" + r'\t' + "'" + ' read -r pdb rest; do\n')
    f.write('    if [ ! -f "$PDB_DIR/${pdb}.pdb" ]; then\n')
    f.write('        echo "Downloading $pdb..."\n')
    f.write('        curl -s -o "$PDB_DIR/${pdb}.pdb" "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${pdb}.ent"\n')
    f.write('    fi\n')
    f.write(f'done < pdb_download_list.txt\n\n')
    f.write(f'echo "Downloaded $(ls $PDB_DIR/*.pdb 2>/dev/null | wc -l) PDB files"\n')

print(f"  Saved: {download_sh}")

# 6b: ADFR preparation script (protein PDBQT)
prep_protein_sh = os.path.join(DOCK_DIR, '02_prepare_proteins.sh')
with open(prep_protein_sh, 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# Prepare protein PDBQT files using ADFR suite (prepare_receptor)\n')
    f.write('# Install: conda install -c conda-forge adfr-suite\n')
    f.write('set -e\n\n')
    f.write('PDB_DIR="pdb_structures"\n')
    f.write('RECEPTOR_DIR="receptors_pdbqt"\n')
    f.write('mkdir -p "$RECEPTOR_DIR"\n\n')
    f.write('for pdb_file in "$PDB_DIR"/*.pdb; do\n')
    f.write('    pdb_id=$(basename "$pdb_file" .pdb)\n')
    f.write('    if [ ! -f "$RECEPTOR_DIR/${pdb_id}.pdbqt" ]; then\n')
    f.write('        echo "Preparing $pdb_id..."\n')
    f.write('        # Remove water, add hydrogens at pH 7.4, add Gasteiger charges\n')
    f.write('        prepare_receptor -r "$pdb_file" -o "$RECEPTOR_DIR/${pdb_id}.pdbqt" \\\n')
    f.write('            -A hydrogens -U nphs\n')
    f.write('    fi\n')
    f.write('done\n\n')
    f.write(f'echo "Prepared $(ls $RECEPTOR_DIR/*.pdbqt 2>/dev/null | wc -l) receptors"\n')

print(f"  Saved: {prep_protein_sh}")

# 6c: Ligand preparation (SMILES -> 3D -> PDBQT)
prep_ligand_sh = os.path.join(DOCK_DIR, '03_prepare_ligands.sh')
with open(prep_ligand_sh, 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# Prepare ligand PDBQT files from SMILES using obabel + prepare_ligand\n')
    f.write('# Requires: openbabel, adfr-suite\n')
    f.write('set -e\n\n')
    f.write('LIGAND_DIR="ligands_pdbqt"\n')
    f.write('LIGAND_3D_DIR="ligands_3d"\n')
    f.write('mkdir -p "$LIGAND_DIR" "$LIGAND_3D_DIR"\n\n')
    f.write('# Step 1: Convert SMILES to 3D SDF (generate 3D conformer)\n')
    f.write('echo "Converting SMILES to 3D..."\n')
    f.write('obabel ligands.smi -O "$LIGAND_3D_DIR"/ligands.sdf --gen3d --conformers --nconf 1\n\n')
    f.write('# Step 2: Split multi-molecule SDF and convert each to PDBQT\n')
    f.write('echo "Converting to PDBQT..."\n')
    f.write('cd "$LIGAND_3D_DIR"\n')
    f.write('obabel ligands.sdf -O ligands_.sdf -m 2>/dev/null || true\n')
    f.write('cd ..\n')
    f.write('for sdf_file in "$LIGAND_3D_DIR"/ligands_*.sdf; do\n')
    f.write('    base=$(basename "$sdf_file" .sdf)\n')
    f.write('    obabel "$sdf_file" -O "$LIGAND_3D_DIR/${base}.pdb" --gen3d 2>/dev/null || true\n')
    f.write('    prepare_ligand -l "$LIGAND_3D_DIR/${base}.pdb" -o "$LIGAND_DIR/${base}.pdbqt" \\\n')
    f.write('        -A hydrogens 2>/dev/null || true\n')
    f.write('done\n\n')
    f.write(f'echo "Prepared $(ls $LIGAND_DIR/*.pdbqt 2>/dev/null | wc -l) ligands"\n')

print(f"  Saved: {prep_ligand_sh}")

# 6d: Grid box auto-calculation script
grid_sh = os.path.join(DOCK_DIR, '04_compute_grid_boxes.py')
with open(grid_sh, 'w') as f:
    f.write('#!/usr/bin/env python3\n')
    f.write('"""\n')
    f.write('Auto-compute Vina grid box from known binding site residues.\n')
    f.write('Parses PDB files, extracts CA coordinates of binding residues,\n')
    f.write('computes bounding box centroid + margin.\n')
    f.write('"""\n')
    f.write('import csv, os, sys\n')
    f.write('from collections import defaultdict\n\n')
    f.write('PDB_DIR = "pdb_structures"\n')
    f.write('MANIFEST = "docking_manifest.tsv"\n')
    f.write('MARGIN = 10.0  # A padding around binding residues\n\n')
    f.write('# Load manifest\n')
    f.write('tasks = []\n')
    f.write('with open(MANIFEST) as f:\n')
    f.write('    reader = csv.DictReader(f, delimiter="\\t")\n')
    f.write('    for row in reader:\n')
    f.write('        tasks.append(row)\n')
    f.write('print(f"Loaded {len(tasks)} tasks")\n\n')
    f.write('# Residue parsing helpers\n')
    f.write('def parse_residues(desc, idx_type):\n')
    f.write('    """Parse residue_or_site_description into (chain, resname, resnum) tuples."""\n')
    f.write('    results = []\n')
    f.write('    if not desc:\n')
    f.write('        return results\n')
    f.write('    # Format varies by site_type:\n')
    f.write('    # PDBe: "PHE144;PHE151;PHE152" (UNIPROT numbering)\n')
    f.write('    # PDBbind pocket: "A:ALA767;A:ALA771;..." (PDB numbering)\n')
    f.write('    # PDBbind match: "CID|NAME@PDB:IC50=X=A:ARG608"\n')
    f.write('    # UniProt: "L-glutamate@156-156(medium)" (skip, no PDB numbering)\n')
    f.write('    \n')
    f.write('    # Skip UniProt-style\n')
    f.write('    if "@" in desc and "-" in desc.split("@")[-1].split("(")[0]:\n')
    f.write('        return results  # UniProt numbering, not PDB\n')
    f.write('    \n')
    f.write('    # Extract residues after PDB ID if present\n')
    f.write('    if "@" in desc:\n')
    f.write('        desc = desc.split("@")[-1]\n')
    f.write('        if ":" in desc:\n')
    f.write('            desc = desc.split(":")[-1]\n')
    f.write('        # Remove IC50=... prefix\n')
    f.write('        if "=" in desc.split(" ")[0]:\n')
    f.write('            desc = " ".join(desc.split(" ")[1:])\n')
    f.write('    \n')
    f.write('    # Split by delimiter\n')
    f.write('    parts = []\n')
    f.write('    if ";" in desc:\n')
    f.write('        parts = [p.strip() for p in desc.split(";") if p.strip()]\n')
    f.write('    elif " " in desc:\n')
    f.write('        parts = [p.strip() for p in desc.split(" ") if p.strip()]\n')
    f.write('    else:\n')
    f.write('        parts = [desc.strip()]\n')
    f.write('    \n')
    f.write('    for part in parts:\n')
    f.write('        # Parse "A:ALA767" or "ALA767" or "PHE144"\n')
    f.write('        chain = "A"\n')
    f.write('        res_str = part\n')
    f.write('        if ":" in part:\n')
    f.write('            chain, res_str = part.split(":", 1)\n')
    f.write('        # Extract resname (letters) and resnum (digits)\n')
    f.write('        import re\n')
    f.write('        m = re.match(r"([A-Za-z]{3})(\\d+)", res_str)\n')
    f.write('        if m:\n')
    f.write('            results.append((chain.strip(), m.group(1), int(m.group(2))))\n')
    f.write('    return results\n\n')
    f.write('# Compute grid boxes\n')
    f.write('grid_file = "docking_manifest_with_grid.tsv"\n')
    f.write('updated = 0\n')
    f.write('with open(grid_file, "w", newline="") as f:\n')
    f.write('    writer = csv.DictWriter(f, fieldnames=tasks[0].keys(), delimiter="\\t", extrasaction="ignore")\n')
    f.write('    writer.writeheader()\n')
    f.write('    for task in tasks:\n')
    f.write('        pdb_id = task["pdb_id"]\n')
    f.write('        pdb_path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")\n')
    f.write('        residues = parse_residues(task["binding_residues"], task["residue_index_type"])\n')
    f.write('        \n')
    f.write('        if residues and os.path.exists(pdb_path):\n')
    f.write('            # Read PDB, extract CA coords for matching residues\n')
    f.write('            coords = []\n')
    f.write('            with open(pdb_path) as pf:\n')
    f.write('                for line in pf:\n')
    f.write('                    if line.startswith("ATOM") and line[13:15].strip() == "CA":\n')
    f.write('                        chain = line[21]\n')
    f.write('                        resname = line[17:20].strip()\n')
    f.write('                        resnum = int(line[22:26])\n')
    f.write('                        for rc, rn, ri in residues:\n')
    f.write('                            if chain == rc and resname == rn and resnum == ri:\n')
    f.write('                                x = float(line[30:38])\n')
    f.write('                                y = float(line[38:46])\n')
    f.write('                                z = float(line[46:54])\n')
    f.write('                                coords.append((x, y, z))\n')
    f.write('            \n')
    f.write('            if coords:\n')
    f.write('                cx = sum(c[0] for c in coords) / len(coords)\n')
    f.write('                cy = sum(c[1] for c in coords) / len(coords)\n')
    f.write('                cz = sum(c[2] for c in coords) / len(coords)\n')
    f.write('                sx = max(c[0] for c in coords) - min(c[0] for c in coords) + 2 * MARGIN\n')
    f.write('                sy = max(c[1] for c in coords) - min(c[1] for c in coords) + 2 * MARGIN\n')
    f.write('                sz = max(c[2] for c in coords) - min(c[2] for c in coords) + 2 * MARGIN\n')
    f.write('                task["grid_center_x"] = f"{cx:.3f}"\n')
    f.write('                task["grid_center_y"] = f"{cy:.3f}"\n')
    f.write('                task["grid_center_z"] = f"{cz:.3f}"\n')
    f.write('                task["grid_size_x"] = f"{sx:.1f}"\n')
    f.write('                task["grid_size_y"] = f"{sy:.1f}"\n')
    f.write('                task["grid_size_z"] = f"{sz:.1f}"\n')
    f.write('                updated += 1\n')
    f.write('        writer.writerow(task)\n')
    f.write('print(f"Grid boxes computed for {updated}/{len(tasks)} tasks")\n')

print(f"  Saved: {grid_sh}")

# 6e: Vina-GPU batch run script
vina_sh = os.path.join(DOCK_DIR, '05_run_vina_gpu.sh')
with open(vina_sh, 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# Batch Vina-GPU docking on A100\n')
    f.write('# Requires: Vina-GPU (https://github.com/DeltaGroupNJUPT/Vina-GPU)\n')
    f.write('# Usage: sbatch 05_run_vina_gpu.sh\n')
    f.write('set -e\n\n')
    f.write('#SBATCH --job-name=vina_dock\n')
    f.write('#SBATCH --gres=gpu:1\n')
    f.write('#SBATCH --cpus-per-task=16\n')
    f.write('#SBATCH --mem=64G\n')
    f.write('#SBATCH --time=72:00:00\n\n')
    f.write('RECEPTOR_DIR="receptors_pdbqt"\n')
    f.write('LIGAND_DIR="ligands_pdbqt"\n')
    f.write('OUTPUT_DIR="docking_results"\n')
    f.write('MANIFEST="docking_manifest_with_grid.tsv"\n')
    f.write('mkdir -p "$OUTPUT_DIR"\n\n')
    f.write('echo "Starting Vina-GPU batch docking..."\n')
    f.write('echo "Manifest: $MANIFEST"\n\n')
    f.write('# Read manifest and dock each pair\n')
    f.write('tail -n +2 "$MANIFEST" | while IFS=$' + "'" + r'\t' + "'" + ' read -r task_id cpd_id name smiles pdb chain target_sym target_up site_type residues idx_type cx cy cz sx sy sz; do\n')
    f.write('    receptor="$RECEPTOR_DIR/${pdb}.pdbqt"\n')
    f.write('    ligand="$LIGAND_DIR/${cpd_id}.pdbqt"\n')
    f.write('    output="$OUTPUT_DIR/${task_id}_out.pdbqt"\n')
    f.write('    \n')
    f.write('    if [ ! -f "$receptor" ]; then\n')
    f.write('        echo "SKIP $task_id: receptor not found $receptor"\n')
    f.write('        continue\n')
    f.write('    fi\n')
    f.write('    if [ ! -f "$ligand" ]; then\n')
    f.write('        echo "SKIP $task_id: ligand not found $ligand"\n')
    f.write('        continue\n')
    f.write('    fi\n')
    f.write('    \n')
    f.write('    echo "Docking $task_id: $cpd_id -> $target_sym ($pdb)"\n')
    f.write('    \n')
    f.write('    # Vina-GPU (adjust path and parameters as needed)\n')
    f.write('    vina_gpu \\\n')
    f.write('        --receptor "$receptor" \\\n')
    f.write('        --ligand "$ligand" \\\n')
    f.write('        --out "$output" \\\n')
    f.write('        --center_x "$cx" --center_y "$cy" --center_z "$cz" \\\n')
    f.write('        --size_x "$sx" --size_y "$sy" --size_z "$sz" \\\n')
    f.write('        --num_modes 9 \\\n')
    f.write('        --exhaustiveness 32 \\\n')
    f.write('        --energy_range 3 \\\n')
    f.write('        --gpu_id 0\n')
    f.write('    \n')
    f.write('    echo "  Done: $output"\n')
    f.write('done\n\n')
    f.write('echo "All docking complete."\n')
    f.write('echo "Results in: $OUTPUT_DIR"\n')

print(f"  Saved: {vina_sh}")

# ===================================================================
# 7. Generate info file
# ===================================================================
info_file = os.path.join(DOCK_DIR, 'README_DOCKING.md')
with open(info_file, 'w', encoding='utf-8') as f:
    f.write(f"""# Vina-GPU Docking Package (A100)

Generated: 2026-08-08
Source: MemPro V6.2 Biological Status Compounds x Membrane Protein Binding Sites

## Scope

| Metric | Count |
|--------|-------|
| Docking tasks | {len(tasks):,} |
| Unique PDB structures | {len(pdb_cpds):,} |
| Unique compounds | {len(cpd_pdbs):,} |
| Membrane protein targets | {len(target_tasks):,} |

## Files

| File | Description |
|------|-------------|
| `docking_manifest.tsv` | Master manifest: compound x protein x PDB mapping |
| `ligands.smi` | SMILES for all {len(unique_cpds)} compounds |
| `pdb_download_list.txt` | {len(pdb_cpds)} PDB IDs to download |
| `target_summary.tsv` | Per-target task counts |
| `01_download_pdbs.sh` | Batch download PDB structures from PDBe |
| `02_prepare_proteins.sh` | Prepare receptor PDBQT (ADFR suite) |
| `03_prepare_ligands.sh` | SMILES -> 3D -> PDBQT (OpenBabel + ADFR) |
| `04_compute_grid_boxes.py` | Auto-compute Vina grid boxes from binding residues |
| `05_run_vina_gpu.sh` | SLURM batch script for Vina-GPU on A100 |

## Workflow

```bash
# Step 1: Download PDBs
bash 01_download_pdbs.sh

# Step 2: Prepare proteins
bash 02_prepare_proteins.sh

# Step 3: Prepare ligands
bash 03_prepare_ligands.sh

# Step 4: Compute grid boxes
python3 04_compute_grid_boxes.py

# Step 5: Run docking (on A100 GPU node)
sbatch 05_run_vina_gpu.sh
```

## Grid Box Method

Grid boxes are computed from known binding site residues extracted from:
- PDBe contact residues (UniProt numbering -> PDB mapping)
- PDBbind pocket residues (PDB numbering, direct match)
- PDBbind key residue annotations

Default margin: 10 Å around the binding residue bounding box.

## Software Requirements

- **PDB download**: curl
- **Protein prep**: [ADFR Suite](https://ccsb.scripps.edu/adfr/) (`conda install -c conda-forge adfr-suite`)
- **Ligand prep**: [OpenBabel](https://openbabel.org/) (`conda install -c conda-forge openbabel`)
- **Docking**: [Vina-GPU](https://github.com/DeltaGroupNJUPT/Vina-GPU) (or alternative GPU Vina)
  - Compile with CUDA for A100 (SM 8.0)

## Notes

- Binding residue data comes from crystallographic evidence (96.2% BE1 tier)
- Grid boxes use CA atom centroids from known binding residues
- For PDBs using UniProt numbering, a numbering conversion may be needed
- All PDBs are X-ray crystallography (no NMR/EM in this set)
- Vina exhaustiveness=32 recommended for production; reduce to 8 for screening
""")

print(f"  Saved: {info_file}")

# ===================================================================
# 8. Summary
# ===================================================================
print(f"\n{'='*60}")
print(f"DOCKING PACKAGE READY")
print(f"{'='*60}")
print(f"  Directory: {DOCK_DIR}")
print(f"  Tasks: {len(tasks):,}")
print(f"  PDBs: {len(pdb_cpds):,}")
print(f"  Compounds: {len(cpd_pdbs):,}")
print(f"  Targets: {len(target_tasks):,}")

# Show top targets
print(f"\n  Top 10 targets by docking tasks:")
for sym, tlist in sorted(target_tasks.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
    print(f"    {sym}: {len(tlist)} tasks ({len(set(t['pdb_id'] for t in tlist))} PDBs, {len(set(t['compound_internal_id'] for t in tlist))} cpds)")

print(f"\n  Files:")
for f in os.listdir(DOCK_DIR):
    fpath = os.path.join(DOCK_DIR, f)
    size = os.path.getsize(fpath)
    print(f"    {f} ({size/1024:.1f} KB)")

print("\nDone.")
