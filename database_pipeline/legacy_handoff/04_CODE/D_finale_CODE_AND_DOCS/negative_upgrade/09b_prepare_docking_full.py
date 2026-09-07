"""
Rebuild complete docking package: ALL bio-status compound × membrane protein pairs.
~152K unique pairs, ~124K with PDB, ~28K need AlphaFold.
"""
import csv, gzip, os, json
from collections import defaultdict, Counter

OUT_DIR = r'D:\finale\negative_upgrade'
DOCK_DIR = os.path.join(OUT_DIR, 'docking_package')
os.makedirs(DOCK_DIR, exist_ok=True)

MASTER = r'D:\finale\01_正式数据_V6.2\small_molecule_master_v1_3.tsv'
FORM = r'D:\finale\01_正式数据_V6.2\compound_form_hierarchy_v1_3.tsv'
EVIDENCE = r'D:\finale\01_正式数据_V6.2\binding_evidence_master_v6_2.tsv.gz'
SITES = r'D:\finale\01_正式数据_V6.2\binding_site_instances_v6_2.tsv.gz'
PROTEIN = r'D:\finale\01_正式数据_V6.2\human_membrane_protein_master_v6_2.tsv'

FLAGS = ['is_approved_drug', 'is_clinical_candidate', 'is_endogenous_ligand',
         'is_natural_product', 'is_chemical_probe']

# ===================================================================
# 1. Load bio-status compound data
# ===================================================================
print("Loading bio-status compounds...")
bio_cpds = {}  # cpd_id -> {smiles, name, flags, preferred_name}
with open(MASTER, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        for flag in FLAGS:
            if row.get(flag, '0') == '1':
                bio_cpds[row['compound_internal_id']] = {
                    'smiles': row.get('standard_smiles', ''),
                    'name': row.get('preferred_name', '')[:80],
                    'mw': row.get('molecular_weight', ''),
                    'xlogp': row.get('xlogp', ''),
                    'flags': [f for f in FLAGS if row.get(f, '0') == '1'],
                }
                break

# Also load compound form SMILES for any missing from master
cpd_smiles_from_form = {}
with open(FORM, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        cpd_id = row.get('compound_internal_id', '')
        if cpd_id in bio_cpds and not bio_cpds[cpd_id]['smiles']:
            smi = row.get('exact_smiles', '') or row.get('canonical_smiles', '')
            if smi:
                cpd_smiles_from_form[cpd_id] = smi

for cpd_id, smi in cpd_smiles_from_form.items():
    bio_cpds[cpd_id]['smiles'] = smi

cpds_with_smiles = sum(1 for v in bio_cpds.values() if v['smiles'])
cpds_without = len(bio_cpds) - cpds_with_smiles
print(f"  Bio-status compounds: {len(bio_cpds):,}")
print(f"    With SMILES: {cpds_with_smiles:,}")
print(f"    Without SMILES: {cpds_without:,}")

# ===================================================================
# 2. Load membrane protein data
# ===================================================================
print("Loading membrane proteins...")
mp_data = {}  # uniprot -> {symbol, name, class}
with open(PROTEIN, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        up = row.get('target_uniprot_id', '')
        mp_data[up] = {
            'symbol': row.get('approved_symbol', ''),
            'name': row.get('preferred_name', ''),
            'mp_class': row.get('membrane_protein_class', ''),
        }
mp_uniprots = set(mp_data.keys())
print(f"  Membrane proteins: {len(mp_uniprots):,}")

# ===================================================================
# 3. Collect all PDBs + binding site residues per (protein, PDB)
# ===================================================================
print("Collecting PDB/residue data from ALL evidence...")
protein_pdbs = defaultdict(set)
pdb_best_info = {}  # (uniprot, pdb) -> best residue data

with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        up = row.get('target_uniprot_id', '')
        pdb_str = row.get('pdb_ids', '').strip()
        residues = row.get('binding_site_residues', '').strip()
        chain = row.get('pdb_chain_ids_v60', '').strip()

        if not up or not pdb_str:
            continue

        for pdb in pdb_str.split(';'):
            pdb = pdb.strip().lower()
            if len(pdb) != 4 or not pdb[0].isdigit():
                continue
            protein_pdbs[up].add(pdb)

            key = (up, pdb)
            if residues and key not in pdb_best_info:
                ch = chain.split(';')[0].strip() if chain else 'A'
                pdb_best_info[key] = {
                    'residues': residues,
                    'chain': ch if ch else 'A',
                    'has_residues': True,
                }

# Fill PDBs without residues
for up, pdbs in protein_pdbs.items():
    for pdb in pdbs:
        key = (up, pdb)
        if key not in pdb_best_info:
            pdb_best_info[key] = {
                'residues': '',
                'chain': 'A',
                'has_residues': False,
            }

prots_with_pdb = len(protein_pdbs)
prots_without = len(mp_uniprots - set(protein_pdbs.keys()))
print(f"  Proteins with PDB: {prots_with_pdb:,}")
print(f"  Proteins without PDB: {prots_without:,}")
print(f"  PDB+residue entries: {sum(1 for v in pdb_best_info.values() if v['has_residues']):,}")
print(f"  PDB residue-less entries: {sum(1 for v in pdb_best_info.values() if not v['has_residues']):,}")

# ===================================================================
# 4. Load binding site residues (site-level, more detailed)
# ===================================================================
print("Loading binding site instances...")
site_residues = {}  # (uniprot, pdb) -> most detailed residue desc
site_types = {}  # (uniprot, pdb) -> site_type

with gzip.open(SITES, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        up = row.get('target_uniprot_id', '')
        pdb_str = row.get('pdb_ids_site', '').strip()
        desc = row.get('residue_or_site_description', '').strip()
        st = row.get('site_type', '')
        idx_type = row.get('residue_index_type_v60', '')

        if not up or not pdb_str:
            continue

        for pdb in pdb_str.split(';'):
            pdb = pdb.strip().lower()
            if len(pdb) != 4:
                continue
            key = (up, pdb)
            if desc and (key not in site_residues or len(desc) > len(site_residues.get(key, ''))):
                site_residues[key] = desc
                site_types[key] = f"{st}|{idx_type}"

# Merge into pdb_best_info
for key, info in pdb_best_info.items():
    if key in site_residues:
        info['site_residues'] = site_residues[key]
        info['site_type'] = site_types.get(key, '')

print(f"  Site-level residue entries: {len(site_residues):,}")

# ===================================================================
# 5. Build complete docking pairs
# ===================================================================
print("Building all docking pairs...")
pairs = []
seen = set()

with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        cpd = row.get('compound_internal_id', '')
        up = row.get('target_uniprot_id', '')
        if cpd not in bio_cpds or up not in mp_uniprots:
            continue

        key = (cpd, up)
        if key in seen:
            continue
        seen.add(key)

        # Only proceed if compound has SMILES
        if not bio_cpds[cpd]['smiles']:
            continue

        # Determine best PDB for this pair
        direct_pdb_str = row.get('pdb_ids', '').strip()
        direct_pdbs = [p.strip().lower() for p in direct_pdb_str.split(';')
                       if len(p.strip()) == 4 and p.strip()[0].isdigit()]

        # Choose best PDB: prefer one with residues
        best_pdb = None
        best_has_res = False
        best_chain = 'A'

        if direct_pdbs:
            # Check direct PDBs for residues
            for pdb in direct_pdbs:
                info = pdb_best_info.get((up, pdb), {})
                if info.get('has_residues'):
                    best_pdb = pdb
                    best_has_res = True
                    best_chain = info.get('chain', 'A')
                    break
            if not best_pdb:
                best_pdb = direct_pdbs[0]
                info = pdb_best_info.get((up, best_pdb), {})
                best_chain = info.get('chain', 'A')
        else:
            # Use best PDB for this protein
            available = protein_pdbs.get(up, set())
            if available:
                # Prefer one with residues
                for pdb in sorted(available):
                    info = pdb_best_info.get((up, pdb), {})
                    if info.get('has_residues'):
                        best_pdb = pdb
                        best_has_res = True
                        best_chain = info.get('chain', 'A')
                        break
                if not best_pdb:
                    best_pdb = sorted(available)[0]
                    info = pdb_best_info.get((up, best_pdb), {})
                    best_chain = info.get('chain', 'A')

        # Get site-level residue details
        site_data = pdb_best_info.get((up, best_pdb), {}) if best_pdb else {}
        site_detail = site_residues.get((up, best_pdb), '') if best_pdb else ''
        site_type = site_types.get((up, best_pdb), '')

        pairs.append({
            'compound_internal_id': cpd,
            'smiles': bio_cpds[cpd]['smiles'],
            'preferred_name': bio_cpds[cpd]['name'],
            'target_uniprot': up,
            'target_symbol': mp_data[up]['symbol'],
            'target_name': mp_data[up]['name'],
            'mp_class': mp_data[up]['mp_class'],
            'pdb_id': best_pdb or '',
            'chain': best_chain,
            'pdb_source': 'direct' if direct_pdbs else ('fallback' if best_pdb else 'none'),
            'has_binding_residues': best_has_res or bool(site_detail),
            'evidence_residues': site_data.get('residues', ''),
            'site_residues': site_detail,
            'site_type': site_type or site_data.get('site_type', ''),
            'bio_flags': ';'.join(bio_cpds[cpd]['flags']),
            'evidence_type': row.get('evidence_type', ''),
            'activity_nM': row.get('standard_value_nM', ''),
        })

print(f"  Total pairs: {len(pairs):,}")
with_pdb = sum(1 for p in pairs if p['pdb_id'])
with_res = sum(1 for p in pairs if p['has_binding_residues'])
without_pdb = sum(1 for p in pairs if not p['pdb_id'])
src_breakdown = Counter(p['pdb_source'] for p in pairs)
src_str = ', '.join(f'{k}:{v:,}' for k, v in src_breakdown.most_common())
print(f"  With PDB: {with_pdb:,} ({src_str})")
print(f"  With binding residues: {with_res:,}")
print(f"  Without PDB (need AlphaFold): {without_pdb:,}")

# Count unique PDBs
all_pdbs = set(p['pdb_id'] for p in pairs if p['pdb_id'])
print(f"  Unique PDBs to download: {len(all_pdbs):,}")

# ===================================================================
# 6. Write docking manifest
# ===================================================================
print("\n--- Writing docking manifest ---")
manifest_file = os.path.join(DOCK_DIR, 'docking_manifest.tsv')
task_id = 0
with open(manifest_file, 'w', encoding='utf-8', newline='') as f:
    fields = [
        'task_id', 'compound_internal_id', 'preferred_name', 'smiles',
        'target_uniprot', 'target_symbol', 'target_name', 'mp_class',
        'pdb_id', 'chain', 'pdb_source',
        'has_binding_residues', 'site_type',
        'site_residues', 'evidence_residues',  # both for grid box computation
        'grid_center_x', 'grid_center_y', 'grid_center_z',
        'grid_size_x', 'grid_size_y', 'grid_size_z',
        'bio_flags', 'evidence_type', 'activity_nM',
    ]
    writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    writer.writeheader()

    for p in pairs:
        task_id += 1
        p['task_id'] = f'DOCK-{task_id:06d}'
        p['grid_center_x'] = ''
        p['grid_center_y'] = ''
        p['grid_center_z'] = ''
        p['grid_size_x'] = '20'
        p['grid_size_y'] = '20'
        p['grid_size_z'] = '20'
        writer.writerow(p)

print(f"  Saved {task_id:,} tasks to {manifest_file}")

# ===================================================================
# 7. Write ligands SMILES
# ===================================================================
print("\n--- Writing ligands SMILES ---")
unique_cpds = {}
for p in pairs:
    cid = p['compound_internal_id']
    if cid not in unique_cpds:
        unique_cpds[cid] = p

smi_file = os.path.join(DOCK_DIR, 'ligands.smi')
with open(smi_file, 'w', encoding='utf-8') as f:
    for cid in sorted(unique_cpds.keys()):
        cpd = unique_cpds[cid]
        f.write(f"{cpd['smiles']}\t{cid}\t{cpd['preferred_name']}\n")
print(f"  Saved {len(unique_cpds):,} compounds")

# ===================================================================
# 8. Write PDB download list
# ===================================================================
print("\n--- Writing PDB download list ---")
pdb_cpds = defaultdict(set)
for p in pairs:
    if p['pdb_id']:
        pdb_cpds[p['pdb_id']].add(p['compound_internal_id'])

pdb_list_file = os.path.join(DOCK_DIR, 'pdb_download_list.txt')
with open(pdb_list_file, 'w') as f:
    for pdb in sorted(pdb_cpds.keys()):
        f.write(f"{pdb}\t{len(pdb_cpds[pdb])} compounds\n")
print(f"  Saved {len(pdb_cpds):,} PDBs")

# ===================================================================
# 9. Write AlphaFold-needed list
# ===================================================================
af_pairs = [p for p in pairs if not p['pdb_id']]
if af_pairs:
    af_uniprots = sorted(set(p['target_uniprot'] for p in af_pairs))
    af_file = os.path.join(DOCK_DIR, 'alphafold_needed.txt')
    with open(af_file, 'w') as f:
        f.write(f"# {len(af_uniprots)} UniProt IDs need AlphaFold structures\n")
        f.write(f"# {len(af_pairs)} pairs affected\n")
        f.write("# Download from: https://alphafold.ebi.ac.uk/files/AF-{UNIPROT}-F1-model_v4.pdb\n\n")
        for up in af_uniprots:
            n_pairs = sum(1 for p in af_pairs if p['target_uniprot'] == up)
            sym = af_pairs[0]['target_symbol'] if af_pairs else ''
            for p in af_pairs:
                if p['target_uniprot'] == up:
                    sym = p['target_symbol']
                    break
            f.write(f"{up}\t{sym}\t{n_pairs} pairs\n")
    print(f"  AlphaFold needed: {len(af_uniprots)} proteins, {len(af_pairs)} pairs → {af_file}")

# ===================================================================
# 10. Statistics
# ===================================================================
print(f"\n{'='*60}")
print(f"PACKAGE SUMMARY")
print(f"{'='*60}")
print(f"  Total docking tasks: {len(pairs):,}")
print(f"  Unique compounds: {len(unique_cpds):,}")
print(f"  Unique PDBs: {len(pdb_cpds):,}")
print(f"  Unique targets: {len(set(p['target_uniprot'] for p in pairs)):,}")
print(f"  With PDB + residues: {with_res:,}")
print(f"  With PDB, no residues: {with_pdb - with_res:,}")
print(f"  Need AlphaFold: {len(af_pairs):,}")

# PDB source breakdown
src_counts = Counter(p['pdb_source'] for p in pairs)
print(f"\n  PDB source:")
for src, cnt in src_counts.most_common():
    print(f"    {src}: {cnt:,}")

# Top targets
target_counts = Counter(p['target_symbol'] for p in pairs)
print(f"\n  Top 15 targets:")
for sym, cnt in target_counts.most_common(15):
    print(f"    {sym:20s}: {cnt:6d} tasks")

print("\nDone. Run generate_vina_configs.py to rebuild .conf files.")
