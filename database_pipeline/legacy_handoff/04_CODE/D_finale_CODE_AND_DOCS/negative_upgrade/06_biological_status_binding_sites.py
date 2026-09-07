"""
Build a detailed table: biological status drugs with membrane protein binding sites.
Includes: compound identity, protein identity, binding site residue-level detail.
"""
import csv, gzip, json, os
from collections import Counter, defaultdict

OUT_DIR = r'D:\finale\negative_upgrade'
MASTER_DIR = r'D:\finale\01_正式数据_V6.2'

COMPOUND_MASTER = os.path.join(MASTER_DIR, 'small_molecule_master_v1_3.tsv')
PROTEIN_MASTER = os.path.join(MASTER_DIR, 'human_membrane_protein_master_v6_2.tsv')
BINDING_EVIDENCE = os.path.join(MASTER_DIR, 'binding_evidence_master_v6_2.tsv.gz')
BINDING_SITES = os.path.join(MASTER_DIR, 'binding_site_instances_v6_2.tsv.gz')

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("Biological Status Drugs × Membrane Protein Binding Sites")
print("=" * 60)

# ==========================================================
# Step 1: Load biological status compounds from master
# ==========================================================
print("\n--- Step 1: Loading biological status compounds ---")

bio_flags = ['is_approved_drug', 'is_clinical_candidate', 'is_endogenous_ligand',
             'is_natural_product', 'is_chemical_probe']

bio_compounds = {}  # compound_internal_id -> compound_info dict
flag_counts = Counter()
flag_compound_sets = {f: set() for f in bio_flags}

with open(COMPOUND_MASTER, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        flags_set = []
        for flag in bio_flags:
            if row.get(flag, '').strip() == '1':
                flags_set.append(flag)
                flag_counts[flag] += 1
                flag_compound_sets[flag].add(row['compound_internal_id'])

        if flags_set:
            cpd_id = row['compound_internal_id']
            bio_compounds[cpd_id] = {
                'compound_internal_id': cpd_id,
                'preferred_name': row.get('preferred_name', ''),
                'standard_smiles': row.get('standard_smiles', ''),
                'standard_inchi': row.get('standard_inchi', ''),
                'standard_inchikey': row.get('standard_inchikey', ''),
                'molecular_formula': row.get('molecular_formula', ''),
                'molecular_weight': row.get('molecular_weight', ''),
                'xlogp': row.get('xlogp', ''),
                'tpsa': row.get('tpsa', ''),
                'pubchem_cids': row.get('pubchem_cids', ''),
                'chembl_ids': row.get('chembl_ids', ''),
                'drugcentral_ids': row.get('drugcentral_ids', ''),
                'chebi_ids': row.get('chebi_ids', ''),
                'gtopdb_ligand_ids': row.get('gtopdb_ligand_ids', ''),
                'hmdb_ids': row.get('hmdb_ids', ''),
                'cas_numbers': row.get('cas_numbers', ''),
                'is_approved_drug': row.get('is_approved_drug', '0'),
                'is_clinical_candidate': row.get('is_clinical_candidate', '0'),
                'is_endogenous_ligand': row.get('is_endogenous_ligand', '0'),
                'is_natural_product': row.get('is_natural_product', '0'),
                'is_chemical_probe': row.get('is_chemical_probe', '0'),
                'development_status': row.get('development_status', ''),
                'max_chembl_phase': row.get('max_chembl_phase', ''),
                'first_approval_year': row.get('first_approval_year', ''),
                'compound_category_tags': row.get('compound_category_tags', ''),
                'source_databases': row.get('source_databases', ''),
                'synonym_count': row.get('synonym_count', ''),
                'best_binding_evidence_level': row.get('best_binding_evidence_level', ''),
                'binding_evidence_count': row.get('binding_evidence_count', ''),
                'identity_confidence': row.get('identity_confidence', ''),
                'record_qc_status': row.get('record_qc_status', ''),
            }

print(f"  Total biological-status compounds: {len(bio_compounds):,}")
for flag in bio_flags:
    print(f"  {flag}: {flag_counts[flag]}")

# Compute intersection/union stats
total_unique = len(bio_compounds)
print(f"  Unique compounds (union): {total_unique:,}")
# How many have multiple flags?
multi_flag = sum(1 for cpd in bio_compounds.values()
                 if sum(int(cpd[f]) for f in bio_flags) > 1)
print(f"  Compounds with multiple flags: {multi_flag}")

# ==========================================================
# Step 2: Load binding evidence for these compounds
# ==========================================================
print("\n--- Step 2: Finding binding evidence for bio-status compounds ---")

bio_cpd_ids = set(bio_compounds.keys())
evidence_records = []  # list of evidence dicts for our compounds
evidence_by_cpd = defaultdict(list)
source_dbs = Counter()

with gzip.open(BINDING_EVIDENCE, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cpd_id = row.get('compound_internal_id', '').strip()
        if cpd_id in bio_cpd_ids:
            # Only keep membrane protein records (E1, E2 tier)
            tier = row.get('protein_release_tier', '').strip()
            if tier in ('E1', 'E2'):
                evidence_records.append(row)
                evidence_by_cpd[cpd_id].append(row['evidence_id'])
                source_dbs[row.get('source_database', '')] += 1

print(f"  Binding evidence records (E1/E2 membrane proteins): {len(evidence_records):,}")
print(f"  Compounds with binding evidence: {len(evidence_by_cpd):,}")

# ==========================================================
# Step 3: Load binding site instances for these evidence records
# ==========================================================
print("\n--- Step 3: Loading binding site instances ---")

evidence_ids = {r['evidence_id'] for r in evidence_records}
site_records = []
site_by_evidence = defaultdict(list)
site_source_dbs = Counter()
site_types = Counter()

with gzip.open(BINDING_SITES, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        ev_id = row.get('evidence_id', '').strip()
        if ev_id in evidence_ids:
            site_records.append(row)
            site_by_evidence[ev_id].append(row)
            site_source_dbs[row.get('source_database', '')] += 1
            site_types[row.get('site_type', '')] += 1

print(f"  Binding site instances: {len(site_records):,}")
print(f"  Evidence records with site data: {len(site_by_evidence):,}")
print(f"  Site types: {dict(site_types)}")

# ==========================================================
# Step 4: Load protein master for membrane protein details
# ==========================================================
print("\n--- Step 4: Loading membrane protein details ---")

# Get unique target_uniprot_ids from our evidence
target_uniprots = {r['target_uniprot_id'] for r in evidence_records}
protein_info = {}

with open(PROTEIN_MASTER, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        up_id = row.get('target_uniprot_id', '').strip()
        if up_id in target_uniprots:
            protein_info[up_id] = {
                'target_uniprot_id': up_id,
                'approved_symbol': row.get('approved_symbol', ''),
                'protein_name': row.get('protein_name', ''),
                'gene_names': row.get('gene_names', ''),
                'membrane_class': row.get('membrane_class_v52', ''),
                'membrane_topology': row.get('membrane_topology', ''),
                'transmembrane_count': row.get('transmembrane_count', ''),
                'functional_primary_class': row.get('functional_primary_class', ''),
                'functional_subclass': row.get('functional_subclass', ''),
                'release_tier': row.get('release_tier_v5', ''),
                'subcellular_location': row.get('subcellular_location', ''),
                'sequence_length': row.get('sequence_length', ''),
                'membrane_scope': row.get('membrane_scope', ''),
                'uniprot_entry_name': row.get('uniprot_entry_name', ''),
            }

print(f"  Membrane proteins with binding: {len(protein_info):,}")

# ==========================================================
# Step 5: Build the merged table (site-level granularity)
# ==========================================================
print("\n--- Step 5: Building merged table ---")

# Build evidence lookup
ev_lookup = {r['evidence_id']: r for r in evidence_records}

output_rows = []

for site in site_records:
    ev_id = site['evidence_id']
    ev = ev_lookup.get(ev_id)
    if not ev:
        continue

    cpd_id = site.get('compound_internal_id', ev.get('compound_internal_id', '')).strip()
    cpd_info = bio_compounds.get(cpd_id, {})
    up_id = site.get('target_uniprot_id', ev.get('target_uniprot_id', '')).strip()
    prot_info = protein_info.get(up_id, {})

    # Determine which biological status flags apply
    active_flags = ';'.join(f for f in bio_flags if cpd_info.get(f, '0') == '1')

    row = {
        # Compound identity
        'compound_internal_id': cpd_id,
        'preferred_name': cpd_info.get('preferred_name', ''),
        'standard_smiles': cpd_info.get('standard_smiles', ''),
        'standard_inchikey': cpd_info.get('standard_inchikey', ''),
        'molecular_formula': cpd_info.get('molecular_formula', ''),
        'molecular_weight': cpd_info.get('molecular_weight', ''),
        'xlogp': cpd_info.get('xlogp', ''),
        'tpsa': cpd_info.get('tpsa', ''),
        'pubchem_cids': cpd_info.get('pubchem_cids', ''),
        'chembl_ids': cpd_info.get('chembl_ids', ''),
        'drugcentral_ids': cpd_info.get('drugcentral_ids', ''),
        'chebi_ids': cpd_info.get('chebi_ids', ''),
        'gtopdb_ligand_ids': cpd_info.get('gtopdb_ligand_ids', ''),
        'cas_numbers': cpd_info.get('cas_numbers', ''),
        # Biological status
        'biological_status_flags': active_flags,
        'is_approved_drug': cpd_info.get('is_approved_drug', '0'),
        'is_clinical_candidate': cpd_info.get('is_clinical_candidate', '0'),
        'is_endogenous_ligand': cpd_info.get('is_endogenous_ligand', '0'),
        'is_natural_product': cpd_info.get('is_natural_product', '0'),
        'is_chemical_probe': cpd_info.get('is_chemical_probe', '0'),
        'development_status': cpd_info.get('development_status', ''),
        'max_chembl_phase': cpd_info.get('max_chembl_phase', ''),
        'first_approval_year': cpd_info.get('first_approval_year', ''),
        'compound_source_databases': cpd_info.get('source_databases', ''),
        'identity_confidence': cpd_info.get('identity_confidence', ''),
        # Protein identity
        'target_uniprot_id': up_id,
        'approved_symbol': prot_info.get('approved_symbol', ''),
        'protein_name': prot_info.get('protein_name', ''),
        'gene_names': prot_info.get('gene_names', ''),
        'membrane_class': prot_info.get('membrane_class', ''),
        'membrane_topology': prot_info.get('membrane_topology', ''),
        'transmembrane_count': prot_info.get('transmembrane_count', ''),
        'functional_primary_class': prot_info.get('functional_primary_class', ''),
        'functional_subclass': prot_info.get('functional_subclass', ''),
        'protein_release_tier': prot_info.get('release_tier', ''),
        'subcellular_location': prot_info.get('subcellular_location', ''),
        # Binding evidence info
        'evidence_id': ev_id,
        'evidence_tier': ev.get('evidence_tier', ''),
        'evidence_type': ev.get('evidence_type', ''),
        'activity_type': ev.get('activity_type', ''),
        'standard_value_nM': ev.get('standard_value_nM', ''),
        'activity_outcome_v60': ev.get('activity_outcome_v60', ''),
        'assay_or_mechanism': ev.get('assay_or_mechanism', ''),
        'pdb_ids_evidence': ev.get('pdb_ids', ''),
        'ligand_het_id': ev.get('ligand_het_id', ''),
        'pubmed_ids': ev.get('pubmed_ids', ''),
        'doi': ev.get('doi', ''),
        'evidence_source_database': ev.get('source_database', ''),
        # Binding site residue-level detail
        'binding_site_instance_id': site.get('binding_site_instance_id', ''),
        'site_type': site.get('site_type', ''),
        'pdb_ids_site': site.get('pdb_ids', ''),
        'residue_or_site_description': site.get('residue_or_site_description', ''),
        'site_compound_specificity': site.get('site_compound_specificity', ''),
        'residue_index_type_v60': site.get('residue_index_type_v60', ''),
        'pdb_chain_ids_v60': site.get('pdb_chain_ids_v60', ''),
        'site_source_database': site.get('source_database', ''),
        'site_record_qc_status': site.get('record_qc_status', ''),
    }
    output_rows.append(row)

print(f"  Total merged rows (site-level): {len(output_rows):,}")

# ==========================================================
# Step 6: Write main output table
# ==========================================================
print("\n--- Step 6: Writing output ---")

output_file = os.path.join(OUT_DIR, 'biological_status_binding_sites_detail.tsv')
fieldnames = list(output_rows[0].keys()) if output_rows else []

with open(output_file, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(output_rows)

print(f"  Detail table: {output_file}")
print(f"  Rows: {len(output_rows):,} x {len(fieldnames)} columns")

# ==========================================================
# Step 7: Build source statistics table
# ==========================================================
print("\n--- Step 7: Building source statistics ---")

# Per-source breakdown
source_stats = defaultdict(lambda: {
    'total_sites': 0,
    'unique_compounds': set(),
    'unique_proteins': set(),
    'unique_pdbs': set(),
    'site_types': Counter(),
    'evidence_tiers': Counter(),
    'biological_status': Counter(),
})

for row in output_rows:
    src = row['site_source_database'] or row['evidence_source_database'] or 'unknown'
    stats = source_stats[src]
    stats['total_sites'] += 1
    stats['unique_compounds'].add(row['compound_internal_id'])
    stats['unique_proteins'].add(row['target_uniprot_id'])
    if row['pdb_ids_site']:
        for pdb in row['pdb_ids_site'].split(';'):
            pdb = pdb.strip()
            if pdb:
                stats['unique_pdbs'].add(pdb)
    stats['site_types'][row['site_type']] += 1
    stats['evidence_tiers'][row['evidence_tier']] += 1
    stats['biological_status'][row['biological_status_flags']] += 1

# Write summary table
summary_file = os.path.join(OUT_DIR, 'biological_status_binding_sites_source_summary.tsv')
with open(summary_file, 'w', encoding='utf-8', newline='') as f:
    f.write('\t'.join([
        'source_database', 'total_binding_sites', 'unique_compounds', 'unique_proteins',
        'unique_pdb_structures', 'site_types_breakdown', 'evidence_tiers_breakdown',
        'biological_status_breakdown'
    ]) + '\n')

    for src in sorted(source_stats.keys()):
        stats = source_stats[src]
        f.write('\t'.join([
            src,
            str(stats['total_sites']),
            str(len(stats['unique_compounds'])),
            str(len(stats['unique_proteins'])),
            str(len(stats['unique_pdbs'])),
            json.dumps(dict(stats['site_types']), ensure_ascii=False),
            json.dumps(dict(stats['evidence_tiers']), ensure_ascii=False),
            json.dumps(dict(stats['biological_status']), ensure_ascii=False),
        ]) + '\n')

print(f"  Source summary: {summary_file}")

# ==========================================================
# Step 8: Compound-centric summary
# ==========================================================
print("\n--- Step 8: Compound-centric summary ---")

cpd_summary = defaultdict(lambda: {
    'compound_internal_id': '',
    'preferred_name': '',
    'biological_status': '',
    'protein_count': set(),
    'site_count': 0,
    'site_types': Counter(),
    'evidence_tiers': Counter(),
    'source_dbs': set(),
    'pdb_ids': set(),
})

for row in output_rows:
    cpd_id = row['compound_internal_id']
    cs = cpd_summary[cpd_id]
    cs['compound_internal_id'] = cpd_id
    cs['preferred_name'] = row['preferred_name']
    cs['biological_status'] = row['biological_status_flags']
    cs['protein_count'].add(row['target_uniprot_id'])
    cs['site_count'] += 1
    cs['site_types'][row['site_type']] += 1
    cs['evidence_tiers'][row['evidence_tier']] += 1
    cs['source_dbs'].add(row['site_source_database'] or row['evidence_source_database'])
    if row['pdb_ids_site']:
        for pdb in row['pdb_ids_site'].split(';'):
            if pdb.strip():
                cs['pdb_ids'].add(pdb.strip())

cpd_file = os.path.join(OUT_DIR, 'biological_status_compound_summary.tsv')
with open(cpd_file, 'w', encoding='utf-8', newline='') as f:
    f.write('\t'.join([
        'compound_internal_id', 'preferred_name', 'biological_status_flags',
        'membrane_protein_count', 'binding_site_count',
        'site_types', 'evidence_tiers', 'source_databases', 'pdb_count'
    ]) + '\n')

    for cpd_id in sorted(cpd_summary.keys()):
        cs = cpd_summary[cpd_id]
        f.write('\t'.join([
            cpd_id,
            cs['preferred_name'],
            cs['biological_status'],
            str(len(cs['protein_count'])),
            str(cs['site_count']),
            json.dumps(dict(cs['site_types']), ensure_ascii=False),
            json.dumps(dict(cs['evidence_tiers']), ensure_ascii=False),
            ';'.join(sorted(cs['source_dbs'])),
            str(len(cs['pdb_ids'])),
        ]) + '\n')

print(f"  Compound summary: {cpd_file}")
print(f"  Total compounds with sites: {len(cpd_summary):,}")

# ==========================================================
# Step 9: Final stats
# ==========================================================
print("\n" + "=" * 60)
print("FINAL STATISTICS")
print("=" * 60)

# How many biological-status compounds have binding sites?
total_bio = len(bio_compounds)
with_sites = len(set(r['compound_internal_id'] for r in output_rows))
without_sites = total_bio - with_sites

print(f"  Biological status compounds total: {total_bio:,}")
print(f"  With membrane protein binding sites: {with_sites:,}")
print(f"  Without binding sites: {without_sites:,}")

# By biological status
for flag in bio_flags:
    flag_cpds = set(r['compound_internal_id'] for r in output_rows
                     if r.get(flag, '0') == '1')
    print(f"  {flag}: {len(flag_cpds)} with sites (of {flag_counts[flag]} total)")

# By evidence tier
tier_counts = Counter(r['evidence_tier'] for r in output_rows)
print(f"\n  Evidence tiers: {dict(tier_counts)}")

# By site type
st_counts = Counter(r['site_type'] for r in output_rows)
print(f"  Site types: {dict(st_counts)}")

# Top proteins
prot_counts = Counter(r['approved_symbol'] for r in output_rows)
print(f"\n  Top 10 membrane proteins by site count:")
for sym, cnt in prot_counts.most_common(10):
    print(f"    {sym}: {cnt}")

print(f"\n  Output files:")
print(f"    1. {output_file}")
print(f"    2. {summary_file}")
print(f"    3. {cpd_file}")
print("\nDone.")
