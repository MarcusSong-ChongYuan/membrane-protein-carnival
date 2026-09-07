"""
Export ALL experimental binding site pairs for external sharing.
No bio-status filter, no E-tier filter — every small-molecule × membrane-protein
pair with experimental binding site evidence from the MemPro V6.2 database.

Output: full_binding_site_pairs.tsv (~95,000 rows x 57 columns)
"""
import csv, gzip, os
from collections import Counter
from datetime import datetime

MASTER_DIR = r'D:\finale\01_正式数据_V6.2'
OUT_DIR = r'D:\finale\negative_upgrade'
OUT_FILE = os.path.join(OUT_DIR, 'full_binding_site_pairs.tsv')

COMPOUND_MASTER = os.path.join(MASTER_DIR, 'small_molecule_master_v1_3.tsv')
PROTEIN_MASTER  = os.path.join(MASTER_DIR, 'human_membrane_protein_master_v6_2.tsv')
BINDING_EVIDENCE = os.path.join(MASTER_DIR, 'binding_evidence_master_v6_2.tsv.gz')
BINDING_SITES    = os.path.join(MASTER_DIR, 'binding_site_instances_v6_2.tsv.gz')

print("=" * 60)
print("Full Binding Site Pair Export")
print("=" * 60)

# ============================================================
# Step 1: Load all compounds (not just bio-status)
# ============================================================
print("\n[1/5] Loading compound master...")
compounds = {}
with open(COMPOUND_MASTER, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        cid = row['compound_internal_id']
        compounds[cid] = {
            'compound_internal_id': cid,
            'preferred_name': row.get('preferred_name', ''),
            'standard_smiles': row.get('standard_smiles', ''),
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
            'cas_numbers': row.get('cas_numbers', ''),
            'is_approved_drug': row.get('is_approved_drug', '0'),
            'is_clinical_candidate': row.get('is_clinical_candidate', '0'),
            'is_endogenous_ligand': row.get('is_endogenous_ligand', '0'),
            'is_natural_product': row.get('is_natural_product', '0'),
            'is_chemical_probe': row.get('is_chemical_probe', '0'),
            'development_status': row.get('development_status', ''),
            'max_chembl_phase': row.get('max_chembl_phase', ''),
            'first_approval_year': row.get('first_approval_year', ''),
            'source_databases': row.get('source_databases', ''),
            'identity_confidence': row.get('identity_confidence', ''),
            'compound_scope_class': row.get('compound_scope_class', ''),
            'molecule_types': row.get('molecule_types', ''),
        }
print(f"  Loaded {len(compounds):,} compounds")

# ============================================================
# Step 2: Load all membrane proteins
# ============================================================
print("\n[2/5] Loading membrane protein master...")
proteins = {}
with open(PROTEIN_MASTER, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        up = row['target_uniprot_id']
        proteins[up] = {
            'target_uniprot_id': up,
            'approved_symbol': row.get('approved_symbol', ''),
            'protein_name': row.get('protein_name', ''),
            'gene_names': row.get('gene_names', ''),
            'membrane_class_v52': row.get('membrane_class_v52', ''),
            'membrane_topology': row.get('membrane_topology', ''),
            'transmembrane_count': row.get('transmembrane_count', ''),
            'functional_primary_class': row.get('functional_primary_class', ''),
            'functional_subclass': row.get('functional_subclass', ''),
            'evidence_level_v52': row.get('evidence_level_v52', ''),
            'membrane_scope': row.get('membrane_scope', ''),
            'subcellular_location': row.get('subcellular_location', ''),
            'sequence_length': row.get('sequence_length', ''),
        }
print(f"  Loaded {len(proteins):,} membrane proteins")

# ============================================================
# Step 3: Load binding sites + join evidence
# ============================================================
print("\n[3/5] Loading binding site instances + joining evidence...")

# First pass: collect all evidence_ids from site instances
site_ev_ids = set()
with gzip.open(BINDING_SITES, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        site_ev_ids.add(row['evidence_id'])
print(f"  Site instances evidence IDs: {len(site_ev_ids):,}")

# Second pass: scan evidence for matching records (membrane protein filter)
ev_data = {}  # evidence_id -> evidence row
ev_in_scope = set()
with gzip.open(BINDING_EVIDENCE, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        ev_id = row['evidence_id']
        if ev_id in site_ev_ids and row.get('target_uniprot_id', '') in proteins:
            ev_data[ev_id] = row
            ev_in_scope.add(ev_id)
print(f"  Evidence records (membrane protein + has site): {len(ev_data):,}")

# Third pass: load site instances that have matching evidence
site_records = []
site_source_dbs = Counter()
site_types = Counter()
with gzip.open(BINDING_SITES, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['evidence_id'] in ev_in_scope:
            site_records.append(row)
            site_source_dbs[row.get('source_database', '')] += 1
            site_types[row.get('site_type', '')] += 1
print(f"  Site instances matched: {len(site_records):,}")
print(f"  Source databases: {dict(site_source_dbs)}")
print(f"  Site types: {dict(site_types)}")

# ============================================================
# Step 4: Merge into output rows
# ============================================================
print("\n[4/5] Merging compound + protein + evidence + site...")

bio_flags = ['is_approved_drug', 'is_clinical_candidate', 'is_endogenous_ligand',
             'is_natural_product', 'is_chemical_probe']

output_rows = []
skipped_no_cpd = 0
skipped_no_prot = 0

for site in site_records:
    ev_id = site['evidence_id']
    ev = ev_data.get(ev_id)
    if not ev:
        continue

    cpd_id = site.get('compound_internal_id', ev.get('compound_internal_id', '')).strip()
    up_id  = site.get('target_uniprot_id', ev.get('target_uniprot_id', '')).strip()

    cpd = compounds.get(cpd_id)
    prot = proteins.get(up_id)

    if not cpd:
        skipped_no_cpd += 1
        continue
    if not prot:
        skipped_no_prot += 1
        continue

    active_flags = ';'.join(f for f in bio_flags if cpd.get(f, '0') == '1')

    row = {
        # ── Compound identity ──
        'compound_internal_id': cpd_id,
        'preferred_name': cpd.get('preferred_name', ''),
        'standard_smiles': cpd.get('standard_smiles', ''),
        'standard_inchikey': cpd.get('standard_inchikey', ''),
        'molecular_formula': cpd.get('molecular_formula', ''),
        'molecular_weight': cpd.get('molecular_weight', ''),
        'xlogp': cpd.get('xlogp', ''),
        'tpsa': cpd.get('tpsa', ''),
        'pubchem_cids': cpd.get('pubchem_cids', ''),
        'chembl_ids': cpd.get('chembl_ids', ''),
        'drugcentral_ids': cpd.get('drugcentral_ids', ''),
        'chebi_ids': cpd.get('chebi_ids', ''),
        'gtopdb_ligand_ids': cpd.get('gtopdb_ligand_ids', ''),
        'cas_numbers': cpd.get('cas_numbers', ''),
        # ── Compound classification ──
        'compound_scope_class': cpd.get('compound_scope_class', ''),
        'molecule_types': cpd.get('molecule_types', ''),
        'identity_confidence': cpd.get('identity_confidence', ''),
        # ── Biological status flags ──
        'biological_status_flags': active_flags,
        'is_approved_drug': cpd.get('is_approved_drug', '0'),
        'is_clinical_candidate': cpd.get('is_clinical_candidate', '0'),
        'is_endogenous_ligand': cpd.get('is_endogenous_ligand', '0'),
        'is_natural_product': cpd.get('is_natural_product', '0'),
        'is_chemical_probe': cpd.get('is_chemical_probe', '0'),
        'development_status': cpd.get('development_status', ''),
        'max_chembl_phase': cpd.get('max_chembl_phase', ''),
        'first_approval_year': cpd.get('first_approval_year', ''),
        'compound_source_databases': cpd.get('source_databases', ''),
        # ── Protein identity ──
        'target_uniprot_id': up_id,
        'approved_symbol': prot.get('approved_symbol', ''),
        'protein_name': prot.get('protein_name', ''),
        'gene_names': prot.get('gene_names', ''),
        'membrane_class_v52': prot.get('membrane_class_v52', ''),
        'membrane_topology': prot.get('membrane_topology', ''),
        'transmembrane_count': prot.get('transmembrane_count', ''),
        'functional_primary_class': prot.get('functional_primary_class', ''),
        'functional_subclass': prot.get('functional_subclass', ''),
        'evidence_level_v52': prot.get('evidence_level_v52', ''),
        'membrane_scope': prot.get('membrane_scope', ''),
        'subcellular_location': prot.get('subcellular_location', ''),
        # ── Evidence info ──
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
        # ── Binding site residue detail ──
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

print(f"  Output rows: {len(output_rows):,}")
print(f"  Skipped (no compound): {skipped_no_cpd}")
print(f"  Skipped (no protein): {skipped_no_prot}")

# ============================================================
# Step 5: Write output
# ============================================================
print(f"\n[5/5] Writing {OUT_FILE}...")

fieldnames = list(output_rows[0].keys())

with open(OUT_FILE, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(output_rows)

file_size_mb = os.path.getsize(OUT_FILE) / (1024 * 1024)

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print(f"EXPORT COMPLETE")
print(f"{'='*60}")
print(f"  File: {OUT_FILE}")
print(f"  Size: {file_size_mb:.1f} MB")
print(f"  Rows: {len(output_rows):,} x {len(fieldnames)} columns")
print(f"")

# Unique counts
uniq_cpds = set(r['compound_internal_id'] for r in output_rows)
uniq_prots = set(r['target_uniprot_id'] for r in output_rows)
uniq_site_pdbs = set()
for r in output_rows:
    if r['pdb_ids_site']:
        for p in r['pdb_ids_site'].split(','):
            p = p.strip()
            if p:
                uniq_site_pdbs.add(p)
uniq_ev_pdbs = set()
for r in output_rows:
    if r['pdb_ids_evidence']:
        for p in r['pdb_ids_evidence'].split(';'):
            p = p.strip()
            if p:
                uniq_ev_pdbs.add(p)

print(f"  Unique compounds:     {len(uniq_cpds):,}")
print(f"  Unique proteins:      {len(uniq_prots):,}")
print(f"  Unique PDBs (site):   {len(uniq_site_pdbs):,}")
print(f"  Unique PDBs (evidence): {len(uniq_ev_pdbs):,}")

# Bio-status flag breakdown
bio_cpd_set = set()
for r in output_rows:
    if r['biological_status_flags']:
        bio_cpd_set.add(r['compound_internal_id'])
print(f"  With bio-status:      {len(bio_cpd_set):,} compounds")
print(f"  Without bio-status:   {len(uniq_cpds - bio_cpd_set):,} compounds")

# Evidence tier
et = Counter(r['evidence_tier'] for r in output_rows)
print(f"\n  Evidence tiers:")
for k, v in et.most_common():
    print(f"    {k}: {v:,}")

# Site types
st = Counter(r['site_type'] for r in output_rows)
print(f"\n  Site types:")
for k, v in st.most_common():
    print(f"    {k}: {v:,}")

# Residue + PDB breakdown
has_both = sum(1 for r in output_rows if r['residue_or_site_description'].strip() and r['pdb_ids_site'].strip())
has_res  = sum(1 for r in output_rows if r['residue_or_site_description'].strip() and not r['pdb_ids_site'].strip())
has_pdb  = sum(1 for r in output_rows if not r['residue_or_site_description'].strip() and r['pdb_ids_site'].strip())
print(f"\n  Residue + PDB status:")
print(f"    Both: {has_both:,}")
print(f"    Residue only (UniProt): {has_res:,}")
print(f"    PDB only (BindingDB):   {has_pdb:,}")

# Membrane class breakdown
mc = Counter(r['membrane_class_v52'] for r in output_rows)
print(f"\n  Membrane class:")
for k, v in mc.most_common():
    print(f"    Class {k}: {v:,}")

# Top targets
tgt = Counter(r['approved_symbol'] for r in output_rows)
print(f"\n  Top 15 targets:")
for sym, cnt in tgt.most_common(15):
    print(f"    {sym}: {cnt:,}")

# Source databases
src = Counter(r['site_source_database'] for r in output_rows)
print(f"\n  Site source databases:")
for k, v in src.most_common():
    print(f"    {k}: {v:,}")

print(f"\nDone.")
