"""Backfill the 3 CAS-resolved SID->CID records."""
import json, subprocess, csv, gzip, shutil, os

OUT_DIR = r'D:\finale\negative_upgrade'
MASTER_DIR = r'D:\finale\01_正式数据_V6.2'
MASTER = os.path.join(MASTER_DIR, 'small_molecule_master_v1_3.tsv')
FORM_MASTER = os.path.join(MASTER_DIR, 'compound_form_hierarchy_v1_3.tsv')
UNMAPPED = os.path.join(MASTER_DIR, 'negative_binding_evidence_unmapped_review_v1_2.tsv.gz')
NEG_EVIDENCE = os.path.join(MASTER_DIR, 'negative_binding_evidence_v1_2.tsv.gz')

RELEASE_DATE = '2026-08-07'

# 1. Query 3 CIDs from PubChem
new_cids = ['19828004', '73894736', '17388989']
props_endpoint = ("IsomericSMILES,InChI,InChIKey,MolecularWeight,MolecularFormula,"
                  "IUPACName,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
                  "RotatableBondCount,HeavyAtomCount,Charge,ExactMass")

props = {}
print("=== Querying 3 new CIDs ===")
for cid in new_cids:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{props_endpoint}/JSON"
    r = subprocess.run(["curl", "-s", "--connect-timeout", "10", "--max-time", "30", url],
        capture_output=True, text=True, timeout=35)
    try:
        d = json.loads(r.stdout)
        plist = d.get('PropertyTable', {}).get('Properties', [])
        if plist:
            entry = plist[0]
            props[cid] = {
                'SMILES': entry.get('SMILES', '') or '',
                'InChI': entry.get('InChI', '') or '',
                'InChIKey': entry.get('InChIKey', '') or '',
                'MolecularWeight': entry.get('MolecularWeight', ''),
                'MolecularFormula': entry.get('MolecularFormula', '') or '',
                'IUPACName': entry.get('IUPACName', '') or '',
                'XLogP': entry.get('XLogP', ''),
                'TPSA': entry.get('TPSA', ''),
                'HBondDonorCount': entry.get('HBondDonorCount', ''),
                'HBondAcceptorCount': entry.get('HBondAcceptorCount', ''),
                'RotatableBondCount': entry.get('RotatableBondCount', ''),
                'HeavyAtomCount': entry.get('HeavyAtomCount', ''),
                'Charge': entry.get('Charge', ''),
                'ExactMass': entry.get('ExactMass', ''),
            }
            print(f"  CID {cid}: {props[cid]['IUPACName'][:80]}")
        else:
            print(f"  CID {cid}: NOT FOUND")
    except Exception as e:
        print(f"  CID {cid}: ERROR {e}")

# 2. Find next IDs
last_cpd_num = 0
with open(MASTER, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        num = int(row['compound_internal_id'].split('-')[-1])
        if num > last_cpd_num:
            last_cpd_num = num
print(f"\nLast compound: HMPD-CMPD-{last_cpd_num:07d}")

last_form_num = 0
with open(FORM_MASTER, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        num = int(row['compound_form_id'].split('-')[-1])
        if num > last_form_num:
            last_form_num = num
print(f"Last form: HMPD-FORM-{last_form_num:07d}")

# 3. Build new compound and form rows
def safe_str(v):
    if v is None or v == '':
        return ''
    return str(v)

cpd_num = last_cpd_num + 1
form_num = last_form_num + 1
new_compounds = []
new_forms = []
cid_to_cpd = {}

for cid in new_cids:
    if cid not in props:
        continue
    p = props[cid]
    cpd_id = f'HMPD-CMPD-{cpd_num:07d}'
    form_id = f'HMPD-FORM-{form_num:07d}'
    cid_to_cpd[cid] = cpd_id

    inchikey = p.get('InChIKey', '') or ''

    cpd_row = {
        'compound_internal_id': cpd_id,
        'preferred_name': p.get('IUPACName', '') or '',
        'preferred_name_source': 'PubChem CID',
        'compound_scope_status': 'core',
        'compound_scope_class': 'core_small_molecule',
        'identity_confidence': 'low',
        'standard_smiles': p.get('SMILES', '') or '',
        'standard_inchi': p.get('InChI', '') or '',
        'standard_inchikey': inchikey,
        'identity_key': f'IK:{inchikey}' if inchikey else '',
        'connectivity_key': inchikey[:14] if inchikey else '',
        'molecular_formula': p.get('MolecularFormula', '') or '',
        'molecular_weight': safe_str(p.get('MolecularWeight', '')),
        'exact_mass': safe_str(p.get('ExactMass', '')),
        'xlogp': safe_str(p.get('XLogP', '')),
        'tpsa': safe_str(p.get('TPSA', '')),
        'hbond_donor_count': safe_str(p.get('HBondDonorCount', '')),
        'hbond_acceptor_count': safe_str(p.get('HBondAcceptorCount', '')),
        'rotatable_bond_count': safe_str(p.get('RotatableBondCount', '')),
        'formal_charge': safe_str(p.get('Charge', '')),
        'heavy_atom_count': safe_str(p.get('HeavyAtomCount', '')),
        'ring_count': '', 'max_ring_size': '', 'amide_bond_count': '',
        'computed_structural_class': '', 'molecule_types': '',
        'form_count': '1', 'source_record_count': '1',
        'source_database_count': '1', 'source_databases': 'PubChem BioAssay',
        'pubchem_cids': cid,
        'chembl_ids': '', 'chebi_ids': '', 'drugcentral_ids': '',
        'gtopdb_ligand_ids': '', 'hmdb_ids': '', 'cas_numbers': '',
        'synonym_count': '0', 'development_status': '',
        'max_chembl_phase': '', 'first_approval_year': '',
        'is_approved_drug': '0', 'is_clinical_candidate': '0',
        'is_endogenous_ligand': '0', 'is_natural_product': '0',
        'is_chemical_probe': '0',
        'compound_category_tags': 'negative_evidence_compound',
        'chebi_direct_classes': '', 'chebi_roles': '',
        'record_qc_status': 'ok',
        'qc_notes': 'auto_added_via_cas_synonym_resolution',
        'protein_target_count': '0', 'best_binding_evidence_level': '',
        'BE1_evidence_count': '0', 'BE2_evidence_count': '0',
        'BE3_evidence_count': '0', 'binding_evidence_count': '0',
        'binding_source_count': '0', 'binding_sources': '',
        'best_standard_value_nM': '',
        'release_version': 'v1.3', 'release_date': RELEASE_DATE,
        'standardization_policy': 'pubchem_properties_only',
        'incremental_evidence_run_v11': '', 'identity_refresh_v12': '',
        'new_source_record_count_v12': '', 'release_version_v12': '',
        'release_date_v12': '',
        'negative_evidence_count_v13': '1', 'negative_target_count_v13': '',
        'positive_negative_conflict_pair_count_v13': '0',
        'negative_identity_refresh_v13': 'cas_synonym_20260807',
        'release_version_v13': 'v1.3', 'release_date_v13': RELEASE_DATE,
    }
    new_compounds.append(cpd_row)

    form_row = {
        'compound_form_id': form_id,
        'compound_internal_id': cpd_id,
        'form_type': 'parent_form',
        'exact_smiles': p.get('SMILES', '') or '',
        'exact_inchi': p.get('InChI', '') or '',
        'exact_inchikey': inchikey,
        'exact_identity_key': f'IK:{inchikey}' if inchikey else '',
        'molecular_formula': p.get('MolecularFormula', '') or '',
        'molecular_weight': safe_str(p.get('MolecularWeight', '')),
        'formal_charge': safe_str(p.get('Charge', '')),
        'source_record_count': '1', 'source_databases': 'PubChem BioAssay',
        'source_compound_ids': f'PubChem CID:{cid}',
        'chembl_parent_ids': '', 'chembl_parent_alignment_status': '',
        'record_qc_status': 'ok',
        'qc_notes': 'auto_added_via_cas_synonym_resolution',
        'release_version_v11': '',
        'parent_inchikey_v12': inchikey,
        'parent_relation_v12': 'self_parent',
    }
    new_forms.append(form_row)
    cpd_num += 1
    form_num += 1

# 4. Align with master columns and append
print(f"\nAppending {len(new_compounds)} to master tables...")

with open(MASTER, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    master_fields = reader.fieldnames
filtered_cpds = [{k: r.get(k, '') for k in master_fields} for r in new_compounds]
with open(MASTER, 'a', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=master_fields, delimiter='\t', extrasaction='ignore')
    writer.writerows(filtered_cpds)

with open(FORM_MASTER, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    form_fields = reader.fieldnames
filtered_forms = [{k: r.get(k, '') for k in form_fields} for r in new_forms]
with open(FORM_MASTER, 'a', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=form_fields, delimiter='\t', extrasaction='ignore')
    writer.writerows(filtered_forms)

print(f"  Appended {len(filtered_cpds)} compounds, {len(filtered_forms)} forms")

# 5. Backfill unmapped records
print(f"\nBackfilling unmapped records...")
with open(os.path.join(OUT_DIR, 'sid_to_cid_via_synonyms.json')) as f:
    sid_to_cid = json.load(f)

newly_mapped = []
still_unmapped = []
with gzip.open(UNMAPPED, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = reader.fieldnames
    for row in reader:
        sid = row.get('pubchem_sid', '').strip()
        if sid in sid_to_cid:
            cid = sid_to_cid[sid]
            if cid in cid_to_cpd:
                row['pubchem_cid'] = cid
                row['compound_source_id'] = f'PubChem CID:{cid}'
                row['compound_internal_id'] = cid_to_cpd[cid]
                row['compound_form_id'] = ''
                row['compound_mapping_status'] = 'mapped_via_cas_synonym'
                row['release_action_v62'] = 'synonym_resolved_cas'
                row['unmapped_reason_v62'] = ''
                newly_mapped.append(row)
                continue
        still_unmapped.append(row)

print(f"  Newly mapped: {len(newly_mapped)}")
print(f"  Still unmapped: {len(still_unmapped)}")

# Save updated unmapped_review
shutil.copy(UNMAPPED, os.path.join(OUT_DIR, 'backup', 'unmapped_review_final.tsv.gz'))
with gzip.open(UNMAPPED, 'wt', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(still_unmapped)

# 6. Append to negative evidence
print(f"\nAppending to negative evidence...")
tmp = os.path.join(OUT_DIR, 'neg_final_cas.tsv')
with gzip.open(NEG_EVIDENCE, 'rt', encoding='utf-8') as fin:
    reader = csv.DictReader(fin, delimiter='\t')
    neg_fields = reader.fieldnames
    with open(tmp, 'w', encoding='utf-8', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=neg_fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        cnt = 0
        for row in reader:
            writer.writerow(row)
            cnt += 1
        print(f"  Copied {cnt:,} existing")
        for row in newly_mapped:
            clean = {k: row.get(k, '') for k in neg_fields}
            writer.writerow(clean)
        print(f"  Appended {len(newly_mapped)}")

subprocess.run(['gzip', '-f', tmp], check=True)
shutil.move(tmp + '.gz', NEG_EVIDENCE)

with gzip.open(NEG_EVIDENCE, 'rt', encoding='utf-8') as f:
    final = sum(1 for _ in f) - 1

print(f"\n=== FINAL STATE ===")
print(f"  Negative evidence total: {final:,}")
print(f"  Unmapped review remaining: {len(still_unmapped):,}")
print(f"Done.")
