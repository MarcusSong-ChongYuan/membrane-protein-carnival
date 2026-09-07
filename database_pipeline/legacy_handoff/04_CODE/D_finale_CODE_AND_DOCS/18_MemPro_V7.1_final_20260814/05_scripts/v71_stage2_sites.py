from __future__ import annotations

import csv, gzip, json, re
from collections import Counter
from pathlib import Path

BASE=Path(r"D:\finale\16_MemPro_V7.0.2_final_20260814")
ROOT=Path(r"D:\finale\17_MemPro_V7.1_candidate_20260814")
OUT=ROOT/'01_release_tables'; QA=ROOT/'04_QA'

def readgz(p):
    with gzip.open(p,'rt',encoding='utf-8-sig',newline='') as f: yield from csv.DictReader(f,delimiter='\t')
def intval(x):
    try:return int(float(x or 0))
    except:return 0
def floatval(x):
    try:return float(x)
    except:return 0.0
def tokens(s):
    if not s:return []
    vals=[]
    for x in re.split(r'[;,|\s]+',s.strip()):
        x=x.strip()
        if x and x not in vals: vals.append(x)
    return vals
def first(*xs):
    return next((x for x in xs if str(x or '').strip()),'')

FIELDS=['site_record_id','evidence_id','target_uniprot_id','compound_internal_id','compound_form_id','source_database','source_record_status','source_site_type','source_pdb_ids','source_reported_residues','mapped_structure_residues','unmapped_source_residues','requested_residue_count','mapped_residue_count','observed_residue_count','mapping_fraction','site_residue_status','coordinate_mapping_status','coordinate_mapping_method','coordinate_confidence','selected_pdb_id','selected_chain','candidate_chains','chain_assignment_status','chain_selection_basis','structure_ligand_relationship','structure_ligand_basis','membrane_side','site_public_interaction_eligible','coordinate_docking_eligible','site_release_layer','limitations','release_version']

def counts(r):
    req=max(intval(r.get('uniprot_requested_count_v703')),intval(r.get('rescue_requested_count_v702')),intval(r.get('requested_residue_count_v7')))
    mapped=max(intval(r.get('uniprot_mapped_count_v703')),intval(r.get('rescue_mapped_count_v702')),intval(r.get('mapped_residue_count_v7')),intval(r.get('pocket_contact_count_v704')))
    obs=max(intval(r.get('uniprot_observed_count_v703')),intval(r.get('rescue_observed_count_v702')),intval(r.get('observed_residue_count_v7')),intval(r.get('pocket_contact_count_v704')))
    return req,mapped,obs

def classify(r):
    src=first(r.get('requested_uniprot_residue_tokens_v7'),r.get('residue_or_site_description'))
    mapped_text=first(r.get('pocket_contact_uniprot_residues_v706'),r.get('uniprot_mapped_pdb_residues_v703'),r.get('rescue_pdb_residues_v702'),r.get('mapped_pdb_residues_v7'))
    req,mapped,obs=counts(r)
    src_t=tokens(src); map_t=tokens(mapped_text)
    if not req and src_t: req=len(src_t)
    if not mapped and map_t: mapped=len(map_t)
    if not obs and mapped: obs=mapped
    frac=(mapped/req) if req else (1.0 if r.get('pocket_rescue_status_v704b') else 0.0)
    pocket=bool(r.get('pocket_rescue_status_v704b'))
    predicted='predict' in (r.get('site_type','')+r.get('coordinate_mapping_basis','')).lower()
    if predicted: residue_status='PREDICTED_POCKET'
    elif pocket: residue_status='RESIDUE_COMPLETE'
    elif src_t or req:
        if mapped and req and mapped>=req: residue_status='RESIDUE_COMPLETE'
        elif mapped: residue_status='RESIDUE_PARTIAL'
        else: residue_status='RESIDUE_SOURCE_ONLY'
    else:
        residue_status='NO_RESIDUE_REPORTED' if r.get('evidence_id') else 'NO_SITE_INFORMATION'

    selected_pdb=first(r.get('pocket_pdb_v704'),r.get('uniprot_rescue_pdb_v703'),r.get('rescue_pdb_v702'),r.get('actual_pdb_used_v7'))
    chains=[]
    for val in [r.get('pdb_chain_ids_v60',''),r.get('pocket_target_chains_v704',''),r.get('uniprot_rescue_chain_v703',''),r.get('rescue_chain_v702',''),r.get('actual_chain_used_v7','')]:
        for c in tokens(val):
            if c not in chains:chains.append(c)
    mapped_chain=first(r.get('pocket_target_chains_v704'),r.get('uniprot_rescue_chain_v703'),r.get('rescue_chain_v702'),r.get('actual_chain_used_v7'))
    mapped_chain_tokens=tokens(mapped_chain)
    if len(mapped_chain_tokens)==1:
        selected_chain=mapped_chain_tokens[0]; chain_status='UNIQUE_CHAIN'; chain_basis='coordinate_mapping_selected_unique_chain'
    elif len(chains)>1:
        selected_chain=''; chain_status='MULTIPLE_EQUIVALENT_CHAINS'; chain_basis='source_or_mapping_reports_multiple_candidate_chains_no_forced_selection'
    elif len(chains)==1:
        selected_chain=chains[0]; chain_status='UNIQUE_CHAIN'; chain_basis='single_source_candidate_chain'
    else:
        selected_chain=''; chain_status='CHAIN_NOT_REPORTED'; chain_basis='source_did_not_report_chain'

    if pocket:
        ligand_rel='EXACT_COCRYSTAL_LIGAND'; ligand_basis='exact_HET_on_SIFTS_target_chain_with_protein_contacts'
    elif r.get('site_compound_specificity','').lower() in {'compound_specific','specific'} and r.get('pdb_ids'):
        ligand_rel='LIGAND_IDENTITY_UNRESOLVED'; ligand_basis='compound_specific_site_but_coordinate_ligand_identity_not_reverified'
    elif r.get('pdb_ids'):
        ligand_rel='STRUCTURE_ONLY_NO_QUERY_LIGAND'; ligand_basis='PDB_association_does_not_by_itself_prove_query_ligand_is_present'
    else:
        ligand_rel='NOT_APPLICABLE'; ligand_basis='no_structure_reported'

    coord_ok=bool(selected_pdb and selected_chain and mapped>0 and obs>0)
    if coord_ok and frac>=0.999: coord_status='COORDINATE_COMPLETE'; conf='HIGH'
    elif coord_ok: coord_status='COORDINATE_PARTIAL'; conf='MEDIUM'
    elif selected_pdb or r.get('pdb_ids'): coord_status='SOURCE_STRUCTURE_ONLY'; conf='NONE'
    else: coord_status='NO_COORDINATE'; conf='NONE'
    unmapped=[]
    if src_t and map_t:
        ms=set(map_t); unmapped=[x for x in src_t if x not in ms]
    limitations=[]
    if residue_status=='NO_RESIDUE_REPORTED':limitations.append('source_did_not_report_residues')
    if residue_status=='RESIDUE_PARTIAL':limitations.append('only_subset_of_source_residues_mapped')
    if chain_status.startswith('MULTIPLE_'):limitations.append('chain_not_uniquely_assigned')
    if ligand_rel in {'STRUCTURE_ONLY_NO_QUERY_LIGAND','LIGAND_IDENTITY_UNRESOLVED'}:limitations.append('query_ligand_not_confirmed_in_structure')
    return {
      'site_record_id':r.get('binding_site_instance_id',''),'evidence_id':r.get('evidence_id',''),'target_uniprot_id':r.get('target_uniprot_id',''),'compound_internal_id':r.get('compound_internal_id',''),'compound_form_id':r.get('compound_form_id',''),'source_database':r.get('source_database',''),'source_record_status':r.get('record_qc_status',''),'source_site_type':r.get('site_type',''),'source_pdb_ids':first(r.get('pdb_ids_original_v64'),r.get('pdb_ids')),'source_reported_residues':src,'mapped_structure_residues':mapped_text,'unmapped_source_residues':';'.join(unmapped),'requested_residue_count':req,'mapped_residue_count':mapped,'observed_residue_count':obs,'mapping_fraction':f'{frac:.6f}','site_residue_status':residue_status,'coordinate_mapping_status':coord_status,'coordinate_mapping_method':first(r.get('pocket_rescue_rule_v704b'),r.get('uniprot_rescue_rule_v703'),r.get('rescue_rule_v702'),r.get('coordinate_mapping_basis'),'NONE'),'coordinate_confidence':conf,'selected_pdb_id':selected_pdb,'selected_chain':selected_chain,'candidate_chains':';'.join(chains),'chain_assignment_status':chain_status,'chain_selection_basis':chain_basis,'structure_ligand_relationship':ligand_rel,'structure_ligand_basis':ligand_basis,'membrane_side':r.get('membrane_side','') or 'UNKNOWN','site_public_interaction_eligible':'1','coordinate_docking_eligible':'1' if coord_ok else '0','site_release_layer':'PUBLIC_SITE_ASSERTION' if residue_status!='PREDICTED_POCKET' else 'PREDICTION_LAYER','limitations':';'.join(limitations),'release_version':'V7.1-candidate'}

def main():
    OUT.mkdir(parents=True,exist_ok=True);QA.mkdir(parents=True,exist_ok=True)
    public=BASE/'01_data'/'binding_site_coordinate_verified_v702.tsv.gz'
    nonpublic=BASE/'02_nonpublic_resolved_status'/'binding_site_coordinate_nonpublic_v702.tsv.gz'
    out=OUT/'interaction_site_v71.tsv.gz'; status=Counter(); coords=Counter(); chains=Counter(); ligand=Counter(); ids=set(); dup=0;n=0
    with gzip.open(out,'wt',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',extrasaction='ignore');w.writeheader()
        for p in [public,nonpublic]:
            for r in readgz(p):
                x=classify(r); n+=1
                if x['site_record_id'] in ids:dup+=1
                ids.add(x['site_record_id']); w.writerow(x)
                status[x['site_residue_status']]+=1;coords[x['coordinate_mapping_status']]+=1;chains[x['chain_assignment_status']]+=1;ligand[x['structure_ligand_relationship']]+=1
    report={'input_rows':n,'output_rows':n,'unique_site_ids':len(ids),'duplicate_site_ids':dup,'residue_status':dict(status),'coordinate_status':dict(coords),'chain_status':dict(chains),'structure_ligand_relationship':dict(ligand),'site_partition_pass':n==55610 and len(ids)==n}
    (QA/'V71_STAGE2_SITE_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
