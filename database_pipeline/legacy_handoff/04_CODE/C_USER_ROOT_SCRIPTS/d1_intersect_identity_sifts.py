#!/usr/bin/env python3
import csv,json
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
with (ROOT/'D1_identity_whitelist_v7_2_1.tsv').open(encoding='utf-8',newline='') as f:
    allowed={r['target_uniprot'] for r in csv.DictReader(f,delimiter='\t')}
src=ROOT/'docking_manifest_D1_sifts_qc.tsv'
with src.open(encoding='utf-8',newline='') as f:
    rd=csv.DictReader(f,delimiter='\t'); fields=rd.fieldnames; rows=list(rd)
extra=['d1_combined_qc_status','d1_combined_qc_reason']
out=[]
for r in rows:
    reasons=[]
    if r['target_uniprot'] not in allowed: reasons.append('MEMBRANE_IDENTITY_NOT_WHITELISTED')
    if r.get('sifts_target_present')!='1': reasons.append('NO_SIFTS_TARGET_MAPPING')
    if r.get('sifts_chain_status') not in {'UNIQUE','MULTIPLE_EQUIVALENT'}: reasons.append('CHAIN_MAPPING_UNRESOLVED')
    try:
        req=int(r.get('requested_residue_count') or 0); cov=int(r.get('sifts_residue_range_covered_count') or 0)
    except ValueError: req=cov=0
    if req==0: reasons.append('NO_PARSEABLE_BINDING_RESIDUES')
    elif cov<req: reasons.append('RESIDUES_OUTSIDE_SIFTS_RANGE')
    x=dict(r);x['d1_combined_qc_status']='PASS_TO_COORDINATE_QC' if not reasons else 'HOLD';x['d1_combined_qc_reason']=';'.join(reasons);out.append(x)
with (ROOT/'docking_manifest_D1_identity_sifts_pass.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields+extra,delimiter='\t');w.writeheader();w.writerows(r for r in out if r['d1_combined_qc_status']=='PASS_TO_COORDINATE_QC')
with (ROOT/'docking_manifest_D1_identity_sifts_hold.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields+extra,delimiter='\t');w.writeheader();w.writerows(r for r in out if r['d1_combined_qc_status']=='HOLD')
summary={'input':len(out),'pass_to_coordinate_qc':sum(r['d1_combined_qc_status']=='PASS_TO_COORDINATE_QC' for r in out),'hold':sum(r['d1_combined_qc_status']=='HOLD' for r in out),'hold_reasons':dict(Counter(x for r in out for x in r['d1_combined_qc_reason'].split(';') if x))}
(ROOT/'metadata'/'D1_IDENTITY_SIFTS_INTERSECTION_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
