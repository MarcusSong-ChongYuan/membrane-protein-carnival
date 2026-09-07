#!/usr/bin/env python3
import csv, gzip, json, re
from pathlib import Path
from collections import Counter

ROOT=Path(r'C:\Users\Administrator')
SEED=ROOT/'MemPro_Docking_D1_AE1_20260813'/'membrane_identity_audit_seed_v7_2.tsv'
MASTER=Path(r'D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay\human_membrane_protein_master_v6_4_candidate.tsv.gz')
OUT=ROOT/'MemPro_V7_2_1_membrane_identity_audit_20260813'; OUT.mkdir(exist_ok=True)
with SEED.open(encoding='utf-8',newline='') as f: seed=list(csv.DictReader(f,delimiter='\t'))
with gzip.open(MASTER,'rt',encoding='utf-8',newline='') as f: master={r['target_uniprot_id']:r for r in csv.DictReader(f,delimiter='\t')}

fields=list(seed[0])+['auto_conflict_flags','auto_membrane_identity_status','auto_docking_identity_eligibility','auto_decision_basis']
rows=[]
for s in seed:
    p=master[s['target_uniprot']]; flags=[]; cl=p.get('membrane_class_v52',''); ev=p.get('evidence_level_v52','')
    loc=(p.get('subcellular_location') or '').lower(); basis=(p.get('membrane_evidence_basis') or '').lower(); decision=p.get('membrane_decision_v5','')
    tm=int(float(p.get('transmembrane_count_v5') or 0)); hpa=p.get('hpa_plasma_membrane_location_v5')=='1'
    if cl=='A' and tm==0: flags.append('A_CLASS_WITH_ZERO_TM_FEATURES')
    if cl=='A' and any(x in loc for x in ['nucleus','nucleolus']) and 'membrane' not in loc: flags.append('A_CLASS_VS_NONMEMBRANE_UNIPROT_LOCATION')
    if decision in {'new_hpa_location_only','new_hpa_predicted_membrane_only'}: flags.append('HPA_ONLY_MEMBRANE_ASSIGNMENT')
    if cl=='C' and hpa and 'membrane' not in loc: flags.append('HPA_VS_UNIPROT_LOCATION_CONFLICT')
    if (p.get('cross_source_conflict_v5') or '').strip().lower() in {'1','true','yes','y'}: flags.append('EXISTING_CROSS_SOURCE_CONFLICT')
    if ev=='E2': flags.append('E2_REVIEW_PRIORITY')
    if cl in {'B','C'}: flags.append('B_OR_C_SCOPE_REVIEW')
    severe=any(x in flags for x in ['A_CLASS_WITH_ZERO_TM_FEATURES','A_CLASS_VS_NONMEMBRANE_UNIPROT_LOCATION','HPA_ONLY_MEMBRANE_ASSIGNMENT','HPA_VS_UNIPROT_LOCATION_CONFLICT','EXISTING_CROSS_SOURCE_CONFLICT'])
    if severe: status='CONFLICT_REVIEW'; eligible='HOLD'
    elif cl=='A' and ev=='E1' and tm>0: status='AUTO_CONFIRMED_A_E1'; eligible='D1_ELIGIBLE_IDENTITY_ONLY'
    else: status='REVIEW_REQUIRED'; eligible='HOLD'
    x=dict(s); x.update({'auto_conflict_flags':';'.join(flags),'auto_membrane_identity_status':status,'auto_docking_identity_eligibility':eligible,'auto_decision_basis':f'class={cl};evidence={ev};tm={tm};decision={decision}'})
    rows.append(x)
with (OUT/'membrane_identity_audit_auto_v7_2.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
with (OUT/'membrane_identity_conflict_review_v7_2.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(r for r in rows if r['auto_membrane_identity_status']=='CONFLICT_REVIEW')
summary={'targets':len(rows),'status':dict(Counter(r['auto_membrane_identity_status'] for r in rows)),'eligibility':dict(Counter(r['auto_docking_identity_eligibility'] for r in rows)),'flags':dict(Counter(x for r in rows for x in r['auto_conflict_flags'].split(';') if x))}
(OUT/'MEMBRANE_IDENTITY_AUTO_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))

