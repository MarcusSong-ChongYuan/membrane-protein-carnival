#!/usr/bin/env python3
import csv,gzip,json
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
PARSED=ROOT/"02_protein_audit"/"protein_identity_audit_all_v7_uniprot_parsed.tsv"
MASTER=Path(r"D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay\human_membrane_protein_master_v6_4_candidate.tsv.gz")
OUT=ROOT/"02_protein_audit"/"protein_identity_audit_all_v7_cross_source_retriage.tsv"
def yes(v):return str(v).strip().lower() in {'1','true','yes','y'}
def n(v):
 try:return int(float(v or 0))
 except:return 0
with gzip.open(MASTER,'rt',encoding='utf-8',newline='') as f:master={r['target_uniprot_id']:r for r in csv.DictReader(f,delimiter='\t')}
with PARSED.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;rows=list(rd)
extra=['external_integral_supports_v7','external_localization_supports_v7','cross_source_proposed_class_v7','cross_source_proposed_evidence_v7','cross_source_priority_v7','cross_source_flags_v7','cross_source_decision_status_v7']
out=[]
for r in rows:
 m=master[r['target_uniprot']]; integ=[]; loc=[]
 if yes(m.get('pdbtm_present_v5')):integ.append('PDBTM')
 if yes(m.get('opm_integral_structure_v5')):integ.append('OPM_integral')
 if 'beta_barrel' in m.get('membrane_evidence_basis','').lower():integ.append('UniProt_beta_barrel_legacy_curated')
 if yes(m.get('htp_present_v5')) and n(m.get('htp_num_tm_v5'))>0:integ.append('HTP_TM')
 if yes(m.get('membranome_present_v5')) and m.get('membranome_segment_v5','').strip():integ.append('Membranome_segment')
 if yes(m.get('hpa_plasma_membrane_location_v5')):loc.append('HPA_plasma_membrane')
 if yes(m.get('hpa_predicted_membrane_v5')):loc.append('HPA_predicted_membrane')
 if yes(m.get('opm_present_v5')) and not yes(m.get('opm_integral_structure_v5')):loc.append('OPM_peripheral_structure')
 ucls=r['proposed_primary_class_v7']; cls=ucls; flags=[]
 if ucls in {'UNRESOLVED','EXCLUDE_CANDIDATE'} and integ:
  cls='A';flags.append('RESCUED_BY_EXTERNAL_INTEGRAL_EVIDENCE')
 elif ucls in {'UNRESOLVED','EXCLUDE_CANDIDATE'} and 'OPM_peripheral_structure' in loc:
  cls='C';flags.append('RESCUED_AS_C_BY_OPM_PERIPHERAL_STRUCTURE')
 elif ucls in {'UNRESOLVED','EXCLUDE_CANDIDATE'} and loc:
  cls='UNRESOLVED';flags.append('LOCALIZATION_ONLY_NO_BINDING_MODE')
 elif ucls=='EXCLUDE_CANDIDATE':flags.append('NO_UNIPROT_OR_EXTERNAL_MEMBRANE_MODE')
 if r['current_class']!=cls and cls in {'A','B','C'}:flags.append('ABC_CHANGE_PROPOSED')
 if r['cross_source_conflict']=='1':flags.append('LEGACY_CROSS_SOURCE_CONFLICT')
 # Evidence is conservative: structure/topology can be E1; HTP/Membranome prediction alone is E2.
 if r['proposed_evidence_level_v7']=='E1':ev='E1'
 elif any(x in integ for x in ['PDBTM','OPM_integral']):ev='E1'
 elif cls in {'A','B','C'}:ev='E2'
 else:ev='E3'
 if cls=='EXCLUDE_CANDIDATE':prio='P0'
 elif cls=='UNRESOLVED' or 'LEGACY_CROSS_SOURCE_CONFLICT' in flags:prio='P1'
 elif 'ABC_CHANGE_PROPOSED' in flags or ev=='E2':prio='P2'
 elif ev=='E3':prio='P3'
 else:prio='P4_CONTROL'
 x=dict(r);x.update({'external_integral_supports_v7':';'.join(integ),'external_localization_supports_v7':';'.join(loc),'cross_source_proposed_class_v7':cls,'cross_source_proposed_evidence_v7':ev,'cross_source_priority_v7':prio,'cross_source_flags_v7':';'.join(flags),'cross_source_decision_status_v7':'CROSS_SOURCE_PROPOSAL_PENDING_HIGH_RISK_REVIEW'});out.append(x)
with OUT.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=base+extra,delimiter='\t');w.writeheader();w.writerows(out)
summary={'input':len(out),'class':dict(Counter(x['cross_source_proposed_class_v7'] for x in out)),'evidence':dict(Counter(x['cross_source_proposed_evidence_v7'] for x in out)),'priority':dict(Counter(x['cross_source_priority_v7'] for x in out)),'flags':dict(Counter(y for x in out for y in x['cross_source_flags_v7'].split(';') if y)),'status':'CROSS_SOURCE_RETRIAGE_COMPLETE_HIGH_RISK_REVIEW_PENDING'}
(ROOT/'02_protein_audit'/'CROSS_SOURCE_RETRIAGE_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))


