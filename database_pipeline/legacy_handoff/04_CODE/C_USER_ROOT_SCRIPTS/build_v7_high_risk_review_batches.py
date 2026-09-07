#!/usr/bin/env python3
import csv,json
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
SRC=ROOT/"02_protein_audit"/"protein_identity_audit_all_v7_cross_source_retriage.tsv"
OUT=ROOT/"02_protein_audit"/"high_risk_batches";OUT.mkdir(exist_ok=True)
with SRC.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames;rows=list(rd)
def bucket(r):
 cls=r['cross_source_proposed_class_v7'];flags=r['cross_source_flags_v7'];old=r['current_class'];ev=r['current_evidence']
 if cls=='EXCLUDE_CANDIDATE':
  if r['signal_peptide_present']=='1' and old=='C':return 'P0_signal_or_secreted_without_membrane_mode'
  if r['hpa_only_or_single_source']=='1':return 'P0_HPA_only_without_membrane_mode'
  if old=='unknown':return 'P0_legacy_unknown_no_membrane_mode'
  if old=='C':return 'P0_legacy_C_no_explicit_peripheral_mode'
  if old=='B':return 'P0_legacy_B_no_lipid_anchor'
  return 'P0_other_no_membrane_mode'
 if cls=='UNRESOLVED':
  if r['hpa_only_or_single_source']=='1':return 'P1_HPA_localization_only'
  if r['cross_source_conflict']=='1':return 'P1_cross_source_conflict_unresolved'
  return 'P1_location_without_binding_mode'
 if 'LEGACY_CROSS_SOURCE_CONFLICT' in flags:return 'P1_confirmed_mode_but_legacy_conflict'
 if 'ABC_CHANGE_PROPOSED' in flags:return 'P2_ABC_reclassification'
 if ev=='E2':return 'P2_E2_confirmation'
 return 'lower_risk'
selected=[]
for r in rows:
 b=bucket(r)
 if b!='lower_risk':
  x=dict(r);x['review_batch_v7']=b;x['manual_decision']='PENDING';x['manual_decision_reason']='';x['reviewer']='';x['review_date']='';selected.append(x)
ofields=fields+['review_batch_v7','manual_decision','manual_decision_reason','reviewer','review_date']
for b in sorted({x['review_batch_v7'] for x in selected}):
 vals=[x for x in selected if x['review_batch_v7']==b]
 with (OUT/(b+'.tsv')).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=ofields,delimiter='\t');w.writeheader();w.writerows(vals)
with (OUT/'ALL_HIGH_RISK_REVIEW_QUEUE.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=ofields,delimiter='\t');w.writeheader();w.writerows(selected)
summary={'total_high_risk_or_reclassification':len(selected),'batch_counts':dict(Counter(x['review_batch_v7'] for x in selected)),'rules':'No automatic deletion. P0/P1 require source review; P2 changes require validation before merge.'}
(OUT/'HIGH_RISK_BATCH_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
