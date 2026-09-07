#!/usr/bin/env python3
import csv,json
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
SRC=ROOT/"02_protein_audit"/"final_adjudication_v7_2"/"protein_identity_final_adjudication_all_v7_2.tsv"
VAL=ROOT/"02_protein_audit"/"t09_validation"/"T09_MECHANISTIC_C_OLD_E1_VALIDATION.tsv"
OUT=ROOT/"02_protein_audit"/"final_adjudication_v7_3";OUT.mkdir(exist_ok=True)
with VAL.open(encoding='utf-8',newline='') as f:val={r['target_uniprot']:r for r in csv.DictReader(f,delimiter='\t')}
with SRC.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;rows=list(rd)
extra=['t09_applied','final_class_v7_3','final_mode_v7_3','final_decision_v7_3','final_evidence_v7_3','public_inclusion_v7_3','final_validation_status_v7_3']
out=[];changes=[]
for r in rows:
 cls=r['final_class_v7_2'];mode=r['final_mode_v7_2'];dec=r['final_decision_v7_2'];ev=r['final_evidence_v7_2'];applied=0
 v=val.get(r['target_uniprot'])
 if v:
  d=v['t09_final_validation_decision']
  if d in {'REVERT_EXCLUDED_NO_MECHANISM','CONFIRM_EXCLUSION'}:
   ncls='EXCLUDED';nmode='no_validated_membrane_binding_mechanism';ndec='EXCLUDED_NO_MEMBRANE_MECHANISM';nev=v['t09_final_validation_reason']
  elif d in {'CONFIRM_C','RESCUE_C'}:
   ncls='C';nmode='validated_direct_or_state_dependent_peripheral';ndec='CONFIRMED_C_VALIDATED';nev=v['t09_final_validation_reason']+' '+v['t09_go_lipid_direct']+' '+v['t09_pubmed_membrane_mechanism_text']
  elif d in {'CONFIRM_C_COMPLEX_MEDIATED','RESCUE_C_COMPLEX_MEDIATED'}:
   ncls='C';nmode='validated_complex_mediated_membrane_association';ndec='CONFIRMED_C_COMPLEX_MEDIATED_VALIDATED';nev=v['t09_final_validation_reason']+' '+v['complexportal_membrane_complex_ids_t09']
  else:ncls,nmode,ndec,nev=cls,mode,dec,ev
  if (ncls,nmode,ndec)!=(cls,mode,dec):
   applied=1;changes.append({'target_uniprot':r['target_uniprot'],'approved_symbol':r['approved_symbol'],'old_class':cls,'new_class':ncls,'old_decision':dec,'new_decision':ndec,'validation_decision':d,'reason':v['t09_final_validation_reason']})
  cls,mode,dec,ev=ncls,nmode,ndec,nev
 inc=int(cls in {'A','B','C'} and r['entity_type_v7_2']=='canonical_protein')
 x=dict(r);x.update({'t09_applied':applied,'final_class_v7_3':cls,'final_mode_v7_3':mode,'final_decision_v7_3':dec,'final_evidence_v7_3':ev,'public_inclusion_v7_3':inc,'final_validation_status_v7_3':'T09_APPLIED_NO_UNRESOLVED_STATUS'});out.append(x)
fields=base+extra
with (OUT/'protein_identity_final_adjudication_all_v7_3.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
for name,vals in [('protein_ABC_release_candidate_v7_3.tsv',[x for x in out if x['public_inclusion_v7_3']==1]),('protein_excluded_audit_v7_3.tsv',[x for x in out if x['public_inclusion_v7_3']==0])]:
 with (OUT/name).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(vals)
cf=list(changes[0]) if changes else []
with (OUT/'T09_CLASSIFICATION_CHANGELOG.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=cf,delimiter='\t');w.writeheader();w.writerows(changes)
summary={'input':len(out),'class':dict(Counter(x['final_class_v7_3'] for x in out)),'decision':dict(Counter(x['final_decision_v7_3'] for x in out)),'ABC_release_candidate':sum(x['public_inclusion_v7_3'] for x in out),'excluded_or_other_entity':sum(1-x['public_inclusion_v7_3'] for x in out),'t09_changes':len(changes),'unresolved':sum(x['final_class_v7_3'] in {'','unknown','UNRESOLVED'} for x in out),'status':'V7_3_T09_APPLIED_PENDING_FINAL_STRATIFIED_SAMPLE'}
(OUT/'V7_3_T09_APPLIED_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
