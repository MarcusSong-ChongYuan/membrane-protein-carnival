#!/usr/bin/env python3
import csv,json,re
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
SRC=ROOT/"02_protein_audit"/"final_adjudication_v7"/"protein_identity_final_adjudication_all_v7.tsv"
OUT=ROOT/"02_protein_audit"/"final_adjudication_v7_2";OUT.mkdir(exist_ok=True)
with SRC.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;rows=list(rd)
seg_re=re.compile(r'^(TRA[VDJ]|TRB[VDJ]|TRG[VDJ]|TRD[VDJ]|IGH[VDJ]|IGK[VDJ]|IGL[VDJ])')
light_constant_re=re.compile(r'^(IGKC|IGLC\d*)$')
heavy_constant_re=re.compile(r'^IGH[AGDEMC]\d*$')
extra=['entity_type_v7_2','form_specificity_v7_2','final_class_v7_2','final_mode_v7_2','final_decision_v7_2','final_evidence_v7_2','public_inclusion_v7_2','validation_status_v7_2']
out=[];segments=[];forms=[]
for r in rows:
 sym=r['approved_symbol'];loc=r['uniprot_subcellular_locations_current'];low=loc.lower()
 cls=r['final_membrane_class_v7'];mode=r['final_membrane_mode_v7'];dec=r['final_decision_v7'];evi=r['final_decision_evidence_v7'];entity='canonical_protein';specific='canonical_or_state'
 if seg_re.match(sym):
  cls='EXCLUDED';mode='immune_receptor_gene_segment_not_complete_chain';dec='MOVED_TO_IMMUNE_GENE_SEGMENT_ENTITY';evi='V/D/J germline segment does not independently encode a complete membrane-anchored receptor chain';entity='immune_gene_segment';specific='not_applicable'
 elif light_constant_re.match(sym):
  cls='EXCLUDED';mode='immune_light_constant_segment_not_complete_chain';dec='MOVED_TO_IMMUNE_GENE_SEGMENT_ENTITY';evi='Light-chain constant-region entry is not an independently membrane-anchored complete receptor chain';entity='immune_constant_segment';specific='not_applicable'
 elif 'lipid-anchor' in low or 'lipid anchor' in low or 'gpi-anchor' in low or 'gpi anchor' in low:
  cls='B';mode='explicit_uniprot_covalent_lipid_anchor';dec='CONFIRMED_B';evi=loc;entity='canonical_protein';specific='membrane_form_specific' if any(x in low for x in ['secreted','cytoplasm','cytosol']) else 'membrane_bound_form'
 elif heavy_constant_re.match(sym) and int(r['uniprot_transmembrane_count_current'] or 0)>0:
  cls='A';mode='immunoglobulin_heavy_chain_membrane_isoform';dec='CONFIRMED_A_FORM_SPECIFIC';evi=r['uniprot_transmembrane_features_current']+'; '+loc;entity='canonical_protein';specific='membrane_isoform_and_secreted_isoform'
 elif cls in {'A','B','C'} and any(x in low for x in ['secreted','cytoplasm','cytosol']) and ('membrane' in low or r['uniprot_transmembrane_count_current']!='0' or r['uniprot_lipidation_count_current']!='0'):
  specific='membrane_state_or_form_specific'
 include=int(cls in {'A','B','C'} and entity=='canonical_protein')
 x=dict(r);x.update({'entity_type_v7_2':entity,'form_specificity_v7_2':specific,'final_class_v7_2':cls,'final_mode_v7_2':mode,'final_decision_v7_2':dec,'final_evidence_v7_2':evi,'public_inclusion_v7_2':include,'validation_status_v7_2':'RULE_ADJUDICATED_PENDING_STRATIFIED_VALIDATION'});out.append(x)
 if entity.startswith('immune_'):segments.append(x)
 if specific not in {'canonical_or_state','not_applicable','membrane_bound_form'}:
  forms.append({'target_uniprot':r['target_uniprot'],'approved_symbol':sym,'entity_type':entity,'form_specificity':specific,'membrane_class':cls,'membrane_mode':mode,'uniprot_isoform_ids':r['uniprot_isoform_ids_current'],'location_and_form_evidence':loc,'decision':dec})
fields=base+extra
with (OUT/'protein_identity_final_adjudication_all_v7_2.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
for name,vals in [('protein_ABC_release_candidate_v7_2.tsv',[x for x in out if x['public_inclusion_v7_2']==1]),('protein_excluded_audit_v7_2.tsv',[x for x in out if x['public_inclusion_v7_2']==0]),('immune_gene_segment_entity_v7_2.tsv',segments)]:
 with (OUT/name).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(vals)
ff=list(forms[0]) if forms else []
with (OUT/'protein_membrane_form_adjudication_v7_2.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=ff,delimiter='\t');w.writeheader();w.writerows(forms)
summary={'input':len(out),'final_class':dict(Counter(x['final_class_v7_2'] for x in out)),'decision':dict(Counter(x['final_decision_v7_2'] for x in out)),'entity_type':dict(Counter(x['entity_type_v7_2'] for x in out)),'form_specificity':dict(Counter(x['form_specificity_v7_2'] for x in out)),'ABC_release_candidate':sum(x['public_inclusion_v7_2'] for x in out),'excluded_or_moved':sum(1-x['public_inclusion_v7_2'] for x in out),'unresolved_count':sum(x['final_class_v7_2'] in {'','unknown','UNRESOLVED'} for x in out),'status':'V7_2_ADJUDICATION_COMPLETE_PENDING_VALIDATION'}
(OUT/'FINAL_ABC_ADJUDICATION_V7_2_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
