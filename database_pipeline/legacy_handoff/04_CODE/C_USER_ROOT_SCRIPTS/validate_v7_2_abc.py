#!/usr/bin/env python3
import csv,json,hashlib
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813\02_protein_audit\final_adjudication_v7_2")
SRC=ROOT/'protein_identity_final_adjudication_all_v7_2.tsv'
with SRC.open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
errors=[];warnings=[]
ids=[r['target_uniprot'] for r in rows]
if len(ids)!=10997:errors.append('ROW_COUNT_NOT_10997')
if len(set(ids))!=len(ids):errors.append('DUPLICATE_UNIPROT')
for r in rows:
 c=r['final_class_v7_2'];inc=r['public_inclusion_v7_2']=='1';mode=r['final_mode_v7_2'];entity=r['entity_type_v7_2']
 if c in {'','unknown','UNRESOLVED'}:errors.append('UNRESOLVED:'+r['target_uniprot'])
 if inc!=(c in {'A','B','C'} and entity=='canonical_protein'):errors.append('INCLUSION_MISMATCH:'+r['target_uniprot'])
 if entity.startswith('immune_') and inc:errors.append('IMMUNE_SEGMENT_IN_ABC:'+r['target_uniprot'])
 if c=='A' and not any(k in mode for k in ['transmembrane','intramembrane','integral','beta_barrel','immunoglobulin_heavy_chain_membrane_isoform']):errors.append('A_MODE_INVALID:'+r['target_uniprot'])
 if c=='B' and not any(k in (mode+' '+r['final_evidence_v7_2']).lower() for k in ['lipid','gpi','myrist','palmit','prenyl']):errors.append('B_MODE_INVALID:'+r['target_uniprot'])
 if c=='C' and not any(k in mode for k in ['peripheral','direct_lipid','surface_binding','state_dependent','complex_mediated']):errors.append('C_MODE_INVALID:'+r['target_uniprot'])
 if c=='EXCLUDED' and inc:warnings.append('EXCLUDED_INCLUDED:'+r['target_uniprot'])
summary={'rows':len(rows),'unique_uniprot':len(set(ids)),'class':dict(Counter(r['final_class_v7_2'] for r in rows)),'entity_type':dict(Counter(r['entity_type_v7_2'] for r in rows)),'blocking_error_count':len(errors),'warning_count':len(warnings),'blocking_errors':errors[:200],'warnings':warnings[:200],'validation_status':'PASS' if not errors else 'FAIL'}
(ROOT/'V7_2_ABC_VALIDATION.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
