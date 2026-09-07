#!/usr/bin/env python3
import csv,json,re
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
SRC=ROOT/"02_protein_audit"/"final_adjudication_v7_2"/"protein_identity_final_adjudication_all_v7_2.tsv"
CP=Path(r"D:\finale\02_V6.3_candidate_20260803\11_source_snapshots\complex\complexportal_9606_20260114.tsv")
OUT=ROOT/"02_protein_audit"/"t09_validation";OUT.mkdir(exist_ok=True)

membrane_terms=re.compile(r'(plasma membrane|cell membrane|endosome membrane|lysosom(?:e|al) membrane|golgi membrane|endoplasmic reticulum membrane|mitochondri(?:on|al) (?:inner |outer )?membrane|nuclear membrane|membrane raft|phagosome membrane|vesicle membrane)',re.I)
cp_by_protein=defaultdict(list)
with CP.open(encoding='utf-8-sig',newline='') as f:
 for r in csv.DictReader((line for line in f if not line.startswith('#')),fieldnames=None,delimiter='\t'):
  pass
# Header begins with #; parse explicitly.
with CP.open(encoding='utf-8-sig',newline='') as f:
 reader=csv.reader(f,delimiter='\t');hdr=next(reader);hdr[0]=hdr[0].lstrip('#')
 for vals in reader:
  if len(vals)!=len(hdr):continue
  r=dict(zip(hdr,vals)); text=' '.join([r.get('Go Annotations',''),r.get('Description',''),r.get('Complex properties','')])
  if not membrane_terms.search(text):continue
  participants=r.get('Expanded participant list','') or r.get('Identifiers (and stoichiometry) of molecules in complex','')
  for acc in re.findall(r'\b[A-Z0-9]{6,10}\b',participants):
   cp_by_protein[acc].append({'complex_ac':r.get('Complex ac',''),'name':r.get('Recommended name',''),'assembly':r.get('Complex assembly',''),'evidence_code':r.get('Evidence Code',''),'experimental_evidence':r.get('Experimental evidence',''),'membrane_context':membrane_terms.search(text).group(0) if membrane_terms.search(text) else ''})

with SRC.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;rows=list(rd)
extra=['complexportal_membrane_complex_count_t09','complexportal_membrane_complex_ids_t09','complexportal_membrane_complex_names_t09','t09_validation_group','t09_validation_flags','t09_recommended_action','t09_evidence_summary']
out=[]
for r in rows:
 acc=r['target_uniprot']; cps=cp_by_protein.get(acc,[]);flags=[]
 if cps:flags.append('COMPLEXPORTAL_MEMBRANE_COMPLEX_MEMBER')
 if r['current_evidence']=='E1' and r['final_class_v7_2']=='EXCLUDED':group='OLD_E1_EXCLUDED'
 elif r['final_decision_v7_2'] in {'CONFIRMED_C_DIRECT','CONFIRMED_C_STATE_DEPENDENT','CONFIRMED_C_COMPLEX_MEDIATED'}:group='NEW_MECHANISTIC_C'
 elif r['final_decision_v7_2']=='EXCLUDED_NO_MEMBRANE_MECHANISM':group='NO_MECHANISM_EXCLUDED'
 else:group='CONTROL_OR_ALREADY_DIRECT'
 action='KEEP_CURRENT_V7_2'
 evidence=r['final_evidence_v7_2']
 if r['entity_type_v7_2'].startswith('immune_'):
  action='KEEP_NON_PROTEIN_SEGMENT_ENTITY';flags.append('IMMUNE_SEGMENT_ENTITY')
 elif r['final_class_v7_2']=='EXCLUDED' and cps:
  action='RESCUE_C_COMPLEX_MEDIATED_CANDIDATE';evidence=';'.join(f"{x['complex_ac']}|{x['name']}|{x['membrane_context']}|{x['evidence_code']}" for x in cps)
 if group=='NEW_MECHANISTIC_C' and r['final_decision_v7_2']=='CONFIRMED_C_DIRECT' and not any(k in r['final_evidence_v7_2'].lower() for k in ['ida:','imp:','hda:','phosphatidyl','phosphoinosit','lipid','opm']):
  flags.append('C_DIRECT_TEXT_ONLY_REVIEW')
 if group=='NEW_MECHANISTIC_C' and r['final_decision_v7_2']=='CONFIRMED_C_COMPLEX_MEDIATED' and not cps:
  flags.append('C_COMPLEX_NOT_CONFIRMED_BY_COMPLEXPORTAL')
 x=dict(r);x.update({'complexportal_membrane_complex_count_t09':len(cps),'complexportal_membrane_complex_ids_t09':';'.join(sorted({z['complex_ac'] for z in cps})),'complexportal_membrane_complex_names_t09':';'.join(sorted({z['name'] for z in cps})),'t09_validation_group':group,'t09_validation_flags':';'.join(flags),'t09_recommended_action':action,'t09_evidence_summary':evidence});out.append(x)
with (OUT/'T09_PROTEIN_VALIDATION_MATRIX.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=base+extra,delimiter='\t');w.writeheader();w.writerows(out)
summary={'input':len(out),'groups':dict(Counter(x['t09_validation_group'] for x in out)),'actions':dict(Counter(x['t09_recommended_action'] for x in out)),'flags':dict(Counter(y for x in out for y in x['t09_validation_flags'].split(';') if y)),'complexportal_membrane_members':sum(bool(x['complexportal_membrane_complex_count_t09']) for x in out),'status':'T09_OFFLINE_EVIDENCE_MATRIX_COMPLETE'}
(OUT/'T09_VALIDATION_MATRIX_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
