#!/usr/bin/env python3
import csv,glob,json,re
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
SRC=ROOT/"02_protein_audit"/"t09_validation"/"T09_PROTEIN_VALIDATION_MATRIX.tsv"
OUT=ROOT/"02_protein_audit"/"t09_validation"
entries={}
for p in glob.glob(str(ROOT/"02_protein_audit"/"uniprot_current_jsonl"/"batch_*.json")):
 d=json.loads(Path(p).read_text(encoding='utf-8'))
 for e in d.get('results',[]):entries[e['primaryAccession']]=e
exp_codes={'IDA','IMP','IGI','IPI','IEP','EXP','HDA','HMP','HGI','HEP'}
lipid_words=('phosphatidyl','phosphoinosit','phospholipid','lipid binding','membrane binding','sulfatide binding')
def evidence(e):
 go_lipid=[];go_mem=[];pubmed_mem_text=[]
 for x in e.get('uniProtKBCrossReferences',[]):
  if x.get('database')!='GO':continue
  p={z.get('key'):z.get('value','') for z in x.get('properties',[])};term=p.get('GoTerm','');ev=p.get('GoEvidenceType','');code=ev.split(':',1)[0]
  if code in exp_codes and term.startswith('F:') and any(k in term.lower() for k in lipid_words):go_lipid.append(term+'|'+ev)
  if code in exp_codes and term.startswith('C:') and 'membrane' in term.lower():go_mem.append(term+'|'+ev)
 for c in e.get('comments',[]):
  if c.get('commentType') not in {'FUNCTION','SUBCELLULAR LOCATION'}:continue
  texts=list(c.get('texts',[]) or [])+list((c.get('note') or {}).get('texts',[]) or [])
  for t in texts:
   val=t.get('value','');evs=t.get('evidences',[]) or []
   if 'membran' in val.lower() and any(z.get('evidenceCode')=='ECO:0000269' and z.get('source')=='PubMed' for z in evs):pubmed_mem_text.append(val)
 return go_lipid,go_mem,pubmed_mem_text
with SRC.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;rows=list(rd)
extra=['t09_go_lipid_direct','t09_go_membrane_experimental','t09_pubmed_membrane_mechanism_text','t09_final_validation_grade','t09_final_validation_decision','t09_final_validation_reason']
out=[]
for r in rows:
 e=entries[r['target_uniprot']];gl,gm,pt=evidence(e);group=r['t09_validation_group'];grade='NOT_APPLICABLE';decision='KEEP_V7_2';reason='Not in targeted T09 validation groups.'
 if group=='NEW_MECHANISTIC_C':
  d=r['final_decision_v7_2']
  if 'OPM' in r['final_evidence_v7_2'] or gl:
   grade='DIRECT';decision='CONFIRM_C';reason='OPM peripheral structure or experimental GO molecular-function lipid/membrane binding.'
  elif pt:
   grade='EXPERIMENTAL_CONTEXT';decision='CONFIRM_C';reason='UniProt function/location text contains PubMed-backed membrane recruitment/association context.'
  elif d=='CONFIRMED_C_COMPLEX_MEDIATED' and gm and int(r['complexportal_membrane_complex_count_t09'] or 0)>0:
   grade='COMPLEX_SUPPORTED';decision='CONFIRM_C_COMPLEX_MEDIATED';reason='Experimental membrane GO plus ComplexPortal membrane-complex membership.'
  else:
   grade='WEAK';decision='REVERT_EXCLUDED_NO_MECHANISM';reason='No direct structure, experimental lipid-binding, PubMed-backed membrane recruitment, or jointly supported complex-mediated evidence.'
 elif group=='OLD_E1_EXCLUDED':
  if gl or pt:
   grade='RECOVERED_EXPERIMENTAL';decision='RESCUE_C';reason='Direct experimental lipid-binding or PubMed-backed membrane mechanism recovered.'
  elif gm and int(r['complexportal_membrane_complex_count_t09'] or 0)>0 and 'secreted' not in r['uniprot_subcellular_locations_current'].lower():
   grade='RECOVERED_COMPLEX';decision='RESCUE_C_COMPLEX_MEDIATED';reason='Experimental membrane localization plus membrane-complex membership; not a soluble secreted ligand.'
  else:
   grade='EXCLUSION_CONFIRMED';decision='CONFIRM_EXCLUSION';reason='Old E1 reflected localization, not a demonstrated membrane-binding mechanism.'
 x=dict(r);x.update({'t09_go_lipid_direct':';'.join(gl),'t09_go_membrane_experimental':';'.join(gm),'t09_pubmed_membrane_mechanism_text':' || '.join(pt),'t09_final_validation_grade':grade,'t09_final_validation_decision':decision,'t09_final_validation_reason':reason});out.append(x)
target=[x for x in out if x['t09_validation_group'] in {'NEW_MECHANISTIC_C','OLD_E1_EXCLUDED'}]
with (OUT/'T09_MECHANISTIC_C_OLD_E1_VALIDATION.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=base+extra,delimiter='\t');w.writeheader();w.writerows(target)
summary={'target_records':len(target),'groups':dict(Counter(x['t09_validation_group'] for x in target)),'decisions':dict(Counter(x['t09_final_validation_decision'] for x in target)),'grades':dict(Counter(x['t09_final_validation_grade'] for x in target)),'status':'T09_TARGETED_VALIDATION_COMPLETE'}
(OUT/'T09_MECHANISTIC_C_OLD_E1_VALIDATION_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
