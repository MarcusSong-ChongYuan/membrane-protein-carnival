#!/usr/bin/env python3
import csv,json,time,urllib.request
from pathlib import Path
ROOT=Path(r"C:\Users\Administrator\MemPro_Docking_D1_AE1_20260813")
SRC=ROOT/'hetatm_entities_in_grid_or_5A.tsv';OUT=ROOT/'hetatm_pdbe_enrichment';OUT.mkdir(exist_ok=True)
with SRC.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;rows=list(rd)
codes=sorted({r['het_resname'].upper() for r in rows}); cache={};fail=[]
for i,code in enumerate(codes,1):
 p=OUT/(code+'.json')
 if p.exists():
  try:cache[code]=json.loads(p.read_text(encoding='utf-8'));continue
  except:pass
 url='https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/'+code
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'MemPro-HETATM-audit/1.0'})
  with urllib.request.urlopen(req,timeout=30) as z:d=json.load(z)
  p.write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8');cache[code]=d
 except Exception as e:fail.append({'code':code,'error':type(e).__name__+':'+str(e)[:120]})
 if i%25==0:print(f'{i}/{len(codes)} failed={len(fail)}',flush=True)
def summary(code):
 d=cache.get(code,{})
 vals=d.get(code.lower()) or d.get(code.upper()) or []
 x=vals[0] if vals else {}
 name=x.get('name') or x.get('compound_name') or ''
 ctype=x.get('compound_type') or x.get('type') or ''
 formula=x.get('formula') or ''
 return name,ctype,formula
extra=['pdbe_name','pdbe_compound_type','pdbe_formula','pdbe_lookup_status','review_action_v2']
out=[]
for r in rows:
 name,ctype,formula=summary(r['het_resname'].upper())
 cat=r['category_candidate']
 if cat=='buffer_or_crystallization_additive_candidate':action='REMOVE_CANDIDATE_MANUAL_CONFIRM'
 elif cat in {'common_ion_or_metal','cofactor_candidate'}:action='HOLD_SPECIAL_CHEMISTRY_REVIEW'
 elif cat=='glycan_or_sugar':action='HOLD_GLYCAN_CONTEXT_REVIEW'
 else:action='CHECK_COGNATE_LIGAND_VS_COFACTOR'
 x=dict(r);x.update({'pdbe_name':name,'pdbe_compound_type':ctype,'pdbe_formula':formula,'pdbe_lookup_status':'PASS' if name else 'NOT_RESOLVED','review_action_v2':action});out.append(x)
with (OUT/'hetatm_entities_pdbe_enriched.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=base+extra,delimiter='\t');w.writeheader();w.writerows(out)
(OUT/'HETATM_PDBE_ENRICHMENT_SUMMARY.json').write_text(json.dumps({'entities':len(out),'unique_codes':len(codes),'resolved_entities':sum(x['pdbe_lookup_status']=='PASS' for x in out),'failed_codes':fail,'rule':'No automatic deletion.'},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'entities':len(out),'unique_codes':len(codes),'resolved_entities':sum(x['pdbe_lookup_status']=='PASS' for x in out),'failed_code_count':len(fail)},indent=2))
