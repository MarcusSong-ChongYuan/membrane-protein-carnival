#!/usr/bin/env python3
import csv,json,time,urllib.request,urllib.error
from pathlib import Path
ROOT=Path(r"C:\Users\Administrator\MemPro_Docking_D1_AE1_20260813")
SRC=ROOT/'hetatm_entities_in_grid_or_5A.tsv';OUT=ROOT/'hetatm_rcsb_enrichment';OUT.mkdir(exist_ok=True)
with SRC.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;rows=list(rd)
codes=sorted({r['het_resname'].upper() for r in rows});cache={};fails=[]
for i,code in enumerate(codes,1):
 p=OUT/(code+'.json')
 if p.exists():
  try:cache[code]=json.loads(p.read_text(encoding='utf-8'));continue
  except:pass
 url='https://data.rcsb.org/rest/v1/core/chemcomp/'+code
 last=''
 for attempt in range(4):
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'MemPro-HETATM-audit/1.0'})
   with urllib.request.urlopen(req,timeout=30) as z:d=json.load(z)
   p.write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8');cache[code]=d;last='';break
  except urllib.error.HTTPError as e:
   last=f'HTTP {e.code}';break
  except Exception as e:last=type(e).__name__+':'+str(e)[:120];time.sleep(2**attempt)
 if last:fails.append({'code':code,'error':last})
 if i%25==0:print(f'{i}/{len(codes)} resolved={len(cache)} failed={len(fails)}',flush=True)
def info(code):
 d=cache.get(code,{})
 c=d.get('chem_comp',{});desc=d.get('rcsb_chem_comp_descriptor',{})
 return c.get('name',''),c.get('type',''),c.get('formula',''),c.get('formula_weight',''),desc.get('InChIKey','') or desc.get('inchi_key','')
extra=['rcsb_name','rcsb_type','rcsb_formula','rcsb_formula_weight','rcsb_inchikey','rcsb_lookup_status','review_action_v3']
out=[]
for r in rows:
 name,ctype,formula,mw,ik=info(r['het_resname'].upper());cat=r['category_candidate']
 if cat=='buffer_or_crystallization_additive_candidate':action='REMOVE_CANDIDATE_MANUAL_CONFIRM'
 elif cat in {'common_ion_or_metal','cofactor_candidate'}:action='HOLD_SPECIAL_CHEMISTRY_REVIEW'
 elif cat=='glycan_or_sugar':action='HOLD_GLYCAN_CONTEXT_REVIEW'
 else:action='CHECK_COGNATE_LIGAND_VS_COFACTOR'
 x=dict(r);x.update({'rcsb_name':name,'rcsb_type':ctype,'rcsb_formula':formula,'rcsb_formula_weight':mw,'rcsb_inchikey':ik,'rcsb_lookup_status':'PASS' if name else 'NOT_RESOLVED','review_action_v3':action});out.append(x)
with (OUT/'hetatm_entities_rcsb_enriched.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=base+extra,delimiter='\t');w.writeheader();w.writerows(out)
summary={'entities':len(out),'unique_codes':len(codes),'resolved_codes':len(cache),'resolved_entities':sum(x['rcsb_lookup_status']=='PASS' for x in out),'failed_codes':fails,'rule':'No automatic deletion.'}
(OUT/'HETATM_RCSB_ENRICHMENT_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in summary.items() if k!='failed_codes'},indent=2))
