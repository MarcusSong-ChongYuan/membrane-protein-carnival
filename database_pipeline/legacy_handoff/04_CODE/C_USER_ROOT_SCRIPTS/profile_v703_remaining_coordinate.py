import csv,gzip,json
from collections import Counter,defaultdict
from pathlib import Path
p=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites\binding_site_coordinate_remaining_v703.tsv.gz")
status=Counter();source=Counter();fractions=Counter();examples=defaultdict(list);total=0
with gzip.open(p,'rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  total+=1;s=r['uniprot_rescue_status_v703'];status[s]+=1;source[(s,r['source_database'])]+=1
  try:req=int(r.get('uniprot_requested_count_v703') or 0);mat=int(r.get('uniprot_mapped_count_v703') or 0);frac=mat/req if req else 0
  except:req=mat=0;frac=0
  if s=='REMAIN_UNIPROT_PARTIAL_MATCH':
   band='>=0.8' if frac>=.8 else ('>=0.5' if frac>=.5 else '<0.5');fractions[band]+=1
  if len(examples[s])<5:examples[s].append({'site':r['binding_site_instance_id'],'source':r['source_database'],'pdb':r['pdb_ids'],'target':r['target_uniprot_id'],'description':r['residue_or_site_description'],'requested':req,'mapped':mat,'fraction':round(frac,3)})
report={'total':total,'status':dict(status),'status_source':{f'{a}|{b}':n for (a,b),n in source.items()},'partial_fraction_bands':dict(fractions),'examples':dict(examples)};Path(r"D:\finale\12_V7_final_completion_20260813\07_qa\V703_REMAINING_COORDINATE_PROFILE.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='examples'},ensure_ascii=False,indent=2))
