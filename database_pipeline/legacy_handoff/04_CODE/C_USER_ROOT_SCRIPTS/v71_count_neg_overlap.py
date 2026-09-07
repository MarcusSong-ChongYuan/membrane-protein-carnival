import csv,gzip,json
from pathlib import Path
B=Path(r'D:\finale\16_MemPro_V7.0.2_final_20260814')
def rows(p):
 with gzip.open(p,'rt',encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f,delimiter='\t')
pairs={(r['target_uniprot_id'],r['compound_internal_id']) for r in rows(B/'01_data'/'protein_compound_pair_v7.tsv.gz')}
n=m=a=0
for p in [B/'01_data'/'negative_evidence_v7.tsv.gz',B/'02_nonpublic_resolved_status'/'negative_evidence_frozen_v7.tsv.gz']:
 for r in rows(p):
  n+=1
  if (r.get('target_uniprot_id',''),r.get('compound_internal_id','')) in pairs:m+=1
  else:a+=1
print(json.dumps({'all':n,'contextual':m,'archive':a}))
