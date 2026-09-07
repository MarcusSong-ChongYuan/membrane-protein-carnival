import csv,gzip,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path
root=Path(r"D:\finale\12_V7_final_completion_20260813\05_hpa_expression");src=root/'expression_measurement_v7_complete.tsv.gz';pub=root/'expression_measurement_v7_deduplicated.tsv.gz';rev=root/'expression_id_conflicts_frozen_v7.tsv.gz'
groups={}; conflict_ids=set(); exact=0
with gzip.open(src,'rt',encoding='utf-8',newline='') as f:
 rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames
 for r in rd:
  eid=r['expression_record_id'];sig=tuple(r[k] for k in fields if k!='expression_record_id')
  if eid not in groups:groups[eid]=(sig,r)
  elif groups[eid][0]==sig:exact+=1
  else:conflict_ids.add(eid)
with gzip.open(src,'rt',encoding='utf-8',newline='') as fi,gzip.open(pub,'wt',encoding='utf-8',newline='') as fp,gzip.open(rev,'wt',encoding='utf-8',newline='') as fr:
 rd=csv.DictReader(fi,delimiter='\t');w=csv.DictWriter(fp,fieldnames=rd.fieldnames,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=rd.fieldnames,delimiter='\t');w.writeheader();wr.writeheader();written=set();np=nr=0
 for r in rd:
  eid=r['expression_record_id']
  if eid in conflict_ids:wr.writerow(r);nr+=1
  elif eid not in written:w.writerow(r);written.add(eid);np+=1
report={'input_rows':np+nr+exact,'public_unique_rows':np,'exact_duplicate_extra_rows_removed':exact,'conflicting_ids_frozen':len(conflict_ids),'conflicting_rows_frozen':nr,'public_duplicate_primary_keys':0}
(root/'T36_T38_HPA_DEDUP_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
