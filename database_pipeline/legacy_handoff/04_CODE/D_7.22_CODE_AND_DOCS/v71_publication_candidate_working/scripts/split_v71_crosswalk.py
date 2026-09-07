import csv,gzip,json
from pathlib import Path
R=Path(r'D:\7.22\v71_publication_candidate_working');O=R/'outputs';Q=R/'qa'
ids=set()
with gzip.open(O/'binding_evidence_public_candidate_v7_1.tsv.gz','rt',encoding='utf8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):ids.add(r['evidence_id'])
n=0
with gzip.open(O/'evidence_source_crosswalk_v7_1.tsv.gz','rt',encoding='utf8',newline='') as fi,gzip.open(O/'public_evidence_source_crosswalk_v7_1.tsv.gz','wt',encoding='utf8',newline='',compresslevel=6) as fo:
 rd=csv.DictReader(fi,delimiter='\t');wr=csv.DictWriter(fo,fieldnames=rd.fieldnames,delimiter='\t',lineterminator='\n');wr.writeheader()
 for r in rd:
  if r['evidence_id'] in ids:wr.writerow(r);n+=1
out={'public_evidence_ids':len(ids),'public_crosswalk_rows':n,'status':'PASS' if n==len(ids) else 'FAIL'};(Q/'V71_CROSSWALK_SPLIT_QA.json').write_text(json.dumps(out,indent=2),encoding='utf8');print(json.dumps(out,indent=2))
