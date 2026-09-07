import asyncio,csv,gzip,json,time,urllib.request,urllib.error
from pathlib import Path
import pandas as pd
ROOT=Path(r'D:\7.22\v71_publication_candidate_working');W=ROOT/'work';Q=ROOT/'qa'
IN=W/'chembl_activity_document_crosswalk_v71.tsv.gz';OUT=W/'chembl_document_citations_v71.tsv';CHK=W/'chembl_document_citations_v71.checkpoint.jsonl'
docs=sorted(set(pd.read_csv(IN,sep='\t',usecols=['chembl_document_id'],dtype=str).chembl_document_id.dropna())-{'nan',''})
done={}
if CHK.exists():
 for line in CHK.read_text(encoding='utf8').splitlines():
  try:r=json.loads(line);done[r['document_chembl_id']]=r
  except:pass
pending=[d for d in docs if d not in done]
def one(d):
 url=f'https://www.ebi.ac.uk/chembl/api/data/document/{d}.json'
 for k in range(5):
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'MemPro-V7.1-publication-citation-recovery/1.0'})
   with urllib.request.urlopen(req,timeout=35) as resp:j=json.load(resp)
   return {'document_chembl_id':d,'pubmed_id':j.get('pubmed_id') or '','doi':j.get('doi') or '','doc_type':j.get('doc_type') or '','patent_id':j.get('patent_id') or '','title':j.get('title') or '','year':j.get('year') or '','journal':j.get('journal') or '','status':'ok','source_url':url}
  except urllib.error.HTTPError as e:
   if e.code==404:return {'document_chembl_id':d,'status':'not_found','source_url':url}
   last='HTTP'+str(e.code)
  except Exception as e:last=type(e).__name__
  time.sleep(min(2**k,15))
 return {'document_chembl_id':d,'status':'failed_'+last,'source_url':url}
from concurrent.futures import ThreadPoolExecutor
with ThreadPoolExecutor(max_workers=12) as ex:
 for i in range(0,len(pending),200):
  rs=list(ex.map(one,pending[i:i+200]))
  with CHK.open('a',encoding='utf8') as f:
   for r in rs:f.write(json.dumps(r,ensure_ascii=False)+'\n');done[r['document_chembl_id']]=r
  print(i+len(rs),len(pending),flush=True)
fields=['document_chembl_id','pubmed_id','doi','doc_type','patent_id','title','year','journal','status','source_url']
with OUT.open('w',encoding='utf8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for d in docs:w.writerow({k:done.get(d,{}).get(k,'') for k in fields})
rep={'documents':len(docs),'ok':sum(r.get('status')=='ok' for r in done.values()),'with_pubmed':sum(bool(r.get('pubmed_id')) for r in done.values()),'with_doi':sum(bool(r.get('doi')) for r in done.values()),'failed':sum(not str(r.get('status','')).startswith('ok') for r in done.values())};(Q/'V71_CHEMBL_CITATION_RECOVERY.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2))
