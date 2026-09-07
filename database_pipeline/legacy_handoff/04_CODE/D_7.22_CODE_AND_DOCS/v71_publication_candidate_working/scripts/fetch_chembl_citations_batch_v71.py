import csv,gzip,json,time,urllib.request,urllib.parse,urllib.error
from pathlib import Path
import pandas as pd
ROOT=Path(r'D:\7.22\v71_publication_candidate_working');W=ROOT/'work';Q=ROOT/'qa'
docs=sorted(set(pd.read_csv(W/'chembl_activity_document_crosswalk_v71.tsv.gz',sep='\t',usecols=['chembl_document_id'],dtype=str).chembl_document_id.dropna())-{'nan',''})
CHK=W/'chembl_document_citations_v71.checkpoint.jsonl';done={}
if CHK.exists():
 for line in CHK.read_text(encoding='utf8').splitlines():
  try:r=json.loads(line);done[r['document_chembl_id']]=r
  except:pass
pending=[d for d in docs if d not in done]
def fetch_batch(batch):
 q=urllib.parse.quote(','.join(batch),safe=',');url=f'https://www.ebi.ac.uk/chembl/api/data/document.json?document_chembl_id__in={q}&limit=500';got={}
 for k in range(5):
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'MemPro-V7.1-publication-citation-recovery/1.0'});j=json.load(urllib.request.urlopen(req,timeout=60));got={x['document_chembl_id']:x for x in j.get('documents',[])};break
  except Exception:time.sleep(min(2**k,15))
 return [{'document_chembl_id':d,'pubmed_id':got.get(d,{}).get('pubmed_id') or '','doi':got.get(d,{}).get('doi') or '','doc_type':got.get(d,{}).get('doc_type') or '','patent_id':got.get(d,{}).get('patent_id') or '','title':got.get(d,{}).get('title') or '','year':got.get(d,{}).get('year') or '','journal':got.get(d,{}).get('journal') or '','status':'ok' if d in got else 'batch_not_found','source_url':f'https://www.ebi.ac.uk/chembl/api/data/document/{d}.json'} for d in batch]
from concurrent.futures import ThreadPoolExecutor
batches=[pending[i:i+180] for i in range(0,len(pending),180)]
with ThreadPoolExecutor(max_workers=8) as ex:
 for i,rs in enumerate(ex.map(fetch_batch,batches),1):
  with CHK.open('a',encoding='utf8') as f:
   for r in rs:f.write(json.dumps(r,ensure_ascii=False)+'\n');done[r['document_chembl_id']]=r
  if i%8==0:print(min(i*180,len(pending)),len(pending),flush=True)
fields=['document_chembl_id','pubmed_id','doi','doc_type','patent_id','title','year','journal','status','source_url']
with (W/'chembl_document_citations_v71.tsv').open('w',encoding='utf8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for d in docs:w.writerow({k:done.get(d,{}).get(k,'') for k in fields})
rep={'documents':len(docs),'ok':sum(r.get('status')=='ok' for r in done.values()),'with_pubmed':sum(bool(r.get('pubmed_id')) for r in done.values()),'with_doi':sum(bool(r.get('doi')) for r in done.values()),'publication':sum(r.get('doc_type')=='PUBLICATION' for r in done.values()),'patent':sum(r.get('doc_type')=='PATENT' for r in done.values()),'failed_or_not_found':sum(r.get('status')!='ok' for r in done.values())};(Q/'V71_CHEMBL_CITATION_RECOVERY.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2))
