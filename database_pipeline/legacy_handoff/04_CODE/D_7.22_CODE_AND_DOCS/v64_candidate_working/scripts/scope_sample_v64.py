from pathlib import Path
import csv,gzip,json,collections,random,hashlib
O=Path(r'D:\7.22\v64_candidate_working');E=O/'candidate_tables'/'binding_evidence_master_v6_4_candidate.tsv.gz';S=O/'stats';Q=O/'qa';R=O/'review_queues'
def truth(v):return str(v).lower() in {'1','true','yes','y'}
allpairs=set();defaultpairs=set();tier=collections.Counter();dtier=collections.Counter();sources=collections.Counter();dsources=collections.Counter();rows=0;drows=0;samples=collections.defaultdict(list)
with gzip.open(E,'rt',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  rows+=1;k=(r['target_uniprot_id'],r['compound_internal_id']);allpairs.add(k);tier[r['evidence_tier']]+=1;sources[r['source_database']]+=1
  inc=truth(r.get('default_release_inclusion')) and truth(r.get('default_release_inclusion_v4_2')) and r.get('record_qc_status','').lower() not in {'review','fail','excluded'}
  if inc:drows+=1;defaultpairs.add(k);dtier[r['evidence_tier']]+=1;dsources[r['source_database']]+=1
  key=(r.get('source_database',''),r.get('evidence_tier',''),r.get('pdb_identity_status_v64',''))
  h=int(hashlib.sha256(r['evidence_id'].encode()).hexdigest()[:16],16)
  a=samples[key]
  if len(a)<5:a.append((h,r.copy()))
  else:
   worst=max(range(len(a)),key=lambda i:a[i][0])
   if h<a[worst][0]:a[worst]=(h,r.copy())
summary={'all_evidence_rows':rows,'all_unique_pairs':len(allpairs),'default_evidence_rows':drows,'default_unique_pairs':len(defaultpairs),'all_tiers':dict(tier),'default_tiers':dict(dtier),'all_sources':dict(sources),'default_sources':dict(dsources)}
(S/'V64_RELEASE_SCOPE_COUNTS.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
with open(R/'V64_STRATIFIED_VALIDATION_SAMPLE.tsv','w',encoding='utf-8',newline='') as f:
 fields=['sample_stratum']+list(next(iter(next(iter(samples.values()))))[1].keys()) if samples else ['sample_stratum'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for key,a in sorted(samples.items()):
  for h,r in sorted(a):w.writerow({'sample_stratum':'|'.join(key),**r})
print(json.dumps(summary,ensure_ascii=False,indent=2));print('sample_rows',sum(len(x) for x in samples.values()))

