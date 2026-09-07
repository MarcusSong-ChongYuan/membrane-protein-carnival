#!/usr/bin/env python3
import csv, gzip, json, time, urllib.parse, urllib.request
from pathlib import Path

SRC=Path(r"D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay\human_membrane_protein_master_v6_4_candidate.tsv.gz")
ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
RAW=ROOT/"02_protein_audit"/"uniprot_current_jsonl"
RAW.mkdir(parents=True,exist_ok=True)
PROGRESS=ROOT/"logs"/"UNIPROT_FULL_AUDIT_PROGRESS.json"
with gzip.open(SRC,"rt",encoding="utf-8",newline="") as f:
    accessions=[r["target_uniprot_id"] for r in csv.DictReader(f,delimiter="\t")]

batch_size=40
batches=[accessions[i:i+batch_size] for i in range(0,len(accessions),batch_size)]
failed=[]; fetched=0; missing=[]
for bi,batch in enumerate(batches):
    out=RAW/f"batch_{bi:04d}.json"
    if out.exists() and out.stat().st_size>10:
        try:
            d=json.loads(out.read_text(encoding="utf-8")); fetched+=len(d.get("results",[])); missing.extend(d.get("missing_accessions",[])); continue
        except Exception: pass
    query=" OR ".join(f"accession:{x}" for x in batch)
    url="https://rest.uniprot.org/uniprotkb/search?"+urllib.parse.urlencode({"query":f"({query})","format":"json","size":500})
    last=""
    for attempt in range(6):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"MemPro-V7-full-audit/1.0"})
            with urllib.request.urlopen(req,timeout=90) as z:data=json.load(z)
            got={x.get("primaryAccession") for x in data.get("results",[])}
            data["requested_accessions"]=batch;data["missing_accessions"]=[x for x in batch if x not in got]
            out.write_text(json.dumps(data,ensure_ascii=False),encoding="utf-8")
            fetched+=len(data.get("results",[]));missing.extend(data["missing_accessions"]);last="";break
        except Exception as e:
            last=type(e).__name__+":"+str(e)[:200];time.sleep(min(60,2**attempt))
    if last:failed.append({"batch":bi,"accessions":batch,"error":last})
    progress={"total_accessions":len(accessions),"total_batches":len(batches),"completed_batches":bi+1,"fetched_entries":fetched,"missing_count":len(set(missing)),"failed_batches":len(failed),"status":"RUNNING"}
    PROGRESS.write_text(json.dumps(progress,indent=2),encoding="utf-8")
    print(f"{bi+1}/{len(batches)} fetched={fetched} missing={len(set(missing))} failed={len(failed)}",flush=True)

final={"total_accessions":len(accessions),"total_batches":len(batches),"fetched_entries":fetched,"unique_missing_accessions":sorted(set(missing)),"failed_batches":failed,"status":"COMPLETE" if not failed else "COMPLETE_WITH_RETRY_QUEUE","retrieval_date":"2026-08-13"}
PROGRESS.write_text(json.dumps(final,ensure_ascii=False,indent=2),encoding="utf-8")
(ROOT/"02_protein_audit"/"UNIPROT_FULL_AUDIT_FETCH_SUMMARY.json").write_text(json.dumps(final,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(final,ensure_ascii=False,indent=2))
