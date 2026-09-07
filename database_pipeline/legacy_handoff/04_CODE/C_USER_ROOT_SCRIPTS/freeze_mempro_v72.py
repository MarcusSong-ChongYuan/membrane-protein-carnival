from __future__ import annotations
import hashlib, json, shutil
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

SRC=Path(r"D:\finale\21_MemPro_V7.2_candidate_20260817")
OUT=Path(r"D:\finale\22_MemPro_V7.2_final_20260817")
qa=json.loads((SRC/'03_QA'/'V72_VALIDATION_REPORT.json').read_text(encoding='utf-8'))
if qa.get('status')!='PASS' or qa.get('blocking_errors'):
    raise SystemExit('Refusing freeze: validation is not PASS')
if OUT.exists():
    for p in OUT.rglob('*'):
        if p.is_file(): p.chmod(0o666)
    shutil.rmtree(OUT)
shutil.copytree(SRC,OUT,copy_function=shutil.copyfile)
qa['version']='MemPro_V7.2'
qa['freeze_status']='FINAL_FROZEN'
qa['frozen_at_utc']=datetime.now(timezone.utc).isoformat()
(OUT/'03_QA'/'V72_VALIDATION_REPORT.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
info=json.loads((OUT/'RELEASE_INFO.json').read_text(encoding='utf-8'))
info.update({'release':'MemPro_V7.2','status':'final_frozen','validation':'PASS','frozen_at_utc':qa['frozen_at_utc']})
(OUT/'RELEASE_INFO.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8')
for stale in [OUT/'02_metadata'/'V72_MANIFEST.tsv',OUT/'02_metadata'/'SHA256SUMS.tsv']:
    if stale.exists(): stale.unlink()
rows=[]
for p in sorted(OUT.rglob('*')):
    if p.is_file() and p.name not in {'V72_MANIFEST.tsv','SHA256SUMS.tsv','FROZEN'}:
        rows.append((str(p.relative_to(OUT)),p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()))
mf=pd.DataFrame(rows,columns=['file','bytes','sha256'])
mf.to_csv(OUT/'02_metadata'/'V72_MANIFEST.tsv',sep='\t',index=False)
mf[['sha256','file']].to_csv(OUT/'02_metadata'/'SHA256SUMS.tsv',sep='\t',index=False,header=False)
validation_hash=hashlib.sha256((OUT/'03_QA'/'V72_VALIDATION_REPORT.json').read_bytes()).hexdigest()
(OUT/'FROZEN').write_text(f"MemPro_V7.2\nstatus=FINAL_FROZEN\nfrozen_at_utc={qa['frozen_at_utc']}\nvalidation_sha256={validation_hash}\n",encoding='utf-8')
for p in OUT.rglob('*'):
    if p.is_file(): p.chmod(0o444)
print(json.dumps({'output':str(OUT),'files':len(list(OUT.rglob('*'))),'manifest_rows':len(mf),'validation_sha256':validation_hash,'status':'FINAL_FROZEN'},indent=2))
