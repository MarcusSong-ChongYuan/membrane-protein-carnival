import csv,hashlib,json
from pathlib import Path
root=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813")
rows=list(csv.DictReader((root/'MANIFEST.tsv').open(encoding='utf-8'),delimiter='\t'));bad=[]
for r in rows:
 p=root/r['relative_path'];h=hashlib.sha256()
 if not p.exists():bad.append({'path':r['relative_path'],'issue':'missing'});continue
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 if str(p.stat().st_size)!=r['bytes'] or h.hexdigest()!=r['sha256']:bad.append({'path':r['relative_path'],'issue':'size_or_sha256_mismatch'})
report={'manifest_rows':len(rows),'verified_files':len(rows)-len(bad),'failures':bad,'status':'PASS_MANIFEST_REVERIFICATION' if not bad else 'FAIL_MANIFEST_REVERIFICATION'}
(root/'02_QA'/'V7_MANIFEST_REVERIFICATION.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
