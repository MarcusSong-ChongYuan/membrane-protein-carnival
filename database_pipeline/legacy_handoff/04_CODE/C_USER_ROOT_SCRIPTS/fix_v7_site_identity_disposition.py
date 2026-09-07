import csv,gzip,json
from collections import Counter
from pathlib import Path
root=Path(r"D:\finale\09_V7_data_freeze_working_20260813");src=root/'03_v7_relation_rebuild'/'binding_site_instances_V7.tsv.gz';out=root/'05_v7_automatic_release_candidate'/'binding_site_identity_qc_V7.tsv.gz'
c=Counter()
with gzip.open(src,'rt',encoding='utf-8',newline='') as fi:
 r=csv.DictReader(fi,delimiter='\t');fields=list(r.fieldnames)+['site_identity_disposition_v7','site_identity_reason_v7']
 with gzip.open(out,'wt',encoding='utf-8',newline='') as fo:
  w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
  for x in r:
   st=x.get('pdb_identity_status_v64','')
   if st=='valid_or_empty':disp='PUBLIC'
   elif st=='cleaned_partial':disp='PUBLIC_WITH_ORIGINAL_PDB_AUDIT_TRAIL'
   else:disp='FROZEN_REVIEW'
   x['site_identity_disposition_v7']=disp;x['site_identity_reason_v7']=st or 'missing_v64_identity_status';c[disp]+=1;w.writerow(x)
report={'binding_site_identity_disposition':dict(c),'policy':{'valid_or_empty':'PUBLIC','cleaned_partial':'PUBLIC_WITH_ORIGINAL_PDB_AUDIT_TRAIL; cleaned current IDs retained and original IDs preserved','other':'FROZEN_REVIEW'}}
(out.parent/'V7_SITE_IDENTITY_DISPOSITION_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
