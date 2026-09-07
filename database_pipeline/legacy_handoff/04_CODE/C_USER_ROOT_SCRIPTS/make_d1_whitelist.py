import csv,json
from pathlib import Path
root=Path(r'C:\Users\Administrator')
audit=root/'MemPro_V7_2_1_membrane_identity_audit_20260813'/'membrane_identity_audit_auto_v7_2.tsv'
out=root/'MemPro_Docking_D1_AE1_20260813'/'D1_identity_whitelist_v7_2_1.tsv'
with audit.open(encoding='utf-8',newline='') as f: rows=[r for r in csv.DictReader(f,delimiter='\t') if r['auto_docking_identity_eligibility']=='D1_ELIGIBLE_IDENTITY_ONLY']
with out.open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['target_uniprot','approved_symbol','auto_membrane_identity_status','auto_decision_basis'],delimiter='\t');w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in rows)
print(out,len(rows))
