import csv,collections,json
from pathlib import Path
root=Path(r'C:\Users\Administrator\MemPro_V7_2_1_membrane_identity_audit_20260813')
src=root/'membrane_identity_audit_auto_v7_2.tsv'
with src.open(encoding='utf-8',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
def rank(r):
 f=set((r['auto_conflict_flags'] or '').split(';'))
 if 'A_CLASS_VS_NONMEMBRANE_UNIPROT_LOCATION' in f: return 'P0_A_vs_location'
 if 'A_CLASS_WITH_ZERO_TM_FEATURES' in f: return 'P0_A_zero_TM'
 if 'HPA_ONLY_MEMBRANE_ASSIGNMENT' in f or 'HPA_VS_UNIPROT_LOCATION_CONFLICT' in f: return 'P1_HPA_conflict'
 if 'EXISTING_CROSS_SOURCE_CONFLICT' in f: return 'P1_existing_conflict'
 if 'E2_REVIEW_PRIORITY' in f: return 'P2_E2'
 if 'B_OR_C_SCOPE_REVIEW' in f: return 'P3_BC_E1'
 return 'P4_control'
for r in rows:r['manual_review_priority_v721']=rank(r)
fields=list(rows[0])+['manual_review_priority_v721']
with (root/'membrane_identity_manual_review_queue_v7_2_1.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(sorted(rows,key=lambda r:(r['manual_review_priority_v721'],r['approved_symbol'])))
print(json.dumps(dict(collections.Counter(r['manual_review_priority_v721'] for r in rows)),indent=2))
