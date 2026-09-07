import csv,gzip,collections,json
from pathlib import Path
m=Path(r'C:\Users\Administrator\Desktop\MemPro_Docking_Dscope_14333_ProductionV3_20260813\docking_manifest.tsv')
p=Path(r'D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay\human_membrane_protein_master_v6_4_candidate.tsv.gz')
with gzip.open(p,'rt',encoding='utf-8',newline='') as f: pm={r['target_uniprot_id']:r for r in csv.DictReader(f,delimiter='\t')}
with m.open(encoding='utf-8',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
cross=collections.Counter((pm[r['target_uniprot']]['membrane_class_v52'],pm[r['target_uniprot']]['evidence_level_v52']) for r in rows)
ucross=collections.Counter((pm[u]['membrane_class_v52'],pm[u]['evidence_level_v52']) for u in {r['target_uniprot'] for r in rows})
# Conservative first-pass flags available locally before PDB/SIFTS QC
flags=collections.Counter()
for r in rows:
 q=pm[r['target_uniprot']]
 if q['membrane_class_v52']=='A' and q['evidence_level_v52']=='E1': flags['A_E1']=flags['A_E1']+1
 if q['membrane_class_v52']=='A' and q['evidence_level_v52']=='E1' and q.get('membrane_decision_v5') in {'retained_v4_integral','external_integral_evidence_promoted','new_external_integral_evidence'}: flags['A_E1_integral_rule']=flags['A_E1_integral_rule']+1
 if q['membrane_class_v52']=='B' and q['evidence_level_v52']=='E1': flags['B_E1']=flags['B_E1']+1
 if q['membrane_class_v52']=='C' and q['evidence_level_v52']=='E1': flags['C_E1']=flags['C_E1']+1
print(json.dumps({'task_cross':{f'{a}/{b}':n for (a,b),n in sorted(cross.items())},'protein_cross':{f'{a}/{b}':n for (a,b),n in sorted(ucross.items())},'provisional':flags},indent=2))
