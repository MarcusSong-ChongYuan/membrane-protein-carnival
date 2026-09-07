import csv,gzip,json,collections
from pathlib import Path
manifest=Path(r'C:\Users\Administrator\Desktop\MemPro_Docking_Dscope_14333_ProductionV3_20260813\docking_manifest.tsv')
protein=Path(r'D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay\human_membrane_protein_master_v6_4_candidate.tsv.gz')
want=['P83881','Q9H0A0','Q9BT17','Q92736','P68104','O60551','P62241']
with gzip.open(protein,'rt',encoding='utf-8',newline='') as f:
    rr=csv.DictReader(f,delimiter='\t')
    pm={r['target_uniprot_id']:r for r in rr}
with manifest.open(encoding='utf-8',newline='') as f: tasks=list(csv.DictReader(f,delimiter='\t'))
ups={r['target_uniprot'] for r in tasks}
missing=sorted(ups-pm.keys())
fields=['membrane_class_v52','evidence_level_v52','release_scope_v52','website_default_v52','membrane_scope_v5','final_membrane_class_v51','membrane_decision_v5']
def ct(key, rows): return dict(collections.Counter((pm.get(r['target_uniprot'],{}).get(key) or 'MISSING') for r in rows))
print(json.dumps({
 'tasks':len(tasks),'unique_proteins':len(ups),'missing_unique':len(missing),'missing_examples':missing[:20],
 'task_counts':{k:ct(k,tasks) for k in fields},
 'unique_counts':{k:dict(collections.Counter((pm.get(u,{}).get(k) or 'MISSING') for u in ups)) for k in fields},
 'strict_tasks_ABC_E1E2':sum(1 for r in tasks if pm.get(r['target_uniprot'],{}).get('membrane_class_v52') in {'A','B','C'} and pm.get(r['target_uniprot'],{}).get('evidence_level_v52') in {'E1','E2'} and pm.get(r['target_uniprot'],{}).get('release_scope_v52')=='confirmed'),
 'strict_unique_ABC_E1E2':len({r['target_uniprot'] for r in tasks if pm.get(r['target_uniprot'],{}).get('membrane_class_v52') in {'A','B','C'} and pm.get(r['target_uniprot'],{}).get('evidence_level_v52') in {'E1','E2'} and pm.get(r['target_uniprot'],{}).get('release_scope_v52')=='confirmed'}),
},ensure_ascii=False,indent=2))
print('EXAMPLES')
for u in want:
 r=pm.get(u,{})
 print('\t'.join(str(x) for x in [u,r.get('approved_symbol'),r.get('protein_name'),r.get('membrane_class_v52'),r.get('evidence_level_v52'),r.get('release_scope_v52'),r.get('website_default_v52'),r.get('membrane_scope_v5'),r.get('membrane_decision_v5'),r.get('membrane_evidence_basis'),r.get('subcellular_location')]))
