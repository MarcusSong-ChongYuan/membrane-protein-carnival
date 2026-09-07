#!/usr/bin/env python3
import csv, gzip, hashlib, json
from pathlib import Path
from collections import Counter

ROOT = Path(r'C:\Users\Administrator')
MANIFEST = ROOT/'Desktop'/'MemPro_Docking_Dscope_14333_ProductionV3_20260813'/'docking_manifest.tsv'
PROTEINS = Path(r'D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay\human_membrane_protein_master_v6_4_candidate.tsv.gz')
OUT = ROOT/'MemPro_Docking_D1_AE1_20260813'
OUT.mkdir(exist_ok=True)

with gzip.open(PROTEINS,'rt',encoding='utf-8',newline='') as f:
    proteins={r['target_uniprot_id']:r for r in csv.DictReader(f,delimiter='\t')}
with MANIFEST.open(encoding='utf-8',newline='') as f:
    reader=csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; tasks=list(reader)

d1=[]
for r in tasks:
    p=proteins[r['target_uniprot']]
    if p.get('membrane_class_v52')=='A' and p.get('evidence_level_v52')=='E1' and p.get('website_default_v52')=='1':
        x=dict(r)
        x.update({'membrane_class_v52':p.get('membrane_class_v52',''),
                  'membrane_evidence_level_v52':p.get('evidence_level_v52',''),
                  'membrane_decision_v5':p.get('membrane_decision_v5',''),
                  'membrane_evidence_basis':p.get('membrane_evidence_basis',''),
                  'd1_stage':'D1_AE1_IDENTITY_PASS_STRUCTURE_PENDING'})
        d1.append(x)

d1_fields=fields+['membrane_class_v52','membrane_evidence_level_v52','membrane_decision_v5','membrane_evidence_basis','d1_stage']
with (OUT/'docking_manifest_D1_AE1_candidate.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=d1_fields,delimiter='\t');w.writeheader();w.writerows(d1)

with (OUT/'pdb_download_list_D1.txt').open('w',encoding='ascii',newline='\n') as f:
    f.write('\n'.join(sorted({r['pdb_id'].lower() for r in d1 if r.get('pdb_id')}))+'\n')
with (OUT/'ligands_D1.smi').open('w',encoding='utf-8',newline='\n') as f:
    seen=set()
    for r in d1:
        k=r['compound_internal_id']
        if k not in seen:
            seen.add(k); f.write(f"{r['smiles']}\t{k}\t{r['preferred_name']}\n")

targets=sorted({r['target_uniprot'] for r in tasks})
audit_fields=['target_uniprot','approved_symbol','protein_name','membrane_class_v52','evidence_level_v52','release_scope_v52','website_default_v52','membrane_scope_v5','membrane_decision_v5','membrane_evidence_basis','subcellular_location','transmembrane_count_v5','opm_present_v5','pdbtm_present_v5','hpa_plasma_membrane_location_v5','cross_source_conflict_v5','audit_priority','audit_status']
audit=[]
for u in targets:
    p=proteins[u]; cl=p.get('membrane_class_v52',''); ev=p.get('evidence_level_v52','')
    priority='P1_E2' if ev=='E2' else ('P2_BC_E1' if cl in {'B','C'} else 'P3_A_E1_CONTROL')
    audit.append({k:(u if k=='target_uniprot' else p.get(k,'')) for k in audit_fields[:-2]} | {'audit_priority':priority,'audit_status':'PENDING'})
with (OUT/'membrane_identity_audit_seed_v7_2.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=audit_fields,delimiter='\t');w.writeheader();w.writerows(audit)

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
report={'source_tasks':len(tasks),'source_unique_proteins':len(targets),'d1_candidate_tasks':len(d1),'d1_unique_proteins':len({r['target_uniprot'] for r in d1}),'d1_unique_pdb':len({r['pdb_id'].lower() for r in d1}),'d1_unique_compounds':len({r['compound_internal_id'] for r in d1}),'audit_priorities':dict(Counter(r['audit_priority'] for r in audit)),'source_manifest_sha256':sha(MANIFEST)}
(OUT/'D1_BUILD_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
with (OUT/'SHA256SUMS.tsv').open('w',encoding='utf-8',newline='') as f:
    f.write('file\tsha256\n')
    for p in sorted(OUT.iterdir()):
        if p.is_file() and p.name!='SHA256SUMS.tsv':f.write(f'{p.name}\t{sha(p)}\n')
print(json.dumps(report,indent=2))
