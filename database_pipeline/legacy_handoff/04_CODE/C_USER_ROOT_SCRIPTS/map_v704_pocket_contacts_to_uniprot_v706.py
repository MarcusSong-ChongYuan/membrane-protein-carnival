import csv,gzip,json,re,time
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites");IN=ROOT/'binding_site_coordinate_rescued_pocket_v704b.tsv.gz';MAP=ROOT/'sifts_residue_mapping_v7.checkpoint.tsv';OUT=ROOT/'binding_site_coordinate_rescued_pocket_mapped_v706.tsv.gz'
sites=[];interest=defaultdict(list)
with gzip.open(IN,'rt',encoding='utf-8',newline='') as f:
 rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames
 for r in rd:
  contacts=[]
  for token in (r.get('pocket_contact_residues_v704') or '').split(';'):
   m=re.fullmatch(r'([^:]+):([A-Z]{3})(-?\d+)([A-Za-z]?)',token.strip())
   if m:contacts.append({'chain':m.group(1),'aa':m.group(2),'num':m.group(3),'icode':m.group(4).upper()})
  idx=len(sites);sites.append({'row':r,'contacts':contacts,'mapped':{}});p=r.get('pocket_pdb_v704','').lower();a=r['target_uniprot_id']
  if p:interest[(p,a)].append(idx)
n=0;start=time.time()
with MAP.open(encoding='utf-8',newline='') as f:
 for m in csv.DictReader(f,delimiter='\t'):
  n+=1;idxs=interest.get((m['pdb_id'].lower(),m['uniprot_acc']))
  if not idxs:continue
  for idx in idxs:
   for i,t in enumerate(sites[idx]['contacts']):
    if t['chain']==m['chain'] and t['num']==m.get('pdb_resnum','') and t['icode']==(m.get('pdb_icode') or '').upper() and t['aa']==m.get('pdb_resname','').upper():sites[idx]['mapped'][i]=m
extra=['pocket_uniprot_mapping_status_v706','pocket_contact_uniprot_residues_v706','pocket_contact_mapping_fraction_v706','pocket_coordinate_confidence_v706'];fields=base+extra;counts=Counter()
with gzip.open(OUT,'wt',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for s in sites:
  total=len(s['contacts']);mapped=len(s['mapped']);frac=mapped/total if total else 0
  if frac==1:status='ALL_CONTACTS_MAPPED_TO_UNIPROT';conf='HIGH'
  elif frac>=.8:status='MOST_CONTACTS_MAPPED_TO_UNIPROT';conf='MEDIUM'
  else:status='AUTHOR_COORDINATE_POCKET_UNIPROT_PARTIAL';conf='MEDIUM_AUTHOR_COORDINATE'
  uni=';'.join(f"{z.get('uniprot_acc','')}:{z.get('uniprot_resnum','')}" for _,z in sorted(s['mapped'].items()));r=s['row'];r.update({'pocket_uniprot_mapping_status_v706':status,'pocket_contact_uniprot_residues_v706':uni,'pocket_contact_mapping_fraction_v706':f'{frac:.4f}','pocket_coordinate_confidence_v706':conf});counts[status]+=1;w.writerow(r)
report={'input':len(sites),'sifts_rows_scanned':n,'status':dict(counts),'runtime_seconds':round(time.time()-start,1),'policy':'exact PDB author chain/resnum/insertion-code/resname to SIFTS UniProt mapping'};(ROOT/'T26_T28_POCKET_UNIPROT_V706_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
