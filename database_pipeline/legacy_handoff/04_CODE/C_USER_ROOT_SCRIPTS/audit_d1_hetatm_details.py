#!/usr/bin/env python3
import csv,math,json
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parent
MAN=ROOT/'docking_manifest_D1_G1_chain_receptor_qc_v2.tsv'
PDB=Path('/home/csong/docking/handoff/MemPro_Docking_Dscope_RCSBclean_14333_20260811/pdb_structures')
OUT=ROOT/'hetatm_audit_chain_v2';OUT.mkdir(exist_ok=True)
waters={'HOH','WAT','DOD'}
ions={'NA','K','CL','CA','MG','MN','ZN','FE','CU','CO','NI','CD','HG','SR','CS','IOD','BR'}
cofactors={'HEM','HEC','FAD','FMN','NAD','NAP','NDP','ADP','ATP','GDP','GTP','PLP','SAM','SAH','COA','ACO','RET'}
sugars={'NAG','BMA','MAN','FUC','GAL','GLC','SIA','NDG'}
buffers={'SO4','PO4','GOL','EDO','PEG','PGE','MPD','ACT','FMT','TRS','MES','HEP','DMS','ACE'}
def category(res):
 if res in waters:return 'water'
 if res in ions:return 'common_ion_or_metal'
 if res in cofactors:return 'cofactor_candidate'
 if res in sugars:return 'glycan_or_sugar'
 if res in buffers:return 'buffer_or_crystallization_additive_candidate'
 return 'ligand_or_unknown_heteroentity'
def inside_shell(x,y,z,r,margin=5):
 vals=[]
 for a,v in zip('xyz',(x,y,z)):
  c=float(r['grid_center_'+a]);half=float(r['grid_size_'+a])/2+margin
  vals.append(abs(v-c)<=half)
 return all(vals)
def center_dist(x,y,z,r):
 return math.sqrt(sum((v-float(r['grid_center_'+a]))**2 for a,v in zip('xyz',(x,y,z))))
with MAN.open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
detail=[]
for r in rows:
 if r['chain_prep_status']!='HOLD_NONWATER_HETATM_NEAR_POCKET':continue
 path=PDB/(r['pdb_id'].lower()+'.pdb');chain=r['actual_chain_used'].strip()
 if not path.exists():continue
 with path.open(encoding='utf-8',errors='ignore') as f:
  for line in f:
   if not line.startswith('HETATM') or line[21:22].strip()!=chain:continue
   res=line[17:20].strip()
   if res in waters:continue
   try:x,y,z=map(float,(line[30:38],line[38:46],line[46:54]))
   except:continue
   if not inside_shell(x,y,z,r):continue
   detail.append({'task_id':r['task_id'],'pdb_id':r['pdb_id'],'author_chain':chain,'target_uniprot':r['target_uniprot'],'compound_internal_id':r['compound_internal_id'],'compound_name':r['preferred_name'],'het_resname':res,'het_resseq':line[22:26].strip(),'insertion_code':line[26:27].strip(),'atom_name':line[12:16].strip(),'element':line[76:78].strip(),'x':x,'y':y,'z':z,'distance_to_box_center_A':round(center_dist(x,y,z,r),3),'category_candidate':category(res),'automatic_action':'REVIEW_DO_NOT_DELETE'})
fields=list(detail[0]) if detail else []
with (OUT/'hetatm_atoms_in_grid_or_5A.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(detail)
entities={}
for x in detail:
 k=(x['task_id'],x['pdb_id'],x['author_chain'],x['het_resname'],x['het_resseq'],x['insertion_code'])
 e=entities.setdefault(k,{q:x[q] for q in ['task_id','pdb_id','author_chain','target_uniprot','compound_internal_id','compound_name','het_resname','het_resseq','insertion_code','category_candidate','automatic_action']})
 e.setdefault('atom_count',0);e['atom_count']+=1
 e.setdefault('min_distance_to_box_center_A',x['distance_to_box_center_A']);e['min_distance_to_box_center_A']=min(e['min_distance_to_box_center_A'],x['distance_to_box_center_A'])
elist=list(entities.values())
ef=list(elist[0]) if elist else []
with (OUT/'hetatm_entities_in_grid_or_5A.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=ef,delimiter='\t');w.writeheader();w.writerows(elist)
summary={'held_tasks':len({x['task_id'] for x in detail}),'hetero_atoms':len(detail),'hetero_entities':len(elist),'resname_counts':dict(Counter(x['het_resname'] for x in elist)),'category_counts':dict(Counter(x['category_candidate'] for x in elist)),'rule':'Candidate classification only. No HETATM is automatically deleted.'}
(OUT/'HETATM_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
