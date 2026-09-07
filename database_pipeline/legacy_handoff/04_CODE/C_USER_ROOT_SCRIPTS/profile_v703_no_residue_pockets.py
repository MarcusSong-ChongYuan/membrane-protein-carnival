import csv,gzip,json,re
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites");REM=ROOT/'binding_site_coordinate_remaining_v703.tsv.gz';EVID=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813\01_data\protein_compound_evidence_v7.tsv.gz")
need=[];eids=set()
with gzip.open(REM,'rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['uniprot_rescue_status_v703']=='REMAIN_NO_RESIDUE_ANNOTATION':need.append(r);eids.add(r['evidence_id'])
ev={}
with gzip.open(EVID,'rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['evidence_id'] in eids:ev[r['evidence_id']]=r
pdbs=set();het=Counter();sources=Counter();pdb_counts=Counter();with_het=0
examples=[]
for r in need:
 sources[r['source_database']]+=1
 ps=[x.lower() for x in re.split(r'[;,| ]+',r.get('pdb_ids','')) if re.fullmatch(r'(?i)[0-9][a-z0-9]{3}',x)];pdbs.update(ps)
 for p in ps:pdb_counts[p]+=1
 e=ev.get(r['evidence_id'],{});h=e.get('ligand_het_id','').strip()
 if h:with_het+=1;het[h]+=1
 if len(examples)<20:examples.append({'site_id':r['binding_site_instance_id'],'evidence_id':r['evidence_id'],'target':r['target_uniprot_id'],'compound_internal_id':r['compound_internal_id'],'pdb_ids':r['pdb_ids'],'source':r['source_database'],'ligand_het_id':h,'compound_name':e.get('compound_name','')})
report={'rows':len(need),'unique_evidence':len(eids),'evidence_joined':len(ev),'unique_pdb':len(pdbs),'rows_with_ligand_het_id':with_het,'source_counts':dict(sources),'top_het_ids':dict(het.most_common(30)),'top_pdbs':dict(pdb_counts.most_common(30)),'examples':examples}
(ROOT/'V703_NO_RESIDUE_POCKET_PROFILE.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');(ROOT/'v703_no_residue_pdb_ids.txt').write_text('\n'.join(sorted(pdbs))+'\n',encoding='ascii');print(json.dumps({k:v for k,v in report.items() if k!='examples'},ensure_ascii=False,indent=2))
