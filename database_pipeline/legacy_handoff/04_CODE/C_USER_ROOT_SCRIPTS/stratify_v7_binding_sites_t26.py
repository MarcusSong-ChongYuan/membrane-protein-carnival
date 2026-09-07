import csv,gzip,json,re
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813");OUT=ROOT/'03_binding_sites';OUT.mkdir(exist_ok=True)
SRC=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data\binding_site_v7.tsv.gz")
VALID=set(json.load(open(r"D:\7.22\v64_candidate_working\qa\rcsb_current_entry_ids_20260812.json")))
VALID={x.lower() for x in VALID}
pdbpat=re.compile(r'^[0-9][a-z0-9]{3}$',re.I);respat=re.compile(r'(?:(?:[A-Z]{3})\s*)?(\d+)([A-Za-z]?)')
def splitp(v):return sorted({x.lower() for x in re.split(r'[;,| ]+',v or '') if pdbpat.fullmatch(x)})
counts=Counter();pdb_need=set();rows=[]
with gzip.open(SRC,'rt',encoding='utf-8',newline='') as fi:
 rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['pdb_snapshot_status_v7_final','coordinate_qc_stratum_v7','requested_uniprot_residue_tokens_v7','sifts_required_v7','direct_docking_box_eligibility_v7']
 with gzip.open(OUT/'binding_site_coordinate_strata_v7.tsv.gz','wt',encoding='utf-8',newline='') as fo:
  w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
  for r in rd:
   pdbs=splitp(r.get('pdb_ids',''));current=[x for x in pdbs if x in VALID];chains=[x for x in re.split(r'[;,| ]+',r.get('pdb_chain_ids_v60','')) if x]
   desc=r.get('residue_or_site_description','');tokens=[''.join(x) for x in respat.findall(desc)] if desc else []
   if not pdbs:st='NO_PDB_CONTEXT';elig='NO'
   elif not current:st='PDB_NOT_CURRENT_IN_FROZEN_RCSB_SNAPSHOT';elig='NO'
   elif not chains:st='CURRENT_PDB_CHAIN_MISSING';elig='NO'
   elif not tokens:st='CURRENT_PDB_CHAIN_PRESENT_RESIDUES_UNPARSEABLE_OR_ABSENT';elig='NO'
   else:st='SIFTS_RESIDUE_MAPPING_REQUIRED';elig='PENDING_SIFTS';pdb_need.update(current)
   r['pdb_snapshot_status_v7_final']='current' if current else ('no_pdb' if not pdbs else 'not_current')
   r['coordinate_qc_stratum_v7']=st;r['requested_uniprot_residue_tokens_v7']=';'.join(tokens);r['sifts_required_v7']='1' if st=='SIFTS_RESIDUE_MAPPING_REQUIRED' else '0';r['direct_docking_box_eligibility_v7']=elig;counts[st]+=1;w.writerow(r)
(OUT/'pdb_ids_requiring_sifts_v7.txt').write_text('\n'.join(sorted(pdb_need))+'\n',encoding='ascii')
report={'rows':sum(counts.values()),'strata':dict(counts),'unique_pdb_requiring_sifts':len(pdb_need),'rcsb_snapshot_entries':len(VALID),'rcsb_snapshot_date':'2026-08-12','policy':'Only current PDB + explicit chain + parseable requested residues proceeds to SIFTS; no coordinate guessing.'}
(OUT/'T26_BINDING_SITE_STRATIFICATION_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
