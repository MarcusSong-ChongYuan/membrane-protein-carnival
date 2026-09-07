import csv,gzip,json,re
from collections import defaultdict,Counter
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites");SITES=ROOT/'binding_site_coordinate_strata_v7.tsv.gz';MAP=ROOT/'sifts_residue_mapping_v7.checkpoint.tsv'
idx=defaultdict(dict)
with MAP.open(encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r.get('mapping_status')!='mapped' or not r.get('uniprot_resnum'):continue
  idx[(r['pdb_id'].lower(),r['chain'],r['uniprot_acc'])][r['uniprot_resnum']] = r
counts=Counter();public=ROOT/'binding_site_coordinate_verified_v7.tsv.gz';review=ROOT/'binding_site_coordinate_frozen_v7.tsv.gz'
with gzip.open(SITES,'rt',encoding='utf-8',newline='') as fi:
 rd=csv.DictReader(fi,delimiter='\t');extra=['actual_pdb_used_v7','actual_chain_used_v7','requested_residue_count_v7','mapped_residue_count_v7','observed_residue_count_v7','matched_uniprot_residues_v7','mapped_pdb_residues_v7','coordinate_mapping_status_v7','coordinate_mapping_confidence_v7','coordinate_mapping_reason_v7'];fields=list(rd.fieldnames)+extra
 with gzip.open(public,'wt',encoding='utf-8',newline='') as fp,gzip.open(review,'wt',encoding='utf-8',newline='') as fr:
  wp=csv.DictWriter(fp,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wp.writeheader();wr.writeheader()
  for r in rd:
   if r['coordinate_qc_stratum_v7']!='SIFTS_RESIDUE_MAPPING_REQUIRED':
    r.update({'coordinate_mapping_status_v7':'FROZEN_REVIEW','coordinate_mapping_confidence_v7':'NONE','coordinate_mapping_reason_v7':r['coordinate_qc_stratum_v7']});wr.writerow(r);counts[r['coordinate_qc_stratum_v7']]+=1;continue
   pdbs=[x.lower() for x in re.split(r'[;,| ]+',r.get('pdb_ids','')) if re.fullmatch(r'[0-9][A-Za-z0-9]{3}',x)];chains=[x for x in re.split(r'[;,| ]+',r.get('pdb_chain_ids_v60','')) if x];req=[x for x in r.get('requested_uniprot_residue_tokens_v7','').split(';') if x];acc=r['target_uniprot_id'];scores=[]
   for p in pdbs:
    for c in chains:
     m=idx.get((p,c,acc),{});hits=[x for x in req if re.match(r'\d+',x) and re.match(r'\d+',x).group() in m];obs=[x for x in hits if m[re.match(r'\d+',x).group()].get('observed')=='observed'];scores.append((len(hits),len(obs),p,c,hits,m))
   scores.sort(reverse=True,key=lambda x:(x[0],x[1]));best=scores[0] if scores else (0,0,'','',[],{});ties=[x for x in scores if (x[0],x[1])==(best[0],best[1])]
   mapped=[]
   for x in best[4]:
    z=best[5][re.match(r'\d+',x).group()];mapped.append(f"{z.get('pdb_resname','')}{z.get('pdb_resnum','')}{z.get('pdb_icode','')}")
   r.update({'actual_pdb_used_v7':best[2],'actual_chain_used_v7':best[3],'requested_residue_count_v7':len(req),'mapped_residue_count_v7':best[0],'observed_residue_count_v7':best[1],'matched_uniprot_residues_v7':';'.join(best[4]),'mapped_pdb_residues_v7':';'.join(mapped)})
   if best[0]==0:status='FROZEN_REVIEW';conf='NONE';reason='NO_SIFTS_RESIDUE_MATCH'
   elif len(ties)>1:status='FROZEN_REVIEW';conf='REVIEW';reason='TIED_BEST_PDB_CHAIN_MAPPING'
   elif best[0]==len(req) and best[1]==len(req):status='PUBLIC_COORDINATE_VERIFIED';conf='HIGH';reason='ALL_REQUESTED_RESIDUES_MAPPED_AND_OBSERVED'
   else:status='PUBLIC_COORDINATE_PARTIAL';conf='MEDIUM';reason='UNIQUE_BEST_MAPPING_PARTIAL_OR_UNOBSERVED'
   r.update({'coordinate_mapping_status_v7':status,'coordinate_mapping_confidence_v7':conf,'coordinate_mapping_reason_v7':reason});counts[status]+=1;(wp if status.startswith('PUBLIC') else wr).writerow(r)
report={'mapping_index_keys':len(idx),'site_status_counts':dict(counts),'policy':'Unique best PDB-chain only; zero-match and tied-best mappings frozen; HIGH requires all requested residues mapped and observed.'}
(ROOT/'T26_T28_SIFTS_COORDINATE_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
