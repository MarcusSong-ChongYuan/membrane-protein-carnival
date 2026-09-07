import csv,gzip,json,re,time
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites");INPUT=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813\02_frozen_review\binding_site_coordinate_frozen_v7.tsv.gz");MAP=ROOT/'sifts_residue_mapping_v7.checkpoint.tsv';OUT=ROOT/'binding_site_coordinate_rescued_v702.tsv.gz';REMAIN=ROOT/'binding_site_coordinate_remaining_v702.tsv.gz'
aa1={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'};aa3=set(aa1.values())
def parse_tokens(text):
 out=[]
 for m in re.finditer(r'(?:(?P<chain>[A-Za-z0-9]+):)?(?P<aa>[A-Z]{1,3})(?P<num>-?\d+)(?P<icode>[A-Za-z]?)\b',text or ''):
  aa=aa1.get(m.group('aa'),m.group('aa'))
  if aa in aa3:out.append({'chain':m.group('chain') or '','aa':aa,'num':m.group('num'),'icode':(m.group('icode') or '').upper()})
 return out
def parse_pdbs(text):return [x.lower() for x in re.split(r'[;,| ]+',text or '') if re.fullmatch(r'(?i)[0-9][a-z0-9]{3}',x)]
sites=[];interest=defaultdict(list)
with gzip.open(INPUT,'rt',encoding='utf-8',newline='') as f:
 rd=csv.DictReader(f,delimiter='\t');base_fields=rd.fieldnames
 for r in rd:
  ts=parse_tokens(r.get('residue_or_site_description',''));ps=parse_pdbs(r.get('pdb_ids',''));explicit={x['chain'] for x in ts if x['chain']};lookup=defaultdict(list)
  for i,t in enumerate(ts):lookup[(t['num'],t['icode'],t['aa'])].append(i)
  idx=len(sites);sites.append({'row':r,'tokens':ts,'pdbs':ps,'explicit':explicit,'lookup':lookup,'chains':set(),'matches':defaultdict(dict)})
  if ts:
   for p in ps:interest[(p,r['target_uniprot_id'])].append(idx)
print(json.dumps({'sites':len(sites),'interest_pdb_uniprot_pairs':len(interest)}),flush=True)
start=time.time();n=0;relevant=0
with MAP.open(encoding='utf-8',newline='') as f:
 for m in csv.DictReader(f,delimiter='\t'):
  n+=1
  if n%5000000==0:print(f'sifts_rows={n} relevant={relevant} elapsed={time.time()-start:.1f}s',flush=True)
  if m.get('mapping_status')!='mapped':continue
  key=(m['pdb_id'].lower(),m['uniprot_acc']);idxs=interest.get(key)
  if not idxs:continue
  relevant+=1;c=m['chain'];num=m.get('pdb_resnum','');icode=(m.get('pdb_icode') or '').upper();aa=m.get('pdb_resname','').upper()
  for idx in idxs:
   s=sites[idx]
   if s['explicit'] and c not in s['explicit']:continue
   s['chains'].add((key[0],c));token_ids=s['lookup'].get((num,icode,aa),[])
   if not token_ids and not icode:token_ids=s['lookup'].get((num,'',aa),[])
   for ti in token_ids:
    t=s['tokens'][ti]
    if t['chain'] and t['chain']!=c:continue
    s['matches'][(key[0],c)][ti]=m
extra=['rescue_status_v702','rescue_rule_v702','rescue_pdb_v702','rescue_chain_v702','rescue_requested_count_v702','rescue_mapped_count_v702','rescue_observed_count_v702','rescue_pdb_residues_v702','rescue_uniprot_residues_v702'];fields=base_fields+extra;counts=Counter();examples=defaultdict(list)
with gzip.open(OUT,'wt',encoding='utf-8',newline='') as fo,gzip.open(REMAIN,'wt',encoding='utf-8',newline='') as fr:
 wo=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wo.writeheader();wr.writeheader()
 for s in sites:
  r=s['row'];req=len(s['tokens']);scores=[]
  for (p,c),matched in s['matches'].items():
   obs=sum(x.get('observed')=='observed' for x in matched.values());scores.append((len(matched),obs,p,c,matched))
  scores.sort(key=lambda x:(x[0],x[1]),reverse=True);best=scores[0] if scores else (0,0,'','',{});ties=[x for x in scores if (x[0],x[1])==(best[0],best[1])]
  if not req:status='REMAIN_NO_RESIDUE_ANNOTATION';rule='PDB context only; co-crystal ligand coordinates required'
  elif best[0]==0:status='REMAIN_NO_AUTHOR_COORDINATE_MATCH';rule='no PDB author-number plus residue-name match on a SIFTS target-UniProt chain'
  elif best[0]<req:status='REMAIN_PARTIAL_AUTHOR_COORDINATE_MATCH';rule='only subset of requested residues matched'
  elif len(ties)==1:status='PUBLIC_RESCUED_AUTHOR_COORDINATE';rule='unique all-residue SIFTS author-coordinate match'
  else:
   sigs={tuple((i,z.get('uniprot_resnum',''),z.get('pdb_resname',''),z.get('observed','')) for i,z in sorted(x[4].items())) for x in ties}
   if len(sigs)==1:best=sorted(ties,key=lambda x:(x[2],x[3]))[0];status='PUBLIC_RESCUED_EQUIVALENT_CHAIN_REPRESENTATIVE';rule='tied chains have identical target-UniProt mapping; deterministic representative'
   else:status='REMAIN_NON_EQUIVALENT_CHAIN_TIE';rule='multiple non-equivalent chains tie'
  mapped=[best[4][i] for i in sorted(best[4])];r.update({'rescue_status_v702':status,'rescue_rule_v702':rule,'rescue_pdb_v702':best[2],'rescue_chain_v702':best[3],'rescue_requested_count_v702':req,'rescue_mapped_count_v702':best[0],'rescue_observed_count_v702':best[1],'rescue_pdb_residues_v702':';'.join(f"{z.get('pdb_resname','')}{z.get('pdb_resnum','')}{z.get('pdb_icode','')}" for z in mapped),'rescue_uniprot_residues_v702':';'.join(z.get('uniprot_resnum','') for z in mapped)});counts[status]+=1
  if len(examples[status])<5:examples[status].append({k:r.get(k,'') for k in ['binding_site_instance_id','target_uniprot_id','source_database','pdb_ids','residue_or_site_description','rescue_pdb_v702','rescue_chain_v702','rescue_mapped_count_v702']})
  (wo if status.startswith('PUBLIC') else wr).writerow(r)
report={'input_rows':len(sites),'sifts_rows_scanned':n,'relevant_sifts_rows':relevant,'status_counts':dict(counts),'rescued_rows':sum(v for k,v in counts.items() if k.startswith('PUBLIC')),'remaining_rows':sum(v for k,v in counts.items() if k.startswith('REMAIN')),'runtime_seconds':round(time.time()-start,1),'algorithm':'streaming inverted site index; no full SIFTS mapping retained in RAM','examples':dict(examples)};(ROOT/'T26_T28_FROZEN_RESCUE_V702_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='examples'},ensure_ascii=False,indent=2))
