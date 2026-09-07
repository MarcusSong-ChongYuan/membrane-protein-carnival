import csv,gzip,json,re,time
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites");INPUT=ROOT/'binding_site_coordinate_remaining_v702.tsv.gz';MAP=ROOT/'sifts_residue_mapping_v7.checkpoint.tsv';OUT=ROOT/'binding_site_coordinate_rescued_uniprot_v703.tsv.gz';REMAIN=ROOT/'binding_site_coordinate_remaining_v703.tsv.gz'
aa1={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'};aa3=set(aa1.values())
def toks(text):
 out=[]
 for m in re.finditer(r'(?:(?P<chain>[A-Za-z0-9]+):)?(?P<aa>[A-Z]{1,3})(?P<num>-?\d+)(?P<icode>[A-Za-z]?)\b',text or ''):
  aa=aa1.get(m.group('aa'),m.group('aa'))
  if aa in aa3:out.append({'chain':m.group('chain') or '','aa':aa,'num':m.group('num')})
 return out
def pdbs(x):return [p.lower() for p in re.split(r'[;,| ]+',x or '') if re.fullmatch(r'(?i)[0-9][a-z0-9]{3}',p)]
sites=[];interest=defaultdict(list)
with gzip.open(INPUT,'rt',encoding='utf-8',newline='') as f:
 rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames
 for r in rd:
  ts=toks(r.get('residue_or_site_description',''));ps=pdbs(r.get('pdb_ids',''));explicit={t['chain'] for t in ts if t['chain']};lookup=defaultdict(list)
  for i,t in enumerate(ts):lookup[(t['num'],t['aa'])].append(i)
  idx=len(sites);sites.append({'row':r,'tokens':ts,'explicit':explicit,'lookup':lookup,'matches':defaultdict(dict)})
  if ts:
   for p in ps:interest[(p,r['target_uniprot_id'])].append(idx)
print(json.dumps({'remaining_sites':len(sites),'interest_pairs':len(interest)}),flush=True)
n=rel=0;start=time.time()
with MAP.open(encoding='utf-8',newline='') as f:
 for m in csv.DictReader(f,delimiter='\t'):
  n+=1
  if n%5000000==0:print(f'sifts_rows={n} relevant={rel} elapsed={time.time()-start:.1f}s',flush=True)
  if m.get('mapping_status')!='mapped' or not m.get('uniprot_resnum'):continue
  key=(m['pdb_id'].lower(),m['uniprot_acc']);idxs=interest.get(key)
  if not idxs:continue
  rel+=1;c=m['chain'];tokenkey=(m['uniprot_resnum'],m.get('pdb_resname','').upper())
  for idx in idxs:
   s=sites[idx]
   if s['explicit'] and c not in s['explicit']:continue
   for ti in s['lookup'].get(tokenkey,[]):
    t=s['tokens'][ti]
    if not t['chain'] or t['chain']==c:s['matches'][(key[0],c)][ti]=m
extra=['uniprot_rescue_status_v703','uniprot_rescue_rule_v703','uniprot_rescue_pdb_v703','uniprot_rescue_chain_v703','uniprot_requested_count_v703','uniprot_mapped_count_v703','uniprot_observed_count_v703','uniprot_mapped_pdb_residues_v703'];fields=base+extra;counts=Counter()
with gzip.open(OUT,'wt',encoding='utf-8',newline='') as fo,gzip.open(REMAIN,'wt',encoding='utf-8',newline='') as fr:
 wo=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wo.writeheader();wr.writeheader()
 for s in sites:
  r=s['row'];req=len(s['tokens']);scores=[]
  for (p,c),matched in s['matches'].items():scores.append((len(matched),sum(x.get('observed')=='observed' for x in matched.values()),p,c,matched))
  scores.sort(key=lambda x:(x[0],x[1]),reverse=True);best=scores[0] if scores else (0,0,'','',{});ties=[x for x in scores if (x[0],x[1])==(best[0],best[1])]
  if not req:status='REMAIN_NO_RESIDUE_ANNOTATION';rule='requires co-crystal ligand pocket extraction'
  elif best[0]==req and len(ties)==1:status='PUBLIC_RESCUED_UNIPROT_COORDINATE';rule='all requested UniProt residues uniquely map to one SIFTS PDB chain'
  elif best[0]==req:
   sigs={tuple((i,z.get('pdb_resnum',''),z.get('pdb_icode',''),z.get('observed','')) for i,z in sorted(x[4].items())) for x in ties}
   if len(sigs)==1:best=sorted(ties,key=lambda x:(x[2],x[3]))[0];status='PUBLIC_RESCUED_UNIPROT_EQUIVALENT_CHAIN';rule='equivalent tied chains; deterministic representative'
   else:status='REMAIN_UNIPROT_NON_EQUIVALENT_CHAIN_TIE';rule='full mapping but non-equivalent chain tie'
  elif best[0]>0:status='REMAIN_UNIPROT_PARTIAL_MATCH';rule='only subset of UniProt residues mapped'
  else:status='REMAIN_NO_UNIPROT_COORDINATE_MATCH';rule='no target-chain UniProt residue match'
  mapped=[best[4][i] for i in sorted(best[4])];r.update({'uniprot_rescue_status_v703':status,'uniprot_rescue_rule_v703':rule,'uniprot_rescue_pdb_v703':best[2],'uniprot_rescue_chain_v703':best[3],'uniprot_requested_count_v703':req,'uniprot_mapped_count_v703':best[0],'uniprot_observed_count_v703':best[1],'uniprot_mapped_pdb_residues_v703':';'.join(f"{z.get('pdb_resname','')}{z.get('pdb_resnum','')}{z.get('pdb_icode','')}" for z in mapped)});counts[status]+=1;(wo if status.startswith('PUBLIC') else wr).writerow(r)
report={'input_rows':len(sites),'sifts_rows_scanned':n,'relevant_rows':rel,'status_counts':dict(counts),'rescued_rows':sum(v for k,v in counts.items() if k.startswith('PUBLIC')),'remaining_rows':sum(v for k,v in counts.items() if k.startswith('REMAIN')),'runtime_seconds':round(time.time()-start,1)};(ROOT/'T26_T28_UNIPROT_RESCUE_V703_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
