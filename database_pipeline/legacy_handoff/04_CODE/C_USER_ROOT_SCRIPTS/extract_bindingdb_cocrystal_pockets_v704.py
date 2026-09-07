import csv,gzip,hashlib,json,os,re,tempfile,time,urllib.request,urllib.error
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import gemmi
ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites");INPUT=ROOT/'binding_site_coordinate_remaining_v703.tsv.gz';MAP=ROOT/'sifts_residue_mapping_v7.checkpoint.tsv';CACHE=ROOT/'v704_cocrystal';CACHE.mkdir(exist_ok=True);STATUS=CACHE/'pdb_status.tsv';CAND=CACHE/'site_pdb_candidates.tsv.gz';REPORT=ROOT/'T26_T28_COCRYSTAL_POCKET_V704_REPORT.json';OUT=ROOT/'binding_site_coordinate_rescued_pocket_v704.tsv.gz';REMAIN=ROOT/'binding_site_coordinate_remaining_v704.tsv.gz'
sites=[];by_pdb=defaultdict(list);interest=set()
with gzip.open(INPUT,'rt',encoding='utf-8',newline='') as f:
 rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames
 for r in rd:
  if r['uniprot_rescue_status_v703']!='REMAIN_NO_RESIDUE_ANNOTATION':continue
  idx=len(sites);sites.append(r);het=r.get('ligand_het_id','').strip().upper()
  for p in [x.lower() for x in re.split(r'[;,| ]+',r.get('pdb_ids','')) if re.fullmatch(r'(?i)[0-9][a-z0-9]{3}',x)]:by_pdb[p].append((idx,het,r['target_uniprot_id']));interest.add((p,r['target_uniprot_id']))
chains=defaultdict(set);n=0
with MAP.open(encoding='utf-8',newline='') as f:
 for m in csv.DictReader(f,delimiter='\t'):
  n+=1;key=(m['pdb_id'].lower(),m['uniprot_acc'])
  if key in interest and m.get('chain'):chains[key].add(m['chain'])
print(json.dumps({'sites':len(sites),'unique_pdb':len(by_pdb),'sifts_rows_scanned':n,'pdb_target_with_chains':len(chains)}),flush=True)
if not STATUS.exists():
 with STATUS.open('w',encoding='utf-8',newline='') as f:csv.writer(f,delimiter='\t').writerow(['pdb_id','status','http_or_error','cif_sha256','resolution','candidate_rows'])
done={}
with STATUS.open(encoding='utf-8') as f:
 for r in csv.DictReader(f,delimiter='\t'):done[r['pdb_id'].lower()]=r
pending=[p for p in sorted(by_pdb) if p not in done]
def process(pdb):
 url=f'https://files.rcsb.org/download/{pdb.upper()}.cif.gz';tmp=None
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'MemPro-V7.0.1-pocket-rescue'})
  with urllib.request.urlopen(req,timeout=90) as resp:blob=resp.read()
  sha=hashlib.sha256(blob).hexdigest();fd,name=tempfile.mkstemp(suffix='.cif.gz');os.close(fd);tmp=Path(name);tmp.write_bytes(blob);st=gemmi.read_structure(str(tmp));model=st[0];resolution=float(st.resolution) if st.resolution and st.resolution>0 else 999.0;ns=gemmi.NeighborSearch(model,st.cell,5.0).populate();results=[]
  for idx,het,target in by_pdb[pdb]:
   target_chains=chains.get((pdb,target),set());ligands=[]
   for ch in model:
    for res in ch:
     if res.name.upper()==het and res.het_flag!='A':ligands.append((ch.name,res))
   best=None
   for lchain,lig in ligands:
    contacts={}
    for atom in lig:
     for mark in ns.find_atoms(atom.pos,radius=5.0):
      cra=mark.to_cra(model)
      if cra.chain.name not in target_chains or cra.residue.het_flag!='A':continue
      key=(cra.chain.name,cra.residue.name,str(cra.residue.seqid.num),cra.residue.seqid.icode.strip())
      contacts[key]=min(contacts.get(key,99.0),atom.pos.dist(cra.atom.pos))
    if contacts:
     rec={'site_index':idx,'pdb_id':pdb,'het_id':het,'target_uniprot':target,'target_chains':';'.join(sorted(target_chains)),'ligand_chain':lchain,'ligand_resnum':str(lig.seqid.num),'ligand_icode':lig.seqid.icode.strip(),'contact_count':len(contacts),'contact_residues':';'.join(f'{c}:{aa}{num}{ic}' for c,aa,num,ic in sorted(contacts)),'min_contact_A':round(min(contacts.values()),3),'resolution_A':resolution,'cif_sha256':sha,'candidate_status':'EXACT_HET_TARGET_CHAIN_CONTACT'}
     if best is None or (rec['contact_count'],-rec['min_contact_A'])>(best['contact_count'],-best['min_contact_A']):best=rec
   if best:results.append(best)
   else:results.append({'site_index':idx,'pdb_id':pdb,'het_id':het,'target_uniprot':target,'target_chains':';'.join(sorted(target_chains)),'ligand_chain':'','ligand_resnum':'','ligand_icode':'','contact_count':0,'contact_residues':'','min_contact_A':'','resolution_A':resolution,'cif_sha256':sha,'candidate_status':'HET_NOT_FOUND' if not ligands else ('TARGET_CHAIN_NOT_MAPPED' if not target_chains else 'NO_TARGET_CHAIN_CONTACT_WITHIN_5A')})
  return pdb,'ok','',sha,resolution,results
 except urllib.error.HTTPError as e:return pdb,'http_error',str(e.code),'',999.0,[]
 except Exception as e:return pdb,'failed',repr(e)[:300],'',999.0,[]
 finally:
  if tmp:tmp.unlink(missing_ok=True)

candidate_fields=['site_index','pdb_id','het_id','target_uniprot','target_chains','ligand_chain','ligand_resnum','ligand_icode','contact_count','contact_residues','min_contact_A','resolution_A','cif_sha256','candidate_status']
existing=[]
if CAND.exists() and CAND.stat().st_size:
 with gzip.open(CAND,'rt',encoding='utf-8',newline='') as f:existing=list(csv.DictReader(f,delimiter='\t'))
all_candidates=existing;batch=100
for start in range(0,len(pending),batch):
 part=pending[start:start+batch];results=[]
 with ThreadPoolExecutor(max_workers=8) as ex:
  futs=[ex.submit(process,p) for p in part]
  for fut in as_completed(futs):results.append(fut.result())
 with STATUS.open('a',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t')
  for pdb,st,msg,sha,res,cands in sorted(results):w.writerow([pdb.upper(),st,msg,sha,res,len(cands)]);all_candidates.extend(cands)
 with gzip.open(CAND,'wt',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=candidate_fields,delimiter='\t');w.writeheader();w.writerows(all_candidates)
 print(f'checkpoint {min(start+len(part),len(pending))}/{len(pending)}',flush=True)
best=defaultdict(list)
for r in all_candidates:
 if r['candidate_status']=='EXACT_HET_TARGET_CHAIN_CONTACT':best[int(r['site_index'])].append(r)
extra=['pocket_rescue_status_v704','pocket_rescue_rule_v704','pocket_pdb_v704','pocket_target_chains_v704','pocket_ligand_het_v704','pocket_ligand_chain_v704','pocket_ligand_resnum_v704','pocket_contact_count_v704','pocket_contact_residues_v704','pocket_resolution_A_v704','pocket_candidate_structure_count_v704'];fields=base+extra;counts=Counter()
with gzip.open(OUT,'wt',encoding='utf-8',newline='') as fo,gzip.open(REMAIN,'wt',encoding='utf-8',newline='') as fr:
 wo=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wo.writeheader();wr.writeheader()
 for idx,r in enumerate(sites):
  candidates=best.get(idx,[]);candidates.sort(key=lambda x:(float(x['resolution_A']) if float(x['resolution_A'])<900 else 999,-int(x['contact_count']),x['pdb_id']));b=candidates[0] if candidates else None
  if b:status='PUBLIC_RESCUED_EXACT_COCRYSTAL_POCKET';rule='exact BindingDB HET ID present; SIFTS target chain has protein contacts within 5A; best-resolution representative selected'
  else:status='REMAIN_NO_EXACT_COCRYSTAL_TARGET_POCKET';rule='no candidate PDB satisfied exact HET plus target-UniProt-chain contact rule'
  r.update({'pocket_rescue_status_v704':status,'pocket_rescue_rule_v704':rule,'pocket_pdb_v704':b['pdb_id'] if b else '','pocket_target_chains_v704':b['target_chains'] if b else '','pocket_ligand_het_v704':b['het_id'] if b else r.get('ligand_het_id',''),'pocket_ligand_chain_v704':b['ligand_chain'] if b else '','pocket_ligand_resnum_v704':b['ligand_resnum'] if b else '','pocket_contact_count_v704':b['contact_count'] if b else 0,'pocket_contact_residues_v704':b['contact_residues'] if b else '','pocket_resolution_A_v704':b['resolution_A'] if b else '','pocket_candidate_structure_count_v704':len(candidates)});counts[status]+=1;(wo if status.startswith('PUBLIC') else wr).writerow(r)
report={'input_sites':len(sites),'unique_pdb':len(by_pdb),'completed_pdb_status_rows':sum(1 for _ in csv.DictReader(STATUS.open(encoding='utf-8'),delimiter='\t')),'candidate_status':dict(Counter(r['candidate_status'] for r in all_candidates)),'site_status':dict(counts),'rescued_sites':counts['PUBLIC_RESCUED_EXACT_COCRYSTAL_POCKET'],'remaining_sites':counts['REMAIN_NO_EXACT_COCRYSTAL_TARGET_POCKET'],'policy':'exact HET ID + SIFTS target chain + <=5A protein contacts; streaming mmCIF deleted after parse'};REPORT.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
