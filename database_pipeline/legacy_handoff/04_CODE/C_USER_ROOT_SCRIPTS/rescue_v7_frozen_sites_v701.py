import csv,gzip,json,re
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites")
INPUT=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813\02_frozen_review\binding_site_coordinate_frozen_v7.tsv.gz")
MAP=ROOT/'sifts_residue_mapping_v7.checkpoint.tsv'
OUT=ROOT/'binding_site_coordinate_rescued_v701.tsv.gz';REMAIN=ROOT/'binding_site_coordinate_remaining_v701.tsv.gz'

# SIFTS author-coordinate and UniProt-coordinate indexes.
chains=defaultdict(list); author=defaultdict(dict); uniprot=defaultdict(dict)
with MAP.open(encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r.get('mapping_status')!='mapped':continue
  p=r['pdb_id'].lower();c=r['chain'];a=r['uniprot_acc'];key=(p,c,a)
  if c and c not in chains[(p,a)]:chains[(p,a)].append(c)
  pn=(r.get('pdb_resnum',''),(r.get('pdb_icode') or '').upper())
  if pn[0]:author[key][pn]=r
  if r.get('uniprot_resnum'):uniprot[key][r['uniprot_resnum']]=r

aa1={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
def tokens(text):
 out=[]
 # Explicit structural formats: A:ALA767, A:S415, GLY23, S415.
 for m in re.finditer(r'(?:(?P<chain>[A-Za-z0-9]+):)?(?P<aa>[A-Z]{1,3})(?P<num>-?\d+)(?P<icode>[A-Za-z]?)\b',text or ''):
  aa=m.group('aa');aa=aa1.get(aa,aa)
  if aa not in set(aa1.values()):continue
  out.append({'chain':m.group('chain') or '','aa':aa,'num':m.group('num'),'icode':(m.group('icode') or '').upper()})
 return out
def pdbs(text):return [x.lower() for x in re.split(r'[;,| ]+',text or '') if re.fullmatch(r'(?i)[0-9][a-z0-9]{3}',x)]
def evaluate(p,c,acc,toks):
 idx=author.get((p,c,acc),{});hits=[];observed=[]
 for t in toks:
  z=idx.get((t['num'],t['icode']))
  if not z and not t['icode']:z=idx.get((t['num'],''))
  if z and (not t['aa'] or z.get('pdb_resname','').upper()==t['aa']):
   hits.append((t,z));observed.append(z.get('observed')=='observed')
 return hits,observed

counts=Counter();examples=defaultdict(list)
with gzip.open(INPUT,'rt',encoding='utf-8',newline='') as fi:
 rd=csv.DictReader(fi,delimiter='\t');extra=['rescue_status_v701','rescue_rule_v701','rescue_pdb_v701','rescue_chain_v701','rescue_requested_count_v701','rescue_mapped_count_v701','rescue_observed_count_v701','rescue_pdb_residues_v701','rescue_uniprot_residues_v701'];fields=list(rd.fieldnames)+extra
 with gzip.open(OUT,'wt',encoding='utf-8',newline='') as fo,gzip.open(REMAIN,'wt',encoding='utf-8',newline='') as ff:
  wo=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');wf=csv.DictWriter(ff,fieldnames=fields,delimiter='\t');wo.writeheader();wf.writeheader()
  for r in rd:
   acc=r['target_uniprot_id'];ts=tokens(r.get('residue_or_site_description',''));ps=pdbs(r.get('pdb_ids',''));explicit={t['chain'] for t in ts if t['chain']};scores=[]
   for p in ps:
    candidate=list(explicit) if explicit else chains.get((p,acc),[])
    for c in candidate:
     subset=[t for t in ts if not t['chain'] or t['chain']==c];h,o=evaluate(p,c,acc,subset);scores.append((len(h),sum(o),p,c,h,o,len(subset)))
   scores.sort(key=lambda x:(x[0],x[1],x[6]),reverse=True);best=scores[0] if scores else (0,0,'','',[],[],len(ts));ties=[x for x in scores if (x[0],x[1],x[6])==(best[0],best[1],best[6])]
   status='';rule=''
   if not ts:status='REMAIN_NO_RESIDUE_ANNOTATION';rule='cannot create residue-defined site from PDB ID alone'
   elif best[0]==0:status='REMAIN_NO_AUTHOR_COORDINATE_MATCH';rule='PDB author residue number plus residue name did not match target UniProt chain'
   elif best[0]<best[6]:status='REMAIN_PARTIAL_AUTHOR_COORDINATE_MATCH';rule='only a subset of requested structural residues matched'
   elif len(ties)==1:status='PUBLIC_RESCUED_AUTHOR_COORDINATE';rule='unique target-UniProt PDB chain; all author-numbered residues and residue names matched SIFTS'
   else:
    sigs={tuple((z.get('uniprot_resnum',''),z.get('pdb_resname',''),z.get('observed','')) for _,z in x[4]) for x in ties}
    if len(sigs)==1:
     best=sorted(ties,key=lambda x:(x[2],x[3]))[0];status='PUBLIC_RESCUED_EQUIVALENT_CHAIN_REPRESENTATIVE';rule='tied chains have identical target-UniProt residue mapping; deterministic representative selected'
    else:status='REMAIN_NON_EQUIVALENT_CHAIN_TIE';rule='multiple chains tie but residue mappings are not equivalent'
   mapped=';'.join(f"{z.get('pdb_resname','')}{z.get('pdb_resnum','')}{z.get('pdb_icode','')}" for _,z in best[4]);uni=';'.join(z.get('uniprot_resnum','') for _,z in best[4])
   r.update({'rescue_status_v701':status,'rescue_rule_v701':rule,'rescue_pdb_v701':best[2],'rescue_chain_v701':best[3],'rescue_requested_count_v701':best[6],'rescue_mapped_count_v701':best[0],'rescue_observed_count_v701':best[1],'rescue_pdb_residues_v701':mapped,'rescue_uniprot_residues_v701':uni})
   counts[status]+=1
   if len(examples[status])<5:examples[status].append({k:r.get(k,'') for k in ['binding_site_instance_id','target_uniprot_id','source_database','pdb_ids','residue_or_site_description','rescue_pdb_v701','rescue_chain_v701','rescue_requested_count_v701','rescue_mapped_count_v701']})
   (wo if status.startswith('PUBLIC') else wf).writerow(r)
report={'input_rows':sum(counts.values()),'status_counts':dict(counts),'rescued_rows':sum(v for k,v in counts.items() if k.startswith('PUBLIC')),'remaining_rows':sum(v for k,v in counts.items() if k.startswith('REMAIN')),'mapping_index_chain_accession_keys':len(author),'numbering_policy':'PDB author residue number + insertion code + residue-name validation; SIFTS target UniProt chain required','examples':dict(examples)}
(ROOT/'T26_T28_FROZEN_RESCUE_V701_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='examples'},ensure_ascii=False,indent=2))
