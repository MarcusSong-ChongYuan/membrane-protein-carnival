from pathlib import Path
import csv,gzip,json,re,collections,hashlib
V62=Path(r'D:\finale\01_正式数据_V6.2');OUT=Path(r'D:\7.22\v64_candidate_working');valid={x.lower() for x in json.load(open(OUT/'qa'/'rcsb_current_entry_ids_20260812.json'))};pat=re.compile(r'^[0-9][a-z0-9]{3}$');unit=re.compile(r'^\d+(?:\.\d+)?(?:pm|nm|um|µm|mm|m)$',re.I)
def clean(v):
 vals=[z.strip() for z in re.split(r'[;,|]',v or '') if z.strip()];good=[];bad=[]
 for z in vals:
  zl=z.lower()
  if zl in valid:good.append(zl)
  else:bad.append((z,'unit_like' if unit.fullmatch(z) else 'not_current_rcsb' if pat.fullmatch(zl) else 'invalid_format'))
 return sorted(set(good)),bad
# Evidence cleaned candidate
ep=V62/'binding_evidence_master_v6_2.tsv.gz';eo=OUT/'candidate_tables'/'binding_evidence_master_v6_4_candidate.tsv.gz';review=OUT/'review_queues'/'pdb_identity_review_v64.tsv'
counts=collections.Counter();src=collections.Counter();affected=[]
with gzip.open(ep,'rt',encoding='utf-8-sig',newline='') as f,gzip.open(eo,'wt',encoding='utf-8',newline='') as g:
 rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames+['pdb_ids_original_v64','pdb_identity_status_v64','structural_evidence_action_v64','release_version_v64'];w=csv.DictWriter(g,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for r in rd:
  good,bad=clean(r.get('pdb_ids',''));orig=r.get('pdb_ids','');action='unchanged';status='valid_or_empty'
  if bad:
   counts['affected_rows']+=1;counts['invalid_tokens']+=len(bad);src[r.get('source_database','')]+=1;status='cleaned_partial' if good else 'all_pdb_ids_invalid';action='retain_nonstructural_relation_review_structure'
   if r.get('evidence_type') in {'experimental_structure_residue_contact','compound_specific_structure_context','compound_specific_structure','experimental_complex_with_affinity'} and not good:
    action='exclude_from_default_structural_layer_pending_review';r['default_release_inclusion']='0';r['record_qc_status']='review';counts['structural_default_excluded']+=1
   affected.append([r.get('evidence_id',''),r.get('target_uniprot_id',''),r.get('compound_internal_id',''),r.get('source_database',''),orig,';'.join(good),';'.join(f'{a}:{b}' for a,b in bad),r.get('evidence_type',''),action])
  r['pdb_ids_original_v64']=orig if bad else '';r['pdb_ids']=';'.join(good);r['pdb_identity_status_v64']=status;r['structural_evidence_action_v64']=action;r['release_version_v64']='MemPro V6.4 candidate';w.writerow(r);counts['rows']+=1
with open(review,'w',encoding='utf-8',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['evidence_id','target_uniprot_id','compound_internal_id','source_database','pdb_ids_original','pdb_ids_valid','invalid_tokens','evidence_type','action']);w.writerows(affected)
# Sites cleaned, propagate evidence action
sp=V62/'binding_site_instances_v6_2.tsv.gz';so=OUT/'candidate_tables'/'binding_site_instances_v6_4_candidate.tsv.gz';sc=collections.Counter()
with gzip.open(sp,'rt',encoding='utf-8-sig',newline='') as f,gzip.open(so,'wt',encoding='utf-8',newline='') as g:
 rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames+['pdb_ids_original_v64','pdb_identity_status_v64','release_version_v64'];w=csv.DictWriter(g,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for r in rd:
  good,bad=clean(r.get('pdb_ids',''));orig=r.get('pdb_ids','');r['pdb_ids_original_v64']=orig if bad else '';r['pdb_ids']=';'.join(good);r['pdb_identity_status_v64']='cleaned_partial' if bad and good else 'all_pdb_ids_invalid' if bad else 'valid_or_empty';r['release_version_v64']='MemPro V6.4 candidate';w.writerow(r);sc['rows']+=1;sc['affected_rows']+=bool(bad);sc['invalid_tokens']+=len(bad)
report={'rcsb_snapshot':str(OUT/'qa'/'rcsb_current_entry_ids_20260812.json'),'rcsb_entry_count':len(valid),'evidence':dict(counts),'binding_sites':dict(sc),'affected_by_source':dict(src),'policy':'preserve original; keep only current RCSB IDs in formal pdb_ids; structural evidence with zero valid PDB is excluded from default structural layer and retained for review'}
(OUT/'qa'/'V64_PDB_IDENTITY_CLEANING.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
