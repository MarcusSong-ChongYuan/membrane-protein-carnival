import csv,gzip,json,hashlib,shutil,os,stat,itertools
from collections import Counter
from pathlib import Path

BASE=Path(r'D:\finale\16_MemPro_V7.0.2_final_20260814');CAND=Path(r'D:\finale\17_MemPro_V7.1_candidate_20260814');FINAL=Path(r'D:\finale\18_MemPro_V7.1_final_20260814')
R=CAND/'01_release_tables';Q=CAND/'04_QA'
def rows(p):
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f,delimiter='\t')
def header(p):
 try:return next(rows(p)).keys()
 except StopIteration:
  op=gzip.open if p.suffix=='.gz' else open
  with op(p,'rt',encoding='utf-8-sig',newline='') as f:return next(csv.reader(f,delimiter='\t'),[])
def count_unique(p,key):
 n=0;s=set();dup=0;blank=0
 for r in rows(p):
  n+=1;k=r.get(key,'')
  if not k:blank+=1
  elif k in s:dup+=1
  else:s.add(k)
 return n,s,dup,blank
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def inherit():
 mapping={
 'protein_master_v7.tsv.gz':'protein_master_v71.tsv.gz','compound_master_v7.tsv.gz':'compound_master_v71.tsv.gz','compound_form_v7.tsv.gz':'compound_form_v71.tsv.gz','protein_compound_pair_v7.tsv.gz':'protein_compound_pair_v71.tsv.gz','protein_membrane_evidence_and_risk_v7.tsv.gz':'protein_membrane_evidence_and_risk_v71.tsv.gz','hpa_protein_coverage_state_v7.tsv':'hpa_protein_coverage_state_v71.tsv'}
 for src,dst in mapping.items():shutil.copy2(BASE/'01_data'/src,R/dst)
 return mapping
def fkcheck(p,field,valid):
 n=bad=blank=0
 for r in rows(p):
  n+=1;v=r.get(field,'')
  if not v:blank+=1
  elif v not in valid:bad+=1
 return {'rows':n,'bad':bad,'blank':blank}
def dictionary(root):
 out=root/'00_schema'/'DATA_DICTIONARY_V71.tsv'
 with out.open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['table','field','ordinal','storage_type','semantics'])
  for p in sorted((root/'01_release_tables').glob('*')):
   if not p.is_file():continue
   for i,col in enumerate(header(p),1):
    sem='source-retained or derived field; see approved rules'
    if col.endswith('_id'):sem='stable identifier or foreign key as defined by table registry'
    if 'source' in col:sem='source provenance retained without overwriting native identity'
    if 'status' in col or 'disposition' in col:sem='controlled status; see controlled vocabulary and module QA'
    if 'release_version' in col:sem='MemPro release/version provenance'
    w.writerow([p.name,col,i,'string',sem])
def main():
 R.mkdir(parents=True,exist_ok=True);Q.mkdir(parents=True,exist_ok=True);inherit()
 proteins_n,proteins,pdup,pblank=count_unique(R/'protein_master_v71.tsv.gz','canonical_uniprot_accession')
 compounds_n,compounds,cdup,cblank=count_unique(R/'compound_master_v71.tsv.gz','compound_internal_id')
 forms_n,forms,fdup,fblank=count_unique(R/'compound_form_v71.tsv.gz','compound_form_id')
 evidence_n,evidence,edup,eblank=count_unique(R/'positive_interaction_evidence_v71.tsv.gz','evidence_id')
 disease_n,diseases,ddup,dblank=count_unique(R/'disease_entity_master_v71.tsv.gz','canonical_disease_id')
 complexes_n,complexes,xdup,xblank=count_unique(R/'protein_complex_v71.tsv.gz','complex_target_id')
 iso_n,isos,idup,iblank=count_unique(R/'isoform_membrane_class_v71.tsv.gz','protein_isoform_entity_id')
 proc_n,procs,prdup,prblank=count_unique(R/'processed_form_membrane_class_v71.tsv.gz','processed_form_id')
 site_n,sites,sdup,sblank=count_unique(R/'interaction_site_v71.tsv.gz','site_record_id')
 qa={'release':'V7.1','status':'PASS','blocking_errors':[],
 'primary_keys':{'proteins':[proteins_n,pdup,pblank],'compounds':[compounds_n,cdup,cblank],'forms':[forms_n,fdup,fblank],'evidence':[evidence_n,edup,eblank],'diseases':[disease_n,ddup,dblank],'complexes':[complexes_n,xdup,xblank],'isoforms':[iso_n,idup,iblank],'processed_forms':[proc_n,prdup,prblank],'sites':[site_n,sdup,sblank]},'foreign_keys':{}}
 for name,p,field,valid in [
 ('evidence_target',R/'positive_interaction_evidence_v71.tsv.gz','target_uniprot_id',proteins),('evidence_compound',R/'positive_interaction_evidence_v71.tsv.gz','compound_internal_id',compounds),('pair_target',R/'protein_compound_pair_v71.tsv.gz','target_uniprot_id',proteins),('pair_compound',R/'protein_compound_pair_v71.tsv.gz','compound_internal_id',compounds),('site_evidence',R/'interaction_site_v71.tsv.gz','evidence_id',evidence),('site_target',R/'interaction_site_v71.tsv.gz','target_uniprot_id',proteins),('classification_target',R/'protein_cross_classification_v71.tsv.gz','target_uniprot_id',proteins),('disease_target',R/'protein_disease_relation_v71.tsv.gz','target_uniprot_id',proteins),('disease_identity',R/'protein_disease_relation_v71.tsv.gz','canonical_disease_id',diseases),('expression_target',R/'expression_measurement_v71.tsv.gz','target_uniprot_id',proteins),('localization_target',R/'subcellular_localization_v71.tsv.gz','target_uniprot_id',proteins),('complex_component',R/'complex_component_v71.tsv.gz','complex_target_id',complexes)]:qa['foreign_keys'][name]=fkcheck(p,field,valid)
 # Business-rule checks.
 business=Counter();negative_counts={}
 for r in rows(R/'interaction_site_v71.tsv.gz'):
  business['site_rows']+=1
  if r['chain_assignment_status'].startswith('MULTIPLE_') and r['selected_chain']:business['multiple_chain_with_forced_selection']+=1
  if r['site_residue_status']=='NO_RESIDUE_REPORTED':business['no_residue_retained']+=1
  if r['structure_ligand_relationship']=='EXACT_COCRYSTAL_LIGAND':business['exact_cocrystal']+=1
 for p,k in [(R/'contextual_negative_evidence_v71.tsv.gz','contextual'),(CAND/'02_archive'/'archived_isolated_negative_evidence_v71.tsv.gz','archive')]:negative_counts[k]=sum(1 for _ in rows(p))
 qa['business_rules']=dict(business);qa['negative_partition']=negative_counts
 expected={'proteins':7800,'compounds':646670,'forms':649792,'evidence':942455,'sites':55610,'isoforms':9240,'processed_forms':8454,'complexes':4584}
 observed={'proteins':proteins_n,'compounds':compounds_n,'forms':forms_n,'evidence':evidence_n,'sites':site_n,'isoforms':iso_n,'processed_forms':proc_n,'complexes':complexes_n}
 qa['expected_counts']=expected;qa['observed_counts']=observed
 for k,v in expected.items():
  if observed[k]!=v:qa['blocking_errors'].append(f'count:{k}:{observed[k]}!={v}')
 for k,v in qa['primary_keys'].items():
  if v[1] or v[2]:qa['blocking_errors'].append(f'pk:{k}:duplicates={v[1]}:blank={v[2]}')
 for k,v in qa['foreign_keys'].items():
  if v['bad']:qa['blocking_errors'].append(f'fk:{k}:bad={v["bad"]}')
 if business['multiple_chain_with_forced_selection']:qa['blocking_errors'].append('multiple_chain_forced_selection')
 if business['exact_cocrystal']!=663:qa['blocking_errors'].append(f'exact_cocrystal:{business["exact_cocrystal"]}')
 if sum(negative_counts.values())!=1936309:qa['blocking_errors'].append('negative_partition')
 qa['status']='FROZEN_PASS' if not qa['blocking_errors'] else 'FAILED'
 (Q/'V71_FINAL_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
 if qa['blocking_errors']:
  print(json.dumps(qa,ensure_ascii=False));raise SystemExit(2)
 dictionary(CAND)
 # Build clean final directory, excluding disposable sqlite work.
 if FINAL.exists():
  raise SystemExit(f'final directory already exists: {FINAL}')
 FINAL.mkdir(parents=True)
 for d in ['00_baseline','00_docs','00_schema','01_release_tables','02_archive','04_QA']:
  shutil.copytree(CAND/d,FINAL/d)
 scripts=FINAL/'05_scripts';scripts.mkdir()
 for name in ['v71_inventory_baseline.py','v71_profile_modules.py','v71_stage1_schema.py','v71_stage2_sites.py','v71_stage2_correct_sites.py','v71_stage3_entities.py','v71_stage4_expression_negative.py','v71_stage4_run_safe_v2.py','v71_stage5_lineage_class_disease.py','v71_stage6_finalize.py']:
  p=Path(r'C:\Users\Administrator')/name
  if p.exists():shutil.copy2(p,scripts/name)
 # Manifest excludes itself and freeze marker.
 rec=[]
 for p in sorted(FINAL.rglob('*')):
  if p.is_file() and p.name not in {'MANIFEST.tsv','MANIFEST_VERIFICATION.json','FREEZE_V71_FINAL.json'}:rec.append({'relative_path':p.relative_to(FINAL).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
 with (FINAL/'MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['relative_path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(rec)
 ver={'checked':0,'missing':0,'mismatched':0,'status':'PASS'}
 for r in rec:
  p=FINAL/r['relative_path'];ver['checked']+=1
  if not p.exists():ver['missing']+=1
  elif p.stat().st_size!=r['bytes'] or sha(p)!=r['sha256']:ver['mismatched']+=1
 if ver['missing'] or ver['mismatched']:ver['status']='FAIL'
 (FINAL/'MANIFEST_VERIFICATION.json').write_text(json.dumps(ver,indent=2),encoding='utf-8')
 freeze={'release':'MemPro V7.1','status':'FROZEN_PASS','baseline':'MemPro V7.0.2','manual_sampling':'WAIVED_BY_USER','blocking_errors':0,'manifest_entries':len(rec),'manifest_verification':ver['status'],'limitations':'Source-insufficient isoform, processed-form, chain, assembly and ligand identity remain explicit; prediction is not promoted to confirmed evidence.'}
 (FINAL/'FREEZE_V71_FINAL.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2),encoding='utf-8')
 for p in FINAL.rglob('*'):
  if p.is_file():p.chmod(stat.S_IREAD)
 print(json.dumps({'final':str(FINAL),'qa':qa['status'],'manifest':ver,'files':len(rec)},ensure_ascii=False))
if __name__=='__main__':main()
