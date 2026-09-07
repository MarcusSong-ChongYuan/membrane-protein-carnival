import csv,gzip,hashlib,json
from collections import Counter
from pathlib import Path
root=Path(r"D:\finale\09_V7_data_freeze_working_20260813");rel=root/'03_v7_relation_rebuild';mod=root/'04_v7_remaining_modules';cand=root/'05_v7_automatic_release_candidate';qa=cand/'QA';qa.mkdir(exist_ok=True)
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
def load(p,key):
 vals=[]
 with op(p) as f:
  for r in csv.DictReader(f,delimiter='\t'):vals.append(r[key])
 return vals
def fkcheck(p,col,parents):
 n=bad=0
 with op(p) as f:
  for r in csv.DictReader(f,delimiter='\t'):n+=1;bad+=r[col] not in parents
 return {'rows':n,'bad':bad}
checks={};warnings={}
proteins=load(cand/'protein_master_V7.tsv.gz','canonical_uniprot_accession');ps=set(proteins);checks['protein_pk']={'rows':len(proteins),'bad':len(proteins)-len(ps)}
compounds=load(rel/'small_molecule_master_V7.tsv.gz','compound_internal_id');cs=set(compounds);checks['compound_pk']={'rows':len(compounds),'bad':len(compounds)-len(cs)}
eids=load(cand/'protein_compound_evidence_V7.tsv.gz','evidence_id');es=set(eids);checks['evidence_pk']={'rows':len(eids),'bad':len(eids)-len(es)}
pairs=load(cand/'protein_compound_pair_V7.tsv.gz','pair_id');checks['pair_pk']={'rows':len(pairs),'bad':len(pairs)-len(set(pairs))}
diseases=load(rel/'disease_entity_master_V7.tsv.gz','canonical_disease_id');ds=set(diseases);checks['disease_pk']={'rows':len(diseases),'bad':len(diseases)-len(ds)}
complexes=load(rel/'complex_target_master_V7.tsv.gz','complex_target_id');xs=set(complexes);checks['complex_pk']={'rows':len(complexes),'bad':len(complexes)-len(xs)}
forms=load(mod/'compound_form_V7.tsv.gz','compound_form_id');checks['compound_form_pk']={'rows':len(forms),'bad':len(forms)-len(set(forms))}
iso=load(mod/'protein_isoform_membrane_state_V7.tsv.gz','protein_isoform_entity_id');checks['isoform_pk']={'rows':len(iso),'bad':len(iso)-len(set(iso))}
proc=load(mod/'protein_processed_form_exact_V7.tsv.gz','processed_form_id');checks['processed_form_pk']={'rows':len(proc),'bad':len(proc)-len(set(proc))}
checks['evidence_target_fk']=fkcheck(cand/'protein_compound_evidence_V7.tsv.gz','target_uniprot_id',ps)
checks['evidence_compound_fk']=fkcheck(cand/'protein_compound_evidence_V7.tsv.gz','compound_internal_id',cs)
checks['form_compound_fk']=fkcheck(mod/'compound_form_V7.tsv.gz','compound_internal_id',cs)
checks['isoform_protein_fk']=fkcheck(mod/'protein_isoform_membrane_state_V7.tsv.gz','canonical_uniprot_accession',ps)
checks['processed_protein_fk']=fkcheck(mod/'protein_processed_form_exact_V7.tsv.gz','target_uniprot_id',ps)
checks['disease_protein_fk']=fkcheck(rel/'protein_disease_relation_V7.tsv.gz','target_uniprot_id',ps)
checks['disease_entity_fk']=fkcheck(rel/'protein_disease_relation_V7.tsv.gz','canonical_disease_id',ds)
checks['complex_component_fk']=fkcheck(rel/'complex_target_components_V7.tsv.gz','complex_target_id',xs)
checks['expression_protein_fk']=fkcheck(rel/'protein_tissue_cell_ihc_expression_V7.tsv.gz','target_uniprot_id',ps)
checks['site_evidence_fk']=fkcheck(rel/'binding_site_instances_V7.tsv.gz','evidence_id',es)

# Business rules and exact pair reconstruction.
bc=Counter();evpairs=set()
with op(cand/'protein_master_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  c=r['membrane_class_v7'];m=r['primary_membrane_mode_v7'];ev=r['decision_evidence_v7']
  if c=='A' and not any(x in m for x in ['transmembrane','intramembrane','integral','beta_barrel','membrane_isoform']):bc['A_without_integral_mode']+=1
  if c=='B' and not any(x in m for x in ['anchor','lipid']):bc['B_without_anchor_mode']+=1
  if c=='C' and not any(x in m for x in ['peripheral','complex_mediated']):bc['C_without_peripheral_mode']+=1
with op(cand/'protein_compound_evidence_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  evpairs.add((r['target_uniprot_id'],r['compound_internal_id']))
  if r['evidence_tier'] not in {'BE1','BE2','BE3'}:bc['invalid_BE_tier']+=1
  if not r['experiment_lineage_key_v7']:bc['missing_lineage_key']+=1
pairkeys=set()
with op(cand/'protein_compound_pair_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):pairkeys.add((r['target_uniprot_id'],r['compound_internal_id']))
checks['business_rules']={'rows':len(ps)+len(es),'bad':sum(bc.values()),'details':dict(bc)}
checks['pair_exact_reconstruction']={'rows':len(pairkeys),'bad':len(evpairs^pairkeys)}

# Review isolation: unresolved form records are allowed only in explicit frozen state.
bad_iso=0
with op(mod/'protein_isoform_membrane_state_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['isoform_membrane_class_v7']=='UNRESOLVED_ISOFORM' and not r['isoform_disposition_v7'].startswith('FROZEN'):bad_iso+=1
checks['isoform_review_isolation']={'rows':len(iso),'bad':bad_iso}
bad_proc=0
with op(mod/'protein_processed_form_exact_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['membrane_class_assignment']=='UNRESOLVED_PROCESSED_FORM' and not r['assignment_status'].startswith('FROZEN'):bad_proc+=1
checks['processed_review_isolation']={'rows':len(proc),'bad':bad_proc}

blocking=sum(v['bad'] for v in checks.values());report={'status':'PASS_AUTOMATIC_V7_RELEASE_CANDIDATE' if blocking==0 else 'FAIL_AUTOMATIC_V7_RELEASE_CANDIDATE','blocking_error_count':blocking,'checks':checks,'manual_review_policy':'OMITTED_BY_USER_APPROVAL; deterministic PUBLIC/FROZEN_REVIEW/EXCLUDED_AUDIT only','counts':{'proteins':len(ps),'compounds':len(cs),'forms':len(forms),'isoforms':len(iso),'processed_forms':len(proc),'evidence':len(es),'pairs':len(pairkeys),'diseases':len(ds),'complexes':len(xs)}}
(qa/'V7_AUTOMATIC_RELEASE_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
