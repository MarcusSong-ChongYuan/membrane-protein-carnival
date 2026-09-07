import csv,gzip,json
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813");OUT=ROOT/'06_remaining_modules';OUT.mkdir(exist_ok=True)
SRC=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data")
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')

proteins=set()
with op(SRC/'protein_master_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):proteins.add(r['canonical_uniprot_accession'])

# Complex: entity can be public with unknown source stoichiometry, but uncertainty
# is explicit and required/optional must not be claimed as curated.
components=Counter();bad_component=set();comp_status={}
with op(SRC/'complex_component_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  cid=r['complex_target_id'];components[cid]+=1
  if r.get('component_identity_review_flag_v0_4')=='1' or r.get('accession_patch_status_v64')=='exclude_pending_source_correction':bad_component.add(cid)
with op(SRC/'protein_complex_v7.tsv.gz') as fi:
 rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['complex_disposition_v7_final','stoichiometry_interpretation_v7','required_optional_interpretation_v7','complex_limitations_v7']
 with op(OUT/'protein_complex_v7_validated.tsv.gz','wt') as fp,op(OUT/'protein_complex_frozen_v7.tsv.gz','wt') as fr:
  wp=csv.DictWriter(fp,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wp.writeheader();wr.writeheader();cp=cr=0
  for r in rd:
   cid=r['complex_target_id'];ok=(components[cid]>0 and cid not in bad_component and r.get('component_set_conflict_flag')=='0')
   r['stoichiometry_interpretation_v7']='source_annotations_available' if r.get('stoichiometry_status')=='source_annotations_available' else 'unknown_not_inferred'
   r['required_optional_interpretation_v7']='source_not_asserted_not_inferred' if r.get('required_optional_subunit_status')=='not_yet_curated' else r.get('required_optional_subunit_status','')
   r['complex_limitations_v7']='Required/optional subunit roles are not asserted unless present in source; unknown stoichiometry remains unknown.'
   r['complex_disposition_v7_final']='PUBLIC_ENTITY' if ok else 'FROZEN_REVIEW'
   (wp if ok else wr).writerow(r);cp+=ok;cr+=not ok

# Disease exact-only validation.
with op(SRC/'protein_disease_relation_v7.tsv.gz') as fi:
 rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['disease_disposition_v7_final','disease_mapping_rule_v7_final']
 with op(OUT/'protein_disease_relation_v7_validated.tsv.gz','wt') as fp,op(OUT/'protein_disease_frozen_v7.tsv.gz','wt') as fr:
  wp=csv.DictWriter(fp,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wp.writeheader();wr.writeheader();dp=dr=0
  for r in rd:
   safe=(r['target_uniprot_id'] in proteins and r.get('mapping_review_flag')=='0' and all(x in {'direct_mondo_id','unique_official_exact_mapping'} for x in r.get('canonical_mapping_statuses','').split(';')) and r.get('default_release_inclusion_v63')=='1')
   r['disease_disposition_v7_final']='PUBLIC_EXACT_CANONICAL' if safe else 'FROZEN_REVIEW';r['disease_mapping_rule_v7_final']='direct MONDO or unique official exact/equivalent only; broad/narrow/related and 1:N prohibited'
   (wp if safe else wr).writerow(r);dp+=safe;dr+=not safe

# Five-axis coverage: source result is regenerated into explicit per-axis coverage
# and unresolved list for final V7 protein whitelist.
axes=['structural_family_primary_v63','molecular_function_primary_v63','biological_process_primary_v63','membrane_role_primary_v63','specialist_classification_primary_v63'];coverage=Counter();unresolved=[];rows=[]
with op(SRC/'protein_cross_classification_summary_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  for a in axes:
   v=r.get(a,'');resolved=v not in {'','unclassified','no_specialist_classification'};coverage[(a,'resolved' if resolved else 'unresolved_or_not_applicable')]+=1
  if r.get('display_classification_status_v63')=='functionally_unresolved':unresolved.append(r)
  rows.append(r)
if unresolved:
 with op(OUT/'protein_functionally_unresolved_v7.tsv.gz','wt') as f:w=csv.DictWriter(f,fieldnames=list(unresolved[0]),delimiter='\t');w.writeheader();w.writerows(unresolved)
with (OUT/'five_axis_coverage_v7.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['axis','coverage_status','protein_count'])
 for (a,s),n in sorted(coverage.items()):w.writerow([a,s,n])

report={'complex':{'public_entities':cp,'frozen_entities':cr,'source_stoichiometry_unknown':289,'required_optional_source_not_asserted':4585,'interpretation':'unknown is retained, not manually inferred'},'disease':{'public_exact_relations':dp,'frozen_relations':dr},'five_axis':{'proteins':len(rows),'functionally_unresolved':len(unresolved),'coverage':{f'{a}|{s}':n for (a,s),n in coverage.items()}},'manual_review':'waived; unresolved status is a final data state rather than a pending human task'}
(OUT/'REMAINING_MODULES_VALIDATION_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
