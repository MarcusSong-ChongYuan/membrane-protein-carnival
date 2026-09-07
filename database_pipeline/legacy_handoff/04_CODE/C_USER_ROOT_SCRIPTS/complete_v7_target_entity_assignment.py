import csv,gzip,json
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813")
SRC=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data")
OLD=Path(r"D:\finale\02_V6.3_candidate_20260803")
OUT=ROOT/'01_entity_assignment'; OUT.mkdir(exist_ok=True)
def op(p,m='rt'): return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')

proteins=set(); iso_by_id={}; iso_accessions={}
with op(SRC/'protein_master_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):proteins.add(r['canonical_uniprot_accession'])
with op(SRC/'protein_isoform_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  iso_by_id[r['protein_isoform_entity_id']]=r
  iso_accessions[r['isoform_uniprot_accession']]=r['protein_isoform_entity_id']

resolution={}; dup=Counter()
with op(OLD/'02_identity'/'binding_evidence_target_resolution_v0_1.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  eid=r['evidence_id']
  if eid in resolution: dup[eid]+=1
  else: resolution[eid]=r

complex_candidates={}
with op(OLD/'03_complex'/'complex_binding_evidence_candidates_v0_2.tsv') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  sid=r.get('source_entity_id','') or r.get('source_complex_id','')
  if sid: complex_candidates.setdefault((r.get('source_database',''),sid),[]).append(r)

counts=Counter(); unresolved=[]
out=OUT/'evidence_target_entity_assignment_v7_final.tsv.gz'
with op(SRC/'protein_compound_evidence_v7.tsv.gz') as fi:
 rd=csv.DictReader(fi,delimiter='\t')
 fields=['evidence_id','source_database','source_record_id','source_target_identifier','canonical_uniprot_accession','target_entity_type_v7','target_entity_id_v7','protein_isoform_entity_id','complex_target_id','assignment_resolution_v7','assignment_disposition_v7','historical_isoform_specificity','identity_review_flag_source','assignment_limitations_v7']
 with op(out,'wt') as fo:
  w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
  for e in rd:
   eid=e['evidence_id']; old=resolution.get(eid); acc=e['target_uniprot_id']; row={'evidence_id':eid,'source_database':e['source_database'],'source_record_id':e['source_record_id'],'source_target_identifier':'','canonical_uniprot_accession':acc,'target_entity_type_v7':'','target_entity_id_v7':'','protein_isoform_entity_id':'','complex_target_id':'','assignment_resolution_v7':'','assignment_disposition_v7':'','historical_isoform_specificity':'','identity_review_flag_source':'','assignment_limitations_v7':''}
   if old:
    row['source_target_identifier']=old.get('source_target_identifier',''); row['historical_isoform_specificity']=old.get('historical_isoform_specificity',''); row['identity_review_flag_source']=old.get('identity_review_flag','')
    iid=old.get('protein_isoform_entity_id','')
    if iid and iid in iso_by_id:
     row.update({'target_entity_type_v7':'protein_isoform','target_entity_id_v7':iid,'protein_isoform_entity_id':iid,'assignment_resolution_v7':'explicit_isoform_foreign_key_recovered','assignment_disposition_v7':'PUBLIC_EXACT_ISOFORM'});counts['public_exact_isoform']+=1
    elif old.get('target_resolution_level')=='canonical_protein' and acc in proteins:
     row.update({'target_entity_type_v7':'canonical_protein','target_entity_id_v7':'CP-'+acc,'assignment_resolution_v7':'canonical_protein_resolved','assignment_disposition_v7':'PUBLIC_CANONICAL'});counts['public_canonical']+=1
     if old.get('historical_isoform_specificity')=='not_recoverable_from_v62_normalized_record':row['assignment_limitations_v7']='source granularity cannot recover historical isoform; canonical relationship only'
    else:
     row.update({'target_entity_type_v7':'unresolved_target','target_entity_id_v7':'UNRESOLVED:'+eid,'assignment_resolution_v7':'legacy_resolution_not_safe_for_current_entity','assignment_disposition_v7':'FROZEN_REVIEW'});counts['frozen_legacy_resolution']+=1;unresolved.append(row.copy())
   elif acc in proteins:
    row.update({'target_entity_type_v7':'canonical_protein','target_entity_id_v7':'CP-'+acc,'assignment_resolution_v7':'current_public_evidence_single_human_uniprot','assignment_disposition_v7':'PUBLIC_CANONICAL','assignment_limitations_v7':'no older identity-resolution row; no isoform specificity inferred'});counts['public_current_canonical_no_legacy_row']+=1
   else:
    row.update({'target_entity_type_v7':'unresolved_target','target_entity_id_v7':'UNRESOLVED:'+eid,'assignment_resolution_v7':'target_not_in_v7_protein_master','assignment_disposition_v7':'FROZEN_REVIEW'});counts['frozen_target_missing']+=1;unresolved.append(row.copy())
   w.writerow(row)

# Complex candidates remain a separate entity-specific candidate layer. Only old
# default-release complex evidence could enter PUBLIC; the current source has none.
complex_out=OUT/'complex_evidence_assignment_v7_final.tsv.gz'; complex_count=0
with op(OLD/'03_complex'/'complex_binding_evidence_candidates_v0_2.tsv') as fi:
 rd=csv.DictReader(fi,delimiter='\t'); fields=list(rd.fieldnames)+['v7_entity_assignment','v7_disposition','v7_reason']
 with op(complex_out,'wt') as fo:
  w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
  for r in rd:
   r['v7_entity_assignment']=r.get('complex_target_id','') or 'unresolved_target'
   safe=(r.get('default_release_inclusion')=='1' and r.get('complex_specificity_status')=='source_asserted_complex_context' and r.get('compound_mapping_status') in {'exact_source_identifier_unique','mapped_core'})
   r['v7_disposition']='PUBLIC_COMPLEX' if safe else 'FROZEN_REVIEW'
   r['v7_reason']='exact complex and compound relationship' if safe else (r.get('review_reason') or r.get('complex_relation_directness_status') or 'complex specificity/directness not uniquely verified')
   if safe:complex_count+=1
   w.writerow(r)

report={'current_evidence_rows':sum(counts.values()),'legacy_resolution_rows':len(resolution),'legacy_duplicate_evidence_ids':len(dup),'assignment_counts':dict(counts),'unresolved_target_rows':len(unresolved),'public_complex_candidate_rows':complex_count,'policy':'No gene-name inference; exact isoform FK only; canonical when explicit/current single-human-UniProt; unsafe complex and unresolved records frozen.'}
(OUT/'T17_T21_ENTITY_ASSIGNMENT_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
