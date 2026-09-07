import csv,gzip,json,hashlib,collections,re
from pathlib import Path
import pandas as pd
ROOT=Path(r'D:\7.22\v71_publication_candidate_working');IN=ROOT/'outputs'/'binding_evidence_v71_preclassified.tsv.gz';OUT=ROOT/'outputs';Q=ROOT/'qa';W=ROOT/'work'
cit=pd.read_csv(W/'chembl_document_citations_v71.tsv',sep='\t',dtype=str).drop_duplicates('document_chembl_id').set_index('document_chembl_id').to_dict('index')
public=OUT/'binding_evidence_public_candidate_v7_1.tsv.gz';frozen=OUT/'binding_evidence_frozen_nonpublic_v7_1.tsv.gz';cross=OUT/'evidence_source_crosswalk_v7_1.tsv.gz'
add=['final_publication_status_v71','final_publication_reason_v71','publication_id_v71']
counts=collections.Counter();pairs=collections.defaultdict(set);exp=collections.defaultdict(set);pubs=set();srcs=collections.Counter()
with gzip.open(IN,'rt',encoding='utf8',newline='') as fi,gzip.open(public,'wt',encoding='utf8',newline='',compresslevel=6) as fp,gzip.open(frozen,'wt',encoding='utf8',newline='',compresslevel=6) as ff,gzip.open(cross,'wt',encoding='utf8',newline='',compresslevel=6) as fx:
 rd=csv.DictReader(fi,delimiter='\t');fields=rd.fieldnames+add;wp=csv.DictWriter(fp,fieldnames=fields,delimiter='\t',lineterminator='\n');wf=csv.DictWriter(ff,fieldnames=fields,delimiter='\t',lineterminator='\n');wp.writeheader();wf.writeheader();cw=csv.DictWriter(fx,fieldnames=['evidence_id','source_database','source_version','source_record_id','source_url','document_id','publication_id','lineage_role'],delimiter='\t',lineterminator='\n');cw.writeheader()
 for r in rd:
  # ChEMBL document citation recovery
  doc=r.get('document_id_v71','');m=cit.get(doc,{}) if doc else {}
  if not r.get('pubmed_id_v71') and m.get('pubmed_id') and str(m.get('pubmed_id'))!='nan':r['pubmed_id_v71']=str(m['pubmed_id'])
  if not r.get('doi_v71') and m.get('doi') and str(m.get('doi'))!='nan':r['doi_v71']=str(m['doi'])
  if r.get('pubmed_id_v71'):pid='PMID:'+str(r['pubmed_id_v71'])
  elif r.get('doi_v71'):pid='DOI:'+str(r['doi_v71']).lower()
  elif r.get('pdb_ids'):pid='PDB:'+r['pdb_ids'].lower()
  else:pid=''
  r['publication_id_v71']=pid
  reasons=[]
  if r.get('protein_evidence_level_v71') not in {'E1','E2'}:reasons.append('protein_not_E1_E2')
  if not r.get('membrane_class_v71'):reasons.append('membrane_class_missing')
  if r.get('compound_mapping_status')!='mapped_core':reasons.append('compound_not_mapped_core')
  if not r.get('source_record_id_v71'):reasons.append('source_record_missing')
  if r.get('source_database')=='v3.1_mixed_sources' and str(r.get('source_record_id_v71','')).startswith('ASSAY_'):reasons.append('legacy_internal_assay_not_native_source_id')
  if not pid:reasons.append('publication_not_verified')
  if r.get('structure_qc_v71')=='artifact_excluded':reasons.append('structure_artifact')
  if r.get('evidence_class_v71')=='P1_functional_pharmacology' and 'pubchem_assay_type=primary' in r.get('local_enrichment_status_v71','').lower():reasons.append('primary_screen_only')
  if r.get('redistribution_status_v71')=='conditional_transform_or_review':reasons.append('license_or_public_field_review')
  # Database-only curated assertions cannot substitute for an experiment or structure.
  if r.get('evidence_type')=='curated_pharmacological_target_assertion':reasons.append('curated_assertion_not_primary_experiment')
  decision='PUBLIC_CANDIDATE' if not reasons else 'FROZEN_NONPUBLIC'
  r['final_publication_status_v71']=decision;r['final_publication_reason_v71']=';'.join(reasons)
  (wp if decision=='PUBLIC_CANDIDATE' else wf).writerow(r);counts[decision]+=1;counts['class|'+r.get('evidence_class_v71','')]+=1
  k=(r.get('target_uniprot_id',''),r.get('compound_internal_id',''));pairs[decision].add(k);exp[decision].add(r.get('independent_experiment_key_v71',''))
  if decision=='PUBLIC_CANDIDATE' and pid:pubs.add(pid)
  srcs[r.get('source_database_v71','')]+=1
  cw.writerow({'evidence_id':r.get('evidence_id',''),'source_database':r.get('source_database_v71',''),'source_version':r.get('source_version_v71',''),'source_record_id':r.get('source_record_id_v71',''),'source_url':r.get('source_url_v71',''),'document_id':doc,'publication_id':pid,'lineage_role':'primary_or_contributing_source'})
report={'rows':dict(counts),'unique_pairs':{k:len(v) for k,v in pairs.items()},'independent_experiment_keys':{k:len(v-{""}) for k,v in exp.items()},'unique_publications_public_candidate':len(pubs),'source_counts':dict(srcs),'status':'CANDIDATE_NOT_FROZEN_LICENSE_AND_EXTERNAL_REPRO_PENDING'};(Q/'V71_FINAL_CANDIDATE_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(report,ensure_ascii=False,indent=2))
