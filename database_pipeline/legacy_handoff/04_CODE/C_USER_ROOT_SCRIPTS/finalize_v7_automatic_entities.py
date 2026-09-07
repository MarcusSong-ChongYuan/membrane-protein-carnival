import csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813");REL=ROOT/'03_v7_relation_rebuild';MOD=ROOT/'04_v7_remaining_modules';OUT=ROOT/'05_v7_automatic_release_candidate';OUT.mkdir(parents=True,exist_ok=True)
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')

# Rebuild canonical master from final adjudication rather than retaining stale V6.4 classes.
src=ROOT/'02_protein_audit'/'V7_PROTEIN_FREEZE_CANDIDATE'/'human_membrane_protein_master_V7.tsv.gz'
canon=[]
with op(src) as f:
 for r in csv.DictReader(f,delimiter='\t'):
  canon.append({'canonical_protein_entity_id':'CP-'+r['target_uniprot'],'canonical_uniprot_accession':r['target_uniprot'],'approved_symbol':r['approved_symbol'],'protein_name':r['protein_name'],'reviewed_status':r['reviewed'],'organism_id':r['organism_id'],'sequence_length':r['uniprot_sequence_length_current'],'membrane_class_v7':r['final_class_v7_3'],'primary_membrane_mode_v7':r['final_mode_v7_3'],'secondary_membrane_modes_v7':r.get('proposed_secondary_modes_v7',''),'membrane_evidence_level_v7':r.get('cross_source_proposed_evidence_v7') or r.get('proposed_evidence_level_v7'),'form_specificity_v7':r.get('form_specificity_v7_2',''),'decision_v7':r['final_decision_v7_3'],'decision_evidence_v7':r['final_evidence_v7_3'],'public_inclusion_v7':'1','source_release':'MemPro V7 automatic candidate'})
fields=list(canon[0]);
with op(OUT/'protein_master_V7.tsv.gz','wt') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(canon)
proteins={r['canonical_uniprot_accession'] for r in canon}

# Evidence entity mapping: all current public rows are already single-protein assignments.
emap=[];source_versions=defaultdict(set);eids=set();pairs=defaultdict(list)
with op(REL/'public_binding_evidence_master_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  eids.add(r['evidence_id']);pairs[(r['target_uniprot_id'],r['compound_internal_id'])].append(r)
  source_versions[r['source_database']].add(r['source_version'])
  explicit_isoform='-' in r['target_uniprot_id']
  emap.append({'evidence_id':r['evidence_id'],'target_entity_type':'protein_isoform' if explicit_isoform else 'canonical_protein','target_entity_id':r['target_uniprot_id'],'canonical_uniprot_accession':r['target_uniprot_id'].split('-')[0],'isoform_assignment':'explicit_source_isoform' if explicit_isoform else 'gene_product_unspecified_canonical_target','relationship_context_id':r['relationship_context_id'],'target_assignment_status':r['target_assignment_status'],'v7_disposition':'PUBLIC'})
with op(OUT/'evidence_target_entity_mapping_V7.tsv.gz','wt') as f:
 fields=list(emap[0]);w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(emap)

# Enriched pair table including lineage counts and conflict state.
lineage={}
with op(MOD/'protein_compound_evidence_lineage_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):lineage[(r['target_uniprot_id'],r['compound_internal_id'])]=r
conflicts=set()
with op(MOD/'evidence_conflict_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):conflicts.add((r['target_uniprot_id'],r['compound_internal_id']))
pairrows=[]
for key,rows in pairs.items():
 lin=lineage.get(key,{})
 pairrows.append({'pair_id':'PAIR-'+hashlib.sha1((key[0]+'|'+key[1]).encode()).hexdigest()[:20].upper(),'target_uniprot_id':key[0],'compound_internal_id':key[1],'evidence_count':len(rows),'best_evidence_tier':min((r['evidence_tier'] for r in rows),key=lambda x:{'BE1':1,'BE2':2,'BE3':3}.get(x,9)),'distinct_database_count':len({r['source_database'] for r in rows}),'source_databases':';'.join(sorted({r['source_database'] for r in rows})),'independent_experiment_count':lin.get('independent_literature_experiment_count_v63',''),'independent_structure_count':lin.get('independent_structure_count_v63',''),'independent_pubchem_assay_count':lin.get('independent_pubchem_assay_count_v63',''),'independent_evidence_modality_count':lin.get('independent_evidence_modality_count_v63',''),'positive_negative_conflict_status':'CONTEXT_DEPENDENT_CONFLICT' if key in conflicts else 'NO_PUBLISHED_NEGATIVE_CONFLICT','release_version':'MemPro V7 automatic candidate'})
with op(OUT/'protein_compound_pair_V7.tsv.gz','wt') as f:
 fields=list(pairrows[0]);w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(pairrows)

# Expression mapping status explicitly distinguishes mapped data from missing.
exprcount=Counter(); layers=defaultdict(set)
with op(REL/'protein_tissue_cell_ihc_expression_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):exprcount[r['target_uniprot_id']]+=1;layers[r['target_uniprot_id']].add(r['measurement_layer'])
with op(OUT/'expression_localization_coverage_V7.tsv.gz','wt') as f:
 fields=['target_uniprot_id','expression_record_count','measurement_layers','expression_mapping_status'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for p in sorted(proteins):w.writerow({'target_uniprot_id':p,'expression_record_count':exprcount[p],'measurement_layers':';'.join(sorted(layers[p])),'expression_mapping_status':'mapped_records_available' if exprcount[p] else 'missing_or_unmapped_no_record'})

# Site identity disposition using V6.4 structural identity state; unknown is explicit.
site_status=Counter()
with op(REL/'binding_site_instances_V7.tsv.gz') as fi,op(OUT/'binding_site_identity_qc_V7.tsv.gz','wt') as fo:
 rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['site_identity_disposition_v7','site_identity_reason_v7'];w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
 for r in rd:
  st=r.get('pdb_identity_status_v64','')
  if st in {'validated_current','current','exact_current','not_applicable_no_pdb'} or not r.get('pdb_ids',''):disp='PUBLIC'
  else:disp='FROZEN_REVIEW'
  r['site_identity_disposition_v7']=disp;r['site_identity_reason_v7']=st or ('no_pdb_identifier' if not r.get('pdb_ids','') else 'missing_structural_identity_status');site_status[disp]+=1;w.writerow(r)

# Source registry derived from actual public evidence plus frozen annotation registries.
with op(OUT/'source_registry_V7.tsv','wt') as f:
 fields=['module','source_database','source_versions','record_scope','registry_basis'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for s,v in sorted(source_versions.items()):w.writerow({'module':'binding_evidence','source_database':s,'source_versions':';'.join(sorted(v)),'record_scope':'public V7 evidence','registry_basis':'derived from retained evidence records'})
 with op(MOD/'source_registry_V7.tsv') as fi:
  for r in csv.DictReader(fi,delimiter='\t'):w.writerow({'module':r['module'],'source_database':r['source'],'source_versions':r['version'],'record_scope':r['role'],'registry_basis':'frozen module registry; sha256='+r['sha256']})

report={'protein_master':len(canon),'evidence_target_mappings':len(emap),'pair_rows':len(pairrows),'expression_proteins_mapped':sum(bool(exprcount[p]) for p in proteins),'expression_proteins_missing_or_unmapped':sum(not exprcount[p] for p in proteins),'binding_site_identity_disposition':dict(site_status),'source_databases':len(source_versions)}
(OUT/'V7_AUTOMATIC_ENTITY_FINALIZATION_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
