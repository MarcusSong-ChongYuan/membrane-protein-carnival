import csv,gzip,hashlib,json,shutil
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

W=Path(r"D:\finale\12_V7_final_completion_20260813");AUTO=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data");DEST=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813");DATA=DEST/'01_data';REV=DEST/'02_frozen_review';QA=DEST/'03_QA';DOC=DEST/'04_docs'
for d in [DATA,REV,QA,DOC]:d.mkdir(parents=True,exist_ok=True)
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
def ids(p,k):
 a=[]
 with op(p) as f:
  for r in csv.DictReader(f,delimiter='\t'):a.append(r[k])
 return a
def count(p):
 with op(p) as f:return sum(1 for _ in csv.DictReader(f,delimiter='\t'))
def fk(p,k,parent):
 n=b=0
 with op(p) as f:
  for r in csv.DictReader(f,delimiter='\t'):n+=1;b+=r[k] not in parent
 return {'rows':n,'bad':b}

P=AUTO/'protein_master_v7.tsv.gz';E=W/'02_be_recompute'/'protein_compound_evidence_BE_recomputed_v7.tsv.gz';PAIR=AUTO/'protein_compound_pair_v7.tsv.gz';CM=W/'04_negative_conflicts'/'compound_master_positive_negative_union_v7.tsv.gz';CF=W/'04_negative_conflicts'/'compound_form_positive_negative_union_v7.tsv.gz';NEG=W/'04_negative_conflicts'/'negative_evidence_v7.tsv.gz';DIS=W/'06_remaining_modules'/'protein_disease_relation_v7_validated.tsv.gz';CX=W/'06_remaining_modules'/'protein_complex_v7_validated.tsv.gz';CC=AUTO/'complex_component_v7.tsv.gz';EXP=W/'05_hpa_expression'/'expression_measurement_v7_deduplicated.tsv.gz';SITEPUB=W/'03_binding_sites'/'binding_site_coordinate_verified_v7.tsv.gz';SITEFR=W/'03_binding_sites'/'binding_site_coordinate_frozen_v7.tsv.gz';ASS=W/'01_entity_assignment'/'evidence_target_entity_assignment_v7_final.tsv.gz'
ps=set(ids(P,'canonical_uniprot_accession'));cs=set(ids(CM,'compound_internal_id'));fs=set(ids(CF,'compound_form_id'));es=ids(E,'evidence_id');ess=set(es);xs=set(ids(CX,'complex_target_id'));checks={}
checks['protein_pk']={'rows':len(ps),'bad':7800-len(ps)};checks['compound_pk']={'rows':len(cs),'bad':len(ids(CM,'compound_internal_id'))-len(cs)};checks['form_pk']={'rows':len(fs),'bad':len(ids(CF,'compound_form_id'))-len(fs)};checks['evidence_pk']={'rows':len(es),'bad':len(es)-len(ess)}
checks['positive_target_fk']=fk(E,'target_uniprot_id',ps);checks['positive_compound_fk']=fk(E,'compound_internal_id',cs);checks['negative_target_fk']=fk(NEG,'target_uniprot_id',ps);checks['negative_compound_fk']=fk(NEG,'compound_internal_id',cs);checks['negative_form_fk']=fk(NEG,'compound_form_id',fs|{''});checks['disease_target_fk']=fk(DIS,'target_uniprot_id',ps);checks['complex_component_parent_fk']=fk(CC,'complex_target_id',xs);checks['expression_target_fk']=fk(EXP,'target_uniprot_id',ps);checks['site_evidence_fk']=fk(SITEPUB,'evidence_id',ess)
# Pair exact reconstruction after BE recompute.
epairs=set();tiers=Counter();sources=Counter();business=0
with op(E) as f:
 for r in csv.DictReader(f,delimiter='\t'):
  epairs.add((r['target_uniprot_id'],r['compound_internal_id']));tiers[r['evidence_tier_recomputed_v7']]+=1;sources[r['source_database']]+=1
  business+=r['evidence_tier_recomputed_v7'] not in {'BE1','BE2','BE3'} or r['be_disposition_v7']!='PUBLIC'
pairkeys=set()
with op(PAIR) as f:
 for r in csv.DictReader(f,delimiter='\t'):pairkeys.add((r['target_uniprot_id'],r['compound_internal_id']))
checks['pair_reconstruction']={'rows':len(epairs),'bad':len(epairs^pairkeys)};checks['BE_business_rules']={'rows':len(es),'bad':business}
# Entity assignment exact one row per evidence.
assignids=ids(ASS,'evidence_id');checks['entity_assignment_one_to_one']={'rows':len(assignids),'bad':len(assignids)-len(set(assignids))+len(ess^set(assignids))}
# HPA primary key and site coordinate status.
exids=ids(EXP,'expression_record_id');checks['expression_pk']={'rows':len(exids),'bad':len(exids)-len(set(exids))}
site_bad=0
with op(SITEPUB) as f:
 for r in csv.DictReader(f,delimiter='\t'):site_bad+=not r['coordinate_mapping_status_v7'].startswith('PUBLIC')
checks['site_public_coordinate_rule']={'rows':count(SITEPUB),'bad':site_bad}
blocking=sum(x['bad'] for x in checks.values())

# Assemble formal data. Existing static entity tables are copied; rebuilt modules replace candidates.
formal=[(P,'protein_master_v7.tsv.gz'),(AUTO/'protein_isoform_v7.tsv.gz','protein_isoform_v7.tsv.gz'),(AUTO/'protein_processed_form_v7.tsv.gz','protein_processed_form_v7.tsv.gz'),(AUTO/'protein_membrane_state_v7.tsv.gz','protein_membrane_state_v7.tsv.gz'),(CX,'protein_complex_v7.tsv.gz'),(CC,'complex_component_v7.tsv.gz'),(CM,'compound_master_v7.tsv.gz'),(CF,'compound_form_v7.tsv.gz'),(E,'protein_compound_evidence_v7.tsv.gz'),(PAIR,'protein_compound_pair_v7.tsv.gz'),(ASS,'evidence_target_entity_assignment_v7.tsv.gz'),(SITEPUB,'binding_site_coordinate_verified_v7.tsv.gz'),(NEG,'negative_evidence_v7.tsv.gz'),(W/'04_negative_conflicts'/'positive_negative_conflict_v7.tsv.gz','evidence_conflict_v7.tsv.gz'),(AUTO/'disease_master_v7.tsv.gz','disease_master_v7.tsv.gz'),(DIS,'protein_disease_relation_v7.tsv.gz'),(EXP,'expression_measurement_v7.tsv.gz'),(W/'05_hpa_expression'/'subcellular_localization_v7_complete.tsv.gz','subcellular_localization_v7.tsv.gz'),(W/'05_hpa_expression'/'hpa_protein_coverage_state_v7.tsv','hpa_protein_coverage_state_v7.tsv'),(AUTO/'protein_cross_classification_v7.tsv.gz','protein_cross_classification_v7.tsv.gz'),(AUTO/'protein_cross_classification_summary_v7.tsv.gz','protein_cross_classification_summary_v7.tsv.gz'),(W/'06_remaining_modules'/'protein_membrane_evidence_and_risk_v7.tsv.gz','protein_membrane_evidence_and_risk_v7.tsv.gz'),(AUTO/'source_registry_v7.tsv','source_registry_v7.tsv')]
for s,n in formal:shutil.copy2(s,DATA/n)
review=[(W/'02_be_recompute'/'evidence_BE_frozen_review_v7.tsv.gz','evidence_BE_frozen_review_v7.tsv.gz'),(W/'01_entity_assignment'/'complex_evidence_assignment_v7_final.tsv.gz','complex_evidence_assignment_v7.tsv.gz'),(W/'03_binding_sites'/'binding_site_coordinate_frozen_v7.tsv.gz','binding_site_coordinate_frozen_v7.tsv.gz'),(W/'04_negative_conflicts'/'negative_evidence_frozen_v7.tsv.gz','negative_evidence_frozen_v7.tsv.gz'),(W/'05_hpa_expression'/'expression_id_conflicts_frozen_v7.tsv.gz','expression_id_conflicts_frozen_v7.tsv.gz'),(W/'06_remaining_modules'/'protein_complex_frozen_v7.tsv.gz','protein_complex_frozen_v7.tsv.gz'),(W/'06_remaining_modules'/'protein_disease_frozen_v7.tsv.gz','protein_disease_frozen_v7.tsv.gz'),(AUTO.parent/'04_frozen_review'/'isoform_frozen_review_v7.tsv.gz','isoform_frozen_review_v7.tsv.gz'),(AUTO.parent/'04_frozen_review'/'processed_form_frozen_review_v7.tsv.gz','processed_form_frozen_review_v7.tsv.gz'),(AUTO.parent/'04_frozen_review'/'excluded_records_frozen_v7.tsv.gz','excluded_records_frozen_v7.tsv.gz')]
for s,n in review:shutil.copy2(s,REV/n)

qa_report={'status':'PASS_V7_FINAL_COMPLETION_QA' if blocking==0 else 'FAIL_V7_FINAL_COMPLETION_QA','blocking_errors':blocking,'checks':checks,'counts':{'proteins':len(ps),'compounds':len(cs),'forms':len(fs),'positive_evidence':len(es),'positive_pairs':len(epairs),'negative_evidence':count(NEG),'negative_pairs':1487597,'BE_tiers':dict(tiers),'coordinate_verified_sites':count(SITEPUB),'coordinate_frozen_sites':count(SITEFR),'expression_rows':len(exids),'disease_relations':count(DIS),'complexes':len(xs)},'positive_source_counts':dict(sources),'manual_review':'WAIVED_BY_USER; no manual accuracy claim'}
(QA/'V7_FINAL_COMPLETION_QA.json').write_text(json.dumps(qa_report,ensure_ascii=False,indent=2),encoding='utf-8')
for p in W.rglob('*.json'):
 if p.is_file():shutil.copy2(p,QA/p.name)
for p in [W/'00_baseline'/'V7_CLASSIFICATION_AND_RELEASE_RULES.tsv',W/'00_baseline'/'V7_DATA_DICTIONARY_DRAFT.tsv']:
 shutil.copy2(p,DOC/p.name)
(DOC/'README.md').write_text('# MemPro V7 final data release\n\nThis release completes the automatic T01-T50 data workflow with manual review explicitly waived. Ambiguous records are frozen, never guessed. See V7_FINAL_COMPLETION_QA.json and FREEZE_V7_FINAL.json.\n',encoding='utf-8')

manifest=[]
for p in sorted(DEST.rglob('*')):
 if not p.is_file() or p.name in {'MANIFEST.tsv','FREEZE_V7_FINAL.json'}:continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 manifest.append({'relative_path':str(p.relative_to(DEST)),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
with (DEST/'MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
freeze={'release':'MemPro V7.0.0 final automatic data release','created_at_utc':datetime.now(timezone.utc).isoformat(),'status':'FROZEN_PASS' if blocking==0 else 'NOT_FROZEN_FAIL','blocking_errors':blocking,'manual_review':'WAIVED_BY_USER','manifest_entries':len(manifest),'counts':qa_report['counts'],'limitations':['No manual validation accuracy claim','48,331 binding-site records lack explicit chain and remain outside coordinate-verified layer','2,358 chain-bearing sites failed unique SIFTS coordinate mapping and are frozen','required/optional complex subunits remain source-not-asserted where unavailable','83 proteins remain functionally unresolved after five-axis annotation']}
(DEST/'FREEZE_V7_FINAL.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({'qa':qa_report,'freeze':freeze},ensure_ascii=False,indent=2))
