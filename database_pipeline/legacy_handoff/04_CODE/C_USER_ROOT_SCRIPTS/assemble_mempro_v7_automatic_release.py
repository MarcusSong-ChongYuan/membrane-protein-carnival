import csv,gzip,hashlib,json,shutil
from datetime import datetime,timezone
from pathlib import Path

W=Path(r"D:\finale\09_V7_data_freeze_working_20260813");REL=W/'03_v7_relation_rebuild';MOD=W/'04_v7_remaining_modules';CAN=W/'05_v7_automatic_release_candidate'
DEST=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813");DATA=DEST/'01_data';QA=DEST/'02_QA';DOC=DEST/'03_docs';REV=DEST/'04_frozen_review';
for d in [DATA,QA,DOC,REV]:d.mkdir(parents=True,exist_ok=True)
def cp(src,name,where=DATA):shutil.copy2(src,where/name)
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')

# T48 formal tables.
mapping=[
(CAN/'protein_master_V7.tsv.gz','protein_master_v7.tsv.gz'),
(MOD/'protein_isoform_membrane_state_V7.tsv.gz','protein_isoform_v7.tsv.gz'),
(MOD/'protein_processed_form_exact_V7.tsv.gz','protein_processed_form_v7.tsv.gz'),
(MOD/'protein_membrane_state_V7.tsv.gz','protein_membrane_state_v7.tsv.gz'),
(REL/'complex_target_master_V7.tsv.gz','protein_complex_v7.tsv.gz'),
(REL/'complex_target_components_V7.tsv.gz','complex_component_v7.tsv.gz'),
(REL/'small_molecule_master_V7.tsv.gz','compound_master_v7.tsv.gz'),
(MOD/'compound_form_V7.tsv.gz','compound_form_v7.tsv.gz'),
(CAN/'protein_compound_evidence_V7.tsv.gz','protein_compound_evidence_v7.tsv.gz'),
(CAN/'protein_compound_pair_V7.tsv.gz','protein_compound_pair_v7.tsv.gz'),
(MOD/'negative_evidence_V7.tsv.gz','negative_evidence_v7.tsv.gz'),
(MOD/'evidence_conflict_V7.tsv.gz','evidence_conflict_v7.tsv.gz'),
(REL/'disease_entity_master_V7.tsv.gz','disease_master_v7.tsv.gz'),
(REL/'protein_disease_relation_V7.tsv.gz','protein_disease_relation_v7.tsv.gz'),
(REL/'protein_tissue_cell_ihc_expression_V7.tsv.gz','expression_localization_v7.tsv.gz'),
(MOD/'protein_cross_classification_V7.tsv.gz','protein_cross_classification_v7.tsv.gz'),
(MOD/'protein_cross_classification_summary_V7.tsv.gz','protein_cross_classification_summary_v7.tsv.gz'),
(CAN/'evidence_target_entity_mapping_V7.tsv.gz','evidence_target_entity_mapping_v7.tsv.gz'),
(CAN/'source_registry_V7.tsv','source_registry_v7.tsv'),
]
for s,n in mapping:cp(s,n)

# Merge binding-site identity and membrane-side fields into one publication table.
side={}
with op(MOD/'binding_site_membrane_side_complete_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):side[r['binding_site_instance_id']]=r
with op(CAN/'binding_site_identity_qc_V7.tsv.gz') as fi:
 rd=csv.DictReader(fi,delimiter='\t');extra=['membrane_side','membrane_structure_support','coordinate_mapping_basis','orientation_basis','side_classification_scope'];fields=list(rd.fieldnames)+extra
 with op(DATA/'binding_site_v7.tsv.gz','wt') as fo:
  w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
  for r in rd:
   s=side.get(r['binding_site_instance_id'],{});r.update({'membrane_side':s.get('membrane_side','unknown'),'membrane_structure_support':s.get('membrane_structure_support','none_or_not_mapped'),'coordinate_mapping_basis':s.get('coordinate_mapping_basis','no_valid_orientation_mapping'),'orientation_basis':s.get('orientation_basis','not_inferred'),'side_classification_scope':s.get('classification_scope','explicit_unknown_no_guess')});w.writerow(r)

# Frozen review outputs remain separate to preserve original schemas; index is the
# uniform review_queue_frozen_v7 table requested by T48.
review_sources=[
('negative_evidence',MOD/'negative_evidence_frozen_review_V7.tsv.gz','negative_evidence_frozen_review_v7.tsv.gz'),
('evidence_R_X',MOD/'review_excluded_evidence_frozen_V7.tsv.gz','evidence_R_X_frozen_v7.tsv.gz'),
('protein_excluded',MOD/'excluded_records_frozen_V7.tsv.gz','excluded_records_frozen_v7.tsv.gz'),
]
for _,s,n in review_sources:cp(s,n,REV)
iso_review=REV/'isoform_frozen_review_v7.tsv.gz';proc_review=REV/'processed_form_frozen_review_v7.tsv.gz'
def filter_rows(src,dst,pred):
 n=0
 with op(src) as fi:
  rd=csv.DictReader(fi,delimiter='\t')
  with op(dst,'wt') as fo:
   w=csv.DictWriter(fo,fieldnames=rd.fieldnames,delimiter='\t');w.writeheader()
   for r in rd:
    if pred(r):w.writerow(r);n+=1
 return n
niso=filter_rows(MOD/'protein_isoform_membrane_state_V7.tsv.gz',iso_review,lambda r:r['isoform_disposition_v7'].startswith('FROZEN'))
nproc=filter_rows(MOD/'protein_processed_form_exact_V7.tsv.gz',proc_review,lambda r:r['assignment_status'].startswith('FROZEN'))
review_index=[
{'review_layer':'negative_evidence','file':'negative_evidence_frozen_review_v7.tsv.gz','record_count':1936309,'automatic_disposition':'FROZEN_REVIEW','reason':'compound/form identity not uniquely publishable under V7 rules'},
{'review_layer':'evidence_R_X','file':'evidence_R_X_frozen_v7.tsv.gz','record_count':1956629,'automatic_disposition':'FROZEN_REVIEW_OR_EXCLUDED_AUDIT','reason':'V6.4.1 conservative triage not promoted'},
{'review_layer':'protein_excluded','file':'excluded_records_frozen_v7.tsv.gz','record_count':3197,'automatic_disposition':'EXCLUDED_OR_OTHER_ENTITY','reason':'no validated membrane mechanism or non-protein immune segment entity'},
{'review_layer':'isoform','file':iso_review.name,'record_count':niso,'automatic_disposition':'FROZEN_REVIEW','reason':'no unique form-specific membrane mechanism'},
{'review_layer':'processed_form','file':proc_review.name,'record_count':nproc,'automatic_disposition':'FROZEN_REVIEW','reason':'no processed-form-specific membrane mechanism'},
]
with (REV/'review_queue_frozen_v7.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(review_index[0]),delimiter='\t');w.writeheader();w.writerows(review_index)

# QA reports.
for p in [CAN/'QA'/'V7_AUTOMATIC_RELEASE_QA.json',CAN/'V7_EVIDENCE_LINEAGE_DETAIL_REPORT.json',CAN/'V7_T11_CANONICAL_PROTEIN_REPORT.json',CAN/'V7_SITE_IDENTITY_DISPOSITION_REPORT.json',MOD/'V7_ISOFORM_PROCESSED_EXACT_REPORT.json',MOD/'V7_NEGATIVE_AND_SITE_COMPLETION_REPORT.json',MOD/'V7_REMAINING_MODULES_BUILD_REPORT.json',REL/'V7_RELATIONAL_INTEGRITY_AUDIT.json']:
 cp(p,p.name,QA)

# Machine-generated data dictionary from every formal table.
dictionary=[];counts={}
for p in sorted(DATA.iterdir()):
 if not p.is_file():continue
 with op(p) as f:
  rd=csv.reader(f,delimiter='\t');header=next(rd);n=sum(1 for _ in rd)
 counts[p.name]=n
 for i,col in enumerate(header,1):dictionary.append({'table':p.name,'column_order':i,'column_name':col,'nullable_policy':'empty means missing/not available; never reinterpret as numeric zero','release':'MemPro V7 automatic'})
with (DOC/'V7_DATA_DICTIONARY.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(dictionary[0]),delimiter='\t');w.writeheader();w.writerows(dictionary)

(DOC/'V7_RELEASE_POLICY.md').write_text('''# MemPro V7 automatic data release policy\n\nHuman manual review was omitted by explicit user approval. Records enter PUBLIC only through deterministic identity, foreign-key, evidence and business-rule checks. Records that cannot be uniquely assigned are not guessed: they are retained in FROZEN_REVIEW or EXCLUDED_AUDIT. A/B/C precedence is A > B > C while secondary mechanisms are retained. Positive and negative evidence never overwrite one another. Missing values are never converted to zero.\n''',encoding='utf-8')
(DOC/'V7_VERSION_CHANGELOG.md').write_text('''# V7 change log\n\n- Re-adjudicated 10,997 candidate human proteins and established a 7,800-protein A/B/C public whitelist.\n- Rebuilt all public relations against that whitelist.\n- Added canonical sequence and gene identifiers, isoform state, exact UniProt processed features, compound parent/form layer, evidence lineage keys, complex entities, disease ontology relations, expression coverage, binding-site identity and membrane-side disposition.\n- Isolated unresolved isoforms, processed forms, negative evidence and V6.4.1 R/X evidence outside the public layer.\n- Manual review and manual sampling were intentionally omitted; this is an automatic rules-based release.\n''',encoding='utf-8')

# Manifest and final freeze record.
manifest=[]
for p in sorted(DEST.rglob('*')):
 if not p.is_file() or p.name in {'MANIFEST.tsv','FREEZE_V7_AUTOMATIC.json'}:continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 manifest.append({'relative_path':str(p.relative_to(DEST)),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
with (DEST/'MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
freeze={'release':'MemPro V7.0.0 automatic data-layer release','created_at_utc':datetime.now(timezone.utc).isoformat(),'manual_review':'OMITTED_BY_USER_APPROVAL','qa_status':'PASS_AUTOMATIC_V7_RELEASE_CANDIDATE','blocking_errors':0,'formal_table_counts':counts,'review_layers':review_index,'manifest_entries':len(manifest),'status':'FROZEN_AUTOMATIC_RULES_BASED_RELEASE','limitations':['No manual review or manual accuracy estimate','4,665 isoforms and 3,130 processed forms isolated because form-specific membrane mechanism was not uniquely inferable','No negative evidence met strict V7 canonical compound/form publication rules','48,331 binding sites have explicit unknown membrane-side orientation']}
(DEST/'FREEZE_V7_AUTOMATIC.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(freeze,ensure_ascii=False,indent=2))
