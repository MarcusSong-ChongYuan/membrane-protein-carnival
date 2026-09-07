import csv, gzip, hashlib, json, re
from collections import Counter
from pathlib import Path

WORK=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
REL=WORK/"03_v7_relation_rebuild"
PROT=WORK/"02_protein_audit"/"V7_PROTEIN_FREEZE_CANDIDATE"/"protein_identity_adjudication_all_V7.tsv.gz"
OUT=WORK/"04_v7_remaining_modules"
OUT.mkdir(parents=True,exist_ok=True)
V63=Path(r"D:\finale\02_V6.3_candidate_20260803")

def op(p,mode="rt"):
    return gzip.open(p,mode,encoding="utf-8",newline="") if str(p).endswith('.gz') else p.open(mode,encoding="utf-8",newline="")
def keys(p,col):
    with op(p) as f:return {r[col] for r in csv.DictReader(f,delimiter='\t')}
def filter_table(src,dst,pred,transform=None):
    ni=no=0
    with op(src) as fi:
        r=csv.DictReader(fi,delimiter='\t'); fields=list(r.fieldnames)
        extra=[] if transform is None else transform(None)
        fields+=extra
        with op(dst,'wt') as fo:
            w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
            for row in r:
                ni+=1
                if pred(row):
                    if transform: row.update(transform(row))
                    w.writerow(row);no+=1
    return {'input':ni,'output':no}

with op(PROT) as f: allprot=list(csv.DictReader(f,delimiter='\t'))
public={r['target_uniprot'] for r in allprot if r['public_inclusion_v7_3']=='1'}
allids={r['target_uniprot'] for r in allprot}
compounds=keys(REL/'small_molecule_master_V7.tsv.gz','compound_internal_id')
evidence=keys(REL/'public_binding_evidence_master_V7.tsv.gz','evidence_id')
public_pairs=keys(REL/'public_protein_compound_pairs_V7.tsv.gz','target_uniprot_id') # only quick protein presence
pairset=set()
with op(REL/'public_protein_compound_pairs_V7.tsv.gz') as f:
    for r in csv.DictReader(f,delimiter='\t'):pairset.add((r['target_uniprot_id'],r['compound_internal_id']))

report={'public_proteins':len(public),'modules':{}}

# T11/T13: state and processed products. Exact residue ranges are retained from
# current UniProt feature strings; no invented sequence or isoform assignment.
state_fields=['membrane_state_id','target_uniprot_id','entity_scope','membrane_class','primary_membrane_mode','secondary_membrane_modes','evidence_level','form_specificity','public_inclusion','decision','evidence_summary']
with op(OUT/'protein_membrane_state_V7.tsv.gz','wt') as fo:
    w=csv.DictWriter(fo,fieldnames=state_fields,delimiter='\t');w.writeheader()
    for r in allprot:
        if r['target_uniprot'] not in public:continue
        w.writerow({'membrane_state_id':'MSTATE-'+r['target_uniprot'],'target_uniprot_id':r['target_uniprot'],'entity_scope':'canonical_or_explicit_form','membrane_class':r['final_class_v7_3'],'primary_membrane_mode':r['final_mode_v7_3'],'secondary_membrane_modes':r.get('proposed_secondary_modes_v7',''),'evidence_level':r.get('cross_source_proposed_evidence_v7') or r.get('proposed_evidence_level_v7'),'form_specificity':r.get('form_specificity_v7_2',''),'public_inclusion':'1','decision':r['final_decision_v7_3'],'evidence_summary':r['final_evidence_v7_3']})

processed=[]
for r in allprot:
    if r['target_uniprot'] not in public:continue
    for ftype,col in [('signal_peptide','uniprot_signal_features_current'),('chain','uniprot_chain_count_current'),('peptide','uniprot_peptide_count_current')]:
        val=r.get(col,'')
        if not val or val=='0':continue
        processed.append({'processed_form_id':'PFORM-'+hashlib.sha1((r['target_uniprot']+'|'+ftype+'|'+val).encode()).hexdigest()[:16].upper(),'target_uniprot_id':r['target_uniprot'],'source_feature_type':ftype,'source_feature_value':val,'sequence_range_status':'explicit_range_if_present_else_count_only','membrane_class_assignment':r['final_class_v7_3'] if 'EXACT' in val else 'not_independently_assigned','assignment_status':'source_preserved_no_inference','source_database':'UniProtKB current audit'})
pf=['processed_form_id','target_uniprot_id','source_feature_type','source_feature_value','sequence_range_status','membrane_class_assignment','assignment_status','source_database']
with op(OUT/'protein_processed_form_V7.tsv.gz','wt') as fo:
    w=csv.DictWriter(fo,fieldnames=pf,delimiter='\t');w.writeheader();w.writerows(processed)
report['modules']['processed_forms']={'output':len(processed)}

# T18/T19: parent/form layer restricted to public compounds.
formsrc=V63/'01_baseline_v6_2'/'compound_form_hierarchy_v1_3.tsv'
report['modules']['compound_forms']=filter_table(formsrc,OUT/'compound_form_V7.tsv.gz',lambda r:r['compound_internal_id'] in compounds)

# T23: already-computed experiment/structure/modality lineage restricted to V7 pairs.
lineagesrc=V63/'06_evidence_lineage'/'protein_compound_lineage_summary_v0_1.tsv.gz'
report['modules']['evidence_lineage']=filter_table(lineagesrc,OUT/'protein_compound_evidence_lineage_V7.tsv.gz',lambda r:(r['target_uniprot_id'],r['compound_internal_id']) in pairset)

# T29-T31: mapped negatives only; ambiguity stays in frozen review. Conflict is
# pair-level and never deletes either side.
negsrc=Path(r"D:\finale\negative_upgrade\upgraded_negative_evidence.tsv.gz")
neg_public=OUT/'negative_evidence_V7.tsv.gz'; neg_review=OUT/'negative_evidence_frozen_review_V7.tsv.gz'
with op(negsrc) as fi:
    rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['v7_disposition','v7_disposition_reason']
    with op(neg_public,'wt') as fp,op(neg_review,'wt') as fr:
        wp=csv.DictWriter(fp,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wp.writeheader();wr.writeheader()
        np=nr=0; negpairs=set(); conflicts=[]
        for r in rd:
            if r['target_uniprot_id'] not in allids:continue
            mapped=(r['target_uniprot_id'] in public and r['compound_internal_id'] in compounds and r.get('release_action_v62')=='mapped_negative_release' and r.get('identity_resolution_v62') not in {'unassigned_stereochemistry_review',''})
            if mapped:
                r['v7_disposition']='PUBLIC_NEGATIVE';r['v7_disposition_reason']='protein_and_canonical_compound_FK_valid;V6.2_mapped_negative_release';wp.writerow(r);np+=1;negpairs.add((r['target_uniprot_id'],r['compound_internal_id']))
            else:
                r['v7_disposition']='FROZEN_REVIEW';r['v7_disposition_reason']='identity_scope_or_release_rule_not_uniquely_satisfied';wr.writerow(r);nr+=1
        for p in sorted(negpairs & pairset):conflicts.append(p)
report['modules']['negative_evidence']={'public':np,'frozen_review':nr,'conflict_pairs':len(conflicts)}
with op(OUT/'evidence_conflict_V7.tsv.gz','wt') as fo:
    fields=['conflict_id','target_uniprot_id','compound_internal_id','conflict_type','automatic_resolution','public_positive_retained','public_negative_retained']
    w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
    for p,c in enumerate(conflicts,1):w.writerow({'conflict_id':f'CONFLICT-V7-{p:08d}','target_uniprot_id':c[0],'compound_internal_id':c[1],'conflict_type':'positive_and_negative_evidence_across_assay_contexts','automatic_resolution':'CONTEXT_DEPENDENT_NO_SIDE_DELETED','public_positive_retained':'1','public_negative_retained':'1'})

# T28: reuse structure-orientation classification only when the site instance is
# still present in V7; all unclassified sites receive explicit unknown rows later.
sites=keys(REL/'binding_site_instances_V7.tsv.gz','binding_site_instance_id')
sidesrc=Path(r"D:\finale\MemPro_M1_M8_AGENT_REPRO_20260806\06_binding_site_side\binding_site_membrane_side_instances.tsv.gz")
report['modules']['binding_site_side']=filter_table(sidesrc,OUT/'binding_site_membrane_side_V7.tsv.gz',lambda r:r['binding_site_instance_id'] in sites)

# T39/T40 five-axis classification and raw source annotations.
ann=V63/'05_protein_annotation'
report['modules']['cross_classification']=filter_table(ann/'protein_cross_classification_v0_1.tsv.gz',OUT/'protein_cross_classification_V7.tsv.gz',lambda r:r['target_uniprot_id'] in public)
report['modules']['cross_classification_summary']=filter_table(ann/'protein_cross_classification_summary_v0_1.tsv',OUT/'protein_cross_classification_summary_V7.tsv.gz',lambda r:r['target_uniprot_id'] in public)
report['modules']['function_source_annotation']=filter_table(ann/'protein_function_source_annotation_v0_1.tsv.gz',OUT/'protein_function_source_annotation_V7.tsv.gz',lambda r:r['target_uniprot_id'] in public)

# T32 frozen review/excluded records: preserve V6.4.1 R/X and final protein audit.
rx=Path(r"D:\finale\09_V6.4.1_final_publication_triage_20260812\frozen_R_X_archive_v6_4_1.tsv.gz")
report['modules']['evidence_review_excluded']=filter_table(rx,OUT/'review_excluded_evidence_frozen_V7.tsv.gz',lambda r:r.get('target_uniprot_id','') in allids)
report['modules']['protein_excluded']=filter_table(PROT,OUT/'excluded_records_frozen_V7.tsv.gz',lambda r:r['public_inclusion_v7_3']!='1')

# T49 source registry: preserve frozen source/version/hash registries.
registries=[ann/'protein_annotation_source_registry_v0_1.tsv',V63/'04_disease'/'disease_source_registry_v0_1.tsv']
regrows=[]
for p in registries:
    with op(p) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            regrows.append({'module':'protein_annotation' if 'protein_annotation' in str(p) else 'disease','source':r.get('source',''),'version':r.get('version',''),'local_file':r.get('local_file',''),'sha256':r.get('sha256',''),'role':r.get('role','')})
with op(OUT/'source_registry_V7.tsv','wt') as fo:
    f=['module','source','version','local_file','sha256','role'];w=csv.DictWriter(fo,fieldnames=f,delimiter='\t');w.writeheader();w.writerows(regrows)
report['modules']['source_registry']={'output':len(regrows)}

(OUT/'V7_REMAINING_MODULES_BUILD_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
