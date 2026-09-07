import csv,gzip,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path

BASE=Path(r'D:\finale\16_MemPro_V7.0.2_final_20260814')
ROOT=Path(r'D:\finale\17_MemPro_V7.1_candidate_20260814');OUT=ROOT/'01_release_tables';QA=ROOT/'04_QA'
def read(rel):
 p=BASE/rel;op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f,delimiter='\t')
def write(name,fields,rows):
 p=OUT/name;n=0
 with gzip.open(p,'wt',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
  for r in rows:w.writerow(r);n+=1
 return n
def h(*xs):return hashlib.sha256('|'.join(str(x or '') for x in xs).encode()).hexdigest()[:20]

def canonical_rows(masters):
 for r in masters:
  yield {'canonical_uniprot_accession':r['canonical_uniprot_accession'],'canonical_protein_entity_id':r['canonical_protein_entity_id'],'gene_entity_id':r['gene_entity_id'],'approved_symbol':r['approved_symbol'],'membrane_class':r['membrane_class_v7'],'primary_membrane_mode':r['primary_membrane_mode_v7'],'secondary_membrane_modes':r['secondary_membrane_modes_v7'],'evidence_level':r['membrane_evidence_level_v7'],'form_specificity':r['form_specificity_v7'],'public_inclusion':r['public_inclusion_v7'],'classification_scope':'CANONICAL_PROTEIN','source_database':'UniProtKB;MemPro-integrated','source_release':r['source_release'],'release_version':'V7.1-candidate'}
def iso_rows():
 seen=set()
 for rel,pred in [('01_data/protein_isoform_v701.tsv.gz',False),('03_prediction_layer/isoform_topology_prediction_v707.tsv.gz',True)]:
  for r in read(rel):
   k=r['protein_isoform_entity_id']
   if k in seen:continue
   seen.add(k);confirmed=r.get('isoform_membrane_class_v701','') if not pred else ''
   predicted=r.get('predicted_membrane_class_v707','') if pred else ''
   yield {'protein_isoform_entity_id':k,'isoform_uniprot_accession':r['isoform_uniprot_accession'],'canonical_protein_entity_id':r['canonical_protein_entity_id'],'canonical_uniprot_accession':r['canonical_uniprot_accession'],'isoform_number':r['isoform_number'],'sequence_length':r['sequence_length'],'sequence_sha256':r['sequence_sha256'],'confirmed_membrane_class':confirmed or 'UNRESOLVED','predicted_membrane_class':predicted,'primary_membrane_class':confirmed or 'UNRESOLVED','evidence_status':'SOURCE_CONFIRMED_OR_RULE_RESOLVED' if confirmed else 'PREDICTION_ONLY_NOT_CONFIRMED','classification_basis':r.get('isoform_membrane_basis_v701') or r.get('topology_prediction_basis_v707',''),'prediction_method':r.get('topology_prediction_method_v707',''),'prediction_release_scope':r.get('prediction_release_scope_v707',''),'source_database':r.get('source_database','UniProtKB'),'source_release':r.get('source_release',''),'release_version':'V7.1-candidate'}
def processed_rows():
 seen=set()
 for rel,pred in [('01_data/protein_processed_form_v701.tsv.gz',False),('03_prediction_layer/processed_form_topology_prediction_v707.tsv.gz',True)]:
  for r in read(rel):
   k=r['processed_form_id']
   if k in seen:continue
   seen.add(k);confirmed=r.get('membrane_class_assignment_v701','') if not pred else '';predicted=r.get('predicted_membrane_class_v707','') if pred else ''
   yield {'processed_form_id':k,'target_uniprot_id':r['target_uniprot_id'],'feature_id':r['feature_id'],'processed_form_type':r['processed_form_type'],'description':r['description'],'start':r['start'],'end':r['end'],'sequence_range_status':r['sequence_range_status'],'confirmed_membrane_class':confirmed or 'UNRESOLVED','predicted_membrane_class':predicted,'primary_membrane_class':confirmed or 'UNRESOLVED','membrane_retention_status':'CONFIRMED_OR_RULE_RESOLVED' if confirmed else 'UNRESOLVED_PREDICTION_AVAILABLE','classification_basis':r.get('assignment_basis_v701') or r.get('topology_prediction_basis_v707',''),'prediction_method':r.get('topology_prediction_method_v707',''),'source_database':r.get('source_database','UniProtKB'),'release_version':'V7.1-candidate'}

def main():
 OUT.mkdir(parents=True,exist_ok=True);QA.mkdir(parents=True,exist_ok=True)
 masters=list(read('01_data/protein_master_v7.tsv.gz')); iso=list(iso_rows());proc=list(processed_rows())
 cf=['canonical_uniprot_accession','canonical_protein_entity_id','gene_entity_id','approved_symbol','membrane_class','primary_membrane_mode','secondary_membrane_modes','evidence_level','form_specificity','public_inclusion','classification_scope','source_database','source_release','release_version']
 counts={};counts['canonical']=write('canonical_protein_membrane_class_v71.tsv.gz',cf,canonical_rows(masters))
 isof=['protein_isoform_entity_id','isoform_uniprot_accession','canonical_protein_entity_id','canonical_uniprot_accession','isoform_number','sequence_length','sequence_sha256','confirmed_membrane_class','predicted_membrane_class','primary_membrane_class','evidence_status','classification_basis','prediction_method','prediction_release_scope','source_database','source_release','release_version']
 counts['isoform']=write('isoform_membrane_class_v71.tsv.gz',isof,iso)
 pf=['processed_form_id','target_uniprot_id','feature_id','processed_form_type','description','start','end','sequence_range_status','confirmed_membrane_class','predicted_membrane_class','primary_membrane_class','membrane_retention_status','classification_basis','prediction_method','source_database','release_version']
 counts['processed']=write('processed_form_membrane_class_v71.tsv.gz',pf,proc)
 iso_by=defaultdict(list);proc_by=Counter()
 for r in iso:iso_by[r['canonical_uniprot_accession']].append(r)
 for r in proc:proc_by[r['target_uniprot_id']]+=1
 def genes():
  for r in masters:
   xs=iso_by[r['canonical_uniprot_accession']]
   yield {'gene_entity_id':r['gene_entity_id'],'approved_symbol':r['approved_symbol'],'canonical_uniprot_accession':r['canonical_uniprot_accession'],'canonical_membrane_class':r['membrane_class_v7'],'canonical_evidence_level':r['membrane_evidence_level_v7'],'isoform_count':len(xs),'isoform_confirmed_class_count':sum(x['confirmed_membrane_class']!='UNRESOLVED' for x in xs),'isoform_prediction_only_count':sum(x['evidence_status'].startswith('PREDICTION') for x in xs),'processed_form_count':proc_by[r['canonical_uniprot_accession']],'gene_summary_semantics':'summary_only_each_protein_form_retains_own_class','release_version':'V7.1-candidate'}
 counts['gene_summary']=write('gene_membrane_summary_v71.tsv.gz',['gene_entity_id','approved_symbol','canonical_uniprot_accession','canonical_membrane_class','canonical_evidence_level','isoform_count','isoform_confirmed_class_count','isoform_prediction_only_count','processed_form_count','gene_summary_semantics','release_version'],genes())
 def states():
  for r in read('01_data/protein_membrane_state_v7.tsv.gz'):
   yield {**r,'state_condition':'DEFAULT_OR_SOURCE_REPORTED_STATE','state_specificity':'NOT_EXPLICITLY_RESOLVED' if r.get('form_specificity','').lower() in {'','canonical'} else r.get('form_specificity'),'release_version':'V7.1-candidate'}
 sf=['membrane_state_id','target_uniprot_id','entity_scope','membrane_class','primary_membrane_mode','secondary_membrane_modes','evidence_level','form_specificity','public_inclusion','decision','evidence_summary','state_condition','state_specificity','release_version']
 counts['states']=write('state_dependent_membrane_association_v71.tsv.gz',sf,states())

 complexes=list(read('01_data/protein_complex_v7.tsv.gz'));comps=list(read('01_data/complex_component_v7.tsv.gz'))
 cxf=list(complexes[0].keys())+['release_version'];counts['complexes']=write('protein_complex_v71.tsv.gz',cxf,({**r,'release_version':'V7.1-candidate'} for r in complexes))
 ccf=list(comps[0].keys())+['release_version'];counts['components']=write('complex_component_v71.tsv.gz',ccf,({**r,'release_version':'V7.1-candidate'} for r in comps))
 def stoich():
  for r in comps:
   yield {'stoichiometry_assertion_id':'STO-'+h(r['complex_target_id'],r['component_assertion_id']),'complex_target_id':r['complex_target_id'],'component_assertion_id':r['component_assertion_id'],'component_entity_type':r['component_entity_type'],'component_identifier':r['component_identifier'],'copy_count':r['copy_count'],'stoichiometry_status':r['stoichiometry_status'],'source_database':r['source_database'],'source_complex_id':r['source_complex_id'],'release_version':'V7.1-candidate'}
 counts['stoichiometry']=write('complex_stoichiometry_v71.tsv.gz',['stoichiometry_assertion_id','complex_target_id','component_assertion_id','component_entity_type','component_identifier','copy_count','stoichiometry_status','source_database','source_complex_id','release_version'],stoich())
 def assemblies():
  for r in complexes:
   yield {'complex_assembly_id':'ASM-'+h(r['complex_target_id']),'complex_target_id':r['complex_target_id'],'biological_assembly_id':'','assembly_status':'NOT_REPORTED_BY_CURRENT_COMPLEX_SOURCE','component_set_status':'CONFLICT_REPORTED' if r['component_set_conflict_flag']=='1' else 'SOURCE_COMPONENT_SET','contributing_databases':r['contributing_databases'],'source_complex_ids':r['source_complex_ids'],'release_version':'V7.1-candidate'}
 counts['assemblies']=write('complex_assembly_v71.tsv.gz',['complex_assembly_id','complex_target_id','biological_assembly_id','assembly_status','component_set_status','contributing_databases','source_complex_ids','release_version'],assemblies())
 def cxstates():
  for r in complexes:
   yield {'complex_state_id':'CXS-'+h(r['complex_target_id']),'complex_target_id':r['complex_target_id'],'state_name':'source_default','state_condition':'NOT_REPORTED','state_dependent_assembly_status':'UNRESOLVED_SOURCE_NOT_REPORTED','release_version':'V7.1-candidate'}
 counts['complex_states']=write('complex_state_v71.tsv.gz',['complex_state_id','complex_target_id','state_name','state_condition','state_dependent_assembly_status','release_version'],cxstates())
 be=list(read('02_nonpublic_resolved_status/complex_evidence_assignment_v7.tsv.gz'));bef=list(be[0].keys())+['binding_target_scope_v71','binding_component_id_v71','release_layer_v71','release_version']
 def bx():
  for r in be:
   scope='COMPLEX' if r.get('complex_target_id') else 'UNRESOLVED_TARGET';disp=r.get('v7_disposition','')
   yield {**r,'binding_target_scope_v71':scope,'binding_component_id_v71':'','release_layer_v71':'PUBLIC_CONTEXT' if disp.startswith('PUBLIC') else 'NONPUBLIC_SOURCE_LIMITED','release_version':'V7.1-candidate'}
 counts['complex_binding']=write('complex_binding_evidence_v71.tsv.gz',bef,bx())

 af=['evidence_id','source_database','source_record_id','source_target_identifier','canonical_uniprot_accession','target_entity_type_v71','target_entity_id_v71','protein_form_specificity_v71','protein_isoform_entity_id','processed_form_id','complex_target_id','assignment_resolution_v71','assignment_limitations_v71','release_version']
 def assignments():
  for r in read('01_data/evidence_target_entity_assignment_v7.tsv.gz'):
   # Canonical accession is a stable identity anchor, but the historical assay usually did not specify an isoform.
   yield {'evidence_id':r['evidence_id'],'source_database':r['source_database'],'source_record_id':r['source_record_id'],'source_target_identifier':r['source_target_identifier'],'canonical_uniprot_accession':r['canonical_uniprot_accession'],'target_entity_type_v71':'canonical_protein','target_entity_id_v71':r['canonical_uniprot_accession'],'protein_form_specificity_v71':'ISOFORM_UNSPECIFIED','protein_isoform_entity_id':'','processed_form_id':'','complex_target_id':r.get('complex_target_id',''),'assignment_resolution_v71':'CANONICAL_IDENTITY_RESOLVED_FORM_UNSPECIFIED','assignment_limitations_v71':'original_normalized_record_does_not_recover_isoform_specificity','release_version':'V7.1-candidate'}
 counts['assignments']=write('evidence_target_entity_assignment_v71.tsv.gz',af,assignments())
 expected={'canonical':7800,'isoform':9240,'processed':8454,'complexes':4584,'components':14825,'assignments':942455}
 failures={k:{'observed':counts.get(k),'expected':v} for k,v in expected.items() if counts.get(k)!=v}
 report={'counts':counts,'expected_core_counts':expected,'failures':failures,'pass':not failures}
 (QA/'V71_STAGE3_ENTITY_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
