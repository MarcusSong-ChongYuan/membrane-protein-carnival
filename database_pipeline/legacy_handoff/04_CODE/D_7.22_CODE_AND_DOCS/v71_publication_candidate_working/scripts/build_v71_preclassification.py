import csv,gzip,json,hashlib,collections,re
from pathlib import Path
ROOT=Path(r'D:\7.22\v71_publication_candidate_working');IN=ROOT/'outputs'/'binding_evidence_v71_locally_enriched.tsv.gz';OUT=ROOT/'outputs';Q=ROOT/'qa'
P0={'KD','KI','KA','KB','KE'}; P1={'IC50','EC50','AC50','ED50','GI50','DC50','IC90','POTENCY','INHIBITION','ACTIVITY'}
restricted={'PDBbind','DrugCentral','PDSP KiDatabase','BioLiP','sc-PDB'}
extra=['evidence_class_v71','prepublication_decision_v71','prepublication_reason_v71','independent_publication_key_v71','independent_experiment_key_v71','independent_structure_key_v71','redistribution_status_v71']
counts=collections.Counter();pairs=collections.defaultdict(set);keys=collections.defaultdict(set)
with gzip.open(IN,'rt',encoding='utf8',newline='') as fi,gzip.open(OUT/'binding_evidence_v71_preclassified.tsv.gz','wt',encoding='utf8',newline='',compresslevel=6) as fo:
 rd=csv.DictReader(fi,delimiter='\t');wr=csv.DictWriter(fo,fieldnames=rd.fieldnames+extra,delimiter='\t',lineterminator='\n');wr.writeheader()
 for r in rd:
  act=r.get('activity_type','').upper();et=r.get('evidence_type','');src=r.get('source_database_v71','');pdb=r.get('pdb_ids','').lower();pm=r.get('pubmed_id_v71','');doi=r.get('doi_v71','').lower()
  if et in {'experimental_structure_residue_contact','compound_specific_structure','compound_specific_structure_context','experimental_complex_with_affinity','protein_level_curated_site_annotation'}:cl='P2_structural_or_site'
  elif act in P0 and et in {'direct_quantitative_binding','pubchem_direct_quantitative_binding','binding_assay_activity'}:cl='P0_direct_binding'
  else:cl='P1_functional_pharmacology'
  r['evidence_class_v71']=cl;r['independent_publication_key_v71']='PMID:'+pm if pm else ('DOI:'+doi if doi else ('PDB:'+pdb if pdb else ''))
  r['independent_experiment_key_v71']=hashlib.sha256('|'.join([r['independent_publication_key_v71'],r.get('assay_id_v71',''),r.get('target_uniprot_id',''),r.get('compound_internal_id',''),act,r.get('standard_value_nM',''),r.get('assay_or_mechanism','')]).encode()).hexdigest()[:24]
  r['independent_structure_key_v71']=hashlib.sha256('|'.join([pdb,r.get('target_uniprot_id',''),r.get('ligand_het_id',''),r.get('binding_site_residues','')]).encode()).hexdigest()[:24] if cl.startswith('P2') and pdb else ''
  r['redistribution_status_v71']='conditional_transform_or_review' if any(x in src for x in restricted) else 'allowed_with_attribution_or_provenance'
  reasons=[]
  if r.get('protein_evidence_level_v71') not in {'E1','E2'}:reasons.append('protein_not_E1_E2')
  if not r.get('membrane_class_v71'):reasons.append('membrane_class_missing')
  if r.get('compound_mapping_status')!='mapped_core':reasons.append('compound_not_mapped_core')
  if not r.get('source_record_id_v71'):reasons.append('source_record_missing')
  if r.get('publication_status_v71')!='citation_or_structure_anchor':reasons.append('citation_not_recovered')
  if r.get('structure_qc_v71')=='artifact_excluded':reasons.append('structure_artifact')
  if r['redistribution_status_v71'].startswith('conditional'):reasons.append('license_or_public_field_review')
  if cl=='P1_functional_pharmacology' and 'primary' in r.get('local_enrichment_status_v71','').lower():reasons.append('primary_screen')
  r['prepublication_decision_v71']='PUBLIC_CANDIDATE' if not reasons else 'FREEZE_PENDING';r['prepublication_reason_v71']=';'.join(reasons)
  counts[r['prepublication_decision_v71']]+=1;counts['class|'+cl]+=1
  k=(r.get('target_uniprot_id',''),r.get('compound_internal_id',''));pairs[r['prepublication_decision_v71']].add(k);keys['experiment'].add(r['independent_experiment_key_v71']);
  if r['independent_structure_key_v71']:keys['structure'].add(r['independent_structure_key_v71'])
  wr.writerow(r)
report={'row_counts':dict(counts),'pair_counts':{k:len(v) for k,v in pairs.items()},'distinct_experiment_keys':len(keys['experiment']),'distinct_structure_keys':len(keys['structure'])};(Q/'V71_PRECLASSIFICATION_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2))
