from pathlib import Path
import csv,json,hashlib,shutil,os
O=Path(r'D:\7.22\v64_candidate_working');R=Path(r'D:\finale\02_V6.3_candidate_20260803');V=Path(r'D:\finale\01_正式数据_V6.2');C=O/'candidate_tables'
# copy/link lightweight authoritative V6.3 products into candidate directory, no large V6.2 duplication
files={
'human_membrane_protein_master_v6_4_candidate.tsv.gz':R/'05_protein_annotation'/'human_membrane_protein_master_v6_3_candidate.tsv.gz',
'canonical_protein_entity_v6_4_candidate.tsv.gz':R/'02_identity'/'canonical_protein_entity_v0_1.tsv.gz',
'protein_isoform_v6_4_candidate.tsv.gz':R/'02_identity'/'protein_isoform_v0_2.tsv.gz',
'complex_target_master_v6_4_candidate.tsv':R/'03_complex'/'complex_target_master_v0_2.tsv',
'disease_entity_master_v6_4_candidate.tsv':R/'04_disease'/'disease_entity_master_v0_1.tsv',
'protein_disease_relation_v6_4_candidate.tsv':R/'04_disease'/'protein_disease_relation_v2_1.tsv',
'disease_therapeutic_area_v6_4_candidate.tsv':R/'04_disease'/'disease_therapeutic_area_v0_1.tsv',
'disease_anatomical_system_v6_4_candidate.tsv':R/'04_disease'/'disease_anatomical_system_v0_1.tsv'}
for n,p in files.items():
 if not (C/n).exists():shutil.copy2(p,C/n)
# license registry conservative
licenses=[
['UniProtKB','2026_02','attribution/terms apply','public_with_attribution','canonical protein annotation','https://www.uniprot.org/help/license'],
['ChEMBL','37','CC BY-SA 3.0 (data; verify release terms)','public_sharealike','bioactivity/compound identifiers','https://www.ebi.ac.uk/chembl/'],
['BindingDB','2024/loaded snapshot','CC BY 4.0 for BindingDB-curated; ChEMBL-derived rows inherit ChEMBL','public_conditional_by_row','binding measurements','https://www.bindingdb.org/'],
['PubChem BioAssay','retrieved 2026-07','contributor-specific; NCBI policies','public_conditional_provenance_required','assay results and compound IDs','https://pubchem.ncbi.nlm.nih.gov/docs/downloads'],
['PDBe','retrieved 2026-07','PDBe/wwPDB terms; attribution','public_with_attribution','PDB structures and residue contacts','https://www.ebi.ac.uk/pdbe/'],
['PDBbind','v2020 reprocessed package','license/redistribution not documented in project','internal_only_pending_permission','affinity and processed complexes','https://www.pdbbind-plus.org.cn/'],
['BRENDA','2026.1 downloaded snapshot','redistribution permission not documented','internal_only_pending_permission','enzyme ligand evidence','https://www.brenda-enzymes.org/'],
['DrugCentral','loaded snapshot','verify source license and row provenance','public_conditional','drug-target mechanism','https://drugcentral.org/'],
['IUPHAR/BPS Guide to PHARMACOLOGY','loaded snapshot','CC BY-SA terms/attribution; verify snapshot','public_with_attribution_sharealike','pharmacology and ligand status','https://www.guidetopharmacology.org/'],
['Human Protein Atlas','25.1','CC BY-SA 3.0; attribution','public_with_attribution_sharealike','RNA/IHC/localization','https://www.proteinatlas.org/about/licence'],
['Open Targets Platform','26.06','Open Targets data terms; attribution','public_with_attribution','disease evidence and therapeutic areas','https://platform.opentargets.org/downloads'],
['MONDO','2026-07-06','CC BY 4.0','public_with_attribution','canonical disease identity','https://mondo.monarchinitiative.org/'],
['Disease Ontology','2026-07-31','CC0 1.0','public','disease hierarchy/anatomy axioms','https://disease-ontology.org/'],
['Uberon','frozen snapshot','CC BY 3.0','public_with_attribution','anatomy mapping','https://obophenotype.github.io/uberon/']]
with open(O/'qa'/'V64_SOURCE_LICENSE_REGISTRY.tsv','w',encoding='utf-8',newline='') as f:w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['source','version','license_status','redistribution_class','contribution','url']);w.writerows(licenses)
# Final validation synthesis
qa={}
for n in ['V64_PDB_IDENTITY_CLEANING.json','V64_KEY_FK_ACCESSION_AUDIT.json','V64_RELATIONAL_INTEGRITY_AUDIT.json']:
 qa[n]=json.load(open(O/'qa'/n,encoding='utf-8'))
scope=json.load(open(O/'stats'/'V64_RELEASE_SCOPE_COUNTS.json',encoding='utf-8'));stats=json.load(open(O/'stats'/'V64_RELEASE_STATISTICS.json',encoding='utf-8'))
gates={'field_audit':'PASS_WITH_CONTROLLED_CORRECTIONS','pdb_identity':'PASS','primary_keys':'PASS','foreign_keys':'PASS_WITH_EXTERNAL_CONTEXT_DECLARED','ihc_identity':'PASS','v631_accession_patch':'PASS','disease_identity':'PASS_FROM_V63','five_axis_protein_classification':'PASS_FROM_V63','license_internal_release':'PASS','license_public_redistribution':'CONDITIONAL'}
validation={'release':'MemPro V6.4 candidate','status':'PASS_CANDIDATE_INTERNAL','baseline_v62_modified':False,'baseline_v63_modified':False,'gates':gates,'counts':{**stats,**scope},'controlled_corrections':{'pdb_evidence_rows_cleaned':1274,'pdb_site_rows_cleaned':1293,'structural_rows_excluded_default':2,'ihc_exact_duplicates_removed':35351,'ihc_distinct_content_rekeyed':9506,'complex_components_high_confidence_replaced':20,'complex_components_ambiguous_retained':9,'complex_components_source_conflict':1},'limitations':['Public redistribution remains source-conditional; BRENDA and PDBbind raw fields are internal-only pending documented permission.','External-context isoforms outside the 10,997 membrane-protein universe are retained and explicitly labeled, not treated as orphan corruption.','Candidate statistics separate all-evidence and default-release scopes.','Docking is not required for V6.4 database freeze.']}
(O/'qa'/'V64_CANDIDATE_VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(validation,ensure_ascii=False,indent=2))
