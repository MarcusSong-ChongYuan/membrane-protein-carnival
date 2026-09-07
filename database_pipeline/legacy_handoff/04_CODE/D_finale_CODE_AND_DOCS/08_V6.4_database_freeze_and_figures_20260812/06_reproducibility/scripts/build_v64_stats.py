from pathlib import Path
import csv,gzip,json,collections,math,statistics
V=Path(r'D:\finale\01_正式数据_V6.2');R=Path(r'D:\finale\02_V6.3_candidate_20260803');O=Path(r'D:\7.22\v64_candidate_working');S=O/'stats';S.mkdir(exist_ok=True)
def truth(v):return str(v).lower() in {'1','true','yes','y'}
def write(name,header,rows):
 with open(S/name,'w',encoding='utf-8',newline='') as f:w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rows)
# protein
pc=collections.Counter();ec=collections.Counter();func=collections.Counter();role=collections.Counter();sf=collections.Counter();hpa=collections.Counter();struct=collections.Counter();tm=[];length=[]
with gzip.open(R/'05_protein_annotation'/'human_membrane_protein_master_v6_3_candidate.tsv.gz','rt',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  pc[r.get('membrane_class_v52') or 'unknown']+=1;ec[r.get('evidence_level_v52') or 'E0']+=1;func[r.get('molecular_function_primary_v63') or 'unclassified']+=1;role[r.get('membrane_role_primary_v63') or 'unclassified']+=1;sf[r.get('structural_family_primary_v63') or 'unclassified']+=1;hpa[r.get('hpa_mapping_status_v62') or 'missing']+=1
  st='Membrane PDB' if truth(r.get('opm_present_v5')) or truth(r.get('pdbtm_present_v5')) else 'Other PDB' if r.get('pdb_ids') else 'AlphaFold only' if r.get('alphafolddb_ids') else 'No structure';struct[st]+=1
  try:tm.append(float(r.get('transmembrane_count_v5') or 0));length.append(float(r.get('sequence_length_v53') or r.get('sequence_length') or 0))
  except:pass
write('protein_membrane_class_counts_v64.tsv',['membrane_class','protein_count'],pc.most_common());write('protein_evidence_level_counts_v64.tsv',['evidence_level','protein_count'],ec.most_common());write('protein_molecular_function_counts_v64.tsv',['molecular_function','protein_count'],func.most_common());write('protein_membrane_role_counts_v64.tsv',['membrane_role','protein_count'],role.most_common());write('protein_structural_family_counts_v64.tsv',['structural_family','protein_count'],sf.most_common());write('protein_hpa_mapping_counts_v64.tsv',['hpa_mapping_status','protein_count'],hpa.most_common());write('protein_structure_coverage_counts_v64.tsv',['structure_status','protein_count'],struct.most_common())
# compounds
status=collections.Counter();cls=collections.Counter();conf=collections.Counter();des={'molecular_weight':[],'xlogp':[],'tpsa':[],'rotatable_bond_count':[]};comp_rows=0
with open(V/'small_molecule_master_v1_3.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  comp_rows+=1;conf[r.get('identity_confidence') or 'missing']+=1;cls[r.get('computed_structural_class') or 'missing']+=1
  flags=[x for x in ['is_approved_drug','is_clinical_candidate','is_endogenous_ligand','is_natural_product','is_chemical_probe'] if truth(r.get(x))]
  status[';'.join(flags) if flags else 'none']+=1
  for k in des:
   try:
    v=float(r.get(k) or '');des[k].append(v) if math.isfinite(v) else None
   except:pass
write('compound_identity_confidence_counts_v64.tsv',['identity_confidence','compound_count'],conf.most_common());write('compound_structural_class_counts_v64.tsv',['structural_class','compound_count'],cls.most_common());write('compound_biostatus_intersections_v64.tsv',['biological_status_intersection','compound_count'],status.most_common())
rows=[]
for k,v in des.items():
 v.sort();q=lambda p:v[min(len(v)-1,int((len(v)-1)*p))] if v else None;rows.append([k,len(v),q(.1),q(.25),q(.5),q(.75),q(.9)])
write('compound_descriptor_summary_v64.tsv',['descriptor','n','p10','p25','median','p75','p90'],rows)
# evidence cleaned
et=collections.Counter();src=collections.Counter();incl=collections.Counter();pdbstat=collections.Counter();pairs=set();prots=set();comps=set()
with gzip.open(O/'candidate_tables'/'binding_evidence_master_v6_4_candidate.tsv.gz','rt',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  et[r.get('evidence_tier') or 'missing']+=1;src[r.get('source_database') or 'missing']+=1;incl[r.get('default_release_inclusion') or 'missing']+=1;pdbstat[r.get('pdb_identity_status_v64') or 'missing']+=1;prots.add(r['target_uniprot_id']);comps.add(r['compound_internal_id']);pairs.add((r['target_uniprot_id'],r['compound_internal_id']))
write('binding_evidence_tier_counts_v64.tsv',['evidence_tier','evidence_rows'],et.most_common());write('binding_source_counts_v64.tsv',['source_database','evidence_rows'],src.most_common());write('binding_pdb_identity_counts_v64.tsv',['pdb_identity_status','evidence_rows'],pdbstat.most_common())
# disease
levels=collections.Counter();areas=collections.Counter();systems=collections.Counter();drel=0
with open(R/'04_disease'/'protein_disease_relation_v2_1.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):drel+=1;levels[r.get('best_evidence_level') or 'missing']+=1
with open(R/'04_disease'/'disease_therapeutic_area_v0_1.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):areas[r.get('therapeutic_area_name') or r.get('therapeutic_area_id') or 'missing']+=1
with open(R/'04_disease'/'disease_anatomical_system_v0_1.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):systems[r.get('anatomical_system_name') or r.get('anatomical_system_id') or 'missing']+=1
write('disease_evidence_level_counts_v64.tsv',['evidence_level','protein_disease_pairs'],levels.most_common());write('disease_therapeutic_area_counts_v64.tsv',['therapeutic_area','relations'],areas.most_common());write('disease_anatomical_system_counts_v64.tsv',['anatomical_system','relations'],systems.most_common())
summary={'release':'MemPro V6.4 candidate','proteins':sum(pc.values()),'protein_isoforms':20211,'protein_complexes':19791,'canonical_compounds':comp_rows,'binding_evidence_rows':sum(et.values()),'unique_binding_pairs':len(pairs),'binding_proteins':len(prots),'binding_compounds':len(comps),'binding_site_instances':95598,'canonical_disease_relations':drel,'ihc_rows_v64':646118,'pdb_affected_evidence_rows':1274,'pdb_structural_default_excluded':2,'ihc_exact_duplicates_removed':35351,'ihc_content_distinct_rekeyed':9506}
(O/'stats'/'V64_RELEASE_STATISTICS.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))

