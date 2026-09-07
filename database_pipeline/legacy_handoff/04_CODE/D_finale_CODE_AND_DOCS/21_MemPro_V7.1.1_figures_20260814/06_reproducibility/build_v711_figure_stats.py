import csv,gzip,json,math,hashlib
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(r'D:\finale\20_MemPro_V7.1.1_final_20260814')
OUT=Path(r'D:\finale\21_MemPro_V7.1.1_figures_20260814')
ST=OUT/'02_statistics'; QA=OUT/'05_QA'
def rows(name):
 p=ROOT/'01_release_tables'/name;op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f,delimiter='\t')
def split(s):return [x.strip() for x in (s or '').replace('|',';').split(';') if x.strip()]
def num(x):
 try:return float(x)
 except:return None
def save_counter(name,c,a='category',b='count'):
 with (ST/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow([a,b]);w.writerows(c.most_common())
def main():
 ST.mkdir(parents=True,exist_ok=True);QA.mkdir(parents=True,exist_ok=True);S={}
 # M1/M2 protein landscape and source support.
 pc=Counter();pe=Counter();pm=Counter();protein_sources=Counter();classes=Counter();axes=defaultdict(set);axis_labels=defaultdict(Counter)
 for r in rows('protein_master_v71.tsv.gz'):
  pc[r['membrane_class_v7']]+=1;pe[r['membrane_evidence_level_v7']]+=1;pm[r['primary_membrane_mode_v7']]+=1
 for r in rows('protein_membrane_evidence_and_risk_v71.tsv.gz'):
  for x in split(r.get('evidence_source_ids_v7','')):protein_sources[x]+=1
 for r in rows('protein_cross_classification_v71.tsv.gz'):
  ax=r['classification_axis'];axes[ax].add(r['target_uniprot_id']);axis_labels[ax][r['classification_label']]+=1
 S['protein']={'class':dict(pc),'evidence':dict(pe),'mode':dict(pm),'source_support':dict(protein_sources),'axis_coverage':{k:len(v) for k,v in axes.items()},'axis_top':{k:v.most_common(12) for k,v in axis_labels.items()}}
 # Positive evidence and source landscape.
 es=Counter();et=Counter();em=Counter();pdbstat=Counter();evidence_to_pair={};pair_source=Counter()
 for r in rows('positive_interaction_evidence_v71.tsv.gz'):
  es[r['source_database']]+=1;et[r.get('evidence_tier_recomputed_v7') or r.get('evidence_tier','')]+=1;em[r.get('evidence_modality_v7','')]+=1;pdbstat[r.get('pdb_identity_status_v64','')]+=1
  evidence_to_pair[r['evidence_id']]=(r['target_uniprot_id'],r['compound_internal_id'])
 S['evidence']={'source':dict(es),'tier':dict(et),'modality':dict(em),'pdb_status':dict(pdbstat)}
 # M3 expression breadth and status.
 exlayer=Counter();exdet=defaultdict(Counter);exmap=Counter();breadth=defaultdict(lambda:defaultdict(set));measured_proteins=defaultdict(set)
 for r in rows('expression_measurement_v71.tsv.gz'):
  layer=r['expression_layer_v711'];det=r['detection_status_v711'];mp=r['mapping_status_v71'];exlayer[layer]+=1;exdet[layer][det]+=1;exmap[mp]+=1
  if r['mapped_measured_denominator_eligible_v71']=='1':measured_proteins[layer].add(r['target_uniprot_id'])
  if det=='DETECTED':
   loc=r.get('cell_type') if layer=='cell_type_RNA' else r.get('tissue')
   if loc:breadth[layer][r['target_uniprot_id']].add(loc)
 bdist={}
 for layer in ['tissue_RNA','cell_type_RNA','normal_tissue_IHC','tissue_protein_MS']:
  ids=measured_proteins[layer]|set(breadth[layer]);bdist[layer]=[len(breadth[layer].get(x,set())) for x in ids]
 S['expression']={'layer':dict(exlayer),'detection':{k:dict(v) for k,v in exdet.items()},'mapping':dict(exmap),'breadth':bdist,'measured_protein_count':{k:len(v) for k,v in measured_proteins.items()}}
 # M4 disease ontology / evidence / sources.
 disease_systems=defaultdict(set);system_rows=Counter()
 for r in rows('disease_anatomical_system_v71.tsv.gz'):
  disease_systems[r['canonical_disease_id']].add(r['anatomical_system_name']);system_rows[r['anatomical_system_name']]+=1
 de=Counter();ds=Counter();system_evidence=defaultdict(Counter);disease_pairs=0
 for r in rows('protein_disease_relation_v71.tsv.gz'):
  disease_pairs+=1;de[r['best_evidence_level']]+=1
  for x in split(r['source_databases']):ds[x]+=1
  for sys in disease_systems.get(r['canonical_disease_id'],{'Unmapped/no anatomy'}):system_evidence[sys][r['best_evidence_level']]+=1
 S['disease']={'evidence':dict(de),'sources':dict(ds),'system_rows':dict(system_rows),'system_evidence':{k:dict(v) for k,v in system_evidence.items()},'relations':disease_pairs,'canonical_diseases':len(disease_systems)}
 # M5 compound status UpSet, structure classes and descriptors.
 bio=Counter();struct=Counter();identity=Counter();desc=defaultdict(list);compounds=0
 flags=['is_approved_drug','is_clinical_candidate','is_endogenous_ligand','is_natural_product','is_chemical_probe']
 for r in rows('compound_master_v71.tsv.gz'):
  compounds+=1;key=tuple(f for f in flags if str(r.get(f,'')).lower() in {'1','true','yes'});bio[key]+=1;struct[r.get('computed_structural_class') or 'Unclassified']+=1;identity[r.get('identity_confidence') or 'Unspecified']+=1
  # deterministic descriptor sample avoids loading all values
  if int(hashlib.sha256(r['compound_internal_id'].encode()).hexdigest()[:8],16)%97==0:
   for k in ['molecular_weight','xlogp','tpsa','rotatable_bond_count']:
    v=num(r.get(k));
    if v is not None and math.isfinite(v):desc[k].append(v)
 S['compound']={'count':compounds,'bio_intersections':{'&'.join(k) if k else 'None':v for k,v in bio.items()},'structural_class':dict(struct),'identity':dict(identity),'descriptor_sample':desc}
 # M6 site status, chain and coordinate coverage; pair-level funnel.
 sr=Counter();sc=Counter();schain=Counter();slig=Counter();sside=Counter();site_pairs=set();coord_pairs=set();complete_pairs=set();docking_pairs=set()
 for r in rows('interaction_site_v71.tsv.gz'):
  sr[r['site_residue_status']]+=1;sc[r['coordinate_mapping_status']]+=1;schain[r['chain_assignment_status']]+=1;slig[r['structure_ligand_relationship']]+=1;sside[r.get('membrane_side') or 'UNKNOWN']+=1
  pair=evidence_to_pair.get(r['evidence_id'])
  if pair:
   site_pairs.add(pair)
   if r['coordinate_mapping_status'] in {'COORDINATE_COMPLETE','COORDINATE_PARTIAL'}:coord_pairs.add(pair)
   if r['coordinate_mapping_status']=='COORDINATE_COMPLETE':complete_pairs.add(pair)
   if r['coordinate_docking_eligible']=='1':docking_pairs.add(pair)
 S['site']={'residue':dict(sr),'coordinate':dict(sc),'chain':dict(schain),'ligand':dict(slig),'membrane_side':dict(sside),'pair_funnel':{'all_positive_pairs':529168,'pairs_with_site_assertion':len(site_pairs),'pairs_with_complete_or_partial_coordinates':len(coord_pairs),'pairs_with_complete_coordinates':len(complete_pairs),'pairs_coordinate_docking_eligible':len(docking_pairs)}}
 # M7 lineage and release QA.
 lin=Counter();exp_mirror=Counter();str_mirror=Counter()
 for r in rows('evidence_lineage_v71.tsv.gz'):
  exp_mirror[r['cross_database_experiment_mirror_status']]+=1;str_mirror[r['cross_database_structure_mirror_status']]+=1;lin[r['evidence_modality']]+=1
 S['lineage']={'experiment_status':dict(exp_mirror),'structure_status':dict(str_mirror),'modality':dict(lin),'negative':{'contextual':1,'isolated_archive':1936308},'qa':{'blocking_errors':0,'manifest_files':66,'hash_mismatches':0}}
 # M8 release entities / docking readiness.
 S['release']={'entities':{'Membrane proteins':7800,'Canonical compounds':646670,'Compound forms':649792,'Positive evidence':942455,'Unique positive pairs':529168,'Binding-site assertions':55610,'Isoforms':9240,'Processed forms':8454,'Protein complexes':4584,'Canonical diseases':3739},'docking_funnel':S['site']['pair_funnel']}
 (ST/'V711_FIGURE_STATISTICS.json').write_text(json.dumps(S,ensure_ascii=False,indent=2),encoding='utf-8')
 # Panel source and unit registry.
 panels=[
 ('M1','Architecture, entities and source contribution','V7.1.1 release tables; source_registry_v71','entities/evidence records'),('M2','Membrane proteome and five-axis annotation','UniProtKB; HPA; HTP; Membranome; OPM; PDBTM; GO; InterPro/Pfam; Reactome; specialist DBs','unique proteins/annotations'),('M3','Expression and localization','Human Protein Atlas; V7.1.1 semantic correction','measurements, proteins and detected locations'),('M4','Disease ontology landscape','Open Targets; UniProtKB; MONDO; DO; Uberon','unique protein-disease relations and ontology assertions'),('M5','Chemical landscape','ChEMBL; BindingDB; PubChem; BRENDA; standardized compound master','canonical compounds'),('M6','Binding evidence and sites','ChEMBL; BindingDB; PubChem; BRENDA; PDBe/PDB; PDBbind transformed fields','evidence/site assertions and unique pairs'),('M7','Provenance, negatives and QA','All contributing databases; experiment/structure lineage keys','records/lineage groups'),('M8','Release and docking readiness','MemPro V7.1.1 final release','entities and unique pairs')]
 with (ST/'PANEL_SOURCE_VERSION_UNIT_V711.tsv').open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['panel','topic','sources','statistical_unit']);w.writerows(panels)
 report={'release':'V7.1.1','protein_count':sum(pc.values()),'compound_count':compounds,'positive_evidence':sum(es.values()),'site_count':sum(sr.values()),'expression_count':sum(exlayer.values()),'disease_relations':disease_pairs,'pass':sum(pc.values())==7800 and compounds==646670 and sum(es.values())==942455 and sum(sr.values())==55610 and sum(exlayer.values())==3046789 and disease_pairs==6753}
 (QA/'V711_FIGURE_STATS_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
