import csv,gzip,json,re,hashlib,collections
from pathlib import Path
import pandas as pd

ROOT=Path(r'D:\7.22\v71_publication_candidate_working');W=ROOT/'work';Q=ROOT/'qa';OUT=ROOT/'outputs'
V7=Path(r'D:\finale\10_MemPro_V7_public_candidate_20260812\01_public_data\binding_evidence_master_v7.tsv.gz')

prot=pd.read_csv(W/'protein_publication_crosswalk_v71.tsv.gz',sep='\t',dtype=str,low_memory=False).drop_duplicates('target_uniprot_id').set_index('target_uniprot_id').to_dict('index')
chem=pd.read_csv(W/'chembl_activity_document_crosswalk_v71.tsv.gz',sep='\t',dtype=str,low_memory=False).set_index('chembl_activity_id').to_dict('index')
pdbe=pd.read_csv(W/'pdbe_structure_qc_crosswalk_v71.tsv.gz',sep='\t',dtype=str,low_memory=False)
pdmap={}
for _,r in pdbe.iterrows():pdmap[str(r.source_evidence_id)]=r.to_dict();pdmap[str(r.source_record_id)]=r.to_dict()
art=pd.read_csv(W/'pdbe_curated_artifact_exclusions_v71.tsv.gz',sep='\t',dtype=str,low_memory=False);artifact_ids=set(art.source_evidence_id.astype(str))|set(art.source_record_id.astype(str))
v3=pd.read_csv(W/'v3_1_activity_source_crosswalk_v71.tsv.gz',sep='\t',dtype=str,low_memory=False).set_index('activity_measurement_id').to_dict('index')

# Compact PubChem key -> citation/assay map; duplicates are kept only if values agree.
pc={};amb=set()
for d in pd.read_csv(W/'pubchem_raw_assay_citation_crosswalk_v71.tsv.gz',sep='\t',dtype=str,chunksize=300000,low_memory=False):
 for _,r in d.iterrows():
  k=(str(r.target_uniprot_id),str(r.AID),str(r.SID),str(r.CID),str(r['Activity Name']),str(r['Activity Value [uM]']))
  val=(str(r.get('PubMed ID','') or ''),str(r.get('Assay Type','') or ''),str(r.get('Assay Name','') or ''))
  if k in pc and pc[k]!=val:amb.add(k)
  else:pc[k]=val
for k in amb:pc.pop(k,None)

extra=['membrane_class_v71','protein_evidence_level_v71','source_database_v71','source_record_id_v71','source_version_v71','source_url_v71','assay_id_v71','document_id_v71','pubmed_id_v71','doi_v71','publication_status_v71','structure_qc_v71','local_enrichment_status_v71']
stats=collections.Counter();layers=collections.Counter();sources=collections.Counter()
with gzip.open(V7,'rt',encoding='utf8',newline='') as fi,gzip.open(OUT/'binding_evidence_v71_locally_enriched.tsv.gz','wt',encoding='utf8',newline='',compresslevel=6) as fo:
 rd=csv.DictReader(fi,delimiter='\t');wr=csv.DictWriter(fo,fieldnames=rd.fieldnames+extra,delimiter='\t',lineterminator='\n');wr.writeheader()
 for r in rd:
  stats['rows']+=1;pr=prot.get(r['target_uniprot_id'],{})
  r['membrane_class_v71']=pr.get('membrane_class_v52','');r['protein_evidence_level_v71']=pr.get('evidence_level_v52','')
  r['source_database_v71']=r.get('source_database','');r['source_record_id_v71']=r.get('source_record_id','');r['source_version_v71']=r.get('source_version','');r['source_url_v71']=r.get('source_url','');r['pubmed_id_v71']=r.get('pubmed_ids','');r['doi_v71']=r.get('doi','');r['structure_qc_v71']='not_structural';notes=[]
  if r['source_database']=='v3.1_mixed_sources':
   m=v3.get(r['source_record_id'],{});r['source_database_v71']=m.get('source_database','PubChem BioAssay');r['assay_id_v71']=m.get('assay_id','');r['source_record_id_v71']=m.get('source_assay_id','') or m.get('assay_id','');notes.append('mixed_source_named_from_historical_activity_assay')
  if 'ChEMBL' in r['source_database']:
   aid=str(r['source_record_id']).split(':')[-1];m=chem.get(aid,{})
   if m:r['assay_id_v71']=m.get('chembl_assay_id','');r['document_id_v71']=m.get('chembl_document_id','');r['source_url_v71']=m.get('source_url','');notes.append('chembl_activity_document_linked')
  if r['source_database']=='PubChem BioAssay':
   m=re.search(r'AID:(\d+);SID:(\d+);CID:(\d+)',r['source_record_id'])
   if m:
    k=(r['target_uniprot_id'],m.group(1),m.group(2),m.group(3),str(r.get('activity_type','')),str(r.get('activity_value','')))
    q=pc.get(k)
    if q:r['assay_id_v71']='AID:'+m.group(1);r['pubmed_id_v71']=q[0] if q[0] not in {'nan','None'} else '';notes.append('pubchem_exact_raw_row_linked')
  if r['source_database']=='PDBe':
   m=pdmap.get(r['evidence_id']) or pdmap.get(r['source_record_id'])
   if m:r['source_url_v71']=m.get('source_url','');r['structure_qc_v71']='artifact_excluded' if (r['evidence_id'] in artifact_ids or r['source_record_id'] in artifact_ids) else 'curated_nonartifact';notes.append('pdbe_structure_qc_linked')
  anchor=bool(str(r['pubmed_id_v71']).strip() or str(r['doi_v71']).strip() or str(r.get('pdb_ids','')).strip())
  r['publication_status_v71']='citation_or_structure_anchor' if anchor else 'citation_not_recovered'
  r['local_enrichment_status_v71']=';'.join(notes) if notes else 'no_local_change'
  stats['citation_anchor']+=int(anchor);stats['protein_E1_E2']+=int(r['protein_evidence_level_v71'] in {'E1','E2'});stats['protein_mapping_missing']+=int(not r['protein_evidence_level_v71']);stats['specific_source']+=int(r['source_database_v71']!='v3.1_mixed_sources');stats['source_record_present']+=int(bool(str(r['source_record_id_v71']).strip()));stats['structure_artifact_excluded']+=int(r['structure_qc_v71']=='artifact_excluded');sources[r['source_database_v71']]+=1;wr.writerow(r)
report={'counts':dict(stats),'source_counts_v71':dict(sources),'pubchem_exact_key_count':len(pc),'pubchem_ambiguous_key_count':len(amb)};(Q/'V71_LOCAL_ENRICHMENT_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(report,ensure_ascii=False,indent=2))
