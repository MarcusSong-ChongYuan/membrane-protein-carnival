import csv,gzip,json,re,collections
from pathlib import Path
import pandas as pd

ROOT=Path(r'D:\7.22\v71_publication_candidate_working');W=ROOT/'work';Q=ROOT/'qa';OUT=ROOT/'outputs'
V7=Path(r'D:\finale\10_MemPro_V7_public_candidate_20260812\01_public_data\binding_evidence_master_v7.tsv.gz')
PROT=W/'protein_publication_crosswalk_v71.tsv.gz'; CHEM=W/'chembl_activity_document_crosswalk_v71.tsv.gz'; V3=W/'v3_1_activity_source_crosswalk_v71.tsv.gz'
PCE=Path(r'D:\7.22\v4_working\releases\release_v4\compound_protein_evidence_v4.tsv.gz')
PD=W/'pdbe_structure_qc_crosswalk_v71.tsv.gz'; ART=W/'pdbe_curated_artifact_exclusions_v71.tsv.gz'; PCRAW=W/'pubchem_raw_assay_citation_crosswalk_v71.tsv.gz'

prot=pd.read_csv(PROT,sep='\t',dtype=str).drop_duplicates('target_uniprot_id').set_index('target_uniprot_id')
chem=pd.read_csv(CHEM,sep='\t',dtype=str).drop_duplicates('chembl_activity_id').set_index('chembl_activity_id')
v3=pd.read_csv(V3,sep='\t',dtype=str).drop_duplicates('activity_measurement_id').set_index('activity_measurement_id')
pdbe=pd.read_csv(PD,sep='\t',dtype=str);pd_by_e=pdbe.drop_duplicates('source_evidence_id').set_index('source_evidence_id');pd_by_r=pdbe.drop_duplicates('source_record_id').set_index('source_record_id')
art=pd.read_csv(ART,sep='\t',dtype=str);artifact=set(art.source_evidence_id)|set(art.source_record_id)

# Recover historical publication by relationship context.
hist=pd.read_csv(PCE,sep='\t',dtype=str,usecols=['relationship_context_id','exp_pubmed','exp_pdb']).drop_duplicates('relationship_context_id').set_index('relationship_context_id')

# Only keep PubChem AIDs present in V7, then derive unambiguous AID -> PubMed mapping.
need=set()
for d in pd.read_csv(V7,sep='\t',dtype=str,usecols=['source_database','source_record_id'],chunksize=200000):
 z=d[d.source_database.eq('PubChem BioAssay')]
 need.update(m.group(1) for v in z.source_record_id.fillna('') if (m:=re.search(r'AID:(\d+)',v)))
aid_pm=collections.defaultdict(set);aid_type=collections.defaultdict(set)
for d in pd.read_csv(PCRAW,sep='\t',dtype=str,usecols=['AID','PubMed ID','Assay Type'],chunksize=600000):
 d=d[d.AID.isin(need)]
 for aid,pm in d.dropna(subset=['PubMed ID']).groupby('AID')['PubMed ID']:aid_pm[aid].update(x for x in pm.astype(str) if x and x!='nan')
 for aid,t in d.dropna(subset=['Assay Type']).groupby('AID')['Assay Type']:aid_type[aid].update(x for x in t.astype(str) if x and x!='nan')
aid_pm={k:';'.join(sorted(v)) for k,v in aid_pm.items() if len(v)==1};aid_type={k:';'.join(sorted(v)) for k,v in aid_type.items()}

extra=['membrane_class_v71','protein_evidence_level_v71','source_database_v71','source_record_id_v71','source_version_v71','source_url_v71','assay_id_v71','document_id_v71','pubmed_id_v71','doi_v71','publication_status_v71','structure_qc_v71','local_enrichment_status_v71']
stats=collections.Counter();sources=collections.Counter()
with gzip.open(V7,'rt',encoding='utf8',newline='') as fi,gzip.open(OUT/'binding_evidence_v71_locally_enriched.tsv.gz','wt',encoding='utf8',newline='',compresslevel=6) as fo:
 rd=csv.DictReader(fi,delimiter='\t');wr=csv.DictWriter(fo,fieldnames=rd.fieldnames+extra,delimiter='\t',lineterminator='\n');wr.writeheader()
 for r in rd:
  stats['rows']+=1;notes=[]
  pr=prot.loc[r['target_uniprot_id']] if r['target_uniprot_id'] in prot.index else None
  r['membrane_class_v71']='' if pr is None else str(pr.get('membrane_class_v52',''));r['protein_evidence_level_v71']='' if pr is None else str(pr.get('evidence_level_v52',''))
  for a,b in [('source_database_v71','source_database'),('source_record_id_v71','source_record_id'),('source_version_v71','source_version'),('source_url_v71','source_url'),('pubmed_id_v71','pubmed_ids'),('doi_v71','doi')]:r[a]=r.get(b,'')
  r['assay_id_v71']='';r['document_id_v71']='';r['structure_qc_v71']='not_structural'
  if r['source_database']=='v3.1_mixed_sources':
   if r['source_record_id'] in v3.index:
    m=v3.loc[r['source_record_id']];r['source_database_v71']=str(m.get('source_database','PubChem BioAssay'));r['assay_id_v71']=str(m.get('assay_id',''));r['source_record_id_v71']=str(m.get('source_assay_id','')) if pd.notna(m.get('source_assay_id')) else str(m.get('assay_id',''));notes.append('historical_activity_assay_recovered')
   if r.get('relationship_context_id','') in hist.index:
    h=hist.loc[r['relationship_context_id']];pm=h.get('exp_pubmed');pdb=h.get('exp_pdb')
    if pd.notna(pm):r['pubmed_id_v71']=str(pm);notes.append('historical_publication_recovered')
    if not r.get('pdb_ids') and pd.notna(pdb):r['pdb_ids']=str(pdb)
  if 'ChEMBL' in r['source_database']:
   aid=str(r['source_record_id']).split(':')[-1]
   if aid in chem.index:
    m=chem.loc[aid];r['assay_id_v71']=str(m.get('chembl_assay_id',''));r['document_id_v71']=str(m.get('chembl_document_id',''));r['source_url_v71']=str(m.get('source_url',''));notes.append('chembl_document_linked')
  if r['source_database']=='PubChem BioAssay':
   m=re.search(r'AID:(\d+)',r['source_record_id'])
   if m:
    aid=m.group(1);r['assay_id_v71']='AID:'+aid
    if not r['pubmed_id_v71'] and aid in aid_pm:r['pubmed_id_v71']=aid_pm[aid];notes.append('pubchem_AID_publication_recovered')
    if aid in aid_type:notes.append('pubchem_assay_type='+aid_type[aid])
  if r['source_database']=='PDBe':
   m=pd_by_e.loc[r['evidence_id']] if r['evidence_id'] in pd_by_e.index else (pd_by_r.loc[r['source_record_id']] if r['source_record_id'] in pd_by_r.index else None)
   if m is not None:r['source_url_v71']=str(m.get('source_url',''));notes.append('pdbe_qc_linked')
   r['structure_qc_v71']='artifact_excluded' if r['evidence_id'] in artifact or r['source_record_id'] in artifact else 'curated_nonartifact'
  anchor=bool(str(r['pubmed_id_v71']).strip() or str(r['doi_v71']).strip() or str(r.get('pdb_ids','')).strip())
  r['publication_status_v71']='citation_or_structure_anchor' if anchor else 'citation_not_recovered';r['local_enrichment_status_v71']=';'.join(notes) if notes else 'no_local_change'
  stats['citation_anchor']+=int(anchor);stats['protein_E1_E2']+=int(r['protein_evidence_level_v71'] in {'E1','E2'});stats['protein_missing']+=int(not r['protein_evidence_level_v71']);stats['specific_source']+=int(r['source_database_v71']!='v3.1_mixed_sources');stats['source_record_present']+=int(bool(str(r['source_record_id_v71']).strip()));stats['structure_artifact_excluded']+=int(r['structure_qc_v71']=='artifact_excluded');sources[r['source_database_v71']]+=1;wr.writerow(r)
report={'counts':dict(stats),'source_counts_v71':dict(sources),'needed_pubchem_AIDs':len(need),'AIDs_with_unique_publication':len(aid_pm),'AIDs_with_assay_type':len(aid_type)};(Q/'V71_LOCAL_ENRICHMENT_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(report,ensure_ascii=False,indent=2))
