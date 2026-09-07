import csv,gzip,json,re,hashlib
from pathlib import Path
import pandas as pd

ROOT=Path(r'D:\7.22\v71_publication_candidate_working');W=ROOT/'work';Q=ROOT/'qa'
V7=Path(r'D:\finale\10_MemPro_V7_public_candidate_20260812\01_public_data\binding_evidence_master_v7.tsv.gz')
CHEM=Path(r'D:\7.22\membrane_master_v5_working\v54_pgd_binding_working\intermediate\chembl_37_binding_activities_v54.tsv')
RAWPC=Path(r'D:\7.22\evidence_expansion_v2_working\raw\pubchem_202607\protein_concise')
PD=Path(r'D:\7.22\evidence_expansion_v2_working\runs\incremental_20260727\staging\pdbe\pdbe_normalized_evidence_v3.tsv.gz')
PDEX=Path(r'D:\7.22\evidence_expansion_v2_working\runs\incremental_20260727\staging\pdbe\pdbe_excluded_artifacts_v3.tsv.gz')

def normnum(v):
 try:return round(float(v),9)
 except:return None

# ChEMBL activity crosswalk: all fields needed to retrieve document citation in a subsequent API pass.
c=pd.read_csv(CHEM,sep='\t',dtype=str,low_memory=False)
c[['chembl_activity_id','chembl_assay_id','chembl_document_id','chembl_target_id','chembl_molecule_id','source_url','source_version']].to_csv(W/'chembl_activity_document_crosswalk_v71.tsv.gz',sep='\t',index=False,compression='gzip')

# PubChem raw evidence lookup keyed by deterministic source tuple.
rows=[]
for f in RAWPC.glob('*.csv.gz'):
 target=f.stem.split('.')[0]
 try:d=pd.read_csv(f,dtype=str,low_memory=False)
 except:continue
 d['target_uniprot_id']=target
 rows.append(d[['target_uniprot_id','AID','SID','CID','Activity Outcome','Activity Name','Activity Qualifier','Activity Value [uM]','Assay Name','Assay Type','PubMed ID']])
pc=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
pc.to_csv(W/'pubchem_raw_assay_citation_crosswalk_v71.tsv.gz',sep='\t',index=False,compression='gzip')

# PDBe exact source-record metadata and already curated artifact identifiers.
pdbe=pd.read_csv(PD,sep='\t',dtype=str,low_memory=False)
pdbe[['source_evidence_id','source_record_id','target_uniprot_id','pdb_ids','pdb_chain_ids','ligand_het_id','binding_site_residues','artifact_class','exclusion_reason','source_url']].to_csv(W/'pdbe_structure_qc_crosswalk_v71.tsv.gz',sep='\t',index=False,compression='gzip')
art=pd.read_csv(PDEX,sep='\t',dtype=str,usecols=['source_evidence_id','source_record_id','pdb_ids','pdb_chain_ids','ligand_het_id','artifact_class','exclusion_reason'],low_memory=False)
art.to_csv(W/'pdbe_curated_artifact_exclusions_v71.tsv.gz',sep='\t',index=False,compression='gzip')

# Coverage audit against V7 source IDs.
chem_ids=set(c.chembl_activity_id.astype(str));pc_aids=set(pc.AID.dropna().astype(str));pd_ids=set(pdbe.source_evidence_id.astype(str))|set(pdbe.source_record_id.astype(str))
cnt={'chembl_rows':0,'chembl_activity_crosswalk':0,'pubchem_rows':0,'pubchem_AID_parseable':0,'pubchem_AID_in_raw':0,'pdbe_rows':0,'pdbe_exact_crosswalk':0};
for d in pd.read_csv(V7,sep='\t',dtype=str,chunksize=120000,low_memory=False):
 z=d[d.source_database.str.contains('ChEMBL',na=False)];cnt['chembl_rows']+=len(z);cnt['chembl_activity_crosswalk']+=sum(str(v).split(':')[-1] in chem_ids for v in z.source_record_id)
 z=d[d.source_database.eq('PubChem BioAssay')];cnt['pubchem_rows']+=len(z)
 for v in z.source_record_id.fillna(''):
  m=re.search(r'AID:(\d+)',v)
  if m:cnt['pubchem_AID_parseable']+=1;cnt['pubchem_AID_in_raw']+=int(m.group(1) in pc_aids)
 z=d[d.source_database.eq('PDBe')];cnt['pdbe_rows']+=len(z);cnt['pdbe_exact_crosswalk']+=sum(str(a) in pd_ids or str(b) in pd_ids for a,b in zip(z.evidence_id,z.source_record_id))
report={'local_resources':cnt,'pubchem_raw_rows':len(pc),'pubchem_unique_AID':pc.AID.nunique(),'chembl_activity_rows':len(c),'chembl_documents':c.chembl_document_id.nunique(),'pdbe_rows':len(pdbe),'pdbe_curated_artifact_rows':len(art)}
(Q/'V71_LOCAL_RECOVERY_COVERAGE.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2))
