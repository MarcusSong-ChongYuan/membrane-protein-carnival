import csv,gzip,json,hashlib,re
from pathlib import Path
import pandas as pd

ROOT=Path(r'D:\7.22\v71_publication_candidate_working'); W=ROOT/'work';Q=ROOT/'qa';W.mkdir(exist_ok=True);Q.mkdir(exist_ok=True)
V7=Path(r'D:\finale\10_MemPro_V7_public_candidate_20260812\01_public_data\binding_evidence_master_v7.tsv.gz')
PROT=Path(r'D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay\human_membrane_protein_master_v6_4_candidate.tsv.gz')
ACT=Path(r'D:\7.22\v4_working\releases\release_v4\activity_measurements_v4.tsv')
ASSAY=Path(r'D:\7.22\v4_working\releases\release_v4\assays_v4.tsv')
PUB=Path(r'D:\7.22\v4_working\releases\release_v4\publications_v4.tsv')

# Protein release-tier and membrane-class crosswalk
p=pd.read_csv(PROT,sep='\t',usecols=['target_uniprot_id','membrane_class_v52','evidence_level_v52','release_scope_v52','website_default_v52','sequence_resolved_accession_v53'],low_memory=False)
p=p.drop_duplicates('target_uniprot_id');p.to_csv(W/'protein_publication_crosswalk_v71.tsv.gz',sep='\t',index=False,compression='gzip')

# Historical v3.1 activity -> assay -> named source reconstruction
a=pd.read_csv(ACT,sep='\t',dtype=str,low_memory=False)
s=pd.read_csv(ASSAY,sep='\t',dtype=str,low_memory=False)
x=a.merge(s[['assay_id','source_database','source_assay_id','assay_description','assay_type_original','source_record_qc_status']],on='assay_id',how='left')
x=x[['activity_measurement_id','relationship_context_id','assay_id','target_uniprot_id','drug_id','source_database','source_assay_id','assay_description','assay_type_original','source_record_qc_status']]
x.to_csv(W/'v3_1_activity_source_crosswalk_v71.tsv.gz',sep='\t',index=False,compression='gzip')

# Publication registry retained from V4
pub=pd.read_csv(PUB,sep='\t',dtype=str,low_memory=False);pub.to_csv(W/'historical_publication_registry_v71.tsv',sep='\t',index=False)

# Audit how much V7 mixed-source material can be named locally by ACT id
actmap=x.set_index('activity_measurement_id').to_dict('index')
counts={};examples=[];n=matched=0
for d in pd.read_csv(V7,sep='\t',dtype=str,chunksize=120000,low_memory=False):
 z=d[d.source_database.eq('v3.1_mixed_sources')]
 for _,r in z.iterrows():
  n+=1;m=actmap.get(r.source_record_id)
  if m and m.get('source_database'):
   matched+=1;src=m['source_database'];counts[src]=counts.get(src,0)+1
   if len(examples)<30:examples.append({'evidence_id':r.evidence_id,'old_source':'v3.1_mixed_sources','recovered_source':src,'source_record_id':r.source_record_id,'assay_id':m.get('assay_id'),'source_assay_id':m.get('source_assay_id'),'assay_description':m.get('assay_description')})
report={'mixed_rows':n,'local_activity_crosswalk_matches':matched,'unmatched':n-matched,'recovered_source_counts':counts,'examples':examples}
(Q/'V71_MIXED_SOURCE_LOCAL_RECOVERY.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(report,ensure_ascii=False,indent=2))
