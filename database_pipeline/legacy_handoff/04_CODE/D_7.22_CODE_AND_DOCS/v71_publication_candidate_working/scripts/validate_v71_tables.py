import csv,gzip,json,collections,re
from pathlib import Path
R=Path(r'D:\7.22\v71_publication_candidate_working');O=R/'outputs';Q=R/'qa'
files={'evidence':O/'binding_evidence_public_candidate_v7_1.tsv.gz','frozen':O/'binding_evidence_frozen_nonpublic_v7_1.tsv.gz','relation':O/'protein_compound_relation_v7_1.tsv.gz','experiment':O/'independent_experiment_master_v7_1.tsv.gz','publication':O/'publication_entity_v7_1.tsv.gz','crosswalk':O/'public_evidence_source_crosswalk_v7_1.tsv.gz'}
def audit(name,p,key=None):
 n=0;ks=set();dup=0;bad=collections.Counter()
 with gzip.open(p,'rt',encoding='utf8',newline='') as f:
  rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames
  for r in rd:
   n+=1
   if key:
    v=r.get(key,'');dup+=int(v in ks);ks.add(v)
   for v in r.values():
    if v in {'#REF!','#DIV/0!','#VALUE!','#NAME?','#N/A'}:bad[v]+=1
 return {'rows':n,'columns':len(fields),'duplicate_key':dup,'formula_like_errors':dict(bad)}
A={k:audit(k,p,'evidence_id' if k=='evidence' else ('independent_experiment_key' if k=='experiment' else ('publication_id' if k=='publication' else None))) for k,p in files.items()}
rep=json.loads((Q/'V71_RELATIONAL_RELEASE_REPORT.json').read_text());checks={'evidence_count_match':A['evidence']['rows']==rep['evidence_rows'],'relation_count_match':A['relation']['rows']==rep['unique_pairs'],'experiment_count_match':A['experiment']['rows']==rep['independent_experiments'],'publication_count_match':A['publication']['rows']==rep['independent_publications'],'crosswalk_count_match':A['crosswalk']['rows']==rep['evidence_rows'],'evidence_id_unique':A['evidence']['duplicate_key']==0};out={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'files':A,'remaining_release_blockers':['exact license version/terms confirmation for conditional sources before public distribution','external-device reproducibility test','refresh figures and presentation from V7.1 frozen counts','website/API implementation and manuscript not yet started']};(Q/'V71_FULL_TABLE_QA.json').write_text(json.dumps(out,indent=2),encoding='utf8');print(json.dumps(out,indent=2))