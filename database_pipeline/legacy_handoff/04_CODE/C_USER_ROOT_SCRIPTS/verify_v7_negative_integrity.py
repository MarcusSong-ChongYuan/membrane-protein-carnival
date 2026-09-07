import csv,gzip,json
from collections import Counter
from pathlib import Path
root=Path(r"D:\finale\12_V7_final_completion_20260813\04_negative_conflicts")
def load(p,col):
 s=[]
 with gzip.open(p,'rt',encoding='utf-8',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):s.append(r[col])
 return s
comp=load(root/'compound_master_positive_negative_union_v7.tsv.gz','compound_internal_id');forms=load(root/'compound_form_positive_negative_union_v7.tsv.gz','compound_form_id');cs=set(comp);fs=set(forms);negids=[];pairs=set();badc=badf=0
with gzip.open(root/'negative_evidence_v7.tsv.gz','rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  nid=r.get('negative_evidence_id_v62') or r['source_evidence_id'];negids.append(nid);badc+=r['compound_internal_id'] not in cs;badf+=bool(r.get('compound_form_id')) and r['compound_form_id'] not in fs;pairs.add((r['target_uniprot_id'],r['compound_internal_id']))
report={'compound_rows':len(comp),'compound_pk_duplicates':len(comp)-len(cs),'form_rows':len(forms),'form_pk_duplicates':len(forms)-len(fs),'negative_rows':len(negids),'negative_id_duplicates':len(negids)-len(set(negids)),'negative_compound_fk_failures':badc,'negative_form_fk_failures':badf,'unique_negative_pairs':len(pairs),'blocking_errors':len(comp)-len(cs)+len(forms)-len(fs)+len(negids)-len(set(negids))+badc+badf}
report['status']='PASS_NEGATIVE_INTEGRITY' if report['blocking_errors']==0 else 'FAIL_NEGATIVE_INTEGRITY';(root/'T30_T31_NEGATIVE_INTEGRITY_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
