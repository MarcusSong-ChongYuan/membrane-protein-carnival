import csv,gzip,json
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites");IN=ROOT/'binding_site_coordinate_remaining_v703.tsv.gz';OUT=ROOT/'binding_site_coordinate_rescued_partial_v705.tsv.gz';REM=ROOT/'binding_site_coordinate_remaining_after_partial_v705.tsv.gz';counts=Counter()
with gzip.open(IN,'rt',encoding='utf-8',newline='') as f:
 rd=csv.DictReader(f,delimiter='\t');extra=['partial_rescue_status_v705','partial_rescue_rule_v705','partial_rescue_fraction_v705'];fields=rd.fieldnames+extra;rows=list(rd)
with gzip.open(OUT,'wt',encoding='utf-8',newline='') as fo,gzip.open(REM,'wt',encoding='utf-8',newline='') as fr:
 wo=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wo.writeheader();wr.writeheader()
 for r in rows:
  try:req=int(r.get('uniprot_requested_count_v703') or 0);mat=int(r.get('uniprot_mapped_count_v703') or 0);obs=int(r.get('uniprot_observed_count_v703') or 0);frac=mat/req if req else 0
  except:req=mat=obs=0;frac=0
  chain=r.get('uniprot_rescue_chain_v703','').strip();eligible=r.get('uniprot_rescue_status_v703')=='REMAIN_UNIPROT_PARTIAL_MATCH' and req>=2 and mat>=2 and frac>=.8 and obs==mat and bool(chain)
  if eligible:status='PUBLIC_COORDINATE_PARTIAL_MEDIUM';rule='unique recorded target chain; >=80% requested UniProt residues mapped; every mapped residue observed'
  else:status='REMAIN';rule='does not meet conservative partial-coordinate publication threshold'
  r.update({'partial_rescue_status_v705':status,'partial_rescue_rule_v705':rule,'partial_rescue_fraction_v705':f'{frac:.4f}'});counts[status]+=1;(wo if status.startswith('PUBLIC') else wr).writerow(r)
report={'input':len(rows),'status':dict(counts),'threshold':'requested>=2, mapped>=2, mapped/requested>=0.8, all mapped observed, chain nonempty','confidence':'MEDIUM; missing requested residues are not assigned coordinates'};(ROOT/'T26_T28_PARTIAL_RESCUE_V705_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
