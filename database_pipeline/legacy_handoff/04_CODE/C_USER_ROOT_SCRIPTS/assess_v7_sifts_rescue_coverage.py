import csv,gzip,re,json
from collections import Counter
from pathlib import Path

site=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813\02_frozen_review\binding_site_coordinate_frozen_v7.tsv.gz")
sifts=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites\sifts_residue_mapping_v7.checkpoint.tsv")
need=set(); need_rows=Counter()
with gzip.open(site,'rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  for p in re.split(r'[;,|]',r['pdb_ids'] or ''):
   p=p.strip().upper()
   if re.fullmatch(r'[0-9][A-Z0-9]{3}',p): need.add(p);need_rows[p]+=1
have=set()
with sifts.open('r',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):have.add(r['pdb_id'].upper())
report={'need_unique_pdb':len(need),'have_unique_pdb':len(have),'covered_needed_pdb':len(need&have),'missing_needed_pdb':len(need-have),'rows_with_at_least_one_covered_pdb':sum(1 for p,n in need_rows.items() if p in have),'top_missing_by_rows':[(p,need_rows[p]) for p in sorted(need-have,key=lambda p:-need_rows[p])[:30]]}
print(json.dumps(report,indent=2))
out=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites\pdb_ids_for_frozen_rescue.txt")
out.write_text('\n'.join(sorted(need))+'\n',encoding='ascii')
