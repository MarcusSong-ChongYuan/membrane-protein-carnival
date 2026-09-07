import csv,gzip
from pathlib import Path

p=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813\01_data\protein_membrane_state_v7.tsv.gz")
with gzip.open(p,'rt',encoding='utf-8',newline='') as f:
 r=csv.DictReader(f,delimiter='\t');print(r.fieldnames)
 for i,row in enumerate(r):
  print({k:v for k,v in row.items() if v})
  if i==5:break
