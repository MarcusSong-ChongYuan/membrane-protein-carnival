import csv
from pathlib import Path

paths = [
    Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites\sifts_residue_mapping_v7.checkpoint.tsv"),
    Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813\01_data\protein_isoform_v7.tsv.gz"),
    Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813\01_data\protein_master_v7.tsv.gz"),
]
import gzip
for p in paths:
    opener = gzip.open if str(p).endswith(".gz") else open
    with opener(p, "rt", encoding="utf-8", newline="") as h:
        r=csv.DictReader(h,delimiter="\t")
        print("FILE",p)
        print("FIELDS",r.fieldnames)
        for i,row in enumerate(r):
            print({k:v for k,v in row.items() if v})
            if i==1: break
