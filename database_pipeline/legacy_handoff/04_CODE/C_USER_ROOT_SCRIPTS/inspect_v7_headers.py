import csv
import gzip
from pathlib import Path

files = [
    r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data\complex_component_v7.tsv.gz",
    r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data\protein_compound_pair_v7.tsv.gz",
    r"D:\finale\12_V7_final_completion_20260813\02_be_recompute\protein_compound_evidence_BE_recomputed_v7.tsv.gz",
]
for path in files:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        print(Path(path).name, next(csv.reader(handle, delimiter="\t")))
