from collections import Counter
from pathlib import Path

import pandas as pd


SOURCE = Path(r"C:\Users\Administrator\Desktop\FINAL\MemPro_V7.2\01_release_tables\positive_interaction_evidence_v72.tsv.gz")
OUT = Path(r"D:\finale\24_MemPro_V7.2_NAR_analysis_and_figures_20260818\03_source_data\F1C_source_evidence_tier_flow.tsv")


def split_values(value):
    return [x.strip() for x in str(value or "").replace("|", ";").split(";") if x.strip()]


counter = Counter()
usecols = ["source_database", "evidence_type", "evidence_tier_recomputed_v7", "be_disposition_v7"]
for chunk in pd.read_csv(SOURCE, sep="\t", compression="gzip", dtype=str, usecols=usecols, chunksize=200_000):
    chunk = chunk[chunk["be_disposition_v7"].eq("PUBLIC")]
    for row in chunk.itertuples(index=False):
        for source in split_values(row.source_database):
            counter[(source, row.evidence_type, row.evidence_tier_recomputed_v7)] += 1

pd.DataFrame([
    {"source_database": key[0], "evidence_type": key[1], "evidence_tier": key[2], "evidence_record_count": value}
    for key, value in counter.items()
]).sort_values("evidence_record_count", ascending=False).to_csv(OUT, sep="\t", index=False)
print(f"rows={len(counter)}")
