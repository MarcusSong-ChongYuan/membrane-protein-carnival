from pathlib import Path

source = Path(r"C:\Users\Administrator\qa_and_freeze_v7_final_completion.py").read_text(encoding="utf-8")
source = source.replace(
    "PAIR=AUTO/'protein_compound_pair_v7.tsv.gz'",
    "PAIR=W/'02_be_recompute'/'protein_compound_pair_BE_recomputed_v7.tsv.gz'",
)
source = source.replace(
    "CC=AUTO/'complex_component_v7.tsv.gz'",
    "CC=W/'06_remaining_modules'/'complex_component_v7_validated.tsv.gz'",
)
source = source.replace(
    "(W/'06_remaining_modules'/'protein_disease_frozen_v7.tsv.gz','protein_disease_frozen_v7.tsv.gz')",
    "(W/'06_remaining_modules'/'complex_component_frozen_v7.tsv.gz','complex_component_frozen_v7.tsv.gz'),(W/'06_remaining_modules'/'protein_disease_frozen_v7.tsv.gz','protein_disease_frozen_v7.tsv.gz')",
)
exec(compile(source, "qa_and_freeze_v7_final_completion.corrected.py", "exec"))
