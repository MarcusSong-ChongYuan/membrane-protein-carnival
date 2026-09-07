from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd


REL = Path(
    r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730"
)
OUT = Path(
    r"D:\7.22\evidence_expansion_v2_working\figures\v62_main_20260730\qa"
)
OUT.mkdir(parents=True, exist_ok=True)


def top(series: pd.Series, n: int = 20) -> dict[str, int]:
    return {
        str(k): int(v)
        for k, v in series.fillna("<NA>").astype(str).value_counts().head(n).items()
    }


protein_cols = [
    "target_uniprot_id",
    "membrane_class_v52",
    "evidence_level_v52",
    "functional_primary_class_v5",
    "transmembrane_count_v5",
    "sequence_length_v53",
    "opm_present_v5",
    "pdbtm_present_v5",
    "pdb_ids",
    "alphafolddb_ids",
    "hpa_predicted_membrane_v5",
    "htp_present_v5",
    "membranome_present_v5",
    "oligomeric_state_consensus_v62",
    "subunit_evidence_grade_v62",
    "subunit_conflict_flag_v62",
    "hpa_mapping_status_v62",
    "hpa_tissue_detected_count_v62",
    "hpa_cell_type_detected_count_v62",
    "hpa_main_locations_v62",
    "hpa_additional_locations_v62",
    "best_binding_evidence_level",
]
protein = pd.read_csv(
    REL / "human_membrane_protein_master_v6_2.tsv",
    sep="\t",
    usecols=protein_cols,
    low_memory=False,
)

disease = pd.read_csv(
    REL / "protein_gene_disease_relations_v6_2.tsv",
    sep="\t",
    low_memory=False,
)

with gzip.open(REL / "protein_compound_summary_v6_2.tsv.gz", "rt", encoding="utf-8") as fh:
    pair = pd.read_csv(fh, sep="\t", low_memory=False)

with gzip.open(
    REL / "protein_compound_negative_summary_v6_2.tsv.gz",
    "rt",
    encoding="utf-8",
) as fh:
    negative = pd.read_csv(fh, sep="\t", low_memory=False)

with gzip.open(REL / "binding_site_instances_v6_2.tsv.gz", "rt", encoding="utf-8") as fh:
    site = pd.read_csv(fh, sep="\t", low_memory=False)

report = {
    "protein_rows": int(len(protein)),
    "membrane_class": top(protein["membrane_class_v52"]),
    "evidence_level": top(protein["evidence_level_v52"]),
    "functional_class": top(protein["functional_primary_class_v5"], 30),
    "tm_count": top(protein["transmembrane_count_v5"], 30),
    "oligomeric_state": top(protein["oligomeric_state_consensus_v62"], 30),
    "subunit_evidence_grade": top(protein["subunit_evidence_grade_v62"], 20),
    "hpa_mapping": top(protein["hpa_mapping_status_v62"]),
    "main_location": top(protein["hpa_main_locations_v62"], 30),
    "disease_rows": int(len(disease)),
    "disease_ontology": top(disease["disease_ontology"], 30),
    "disease_evidence_level": top(disease["disease_evidence_level"], 20),
    "disease_names": top(disease["disease_name"], 30),
    "disease_sources": top(disease["source_database"], 30),
    "pair_rows": int(len(pair)),
    "pair_best_evidence": top(pair["best_binding_evidence_level"], 10),
    "pair_sources": top(pair["independent_source_count"], 20),
    "negative_pair_rows": int(len(negative)),
    "negative_conflict": top(negative["positive_negative_conflict_flag"], 10),
    "site_rows": int(len(site)),
    "site_type": top(site["site_type"], 30),
    "site_qc": top(site["record_qc_status"], 20),
    "site_sources": top(site["source_database"], 30),
}
path = OUT / "FIGURE_INPUT_PROFILE_V62.json"
path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False, indent=2))
