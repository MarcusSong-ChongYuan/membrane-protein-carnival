import csv
import gzip
import hashlib
import json
from pathlib import Path

AUTO = Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data")
WORK = Path(r"D:\finale\12_V7_final_completion_20260813")
EVIDENCE = WORK / "02_be_recompute" / "protein_compound_evidence_BE_recomputed_v7.tsv.gz"
OLD_PAIR = AUTO / "protein_compound_pair_v7.tsv.gz"
NEW_PAIR = WORK / "02_be_recompute" / "protein_compound_pair_BE_recomputed_v7.tsv.gz"
COMPLEX = WORK / "06_remaining_modules" / "protein_complex_v7_validated.tsv.gz"
OLD_COMPONENT = AUTO / "complex_component_v7.tsv.gz"
NEW_COMPONENT = WORK / "06_remaining_modules" / "complex_component_v7_validated.tsv.gz"
FROZEN_COMPONENT = WORK / "06_remaining_modules" / "complex_component_frozen_v7.tsv.gz"

def reader(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="")

def writer(path):
    return gzip.open(path, "wt", encoding="utf-8", newline="")

with reader(COMPLEX) as handle:
    public_complexes = {r["complex_target_id"] for r in csv.DictReader(handle, delimiter="\t")}

public_n = frozen_n = 0
with reader(OLD_COMPONENT) as source, writer(NEW_COMPONENT) as out1, writer(FROZEN_COMPONENT) as out2:
    rows = csv.DictReader(source, delimiter="\t")
    w1 = csv.DictWriter(out1, fieldnames=rows.fieldnames, delimiter="\t")
    w2 = csv.DictWriter(out2, fieldnames=rows.fieldnames, delimiter="\t")
    w1.writeheader(); w2.writeheader()
    for row in rows:
        if row["complex_target_id"] in public_complexes:
            w1.writerow(row); public_n += 1
        else:
            w2.writerow(row); frozen_n += 1

old = {}
with reader(OLD_PAIR) as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        old[(row["target_uniprot_id"], row["compound_internal_id"])] = row

agg = {}
with reader(EVIDENCE) as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        key = (row["target_uniprot_id"], row["compound_internal_id"])
        a = agg.setdefault(key, {"n": 0, "tiers": set(), "db": set(), "exp": set(), "str": set(), "aid": set(), "mod": set()})
        a["n"] += 1
        a["tiers"].add(row["evidence_tier_recomputed_v7"])
        if row["source_database"]: a["db"].add(row["source_database"])
        if row["experiment_lineage_key_v7"]: a["exp"].add(row["experiment_lineage_key_v7"])
        if row["structure_lineage_key_v7"]: a["str"].add(row["structure_lineage_key_v7"])
        if row["evidence_modality_v7"]: a["mod"].add(row["evidence_modality_v7"])
        if row["source_database"].lower().startswith("pubchem") and row["source_record_id"]: a["aid"].add(row["source_record_id"])

fields = ["pair_id","target_uniprot_id","compound_internal_id","evidence_count","best_evidence_tier","distinct_database_count","source_databases","independent_experiment_count","independent_structure_count","independent_pubchem_assay_count","independent_evidence_modality_count","positive_negative_conflict_status","release_version"]
rank = {"BE1": 1, "BE2": 2, "BE3": 3}
with writer(NEW_PAIR) as handle:
    out = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); out.writeheader()
    for key in sorted(agg):
        target, compound = key; a = agg[key]; prior = old.get(key, {})
        pair_id = prior.get("pair_id") or "PAIR-" + hashlib.sha256((target + "\t" + compound).encode()).hexdigest()[:16].upper()
        out.writerow({"pair_id":pair_id,"target_uniprot_id":target,"compound_internal_id":compound,"evidence_count":a["n"],"best_evidence_tier":min(a["tiers"],key=rank.get),"distinct_database_count":len(a["db"]),"source_databases":";".join(sorted(a["db"])),"independent_experiment_count":len(a["exp"]),"independent_structure_count":len(a["str"]),"independent_pubchem_assay_count":len(a["aid"]),"independent_evidence_modality_count":len(a["mod"]),"positive_negative_conflict_status":prior.get("positive_negative_conflict_status","NO_CONFLICT"),"release_version":"MemPro_V7.0.0"})

report = {"public_complexes":len(public_complexes),"public_complex_components":public_n,"frozen_complex_components":frozen_n,"recomputed_positive_pairs":len(agg),"old_positive_pairs":len(old),"pair_key_symmetric_difference":len(set(agg)^set(old))}
(WORK/"07_qa"/"FINAL_DEPENDENT_RELATIONS_REPORT.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(report,ensure_ascii=False,indent=2))
