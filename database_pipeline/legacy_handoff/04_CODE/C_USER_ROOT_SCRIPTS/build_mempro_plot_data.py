import csv, gzip, json, math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(r"D:\finale\02_V6.3_candidate_20260803")
OUT = Path(r"C:\Users\Administrator\mempro_ppt_final_redraw\data\plot_data.json")

def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else open(path, "r", encoding="utf-8", newline="")

def rows(path):
    with open_text(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")

def num(v, default=0.0):
    try:
        return float(v)
    except Exception:
        return default

def quantiles(values):
    vals = sorted(float(x) for x in values)
    if not vals:
        return {"n": 0, "q1": 0, "median": 0, "q3": 0, "min": 0, "max": 0}
    def q(p):
        pos = (len(vals)-1)*p
        lo, hi = int(math.floor(pos)), int(math.ceil(pos))
        return vals[lo] if lo == hi else vals[lo]*(hi-pos)+vals[hi]*(pos-lo)
    return {"n": len(vals), "q1": q(.25), "median": q(.5), "q3": q(.75), "min": vals[0], "max": vals[-1]}

data = {"release": "MemPro V6.3.1", "generated_from": "V6.3 candidate + V6.2 frozen modules"}

# M2 protein classification.
protein_path = ROOT / "05_protein_annotation" / "human_membrane_protein_master_v6_3_candidate.tsv.gz"
class_evidence = defaultdict(Counter)
role_counts = Counter()
tm_bins = Counter()
protein_class = {}
for r in rows(protein_path):
    cls = r.get("membrane_class_v52") or r.get("final_membrane_class_v51") or "unknown"
    ev = r.get("evidence_level_v52") or "E0"
    class_evidence[cls][ev] += 1
    role = r.get("membrane_role_primary_v63") or "unclassified"
    role_counts[role] += 1
    protein_class[r["target_uniprot_id"]] = cls
    tm = int(num(r.get("transmembrane_count_v5") or r.get("transmembrane_count"), 0))
    if tm == 0: tm_bins["0"] += 1
    elif tm == 1: tm_bins["1"] += 1
    elif tm <= 6: tm_bins["2–6"] += 1
    elif tm <= 12: tm_bins["7–12"] += 1
    else: tm_bins[">12"] += 1

axis_labels = {
    "结构家族": 98.7178321360371,
    "生物过程": 89.55169591706829,
    "分子功能": 83.35000454669455,
    "具体膜角色": 75.96617259252524,
    "Reactome通路": 63.435482404292074,
    "专科层级": 38.04674002000545,
}
data["M2"] = {
    "total": sum(sum(c.values()) for c in class_evidence.values()),
    "class_evidence": {k: dict(v) for k,v in class_evidence.items()},
    "axis_coverage": axis_labels,
    "top_roles": role_counts.most_common(9),
    "tm_bins": dict(tm_bins),
}

# M3 expression and localization.
expr_path = ROOT / "01_baseline_v6_2" / "expression_location_summary_v2.tsv.gz"
mapping = Counter()
breadth = defaultdict(list)
for r in rows(expr_path):
    status = r.get("hpa_mapping_status_v62", "")
    detected = int(num(r.get("hpa_tissue_detected_count_v62"), 0))
    if status == "mapped" and detected > 0: mapping["mapped_detected"] += 1
    elif status == "mapped": mapping["mapped_zero"] += 1
    else: mapping["unmapped_missing"] += 1
    cls = protein_class.get(r.get("target_uniprot_id", ""), "unknown")
    if status == "mapped": breadth[cls].append(detected)
loc_counts = Counter()
for r in rows(ROOT / "01_baseline_v6_2" / "protein_subcellular_localization_v2.tsv.gz"):
    loc_counts[r.get("location_term") or "Unknown"] += 1
zscore = []
zpath = Path(r"D:\finale\05_V6.3.1_双页模块PPT重制_20260803\figures\data\M3_tissue_expression_zscore.tsv")
for r in rows(zpath): zscore.append(r)
data["M3"] = {
    "mapping": dict(mapping),
    "breadth": {k: quantiles(v) for k,v in breadth.items()},
    "top_localizations": loc_counts.most_common(10),
    "zscore": zscore,
}

# M4 disease identity, ontology, therapeutic areas, anatomy and evidence channels.
dval = json.loads((ROOT / "09_qa" / "disease_v63_working" / "DISEASE_MODULE_V0_1_VALIDATION.json").read_text(encoding="utf-8"))
ta = Counter()
for r in rows(ROOT / "08_figures" / "disease_v63_working" / "M4_therapeutic_area_counts.tsv"):
    ta[r["therapeutic_area_name"]] = int(r["protein_disease_relation_count"])
evidence_flags = Counter()
flag_map = {
    "human_genetic_flag":"Human genetics", "clinical_genetic_flag":"Clinical genetics",
    "somatic_mutation_flag":"Somatic mutation", "functional_flag":"Functional experiment",
    "animal_model_flag":"Animal model", "expression_flag":"Disease expression",
    "literature_flag":"Literature", "uniprot_disease_flag":"UniProt curated",
}
for r in rows(ROOT / "01_baseline_v6_2" / "protein_gene_disease_relations_v6_2.tsv"):
    for col,label in flag_map.items():
        if str(r.get(col, "")).lower() in {"true","1","yes"}: evidence_flags[label] += 1
data["M4"] = {
    "relations": dval["counts"]["canonical_protein_disease_relations"],
    "source_diseases": dval["counts"]["unique_source_diseases"],
    "canonicalization": dval["counts"]["canonicalization_status"],
    "therapeutic_areas": ta.most_common(),
    "anatomy_assertion": dval["counts"]["anatomy_assertion_status"],
    "evidence_channels": evidence_flags.most_common(),
}

# M5 compound status UpSet, Lipinski description and scaffold diversity.
status = []
for r in rows(Path(r"D:\finale\05_V6.3.1_双页模块PPT重制_20260803\figures\data\M5_biological_status_intersections.tsv")):
    status.append({"combination": r["combination"], "count": int(r["count"])})
lip = []
for r in rows(Path(r"D:\finale\05_V6.3.1_双页模块PPT重制_20260803\figures\data\M5_lipinski_descriptor_pass.tsv")):
    label = next(k for k in r if k != "percent")
    lip.append({"label": r[label], "percent": float(r["percent"])})
scaf = next(rows(ROOT / "04_disease" / "M5E_scaffold_diversity_summary.tsv"))
scaf_class = []
for r in rows(ROOT / "04_disease" / "scaffold_diversity_by_structural_class_v0_1.tsv"):
    scaf_class.append({"class":r["computed_structural_class"],"compounds":int(r["scaffold_bearing_compounds"]),"scaffolds":int(r["unique_murcko_scaffolds"]),"ratio":float(r["scaffold_diversity_ratio"])})
data["M5"] = {
    "core_qc": int(scaf["core_qc_compounds"]),
    "status_intersections": status,
    "lipinski": lip,
    "scaffold": {k:(float(v) if "." in v else int(v)) for k,v in scaf.items() if k not in {"sampling_used","rdkit_version"}},
    "scaffold_by_class": scaf_class,
}

# M6 binding evidence source/tier, lineage and site readiness.
source_tier = defaultdict(Counter)
for r in rows(ROOT / "01_baseline_v6_2" / "binding_evidence_master_v6_2.tsv.gz"):
    if r.get("default_release_inclusion") not in {"1","True","true"}: continue
    source_tier[r.get("source_database") or "Unknown"][r.get("evidence_tier") or "Unknown"] += 1
site_type = Counter(); site_specificity = Counter(); site_total = 0
for r in rows(ROOT / "01_baseline_v6_2" / "binding_site_instances_v6_2.tsv.gz"):
    if r.get("record_qc_status") != "ok": continue
    site_total += 1
    site_type[r.get("site_type") or "Unknown"] += 1
    site_specificity[r.get("site_compound_specificity") or "Unknown"] += 1
lineage = {"structure_lineage":22916,"pubchem_assay_lineage":1334718,"literature_experiment_proxy":210003,"multi_modality_pairs":324950}
data["M6"] = {
    "positive_pairs": 1502456,
    "source_tier": {k:dict(v) for k,v in sorted(source_tier.items(), key=lambda kv:-sum(kv[1].values()))[:10]},
    "lineage": lineage,
    "site_total": site_total,
    "site_type": site_type.most_common(8),
    "site_specificity": dict(site_specificity),
}

# M7 negative evidence, conflicts and release QA.
nval = json.loads(Path(r"D:\finale\05_质量控制\v62_completion_qa\NEGATIVE_EVIDENCE_RESOLUTION_V62_VALIDATION.json").read_text(encoding="utf-8"))
review_reasons = Counter()
for r in rows(Path(r"D:\finale\03_V6.3.1_readiness_20260803\02_review_queues\negative_unmapped_review_summary_v631.tsv")):
    reason = r.get("unmapped_reason_v62") or "unknown"
    if reason.startswith("Can't kekulize"): reason = "canonicalization / kekulization"
    elif reason.startswith("full_inchikey_mismatch"): reason = "full InChIKey mismatch"
    review_reasons[reason] += int(r.get("record_count") or 0)
data["M7"] = {
    "negative_table": nval["negative_table_counts"],
    "identity_resolution": nval["identity_resolution_counts"],
    "review_reasons": review_reasons.most_common(7),
    "qa": nval["quality_checks"],
}

# M8 docking prioritization and preparation readiness.
tiers = Counter(); structures = Counter(); receptor = Counter(); ligand = Counter(); boxes = Counter(); targets=set()
for r in rows(ROOT / "07_docking" / "docking_priority_v6_3_candidate.tsv.gz"):
    tiers[r.get("docking_tier") or "Unknown"] += 1
    structures[r.get("structure_grade") or "Unknown"] += 1
    receptor[r.get("receptor_preparation_status_v2") or "not_assessed"] += 1
    ligand[r.get("ligand_preparation_status_v2") or "not_assessed"] += 1
    boxes[r.get("box_definition_status_v2") or "not_assessed"] += 1
    targets.add(r.get("target_uniprot_id"))
data["M8"] = {
    "pairs": sum(tiers.values()), "targets": len(targets), "tiers": dict(tiers),
    "structure_grade": structures.most_common(), "receptor_prep": receptor.most_common(),
    "ligand_prep": ligand.most_common(), "box_status": boxes.most_common(),
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
print(OUT)
