import csv, gzip, hashlib, json, shutil
from collections import Counter
from pathlib import Path

ROOT = Path(r"D:\finale\09_V7_data_freeze_working_20260813")
PROTEIN = ROOT / "02_protein_audit" / "V7_PROTEIN_FREEZE_CANDIDATE" / "human_membrane_protein_master_V7.tsv"
V64 = Path(r"D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay")
V641 = Path(r"D:\finale\09_V6.4.1_final_publication_triage_20260812")
COMPOUND = Path(r"D:\finale\02_V6.3_candidate_20260803\01_baseline_v6_2\small_molecule_master_v1_3.tsv")
OUT = ROOT / "03_v7_relation_rebuild"
OUT.mkdir(parents=True, exist_ok=True)

def op(path, mode="rt"):
    return gzip.open(path, mode, encoding=None if "b" in mode else "utf-8", newline="" if "t" in mode else None) if str(path).endswith(".gz") else path.open(mode, encoding=None if "b" in mode else "utf-8", newline="" if "t" in mode else None)

def filter_table(src, dst, predicate, collect=None):
    n_in=n_out=0
    with op(src) as fi, op(dst,"wt") as fo:
        r=csv.DictReader(fi,delimiter="\t"); w=csv.DictWriter(fo,fieldnames=r.fieldnames,delimiter="\t"); w.writeheader()
        for row in r:
            n_in+=1
            if predicate(row):
                w.writerow(row); n_out+=1
                if collect: collect(row)
    return {"input":n_in,"output":n_out}

with PROTEIN.open(encoding="utf-8",newline="") as f:
    proteins={r["target_uniprot"] for r in csv.DictReader(f,delimiter="\t")}

report={"protein_whitelist":len(proteins),"tables":{}}
evidence_ids=set(); compounds=set(); pairs=set(); tiers=Counter(); sources=Counter()
def collect_ev(r):
    evidence_ids.add(r["evidence_id"]); compounds.add(r["compound_internal_id"]); pairs.add((r["target_uniprot_id"],r["compound_internal_id"])); tiers[r["evidence_tier"]]+=1; sources[r["source_database"]]+=1
report["tables"]["public_binding_evidence_V7"] = filter_table(
    V641/"public_binding_evidence_master_v6_4_1.tsv.gz", OUT/"public_binding_evidence_master_V7.tsv.gz",
    lambda r:r["target_uniprot_id"] in proteins, collect_ev)

report["tables"]["public_pairs_V7"] = filter_table(
    V641/"public_protein_compound_pairs_v6_4_1.tsv.gz", OUT/"public_protein_compound_pairs_V7.tsv.gz",
    lambda r:(r["target_uniprot_id"],r["compound_internal_id"]) in pairs)
report["tables"]["binding_sites_V7"] = filter_table(
    V64/"binding_site_instances_v6_4_candidate.tsv.gz", OUT/"binding_site_instances_V7.tsv.gz",
    lambda r:r["evidence_id"] in evidence_ids and r["target_uniprot_id"] in proteins)
report["tables"]["canonical_protein_V7"] = filter_table(
    V64/"canonical_protein_entity_v6_4_candidate.tsv.gz", OUT/"canonical_protein_entity_V7.tsv.gz",
    lambda r:r["canonical_uniprot_accession"] in proteins)
report["tables"]["protein_isoform_V7"] = filter_table(
    V64/"protein_isoform_v6_4_candidate.tsv.gz", OUT/"protein_isoform_V7.tsv.gz",
    lambda r:r["canonical_uniprot_accession"] in proteins)
report["tables"]["expression_V7"] = filter_table(
    V64/"protein_tissue_cell_ihc_expression_v6_4_candidate.tsv.gz", OUT/"protein_tissue_cell_ihc_expression_V7.tsv.gz",
    lambda r:r["target_uniprot_id"] in proteins)

disease_ids=set()
def collect_dis(r): disease_ids.add(r["canonical_disease_id"])
report["tables"]["protein_disease_V7"] = filter_table(
    V64/"protein_disease_relation_v6_4_candidate.tsv", OUT/"protein_disease_relation_V7.tsv.gz",
    lambda r:r["target_uniprot_id"] in proteins and r.get("default_release_inclusion_v63","1") in {"1","True","true"}, collect_dis)

for srcname,dstname,key in [
    ("disease_entity_master_v6_4_candidate.tsv","disease_entity_master_V7.tsv.gz","canonical_disease_id"),
    ("disease_anatomical_system_v6_4_candidate.tsv","disease_anatomical_system_V7.tsv.gz","canonical_disease_id"),
    ("disease_therapeutic_area_v6_4_candidate.tsv","disease_therapeutic_area_V7.tsv.gz","canonical_disease_id")]:
    report["tables"][dstname]=filter_table(V64/srcname,OUT/dstname,lambda r,k=key:r.get(k,"") in disease_ids)

# Select complexes that contain at least one final membrane protein, then retain
# every component of those complexes so stoichiometry is not destroyed.
selected_complexes=set()
with op(V64/"complex_target_components_v6_4_candidate.tsv.gz") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["component_uniprot_id"] in proteins: selected_complexes.add(r["complex_target_id"])
report["tables"]["complex_components_V7"] = filter_table(
    V64/"complex_target_components_v6_4_candidate.tsv.gz",OUT/"complex_target_components_V7.tsv.gz",
    lambda r:r["complex_target_id"] in selected_complexes)
report["tables"]["complex_master_V7"] = filter_table(
    V64/"complex_target_master_v6_4_candidate.tsv",OUT/"complex_target_master_V7.tsv.gz",
    lambda r:r["complex_target_id"] in selected_complexes and r.get("default_release_inclusion","1") in {"1","True","true"})

report["tables"]["small_molecule_master_V7"] = filter_table(
    COMPOUND, OUT/"small_molecule_master_V7.tsv.gz", lambda r:r["compound_internal_id"] in compounds)

report.update({"evidence_ids":len(evidence_ids),"unique_pairs":len(pairs),"unique_compounds":len(compounds),"disease_ids":len(disease_ids),"selected_complexes":len(selected_complexes),"evidence_tiers":dict(tiers),"evidence_sources":dict(sources)})

manifest=[]
for p in sorted(OUT.glob("*")):
    if not p.is_file() or p.name in {"V7_RELATION_REBUILD_REPORT.json","MANIFEST.tsv"}: continue
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    manifest.append({"filename":p.name,"bytes":p.stat().st_size,"sha256":h.hexdigest()})
with (OUT/"MANIFEST.tsv").open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["filename","bytes","sha256"],delimiter="\t");w.writeheader();w.writerows(manifest)
report["manifest_files"]=len(manifest)
(OUT/"V7_RELATION_REBUILD_REPORT.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(report,ensure_ascii=False,indent=2))
