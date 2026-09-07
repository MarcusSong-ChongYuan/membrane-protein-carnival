import csv, hashlib, json, shutil
from datetime import datetime, timezone
from pathlib import Path

WORK=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
DEST=Path(r"D:\finale\10_MemPro_V7_data_freeze_20260813")
DATA=DEST/"01_data"
QA=DEST/"02_QA"
DATA.mkdir(parents=True,exist_ok=True); QA.mkdir(parents=True,exist_ok=True)

protein=WORK/"02_protein_audit"/"V7_PROTEIN_FREEZE_CANDIDATE"
relations=WORK/"03_v7_relation_rebuild"
files=[
 protein/"human_membrane_protein_master_V7.tsv.gz",
 protein/"protein_excluded_or_other_entity_audit_V7.tsv.gz",
 protein/"protein_identity_adjudication_all_V7.tsv.gz",
 relations/"public_binding_evidence_master_V7.tsv.gz",
 relations/"public_protein_compound_pairs_V7.tsv.gz",
 relations/"binding_site_instances_V7.tsv.gz",
 relations/"canonical_protein_entity_V7.tsv.gz",
 relations/"protein_isoform_V7.tsv.gz",
 relations/"protein_tissue_cell_ihc_expression_V7.tsv.gz",
 relations/"protein_disease_relation_V7.tsv.gz",
 relations/"disease_entity_master_V7.tsv.gz",
 relations/"disease_anatomical_system_V7.tsv.gz",
 relations/"disease_therapeutic_area_V7.tsv.gz",
 relations/"complex_target_master_V7.tsv.gz",
 relations/"complex_target_components_V7.tsv.gz",
 relations/"small_molecule_master_V7.tsv.gz",
]
for src in files:
    if not src.exists(): raise FileNotFoundError(src)
    shutil.copy2(src,DATA/src.name)

qa_files=[
 protein/"V7_PROTEIN_FREEZE_REPORT.json",
 WORK/"02_protein_audit"/"final_validation_v7_3_1"/"V7_3_FINAL_VALIDATION_SUMMARY.json",
 WORK/"02_protein_audit"/"final_validation_v7_3_1"/"V7_3_FINAL_VALIDATION_BLOCKING.tsv",
 WORK/"02_protein_audit"/"final_validation_v7_3_1"/"V7_3_FINAL_VALIDATION_WARNINGS.tsv",
 WORK/"02_protein_audit"/"final_adjudication_v7_3_1"/"V7_3_1_BOUNDARY_CORRECTION.json",
 relations/"V7_RELATION_REBUILD_REPORT.json",
 relations/"V7_RELATIONAL_INTEGRITY_AUDIT.json",
]
for src in qa_files: shutil.copy2(src,QA/src.name)

manifest=[]
for p in sorted(DEST.rglob("*")):
    if not p.is_file() or p.name in {"MANIFEST.tsv","FREEZE_V7_DATA.json"}:continue
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    manifest.append({"relative_path":str(p.relative_to(DEST)),"bytes":p.stat().st_size,"sha256":h.hexdigest()})
with (DEST/"MANIFEST.tsv").open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["relative_path","bytes","sha256"],delimiter="\t");w.writeheader();w.writerows(manifest)

report={
 "release":"MemPro V7 data-layer freeze",
 "frozen_at_utc":datetime.now(timezone.utc).isoformat(),
 "protein_records_ABC":7800,"protein_A":5680,"protein_B":613,"protein_C":1507,
 "excluded_or_other_entity_records":3197,
 "public_binding_evidence_rows":942455,"public_unique_protein_compound_pairs":529168,
 "public_compounds":240055,"protein_disease_relations":6753,"canonical_diseases":3739,
 "binding_site_instances":55610,"expression_records":438596,"isoform_records":9240,
 "complexes":4585,"complex_components":14827,
 "protein_blocking_errors":0,"protein_warnings":0,"relational_blocking_errors":0,
 "status":"FROZEN_DATA_LAYER_PASS",
 "scope_note":"Data layer only. Figures and manuscript must be regenerated from this freeze; docking remains a parallel downstream analysis.",
 "manifest_entries":len(manifest)
}
(DEST/"FREEZE_V7_DATA.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(report,ensure_ascii=False,indent=2))
