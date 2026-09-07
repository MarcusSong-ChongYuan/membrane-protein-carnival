import csv, gzip, hashlib, json
from pathlib import Path

ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
OUT=ROOT/"03_v7_relation_rebuild"
V64=Path(r"D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay")

def op(p,mode="rt"):
    return gzip.open(p,mode,encoding="utf-8",newline="") if str(p).endswith(".gz") else p.open(mode,encoding="utf-8",newline="")
def vals(p,col):
    with op(p) as f:return {r[col] for r in csv.DictReader(f,delimiter="\t")}
def count_bad(p,pred):
    n=bad=0
    with op(p) as f:
        for r in csv.DictReader(f,delimiter="\t"):
            n+=1; bad+=not pred(r)
    return n,bad
def filter_table(src,dst,pred):
    ni=no=0
    with op(src) as fi,op(dst,"wt") as fo:
        r=csv.DictReader(fi,delimiter="\t");w=csv.DictWriter(fo,fieldnames=r.fieldnames,delimiter="\t");w.writeheader()
        for row in r:
            ni+=1
            if pred(row):w.writerow(row);no+=1
    return ni,no

proteins=vals(OUT/"canonical_protein_entity_V7.tsv.gz","canonical_uniprot_accession")
evidence=vals(OUT/"public_binding_evidence_master_V7.tsv.gz","evidence_id")
compounds=vals(OUT/"small_molecule_master_V7.tsv.gz","compound_internal_id")
diseases=vals(OUT/"disease_entity_master_V7.tsv.gz","canonical_disease_id")

# Repair complex parent/child scope: select release-included parent complexes that
# contain at least one V7 membrane protein, then retain all their components.
candidate=set()
with op(V64/"complex_target_components_v6_4_candidate.tsv.gz") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["component_uniprot_id"] in proteins:candidate.add(r["complex_target_id"])
master_tmp=OUT/"complex_target_master_V7.tsv.gz"
filter_table(V64/"complex_target_master_v6_4_candidate.tsv",master_tmp,
             lambda r:r["complex_target_id"] in candidate and r.get("default_release_inclusion","1") in {"1","True","true"})
complexes=vals(master_tmp,"complex_target_id")
filter_table(V64/"complex_target_components_v6_4_candidate.tsv.gz",OUT/"complex_target_components_V7.tsv.gz",
             lambda r:r["complex_target_id"] in complexes)

checks={}
checks["evidence_target_fk"]=count_bad(OUT/"public_binding_evidence_master_V7.tsv.gz",lambda r:r["target_uniprot_id"] in proteins)
checks["evidence_compound_fk"]=count_bad(OUT/"public_binding_evidence_master_V7.tsv.gz",lambda r:r["compound_internal_id"] in compounds)
checks["pair_target_compound_fk"]=count_bad(OUT/"public_protein_compound_pairs_V7.tsv.gz",lambda r:r["target_uniprot_id"] in proteins and r["compound_internal_id"] in compounds)
checks["site_evidence_target_fk"]=count_bad(OUT/"binding_site_instances_V7.tsv.gz",lambda r:r["evidence_id"] in evidence and r["target_uniprot_id"] in proteins)
checks["isoform_canonical_fk"]=count_bad(OUT/"protein_isoform_V7.tsv.gz",lambda r:r["canonical_uniprot_accession"] in proteins)
checks["expression_target_fk"]=count_bad(OUT/"protein_tissue_cell_ihc_expression_V7.tsv.gz",lambda r:r["target_uniprot_id"] in proteins)
checks["disease_relation_fk"]=count_bad(OUT/"protein_disease_relation_V7.tsv.gz",lambda r:r["target_uniprot_id"] in proteins and r["canonical_disease_id"] in diseases)
checks["disease_anatomy_fk"]=count_bad(OUT/"disease_anatomical_system_V7.tsv.gz",lambda r:r["canonical_disease_id"] in diseases)
checks["disease_therapeutic_fk"]=count_bad(OUT/"disease_therapeutic_area_V7.tsv.gz",lambda r:r["canonical_disease_id"] in diseases)
checks["complex_component_parent_fk"]=count_bad(OUT/"complex_target_components_V7.tsv.gz",lambda r:r["complex_target_id"] in complexes)

# Exact public pair reconciliation against evidence-derived unique pair keys.
evpairs=set()
with op(OUT/"public_binding_evidence_master_V7.tsv.gz") as f:
    for r in csv.DictReader(f,delimiter="\t"):evpairs.add((r["target_uniprot_id"],r["compound_internal_id"]))
pairkeys=[]
with op(OUT/"public_protein_compound_pairs_V7.tsv.gz") as f:
    pairkeys=[(r["target_uniprot_id"],r["compound_internal_id"]) for r in csv.DictReader(f,delimiter="\t")]
pairset=set(pairkeys)

blocking={k:{"rows":v[0],"bad":v[1]} for k,v in checks.items()}
blocking["pair_exact_reconciliation"]={"rows":len(pairkeys),"bad":len(evpairs^pairset)+max(0,len(pairkeys)-len(pairset))}
blocking_count=sum(v["bad"] for v in blocking.values())
report={"status":"PASS_V7_RELATIONAL_INTEGRITY" if blocking_count==0 else "FAIL_V7_RELATIONAL_INTEGRITY","blocking_error_count":blocking_count,"checks":blocking,"counts":{"proteins":len(proteins),"evidence":len(evidence),"pairs":len(pairset),"compounds":len(compounds),"diseases":len(diseases),"complexes":len(complexes)}}
(OUT/"V7_RELATIONAL_INTEGRITY_AUDIT.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")

# Rebuild manifest after complex repair and QA.
manifest=[]
for p in sorted(OUT.iterdir()):
    if not p.is_file() or p.name=="MANIFEST.tsv":continue
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    manifest.append({"filename":p.name,"bytes":p.stat().st_size,"sha256":h.hexdigest()})
with (OUT/"MANIFEST.tsv").open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["filename","bytes","sha256"],delimiter="\t");w.writeheader();w.writerows(manifest)
print(json.dumps(report,ensure_ascii=False,indent=2))
