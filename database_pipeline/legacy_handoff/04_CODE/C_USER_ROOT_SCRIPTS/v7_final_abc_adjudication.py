#!/usr/bin/env python3
import csv, glob, gzip, json, re
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
SRC=ROOT/"02_protein_audit"/"protein_identity_audit_all_v7_cross_source_retriage.tsv"
MASTER=Path(r"D:\finale\08_V6.4_database_freeze_and_figures_20260812\01_release_overlay\human_membrane_protein_master_v6_4_candidate.tsv.gz")
OUT=ROOT/"02_protein_audit"/"final_adjudication_v7"
OUT.mkdir(exist_ok=True)

with gzip.open(MASTER,"rt",encoding="utf-8",newline="") as f:
    master={r["target_uniprot_id"]:r for r in csv.DictReader(f,delimiter="\t")}
entries={}
for p in glob.glob(str(ROOT/"02_protein_audit"/"uniprot_current_jsonl"/"batch_*.json")):
    d=json.loads(Path(p).read_text(encoding="utf-8"))
    for e in d.get("results",[]):entries[e["primaryAccession"]]=e
with SRC.open(encoding="utf-8",newline="") as f:
    rd=csv.DictReader(f,delimiter="\t");base=rd.fieldnames;rows=list(rd)

experimental={"IDA","IMP","IGI","IPI","IEP","EXP","HDA","HMP","HGI","HEP"}
lipid_terms=("phospholipid binding","phosphatidylinositol","phosphoinositide","lipid binding","sulfatide binding","membrane binding")
membrane_cc=("plasma membrane","cell membrane","endosome membrane","lysosomal membrane","late endosome membrane","early endosome membrane","golgi membrane","endoplasmic reticulum membrane","mitochondrial membrane","nuclear membrane","membrane raft","phagosome membrane","vesicle membrane")
direct_patterns=(r"binds? (?:to )?(?:phosphatidyl|phosphoinositide|phospholipid|lipid)",r"lipid-binding",r"membrane-binding",r"association with (?:the )?membrane",r"associat(?:es|ed|ion) with .*membrane")
recruit_patterns=(r"recruit(?:ed|ment|s) .* membrane",r"translocat(?:es|ed|ion) .* membrane",r"membrane association",r"localiz(?:es|ed|ation) to .* membrane",r"anchored .* membrane")

def go_data(e):
    terms=[]
    for x in e.get("uniProtKBCrossReferences",[]):
        if x.get("database")!="GO":continue
        props={p.get("key"):p.get("value","") for p in x.get("properties",[])}
        term=props.get("GoTerm","");ev=props.get("GoEvidenceType","")
        code=ev.split(":",1)[0]
        terms.append((term,ev,code in experimental))
    return terms

def comments(e):
    vals=[]
    for c in e.get("comments",[]):
        if c.get("commentType") in {"FUNCTION","SUBCELLULAR LOCATION","SUBUNIT"}:
            for t in c.get("texts",[]) or []:vals.append(t.get("value",""))
            note=c.get("note") or {}
            for t in note.get("texts",[]) or []:vals.append(t.get("value",""))
    return " ".join(vals)

def adjudicate(r):
    e=entries[r["target_uniprot"]];m=master[r["target_uniprot"]]
    cross=r["cross_source_proposed_class_v7"]
    if cross=="A":
        mode="transmembrane_or_intramembrane"
        if "PDBTM" in r["external_integral_supports_v7"] or "OPM_integral" in r["external_integral_supports_v7"]:mode="structure_supported_integral_membrane"
        elif "beta_barrel" in r["external_integral_supports_v7"].lower():mode="beta_barrel_integral_membrane"
        return "A",mode,"CONFIRMED_A",r["external_integral_supports_v7"] or r["uniprot_transmembrane_features_current"] or r["uniprot_intramembrane_features_current"]
    if cross=="B":
        return "B","covalent_lipid_anchor","CONFIRMED_B",r["uniprot_lipidation_features_current"]
    if cross=="C":
        mode="structure_supported_peripheral" if "OPM_peripheral_structure" in r["external_localization_supports_v7"] else "uniprot_peripheral_membrane_annotation"
        return "C",mode,"CONFIRMED_C",r["external_localization_supports_v7"] or r["uniprot_subcellular_locations_current"]
    go=go_data(e); text=comments(e).lower(); kws={x.get("name","").lower() for x in e.get("keywords",[])}
    go_lipid=[(t,ev) for t,ev,exp in go if exp and t.startswith("F:") and any(k in t.lower() for k in lipid_terms)]
    go_mem=[(t,ev) for t,ev,exp in go if exp and t.startswith("C:") and any(k in t.lower() for k in membrane_cc)]
    direct=go_lipid or "lipid-binding" in kws or any(re.search(p,text) for p in direct_patterns)
    recruited=any(re.search(p,text) for p in recruit_patterns)
    complex_membrane=bool(go_mem) and any(k in text for k in ("complex","recruit","scaffold","coat","retromer","ragulator","v-atpase","clathrin"))
    if direct:
        evidence=";".join(t+"|"+ev for t,ev in go_lipid) or "UniProt function/keyword direct lipid or membrane binding"
        return "C","direct_lipid_or_membrane_surface_binding","CONFIRMED_C_DIRECT",evidence
    if recruited:
        return "C","state_dependent_peripheral_recruitment","CONFIRMED_C_STATE_DEPENDENT","UniProt function/location note explicitly describes membrane recruitment/association"
    if complex_membrane:
        return "C","complex_mediated_membrane_association","CONFIRMED_C_COMPLEX_MEDIATED",";".join(t+"|"+ev for t,ev in go_mem)
    signal=int(r["uniprot_signal_count_current"] or 0)>0
    loc=r["uniprot_subcellular_locations_current"].lower()
    if signal and any(k in loc for k in ("secreted","extracellular","lysosome","lumen")):
        return "EXCLUDED","soluble_secreted_or_lumenal","EXCLUDED_NON_MEMBRANE","Signal peptide plus soluble secreted/lumenal localization; no transmembrane, lipid anchor, direct peripheral, or complex-mediated membrane mechanism"
    if go_mem:
        return "EXCLUDED","localization_without_binding_mechanism","EXCLUDED_NO_MEMBRANE_MECHANISM","Experimental membrane-location GO term exists, but no direct/anchored/complex-mediated attachment mechanism was found"
    return "EXCLUDED","no_membrane_binding_mechanism","EXCLUDED_NO_MEMBRANE_MECHANISM","No transmembrane/intramembrane segment, covalent lipid anchor, direct peripheral binding, or explicit complex-mediated membrane association in audited sources"

extra=["final_membrane_class_v7","final_membrane_mode_v7","final_decision_v7","final_decision_evidence_v7","final_public_inclusion_v7","final_review_status_v7"]
out=[]
for r in rows:
    cls,mode,dec,evidence=adjudicate(r)
    x=dict(r);x.update({"final_membrane_class_v7":cls,"final_membrane_mode_v7":mode,"final_decision_v7":dec,"final_decision_evidence_v7":evidence,"final_public_inclusion_v7":1 if cls in {"A","B","C"} else 0,"final_review_status_v7":"ADJUDICATED_NO_UNRESOLVED_STATUS"});out.append(x)
with (OUT/"protein_identity_final_adjudication_all_v7.tsv").open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=base+extra,delimiter="\t");w.writeheader();w.writerows(out)
included=[x for x in out if x["final_public_inclusion_v7"]==1];excluded=[x for x in out if x["final_public_inclusion_v7"]==0]
for name,vals in [("protein_ABC_release_candidate_v7.tsv",included),("protein_excluded_audit_v7.tsv",excluded)]:
    with (OUT/name).open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=base+extra,delimiter="\t");w.writeheader();w.writerows(vals)
summary={"input":len(out),"final_class":dict(Counter(x["final_membrane_class_v7"] for x in out)),"final_decision":dict(Counter(x["final_decision_v7"] for x in out)),"final_mode":dict(Counter(x["final_membrane_mode_v7"] for x in out)),"release_candidate_count":len(included),"excluded_count":len(excluded),"unresolved_count":sum(x["final_membrane_class_v7"] in {"","UNRESOLVED","unknown"} for x in out),"status":"DETERMINISTIC_ADJUDICATION_COMPLETE_PENDING_VALIDATION_SAMPLE"}
(OUT/"FINAL_ABC_ADJUDICATION_SUMMARY.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(summary,ensure_ascii=False,indent=2))
