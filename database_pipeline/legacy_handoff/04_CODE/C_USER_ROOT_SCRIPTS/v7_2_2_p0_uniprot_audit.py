#!/usr/bin/env python3
import csv, json, time, urllib.parse, urllib.request
from collections import Counter
from pathlib import Path

SRC=Path(r"C:\Users\Administrator\MemPro_V7_2_1_membrane_identity_audit_20260813\membrane_identity_manual_review_queue_v7_2_1.tsv")
OUT=Path(r"C:\Users\Administrator\MemPro_V7_2_2_membrane_identity_audit_20260813")
OUT.mkdir(exist_ok=True)

with SRC.open(encoding="utf-8",newline="") as f:
    rr=csv.reader(f,delimiter="\t"); hdr=next(rr); rows=list(rr)
idx={x:i for i,x in enumerate(hdr)}; pri=len(hdr)-1
p0=[r for r in rows if r[pri].startswith("P0")]

def get_json(acc):
    url="https://rest.uniprot.org/uniprotkb/"+urllib.parse.quote(acc)+".json"
    req=urllib.request.Request(url,headers={"User-Agent":"MemPro-V7.2.2-audit/1.0"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req,timeout=45) as z: return json.load(z),"PASS"
        except Exception as e:
            if attempt==4:return None,type(e).__name__+":"+str(e)[:160]
            time.sleep(2**attempt)

def locations(entry):
    out=[]
    for c in entry.get("comments",[]):
        if c.get("commentType")!="SUBCELLULAR LOCATION":continue
        for s in c.get("subcellularLocations",[]):
            loc=(s.get("location") or {}).get("value","")
            topo=(s.get("topology") or {}).get("value","")
            ori=(s.get("orientation") or {}).get("value","")
            out.append(" | ".join(x for x in (loc,topo,ori) if x))
    return "; ".join(out)

def classify(entry):
    feats=entry.get("features",[])
    trans=[x for x in feats if x.get("type")=="Transmembrane"]
    intra=[x for x in feats if x.get("type")=="Intramembrane"]
    lipid=[x for x in feats if x.get("type")=="Lipidation"]
    loc=locations(entry).lower()
    evidence=[]
    if trans:evidence.append("UniProt_transmembrane_feature")
    if intra:evidence.append("UniProt_intramembrane_feature")
    if lipid:evidence.append("UniProt_lipidation_feature")
    if trans or intra:return "A","CONFIRMED_OR_RECLASSIFIED",evidence
    if lipid:return "B","CONFIRMED_RECLASSIFY_A_TO_B",evidence
    if "peripheral membrane protein" in loc or "membrane-associated" in loc:
        return "C","CONFIRMED_RECLASSIFY_A_TO_C",["UniProt_peripheral_membrane_annotation"]
    # Mere organelle residence, secretion, or OPM orientation is not membrane identity.
    if any(x in loc for x in ("membrane","cell surface")):
        return "CANDIDATE","MANUAL_REVIEW_LOCATION_WITHOUT_MEMBRANE_MODE",["UniProt_location_without_topology_or_anchor"]
    return "EXCLUDE_CANDIDATE","LIKELY_NONMEMBRANE_NEEDS_FINAL_REVIEW",["No_current_UniProt_membrane_topology_anchor_or_peripheral_mode"]

fields=["target_uniprot","approved_symbol","old_class","old_evidence_level","old_tm_count","old_basis","old_location","uniprot_fetch_status","uniprot_reviewed","current_uniprot_location","current_transmembrane_count","current_intramembrane_count","current_lipidation_count","proposed_membrane_class","proposed_audit_decision","decision_evidence","manual_final_required"]
out=[]
for n,r in enumerate(p0,1):
    acc=r[idx["target_uniprot"]]; entry,status=get_json(acc)
    if entry:
        feats=entry.get("features",[]); cls,dec,ev=classify(entry)
        reviewed="reviewed" if entry.get("entryType","").lower().startswith("uniprotkb reviewed") else entry.get("entryType","")
        loc=locations(entry)
        counts={t:sum(x.get("type")==t for x in feats) for t in ("Transmembrane","Intramembrane","Lipidation")}
    else:
        cls,dec,ev,reviewed,loc="UNRESOLVED","FETCH_FAILED",[],"","";counts={t:0 for t in ("Transmembrane","Intramembrane","Lipidation")}
    out.append({"target_uniprot":acc,"approved_symbol":r[idx["approved_symbol"]],"old_class":r[idx["membrane_class_v52"]],"old_evidence_level":r[idx["evidence_level_v52"]],"old_tm_count":r[idx["transmembrane_count_v5"]],"old_basis":r[idx["membrane_evidence_basis"]],"old_location":r[idx["subcellular_location"]],"uniprot_fetch_status":status,"uniprot_reviewed":reviewed,"current_uniprot_location":loc,"current_transmembrane_count":counts["Transmembrane"],"current_intramembrane_count":counts["Intramembrane"],"current_lipidation_count":counts["Lipidation"],"proposed_membrane_class":cls,"proposed_audit_decision":dec,"decision_evidence":";".join(ev),"manual_final_required":1 if cls in {"CANDIDATE","EXCLUDE_CANDIDATE","UNRESOLVED"} else 0})
    print(f"{n}/{len(p0)} {acc} {cls} {dec}",flush=True)

dest=OUT/"P0_membrane_identity_uniprot_audit_v7_2_2.tsv"
with dest.open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(out)
summary={"input_p0":len(p0),"fetch_status":dict(Counter(x["uniprot_fetch_status"] for x in out)),"proposed_class":dict(Counter(x["proposed_membrane_class"] for x in out)),"decision":dict(Counter(x["proposed_audit_decision"] for x in out)),"manual_final_required":sum(x["manual_final_required"] for x in out),"rule":"Automated proposal only; CANDIDATE/EXCLUDE_CANDIDATE require source/structure/manual review before release change.","source":"UniProt REST fetched 2026-08-13"}
(OUT/"P0_AUDIT_SUMMARY_v7_2_2.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(summary,ensure_ascii=False,indent=2))
