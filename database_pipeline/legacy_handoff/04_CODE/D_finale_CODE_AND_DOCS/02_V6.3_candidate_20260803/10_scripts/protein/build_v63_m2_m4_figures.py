#!/usr/bin/env python3
"""Generate publication-oriented M2 and M4 panels for the V6.3 candidate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


C = {"navy":"#153B55","blue":"#3B73A3","teal":"#2A9D8F","orange":"#F4A261","red":"#C95B5B","purple":"#7A5195","gray":"#AEB7C1","light":"#E8EDF1","dark":"#5F6B75","green":"#6A994E"}


def arguments():
    p=argparse.ArgumentParser()
    p.add_argument("--protein-master",type=Path,required=True)
    p.add_argument("--protein-summary",type=Path,required=True)
    p.add_argument("--protein-review",type=Path,required=True)
    p.add_argument("--disease-source",type=Path,required=True)
    p.add_argument("--disease-relations",type=Path,required=True)
    p.add_argument("--therapeutic-area",type=Path,required=True)
    p.add_argument("--anatomical-system",type=Path,required=True)
    p.add_argument("--protein-figure-dir",type=Path,required=True)
    p.add_argument("--disease-figure-dir",type=Path,required=True)
    p.add_argument("--qa-dir",type=Path,required=True)
    return p.parse_args()


def setup():
    mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.titlesize":9,"axes.titleweight":"bold","axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,"figure.facecolor":"white","savefig.facecolor":"white","pdf.fonttype":42,"svg.fonttype":"none"})


def panel(ax, letter, title):
    ax.text(-.10,1.06,letter,transform=ax.transAxes,fontsize=11,fontweight="bold",color=C["navy"])
    ax.set_title(title,loc="left",pad=7)


def save(fig, base: Path, stem: str):
    paths={}
    for ext in ("png","svg","pdf"):
        d=base/ext; d.mkdir(parents=True,exist_ok=True)
        p=d/f"{stem}.{ext}"; fig.savefig(p,dpi=450 if ext=="png" else None,bbox_inches="tight"); paths[ext]=str(p)
    plt.close(fig); return paths


def make_m2(a):
    master=pd.read_csv(a.protein_master,sep="\t",dtype=str,keep_default_na=False)
    summary=pd.read_csv(a.protein_summary,sep="\t",dtype=str,keep_default_na=False)
    review=pd.read_csv(a.protein_review,sep="\t",dtype=str,keep_default_na=False)
    merged=master[["target_uniprot_id","classification_status_v5","evidence_level_v52"]].merge(summary,on="target_uniprot_id",validate="one_to_one")
    roles=summary.membrane_role_primary_v63.value_counts()
    coverage={
      "Structural family":(~summary.structural_family_primary_v63.isin(["","unclassified"])).mean(),
      "Molecular function":(~summary.molecular_function_primary_v63.isin(["","unclassified"])).mean(),
      "Biological process":(~summary.biological_process_primary_v63.isin(["","unclassified"])).mean(),
      "Specific membrane role":(~summary.membrane_role_primary_v63.isin(["","unclassified","family_defined_membrane_role_unresolved"])).mean(),
      "Specialist hierarchy":(~summary.specialist_classification_primary_v63.isin(["","no_specialist_classification"])).mean(),
      "Reactome pathway":summary.reactome_top_level_v63.ne("").mean(),
    }
    fig=plt.figure(figsize=(11.6,7.6)); gs=fig.add_gridspec(2,2,hspace=.42,wspace=.32)
    axa,axb,axc,axd=[fig.add_subplot(gs[i]) for i in range(4)]
    panel(axa,"A","Legacy single-label status retained for provenance")
    legacy=master.classification_status_v5.replace({"family_level_classified_function_unresolved":"Family defined; function unresolved","unclassified":"Unclassified"})
    legacy=legacy.where(legacy.isin(["Family defined; function unresolved","Unclassified"]),"Function/specialist classified").value_counts()
    labs=["Function/specialist classified","Family defined; function unresolved","Unclassified"]
    vals=[legacy.get(x,0) for x in labs]
    axa.barh(np.arange(3),vals,color=[C["teal"],C["orange"],C["gray"]]); axa.set_yticks(np.arange(3),labs); axa.invert_yaxis(); axa.set_xlabel("Proteins")
    for i,v in enumerate(vals): axa.text(v,i,f" {v:,} ({v/len(master):.1%})",va="center",fontsize=7)
    axa.grid(axis="x",color=C["light"],lw=.6)

    panel(axb,"B","Primary membrane role in the five-axis display")
    role_labels={"membrane_associated_enzyme":"Enzyme","family_defined_membrane_role_unresolved":"Family-defined; role unresolved","receptor":"Receptor","transporter":"Transporter","membrane_scaffold_or_linker":"Scaffold/linker","signaling_regulator":"Signaling regulator","ion_channel":"Ion channel","adhesion_recognition":"Adhesion/recognition","immune_or_cell_recognition":"Immune/cell recognition","membrane_trafficking":"Membrane trafficking","membrane_organizer":"Membrane organizer","junction_or_adhesion":"Junction/adhesion"}
    top=roles.sort_values().tail(12)
    axb.barh(range(len(top)),top.values,color=[C["orange"] if x=="family_defined_membrane_role_unresolved" else C["blue"] for x in top.index])
    axb.set_yticks(range(len(top)),[role_labels.get(x,x) for x in top.index]); axb.set_xlabel("Proteins"); axb.grid(axis="x",color=C["light"],lw=.6)
    for i,v in enumerate(top.values): axb.text(v,i,f" {v:,}",va="center",fontsize=6.8)

    panel(axc,"C","Annotation coverage by independent classification axis")
    names=list(coverage); vals=np.array(list(coverage.values()))*100
    axc.barh(range(len(names)),vals,color=[C["teal"],C["teal"],C["teal"],C["blue"],C["purple"],C["green"]]); axc.set_yticks(range(len(names)),names); axc.invert_yaxis(); axc.set_xlim(0,105); axc.set_xlabel("Proteins with an annotation (%)"); axc.grid(axis="x",color=C["light"],lw=.6)
    for i,v in enumerate(vals): axc.text(v,i,f" {v:.1f}%",va="center",fontsize=7)

    panel(axd,"D","Priority E1/E2 reclassification set (n=2,426)")
    tab=pd.crosstab(review.membrane_evidence_level,review.automatic_reclassification_status).reindex(["E1","E2"]).fillna(0)
    specific=tab.get("specific_role_proposed",pd.Series(0,index=tab.index)).astype(int)
    context=tab.get("multiaxis_context_added_role_unresolved",pd.Series(0,index=tab.index)).astype(int)
    x=np.arange(2); axd.bar(x,specific,color=C["teal"],label="Specific role proposed"); axd.bar(x,context,bottom=specific,color=C["orange"],label="Context added; role unresolved")
    axd.set_xticks(x,["E1","E2"]); axd.set_ylabel("Proteins"); axd.grid(axis="y",color=C["light"],lw=.6); axd.legend(frameon=False,fontsize=7,loc="upper right")
    for i in range(2): axd.text(i,(specific.iloc[i]+context.iloc[i])+10,f"{specific.iloc[i]+context.iloc[i]:,}",ha="center",fontsize=7)
    axd.text(.02,.04,"Only proteins already in the E1/E2 membrane-evidence layer are in this priority set.\nNo family-only record is forced into a specific function.",transform=axd.transAxes,fontsize=6.7,color=C["dark"],va="bottom")
    fig.suptitle("M2 — Orthogonal functional classification of the human membrane-protein catalog",x=.06,ha="left",fontsize=12,fontweight="bold",color=C["navy"])
    fig.text(.06,.93,"GO ontology ancestry, UniProt families, InterPro/Pfam, Reactome and specialist hierarchies are retained as separate axes (n=10,997).",ha="left",fontsize=7.5,color=C["dark"])
    fig.subplots_adjust(left=.13,right=.98,top=.86,bottom=.08)
    paths=save(fig,a.protein_figure_dir,"M2_five_axis_protein_classification_v63")
    data=pd.DataFrame({"axis":names,"annotated_fraction":list(coverage.values())}); data.to_csv(a.protein_figure_dir/"M2_annotation_coverage.tsv",sep="\t",index=False)
    caption="""# M2 caption\n\n**M2. Orthogonal classification of the 10,997-protein membrane catalog.** The V5 legacy label is retained only for provenance (A). V6.3 does not replace family identity with function: each protein receives separate structural-family, molecular-function, biological-process, membrane-role and specialist-classification fields (B–C), with Reactome recorded as an additional pathway layer. Panel D shows the predeclared 2,426-protein E1/E2 priority set. A specific role is proposed only when a specialist annotation, GO ancestry, a prior functional class or a conservative family rule supports it; otherwise the protein remains explicitly family-defined with unresolved membrane role. Counts describe annotation coverage, not experimental proof of function.\n"""
    (a.protein_figure_dir/"CAPTION_M2_FIVE_AXIS_CLASSIFICATION.md").write_text(caption,encoding="utf-8")
    return paths,{"proteins":len(master),"priority_set":len(review),"specific_role_proposed":int((review.automatic_reclassification_status=="specific_role_proposed").sum()),"coverage":coverage}


def make_m4(a):
    source=pd.read_csv(a.disease_source,sep="\t",dtype=str,keep_default_na=False)
    rel=pd.read_csv(a.disease_relations,sep="\t",dtype=str,keep_default_na=False)
    ta=pd.read_csv(a.therapeutic_area,sep="\t",dtype=str,keep_default_na=False)
    anatomy=pd.read_csv(a.anatomical_system,sep="\t",dtype=str,keep_default_na=False)
    default=rel[rel.default_release_inclusion_v63.astype(str).eq("1")]
    fig=plt.figure(figsize=(11.6,7.8)); gs=fig.add_gridspec(2,2,hspace=.48,wspace=.36)
    axa,axb,axc,axd=[fig.add_subplot(gs[i]) for i in range(4)]
    panel(axa,"A","Exact-only disease identity resolution")
    status=source.canonicalization_status.value_counts()
    labels={"unique_official_exact_mapping":"Official 1:1 exact/equivalent","direct_mondo_id":"Direct MONDO ID","non_disease_mondo_entity_review":"Non-disease MONDO review","unmapped_source_only":"Source-only (unmapped)","unique_mapping_to_obsolete_mondo_review":"Obsolete MONDO review"}
    ordered=status.sort_values()
    colors=[C["teal"] if x in {"unique_official_exact_mapping","direct_mondo_id"} else C["orange"] for x in ordered.index]
    axa.barh(range(len(ordered)),ordered.values,color=colors); axa.set_yticks(range(len(ordered)),[labels.get(x,x) for x in ordered.index]); axa.set_xlabel("Source disease entities"); axa.grid(axis="x",color=C["light"],lw=.6)
    for i,v in enumerate(ordered.values): axa.text(v,i,f" {v:,}",va="center",fontsize=7)

    panel(axb,"B","Open Targets 26.06 therapeutic areas (multi-label)")
    pair_ta=default[["target_uniprot_id","canonical_disease_id"]].merge(ta[["canonical_disease_id","therapeutic_area_name"]].drop_duplicates(),on="canonical_disease_id",how="inner")
    ta_counts=pair_ta.drop_duplicates().groupby("therapeutic_area_name").size().sort_values().tail(12)
    axb.barh(range(len(ta_counts)),ta_counts.values,color=C["purple"]); axb.set_yticks(range(len(ta_counts)),ta_counts.index); axb.set_xlabel("Unique protein–disease relations"); axb.grid(axis="x",color=C["light"],lw=.6)
    for i,v in enumerate(ta_counts.values): axb.text(v,i,f" {v:,}",va="center",fontsize=6.7)

    panel(axc,"C","Disease–anatomy mapping provenance")
    ad=anatomy[["canonical_disease_id","assertion_status"]].drop_duplicates().assertion_status.value_counts()
    ad=ad.reindex([x for x in ["asserted","inferred_from_mondo_ancestor","inferred_from_doid_ancestor"] if x in ad.index])
    axc.bar(range(len(ad)),ad.values,color=[C["teal"],C["blue"],C["green"]][:len(ad)]); axc.set_xticks(range(len(ad)),[x.replace("inferred_from_","").replace("_ancestor","") for x in ad.index],rotation=18,ha="right"); axc.set_ylabel("Unique disease–system assignments"); axc.grid(axis="y",color=C["light"],lw=.6)
    for i,v in enumerate(ad.values): axc.text(i,v,f"{v:,}",ha="center",va="bottom",fontsize=7)

    panel(axd,"D","Ontology-derived anatomical systems (multi-label)")
    pair_an=default[["target_uniprot_id","canonical_disease_id"]].merge(anatomy[["canonical_disease_id","anatomical_system_name"]].drop_duplicates(),on="canonical_disease_id",how="inner")
    an_counts=pair_an.drop_duplicates().groupby("anatomical_system_name").size().sort_values().tail(12)
    axd.barh(range(len(an_counts)),an_counts.values,color=C["blue"]); axd.set_yticks(range(len(an_counts)),an_counts.index); axd.set_xlabel("Unique protein–disease relations"); axd.grid(axis="x",color=C["light"],lw=.6)
    for i,v in enumerate(an_counts.values): axd.text(v,i,f" {v:,}",va="center",fontsize=6.7)
    fig.suptitle("M4 — Canonical disease identities and ontology-derived clinical context",x=.06,ha="left",fontsize=12,fontweight="bold",color=C["navy"])
    fig.text(.06,.93,"MONDO exact/equivalent mappings define identity; Open Targets therapeutic areas and DO/MONDO–Uberon anatomy are multi-label annotations, never keyword bins.",ha="left",fontsize=7.5,color=C["dark"])
    fig.subplots_adjust(left=.17,right=.98,top=.86,bottom=.10)
    paths=save(fig,a.disease_figure_dir,"M4_disease_ontology_context_v63")
    ta_counts.rename("protein_disease_relation_count").to_csv(a.disease_figure_dir/"M4_therapeutic_area_counts.tsv",sep="\t")
    an_counts.rename("protein_disease_relation_count").to_csv(a.disease_figure_dir/"M4_anatomical_system_counts.tsv",sep="\t")
    caption="""# M4 caption\n\n**M4. Canonical disease identities and ontology-derived context.** Source disease identifiers are merged only through direct MONDO identifiers or official one-to-one exact/equivalent mappings; non-disease MONDO nodes, obsolete targets and unmapped records remain review/source-only entities (A). Open Targets Platform 26.06 therapeutic-area assignments are shown as multi-label memberships (B). Anatomical associations derive from asserted or logically inherited DO/MONDO axioms pointing to Uberon, with inference provenance retained (C–D). A protein–disease relation can contribute to more than one therapeutic area or anatomical system; bars therefore must not be summed as mutually exclusive totals. Only V6.3 default relations are plotted.\n"""
    (a.disease_figure_dir/"CAPTION_M4_DISEASE_ONTOLOGY_CONTEXT.md").write_text(caption,encoding="utf-8")
    return paths,{"source_diseases":len(source),"canonical_relations":len(rel),"default_relations":len(default),"therapeutic_area_memberships":len(pair_ta.drop_duplicates()),"anatomical_memberships":len(pair_an.drop_duplicates())}


def main():
    a=arguments(); setup(); a.protein_figure_dir.mkdir(parents=True,exist_ok=True); a.disease_figure_dir.mkdir(parents=True,exist_ok=True); a.qa_dir.mkdir(parents=True,exist_ok=True)
    m2p,m2c=make_m2(a); m4p,m4c=make_m4(a)
    validation={"status":"PASS","figures":{"M2":m2p,"M4":m4p},"counts":{"M2":m2c,"M4":m4c}}
    (a.qa_dir/"V63_M2_M4_FIGURE_VALIDATION.json").write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(validation,ensure_ascii=False,indent=2))


if __name__=="__main__": main()
