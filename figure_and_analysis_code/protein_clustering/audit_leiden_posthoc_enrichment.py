#!/usr/bin/env python3
"""Independent, FDR-controlled interpretation audit for Leiden communities."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parent
INPUT=ROOT/"inputs"; OUT=ROOT/"leiden_structural_results"

def read(path,**kw): return pd.read_csv(path,sep="\t",low_memory=False,**kw)

def enrich_binary(meta: pd.DataFrame, subject: pd.DataFrame, category_col: str, output_axis: str) -> pd.DataFrame:
    """Community/category Fisher enrichment where every protein is one subject."""
    subjects=meta[["uniprot_id","macro_community"]].merge(subject[["uniprot_id",category_col]].drop_duplicates(),on="uniprot_id",how="left")
    subjects=subjects.dropna(subset=[category_col]).drop_duplicates(["uniprot_id",category_col])
    values=subjects[category_col].astype(str)
    bad={"","nan","none","unknown","unclassified","family_defined_membrane_role_unresolved","no_specialist_classification","no_molecular_function_classification","no_biological_process_classification"}
    values=values[~values.str.lower().isin(bad)]
    subjects=subjects.loc[values.index].copy(); subjects[category_col]=values
    n_total=len(meta); rows=[]
    for value, vv in subjects.groupby(category_col):
        labelled=set(vv.uniprot_id); total_with=len(labelled)
        if total_with<10 or total_with>n_total-10: continue
        for c, gg in meta.groupby("macro_community"):
            members=set(gg.uniprot_id); a=len(members & labelled); b=len(members)-a; cc=total_with-a; d=n_total-a-b-cc
            if a<5: continue
            _,p=fisher_exact([[a,b],[cc,d]],alternative="greater")
            # Haldane–Anscombe correction only for effect-size stability.
            log2or=float(np.log2(((a+.5)*(d+.5))/((b+.5)*(cc+.5))))
            rows.append({"community":int(c),"annotation_axis":output_axis,"annotation_label":str(value),"community_support":int(len(members)),"label_support":int(total_with),"observed":int(a),"expected":len(members)*total_with/n_total,"log2_odds_ratio":log2or,"fisher_p":float(p)})
    return pd.DataFrame(rows)

def main():
    meta=read(OUT/"protein_leiden_structural_umap3d.tsv",dtype={"uniprot_id":str})
    ann=read(INPUT/"protein_cross_classification_v72.tsv.gz",usecols=["target_uniprot_id","classification_axis","classification_label"],dtype=str).rename(columns={"target_uniprot_id":"uniprot_id"})
    # These axes were excluded from the Leiden graph. Structural family is not
    # tested because it was an input feature.
    pieces=[]
    for axis in ("molecular_function","biological_process","specialist_classification"):
        a=ann[ann.classification_axis.eq(axis)][["uniprot_id","classification_label"]].rename(columns={"classification_label":"label"})
        pieces.append(enrich_binary(meta,a.rename(columns={"label":"value"}),"value",axis))
    role=meta[["uniprot_id","posthoc_membrane_role"]].rename(columns={"posthoc_membrane_role":"value"})
    pieces.append(enrich_binary(meta,role,"value","posthoc_membrane_role"))
    all_enrich=pd.concat([x for x in pieces if not x.empty],ignore_index=True)
    all_enrich["bh_fdr_q"]=multipletests(all_enrich.fisher_p,method="fdr_bh")[1]
    all_enrich["within_community_fraction"]=all_enrich.observed/all_enrich.community_support
    all_enrich["passes_interpretation_threshold"]=(all_enrich.bh_fdr_q<=.05)&(all_enrich.log2_odds_ratio>=1)&(all_enrich.observed>=10)
    all_enrich.sort_values(["community","passes_interpretation_threshold","bh_fdr_q","log2_odds_ratio","observed"],ascending=[True,False,True,False,False]).to_csv(OUT/"leiden_posthoc_enrichment_all.tsv",sep="\t",index=False)

    # Ligand status enrichment is protein-level: a protein is counted once per
    # status regardless of its number of compound pairs.
    pairs=read(INPUT/"protein_compound_pair_v72.tsv.gz",usecols=["target_uniprot_id","compound_internal_id"],dtype=str).rename(columns={"target_uniprot_id":"uniprot_id"})
    compounds=read(INPUT/"compound_master_v72.tsv.gz",usecols=["compound_internal_id","is_approved_drug","is_clinical_candidate","is_endogenous_ligand","is_natural_product","is_chemical_probe"],dtype=str)
    pc=pairs.merge(compounds,on="compound_internal_id",how="left")
    lig=[]
    for col in ["is_approved_drug","is_clinical_candidate","is_endogenous_ligand","is_natural_product","is_chemical_probe"]:
        yes=pc.loc[pc[col].astype(str).str.lower().isin({"true","1","yes"}),"uniprot_id"].drop_duplicates()
        lig.append(pd.DataFrame({"uniprot_id":yes,"value":f"has_{col.removeprefix('is_')}"}))
    ligand=enrich_binary(meta,pd.concat(lig,ignore_index=True),"value","protein_level_ligand_status")
    ligand["bh_fdr_q"]=multipletests(ligand.fisher_p,method="fdr_bh")[1]
    ligand["passes_interpretation_threshold"]=(ligand.bh_fdr_q<=.05)&(ligand.log2_odds_ratio>=1)&(ligand.observed>=10)
    ligand.sort_values(["community","passes_interpretation_threshold","bh_fdr_q","log2_odds_ratio"],ascending=[True,False,True,False]).to_csv(OUT/"leiden_ligand_status_enrichment.tsv",sep="\t",index=False)

    # Candidate names require independent evidence.  Keep all annotations even
    # when no candidate passes; no community is forced to receive a name.
    sig=all_enrich[all_enrich.passes_interpretation_threshold].copy()
    # Generic ontology ancestors are evidence of broad enrichment but are not
    # acceptable community names.  They remain in the complete audit table.
    generic_name_terms={
        "binding_activity", "catalytic_activity", "structural_molecule_activity",
        "organelle_organization", "protein_metabolic_process", "signal_transduction",
        "transport_and_localization", "metabolic_process", "cellular_process",
    }
    sig["label_specificity_eligible"]=(~sig.annotation_label.str.lower().isin(generic_name_terms)) & (sig.within_community_fraction>=0.10)
    sig["axis_priority"]=sig.annotation_axis.map({
        "specialist_classification":0, "posthoc_membrane_role":1,
        "molecular_function":2, "biological_process":3,
    }).fillna(9)
    labels=[]
    for c,g in sig.groupby("community"):
        g=g[g.label_specificity_eligible].sort_values(["axis_priority","bh_fdr_q","log2_odds_ratio","observed"],ascending=[True,True,False,False])
        if g.empty:
            continue
        top=g.iloc[0]
        labels.append({"community":int(c),"candidate_annotation_axis":top.annotation_axis,"candidate_label":top.annotation_label,"observed":int(top.observed),"within_community_fraction":float(top.within_community_fraction),"log2_odds_ratio":float(top.log2_odds_ratio),"bh_fdr_q":float(top.bh_fdr_q),"label_status":"CANDIDATE_REQUIRES_MANUAL_SCIENTIFIC_REVIEW"})
    label_df=pd.DataFrame(labels)
    all_communities=pd.DataFrame({"community":sorted(meta.macro_community.unique())})
    label_df=all_communities.merge(label_df,on="community",how="left")
    label_df["label_status"]=label_df.label_status.fillna("UNRESOLVED_NO_INDEPENDENT_ENRICHMENT")
    label_df.to_csv(OUT/"leiden_macro_module_candidate_labels.tsv",sep="\t",index=False)

    summary={"n_proteins":int(len(meta)),"n_macro_communities":int(meta.macro_community.nunique()),"input_excluded_from_interpretation_tests":["structural_family"],"tested_independent_axes":["molecular_function","biological_process","specialist_classification","posthoc_membrane_role","protein_level_ligand_status"],"annotation_tests":int(len(all_enrich)),"annotation_significant_after_BH":int(all_enrich.passes_interpretation_threshold.sum()),"ligand_tests":int(len(ligand)),"ligand_significant_after_BH":int(ligand.passes_interpretation_threshold.sum()),"communities_with_candidate_label":int((label_df.label_status=="CANDIDATE_REQUIRES_MANUAL_SCIENTIFIC_REVIEW").sum()),"interpretation_rule":"Fisher one-sided enrichment; BH across all tests in each output family; observed >=10; log2 OR >=1; q <=0.05. A community-name candidate additionally requires one specific held-out label in >=10% of that community. Candidate labels are not automatic final biological names."}
    (OUT/"leiden_posthoc_enrichment_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    (OUT/"leiden_posthoc_enrichment_summary.md").write_text("# Leiden post-hoc interpretation audit\n\n"+"\n".join(f"- **{k}**: {v}" for k,v in summary.items()),encoding="utf-8")
    print(json.dumps(summary,indent=2))

if __name__=="__main__": main()
