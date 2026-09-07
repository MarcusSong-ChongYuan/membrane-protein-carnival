"""C11: Exact scaffold × membrane-role enrichment on canonical pair units."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compound_figure_pool"
FORMAL = Path(r"D:\finale\FORMAL")
PAIR = FORMAL / "01_core_v72" / "01_release_tables" / "protein_compound_pair_v72.tsv.gz"
ROLE = FORMAL / "03_publication_repairs" / "protein_membrane_role_FORMAL.tsv"
CATEGORY = FORMAL / "02_companion_v721" / "02_companion_tables" / "compound_chemical_classification_v721.tsv.gz"
MASTER = FORMAL / "01_core_v72" / "01_release_tables" / "compound_master_v72.tsv.gz"
QUARANTINE = FORMAL / "03_publication_repairs" / "COMPOUND_STRUCTURE_REPAIR_AND_QUARANTINE.tsv"
C3 = ROOT / "results" / "compound_figure_pool" / "C3_scaffold_rank_frequency.tsv"
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans","sans-serif"],"svg.fonttype":"none","pdf.fonttype":42,"figure.facecolor":"white","savefig.facecolor":"white","axes.linewidth":0.8})
INK,SLATE="#24364B","#5E7182"

def bh(p: np.ndarray) -> np.ndarray:
    p=np.asarray(p,float); n=len(p); idx=np.argsort(p); ranks=np.arange(1,n+1); q=np.empty(n); q[idx]=np.minimum.accumulate((p[idx]*n/ranks)[::-1])[::-1]; return np.clip(q,0,1)

def main() -> None:
    OUT.mkdir(parents=True,exist_ok=True)
    scaff=pd.read_csv(C3,sep="\t")
    # Predeclared minimum supports are based on the observed scaffold-count
    # distribution; every passing scaffold is tested, not only visually chosen scaffolds.
    # High-support display stratum. The initial 30/20 support screen yielded
    # 185 scaffold rows, which is not legible in a single publication figure.
    # These stricter support minima retain the high-frequency joint tail while
    # preserving every lower-support scaffold in the C3 source table.
    candidates=scaff.loc[(scaff.unique_compound_count>=150)&(scaff.unique_protein_count>=70)].copy()
    candidates["scaffold_short_id"]="S"+candidates.scaffold_rank.astype(str)
    candidates.to_csv(OUT/"C11_scaffold_filtering_candidates.tsv",sep="\t",index=False)
    pairs=pd.read_csv(PAIR,sep="\t",compression="gzip",usecols=["pair_id","target_uniprot_id","compound_internal_id"],low_memory=False)
    roles=pd.read_csv(ROLE,sep="\t",usecols=["canonical_uniprot_accession","formal_primary_membrane_role"]).drop_duplicates("canonical_uniprot_accession")
    master=pd.read_csv(MASTER,sep="\t",compression="gzip",usecols=["compound_internal_id","compound_scope_class"],low_memory=False)
    master=master.loc[master.compound_scope_class.isin(["core_small_molecule","structure_resolved_small_molecule"])]
    quarantine=pd.read_csv(QUARANTINE,sep="\t",usecols=["compound_internal_id","eligible_for_scaffold_fingerprint_pagtn_similarity"],low_memory=False).drop_duplicates("compound_internal_id")
    master=master.merge(quarantine,on="compound_internal_id",how="left")
    master=master.loc[master.eligible_for_scaffold_fingerprint_pagtn_similarity.fillna(True).astype(str).str.lower().isin(["true","1","yes"])]
    compound_scaffold=pd.read_csv(CATEGORY,sep="\t",compression="gzip",usecols=["compound_internal_id","exact_murcko_scaffold_smiles","scaffold_status"],low_memory=False)
    compound_scaffold=compound_scaffold.merge(master[["compound_internal_id"]],on="compound_internal_id",how="inner",validate="one_to_one")
    compound_scaffold=compound_scaffold.loc[compound_scaffold.scaffold_status.eq("CYCLIC_SCAFFOLD") & compound_scaffold.exact_murcko_scaffold_smiles.notna()]
    compound_scaffold=compound_scaffold.merge(scaff[["exact_murcko_scaffold_smiles","scaffold_rank"]],on="exact_murcko_scaffold_smiles",how="inner",validate="many_to_one")
    assoc=pairs.merge(compound_scaffold[["compound_internal_id","scaffold_rank"]],on="compound_internal_id",how="inner",validate="many_to_one").merge(roles,left_on="target_uniprot_id",right_on="canonical_uniprot_accession",how="left",validate="many_to_one")
    assoc["formal_primary_membrane_role"]=assoc["formal_primary_membrane_role"].fillna("unresolved_target_role_mapping")
    assoc=assoc.drop_duplicates("pair_id")
    test=assoc.loc[assoc.scaffold_rank.isin(candidates.scaffold_rank)].copy()
    role_order=assoc.formal_primary_membrane_role.value_counts().index.tolist()
    total_by_role=assoc.formal_primary_membrane_role.value_counts().reindex(role_order,fill_value=0)
    total_pairs=len(assoc)
    scaffold_totals=test.scaffold_rank.value_counts()
    rows=[]
    for s,total_s in scaffold_totals.items():
        sub=test.loc[test.scaffold_rank.eq(s),"formal_primary_membrane_role"].value_counts()
        for role in role_order:
            a=int(sub.get(role,0)); b=int(total_s-a); c=int(total_by_role[role]-a); d=int(total_pairs-a-b-c)
            _,p=fisher_exact([[a,b],[c,d]],alternative="two-sided")
            log2or=float(np.log2(((a+.5)*(d+.5))/((b+.5)*(c+.5))))
            rows.append({"scaffold_rank":int(s),"formal_primary_membrane_role":role,"observed_unique_protein_compound_pairs":a,"scaffold_pair_total":int(total_s),"role_pair_total":int(total_by_role[role]),"background_pair_total":int(total_pairs),"log2_odds_ratio_haldane_anscombe":log2or,"raw_p":float(p)})
    stats=pd.DataFrame(rows);stats["bh_fdr_q"]=bh(stats.raw_p.to_numpy()); stats=stats.merge(candidates[["scaffold_rank","scaffold_short_id","exact_murcko_scaffold_smiles","unique_compound_count","unique_protein_count","unique_protein_compound_pair_count"]],on="scaffold_rank",how="left")
    stats["significant_with_support"]=stats.bh_fdr_q.lt(.05)&stats.observed_unique_protein_compound_pairs.ge(20)
    stats.to_csv(OUT/"C11_scaffold_role_enrichment_statistics.tsv",sep="\t",index=False)
    counts=stats[["scaffold_rank","scaffold_short_id","exact_murcko_scaffold_smiles","unique_compound_count","unique_protein_count","unique_protein_compound_pair_count"]].drop_duplicates().sort_values("scaffold_rank")
    counts.to_csv(OUT/"C11_scaffold_role_counts.tsv",sep="\t",index=False)
    shown=stats.groupby("scaffold_rank")["significant_with_support"].any(); shown_ids=shown.index[shown].tolist()
    # If no association meets FDR/support, a figure of all testable scaffolds is
    # still produced to make the null result transparent.
    display=stats.loc[stats.scaffold_rank.isin(shown_ids if shown_ids else candidates.scaffold_rank)].copy()
    display_order=display.groupby(["scaffold_rank","scaffold_short_id"])["observed_unique_protein_compound_pairs"].sum().sort_values(ascending=True).index.tolist()
    display["scaffold_short_id"]=pd.Categorical(display.scaffold_short_id,categories=[x[1] for x in display_order],ordered=True)
    display["formal_primary_membrane_role"]=pd.Categorical(display.formal_primary_membrane_role,categories=role_order,ordered=True)
    vmax=max(1.5,float(np.quantile(np.abs(display.log2_odds_ratio_haldane_anscombe),.98)))
    fig_h=max(102,5.0*len(display_order)+40)
    fig=plt.figure(figsize=(178/25.4,fig_h/25.4),dpi=600);ax=fig.add_axes([.26,.12,.62,.80])
    x=display.formal_primary_membrane_role.cat.codes.to_numpy(); y=display.scaffold_short_id.cat.codes.to_numpy(); sizes=8+74*np.sqrt(display.observed_unique_protein_compound_pairs.to_numpy()/max(1,display.observed_unique_protein_compound_pairs.max()))
    sc=ax.scatter(x,y,s=sizes,c=display.log2_odds_ratio_haldane_anscombe,cmap="coolwarm",norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax),alpha=.88,edgecolors=np.where(display.significant_with_support,INK,"none"),linewidths=np.where(display.significant_with_support,.55,0),zorder=3)
    ax.set(xticks=np.arange(len(role_order)),xticklabels=[r.replace("_"," ").title() for r in role_order],yticks=np.arange(len(display_order)),yticklabels=[x[1] for x in display_order],xlabel="Primary membrane role",ylabel="Exact scaffold (rank-derived ID)")
    ax.tick_params(axis="x",labelrotation=40,labelsize=6.7);plt.setp(ax.get_xticklabels(),ha="right");ax.tick_params(axis="y",labelsize=6.7)
    ax.grid(color="#EEF2F4",lw=.55);ax.set_axisbelow(True);ax.spines[["top","right"]].set_visible(False)
    cbar=fig.colorbar(sc,ax=ax,pad=.015,fraction=.035);cbar.set_label("log2 odds ratio",fontsize=7);cbar.ax.tick_params(labelsize=6.5)
    ax.text(0,1.025,f"Tested {len(candidates):,} high-support scaffolds (≥150 compounds and ≥70 proteins); displayed {len(display_order):,} with BH-FDR < 0.05 and observed support ≥20",transform=ax.transAxes,fontsize=6.7,color=SLATE)
    ax.text(0,-.19,"Bubble area: observed unique protein–compound pairs. Black outline: BH-FDR < 0.05 with support ≥20. All tested scaffold–role statistics are in Source Data.",transform=ax.transAxes,fontsize=6.4,color=SLATE)
    for sfx,k in {".svg":{},".pdf":{},".png":{"dpi":600}}.items():fig.savefig(OUT/f"C11_scaffold_role_enrichment{sfx}",bbox_inches="tight",**k)
    plt.close(fig)
    (OUT/"C11_scaffold_filtering_rules.txt").write_text(f"Analysis unit\tunique canonical protein–compound pair\nExact scaffold\texisting exact_murcko_scaffold_smiles from C3\nCandidate rule\t≥150 unique canonical compounds AND ≥70 unique proteins (high-support joint tail, selected before Fisher testing for final-size readability)\nCandidate scaffolds\t{len(candidates):,}\nHypotheses\t{len(stats):,} scaffold × role Fisher tests\nDisplay rule\tall candidate scaffolds with at least one BH-FDR <0.05 cell and observed pair support ≥20; no manual scaffold selection\n",encoding="utf-8")
    (OUT/"C11_reproducibility_manifest.json").write_text(json.dumps({"figure":"C11_scaffold_role_enrichment","inputs":[str(C3),str(PAIR),str(ROLE)],"test":"two-sided Fisher exact; Haldane–Anscombe log2 odds ratio; BH-FDR across all tested cells"},indent=2),encoding="utf-8")
if __name__=="__main__":main()
