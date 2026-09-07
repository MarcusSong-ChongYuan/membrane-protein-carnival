"""C10: Cross-role breadth of canonical MemPro compounds."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compound_figure_pool"
FORMAL = Path(r"D:\finale\FORMAL")
PAIR = FORMAL / "01_core_v72" / "01_release_tables" / "protein_compound_pair_v72.tsv.gz"
ROLE = FORMAL / "03_publication_repairs" / "protein_membrane_role_FORMAL.tsv"
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans","sans-serif"],"svg.fonttype":"none","pdf.fonttype":42,"figure.facecolor":"white","savefig.facecolor":"white","axes.linewidth":0.8})
INK, SLATE, NAVY, ACCENT = "#24364B", "#5E7182", "#173B6C", "#D6A23A"

def entropy(counts: np.ndarray) -> float:
    p = counts / counts.sum()
    return float(-(p * np.log2(p)).sum())

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pair = pd.read_csv(PAIR, sep="\t", compression="gzip", usecols=["pair_id","target_uniprot_id","compound_internal_id"], low_memory=False)
    role = pd.read_csv(ROLE, sep="\t", usecols=["canonical_uniprot_accession","formal_primary_membrane_role"], low_memory=False).drop_duplicates("canonical_uniprot_accession")
    d = pair.merge(role, left_on="target_uniprot_id", right_on="canonical_uniprot_accession", how="left", validate="many_to_one")
    d = d.drop_duplicates(["compound_internal_id","target_uniprot_id"])
    d["formal_primary_membrane_role"] = d["formal_primary_membrane_role"].fillna("unresolved_target_role_mapping")
    counts = d.groupby("compound_internal_id", as_index=False).agg(
        number_of_distinct_target_roles=("formal_primary_membrane_role","nunique"),
        number_of_unique_membrane_targets=("target_uniprot_id","nunique"),
        number_of_canonical_pairs=("pair_id","nunique"),
    )
    role_counts = d.groupby(["compound_internal_id","formal_primary_membrane_role"]).size().rename("pair_count").reset_index()
    ent = role_counts.groupby("compound_internal_id")["pair_count"].apply(lambda x: entropy(x.to_numpy())).rename("target_role_shannon_entropy").reset_index()
    counts = counts.merge(ent, on="compound_internal_id", how="left")
    counts["role_breadth_display"] = counts["number_of_distinct_target_roles"].clip(upper=5).astype(str).replace({"5":"5+"})
    counts["role_breadth_interpretation"] = np.where(counts["number_of_distinct_target_roles"].eq(1), "role-restricted", "cross-role")
    counts.to_csv(OUT / "C10_compound_role_breadth.tsv", sep="\t", index=False)
    ent.to_csv(OUT / "C10_role_entropy.tsv", sep="\t", index=False)
    order = ["1","2","3","4","5+"]
    freq = counts["role_breadth_display"].value_counts().reindex(order, fill_value=0).rename_axis("role_breadth").reset_index(name="canonical_compound_count")
    freq["fraction"] = freq["canonical_compound_count"] / len(counts)
    freq.to_csv(OUT / "C10_compound_role_breadth_statistics.tsv", sep="\t", index=False)
    fig = plt.figure(figsize=(178/25.4, 76/25.4), dpi=600); ax = fig.add_axes([.16,.22,.70,.63])
    x=np.arange(len(freq)); ax.vlines(x,0,freq["fraction"],color="#B7C4CE",lw=2.0,zorder=1)
    cols=[NAVY if v=="1" else "#168C8C" for v in freq.role_breadth]
    ax.scatter(x,freq["fraction"],s=72,color=cols,edgecolor="white",linewidth=.8,zorder=3)
    for i,row in freq.iterrows(): ax.text(i,row.fraction+.018,f"{row.canonical_compound_count:,}\n({row.fraction:.1%})",ha="center",va="bottom",fontsize=7,color=INK)
    ax.set(xticks=x,xticklabels=["1\nrole-restricted","2\ncross-role","3\ncross-role","4\ncross-role","5+\ncross-role"],ylim=(0,max(.06,float(freq.fraction.max())*1.31)),ylabel="Fraction of canonical compounds",xlabel="Number of distinct membrane-protein roles")
    ax.tick_params(axis="both",labelsize=7,colors=INK); ax.grid(axis="y",color="#E4EAED",lw=.65); ax.set_axisbelow(True); ax.spines[["top","right"]].set_visible(False)
    ax.text(0,1.04,f"Analysis N = {len(counts):,} canonical compounds; role breadth is a multi-label compound property",transform=ax.transAxes,fontsize=7.2,color=SLATE)
    ax.text(0,-.27,"A compound may connect to multiple membrane-protein roles. ‘Role-restricted’ describes database breadth, not pharmacological selectivity.",transform=ax.transAxes,fontsize=6.8,color=SLATE)
    for s,k in {".svg":{},".pdf":{},".png":{"dpi":600}}.items(): fig.savefig(OUT/f"C10_compound_role_breadth{s}",bbox_inches="tight",**k)
    plt.close(fig)
    (OUT/"C10_role_breadth_statistics.txt").write_text(f"Analysis unit\tone canonical compound\nN\t{len(counts):,}\nInput canonical protein–compound pairs\t{len(pair):,}\nRole mapping missing\t{int(d.formal_primary_membrane_role.eq('unresolved_target_role_mapping').sum()):,} pair rows\n",encoding="utf-8")
    (OUT/"C10_reproducibility_manifest.json").write_text(json.dumps({"figure":"C10_compound_role_breadth","analysis_unit":"one canonical compound","inputs":[str(PAIR),str(ROLE)],"deduplication_key":["compound_internal_id","target_uniprot_id"]},indent=2),encoding="utf-8")
if __name__ == "__main__": main()
