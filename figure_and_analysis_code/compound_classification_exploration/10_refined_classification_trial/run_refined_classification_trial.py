from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Lipinski, rdMolDescriptors
from scipy.stats import chi2_contingency
from sklearn.compose import ColumnTransformer
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MultiLabelBinarizer, OneHotEncoder
from sklearn.linear_model import LogisticRegression

ROOT = Path(r"D:\finale\compound_classification_exploration")
OUT = ROOT / "10_refined_classification_trial"
DATA = OUT / "data"; FIG = OUT / "figures"; QA = OUT / "qa"
for p in (DATA, FIG, QA): p.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT / "02_classification_tables" / "INTERACTION_LINKED_COMPOUND_ANALYSIS.tsv.gz"

SEED = 42
PALETTE = ["#7b95c6", "#49c2d9", "#67a583", "#a2c986", "#fded95", "#ffc1a6", "#f59c7c", "#c85e62", "#8e7db5"]
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,
                     "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.8,
                     "svg.fonttype":"none","pdf.fonttype":42,"legend.frameon":False})

FG_SMARTS = {
    "alcohol": "[OX2H][CX4;!$(C=O)]", "phenol": "[OX2H]c", "ether": "[OD2]([#6])[#6]",
    "ester": "[CX3](=O)[OX2][#6]", "amide": "[NX3][CX3](=[OX1])", "primary_or_secondary_amine": "[NX3;H1,H2;!$(NC=O)]",
    "tertiary_amine": "[NX3;H0;!$(N-C=O);!$(N=*);!$([N+])]", "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ketone": "[#6][CX3](=O)[#6]", "aldehyde": "[CX3H1](=O)[#6]", "sulfonamide": "[SX4](=O)(=O)[NX3]",
    "urea": "[NX3][CX3](=[OX1])[NX3]", "carbamate": "[NX3][CX3](=[OX1])[OX2]", "nitrile": "[CX2]#N",
    "nitro": "[$([NX3](=O)=O),$([N+](=O)[O-])]", "halogen": "[F,Cl,Br,I]", "phosphate": "P(=O)(O)O",
    "sulfate_or_sulfonate": "S(=O)(=O)[O-,$([OH])]",
}
MOTIF_SMARTS = {
    "benzene": "c1ccccc1", "pyridine": "n1ccccc1", "pyrimidine": "n1ccnc(n1)", "purine": "n1cnc2ncnc12",
    "indole": "c1ccc2[nH]ccc2c1", "quinoline": "c1ccc2ncccc2c1", "isoquinoline": "c1ccc2ccncc2c1",
    "imidazole": "c1ncc[nH]1", "benzimidazole": "c1ccc2[nH]cnc2c1", "piperidine": "N1CCCCC1",
    "piperazine": "N1CCNCC1", "morpholine": "O1CCNCC1",
    "steroid_like": "C1CCC2C3CCC4CCCCC4C3CCC12", "flavonoid_like": "O=c1cc(-c2ccccc2)oc2ccccc12",
}
FG = {k: Chem.MolFromSmarts(v) for k,v in FG_SMARTS.items()}
MOTIFS = {k: Chem.MolFromSmarts(v) for k,v in MOTIF_SMARTS.items()}

def export(fig, base: Path):
    fig.savefig(base.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)

def classify_mol(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    rings = list(mol.GetRingInfo().AtomRings())
    nring = len(rings); macro = any(len(r)>=12 for r in rings)
    spiro_n = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    bridge_n = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    fused = any(len(set(rings[i]) & set(rings[j])) >= 2 for i in range(nring) for j in range(i))
    if nring == 0: refined = "Acyclic"
    elif macro: refined = "Macrocyclic"
    elif spiro_n > 0: refined = "Spirocyclic"
    elif bridge_n > 0: refined = "Bridged polycyclic"
    elif fused: refined = "Fused polycyclic"
    elif nring == 1: refined = "Monocyclic"
    else: refined = "Non-fused polycyclic"
    old = "Acyclic" if nring==0 else ("Macrocyclic" if macro else ("Monocyclic" if nring==1 else "Other polycyclic"))
    pos = any(a.GetFormalCharge()>0 for a in mol.GetAtoms()); neg = any(a.GetFormalCharge()<0 for a in mol.GetAtoms())
    charge = "Zwitterionic" if pos and neg else ("Cationic" if pos else ("Anionic" if neg else "Neutral"))
    fgs = [k for k,q in FG.items() if q is not None and mol.HasSubstructMatch(q)] or ["none_of_selected_groups"]
    motifs = [k for k,q in MOTIFS.items() if q is not None and mol.HasSubstructMatch(q)] or ["none_of_selected_motifs"]
    return old, refined, nring, spiro_n, bridge_n, int(fused), charge, ";".join(fgs), ";".join(motifs)

def bins(row):
    mw = "MW<300" if row.molecular_weight<300 else ("MW300-500" if row.molecular_weight<=500 else "MW>500")
    lp = "LogP<1" if row.logp<1 else ("LogP1-3" if row.logp<3 else ("LogP3-5" if row.logp<=5 else "LogP>5"))
    psa = "TPSA<60" if row.tpsa<60 else ("TPSA60-120" if row.tpsa<=120 else "TPSA>120")
    flex = "RB0-3" if row.rotatable_bonds<=3 else ("RB4-8" if row.rotatable_bonds<=8 else "RB>8")
    sp3 = "Fsp3<0.25" if row.fraction_csp3<.25 else ("Fsp3=0.25-0.50" if row.fraction_csp3<=.5 else "Fsp3>0.50")
    return mw,lp,psa,flex,sp3

def cramers_v(x, y):
    tab=pd.crosstab(x,y); chi2=chi2_contingency(tab,correction=False)[0]; n=tab.values.sum(); r,k=tab.shape
    phi2=chi2/n; phi2c=max(0,phi2-((k-1)*(r-1))/(n-1)); rc=r-((r-1)**2)/(n-1); kc=k-((k-1)**2)/(n-1)
    return math.sqrt(phi2c/max(1e-12,min(kc-1,rc-1))), r, k

def onehot_strings(series):
    enc=OneHotEncoder(handle_unknown="ignore", sparse_output=True)
    return enc.fit_transform(series.to_numpy().reshape(-1,1))

def multi_hot(series):
    mlb=MultiLabelBinarizer(sparse_output=True)
    return mlb.fit_transform(series.str.split(";"))

def evaluate(X, y, train_idx, test_idx):
    model=LogisticRegression(max_iter=300,solver="saga",class_weight="balanced",n_jobs=1,random_state=SEED)
    model.fit(X[train_idx],y[train_idx]); pred=model.predict(X[test_idx])
    return balanced_accuracy_score(y[test_idx],pred), f1_score(y[test_idx],pred,average="macro")

def main():
    d=pd.read_csv(SOURCE,sep="\t",low_memory=False)
    records=[]
    for i,smi in enumerate(d.canonical_smiles_recalculated):
        c=classify_mol(smi)
        if c is None: records.append((None,)*9)
        else: records.append(c)
        if (i+1)%25000==0: print(f"classified {i+1}/{len(d)}",flush=True)
    cols=["old_ring_class","refined_ring_topology","rdkit_ring_count","spiro_atom_count","bridgehead_atom_count","has_fused_rings","record_charge_class","functional_groups_multilabel","named_scaffold_motifs_multilabel"]
    c=pd.DataFrame(records,columns=cols)
    b=pd.DataFrame([bins(x) for x in d.itertuples()],columns=["mw_bin","logp_bin","tpsa_bin","rotatable_bond_bin","fraction_csp3_bin"])
    out=pd.concat([d[["compound_internal_id","dominant_membrane_role_descriptive"]].reset_index(drop=True),c,b],axis=1)
    out["multi_axis_profile"] = out[["refined_ring_topology","record_charge_class","mw_bin","logp_bin","tpsa_bin","fraction_csp3_bin"]].astype(str).agg("|".join,axis=1)
    out.to_csv(DATA/"REFINED_COMPOUND_CLASSIFICATION.tsv.gz",sep="\t",index=False,compression="gzip")

    summary=[]
    for col in ["old_ring_class","refined_ring_topology","record_charge_class","mw_bin","logp_bin","tpsa_bin","fraction_csp3_bin","chemical_regime_local"]:
        s = out[col] if col in out else d[col]
        vc=s.value_counts(); p=vc/vc.sum(); summary.append({"classification":col,"compound_n":int(s.notna().sum()),"category_n":len(vc),"effective_category_n":float(np.exp(-(p*np.log(p)).sum())),"largest_category_fraction":float(p.max()),"cramers_v_vs_dominant_role":cramers_v(s,d.dominant_membrane_role_descriptive)[0]})
    summary=pd.DataFrame(summary)

    y=d.dominant_membrane_role_descriptive.fillna("Unresolved").astype(str).to_numpy()
    idx=np.arange(len(d)); tr,te=train_test_split(idx,test_size=.2,random_state=SEED,stratify=y)
    feature_sets={
        "Old ring class": onehot_strings(out.old_ring_class),
        "Refined ring topology": onehot_strings(out.refined_ring_topology),
        "Local chemical regime": onehot_strings(d.chemical_regime_local),
        "Functional groups": multi_hot(out.functional_groups_multilabel),
        "Named scaffold motifs": multi_hot(out.named_scaffold_motifs_multilabel),
        "Multi-axis physicochemical": pd.concat([out[["refined_ring_topology","record_charge_class","mw_bin","logp_bin","tpsa_bin","rotatable_bond_bin","fraction_csp3_bin"]]],axis=1),
    }
    # encode combined categorical frame separately
    ct=ColumnTransformer([("cat",OneHotEncoder(handle_unknown="ignore"),feature_sets["Multi-axis physicochemical"].columns)])
    feature_sets["Multi-axis physicochemical"]=ct.fit_transform(feature_sets["Multi-axis physicochemical"])
    metrics=[]
    for name,X in feature_sets.items():
        ba,mf=evaluate(X,y,tr,te); metrics.append({"feature_set":name,"balanced_accuracy":ba,"macro_f1":mf,"feature_count":X.shape[1],"train_n":len(tr),"test_n":len(te)})
        print(name,ba,mf,flush=True)
    metrics=pd.DataFrame(metrics)
    summary.to_csv(DATA/"CLASSIFICATION_RESOLUTION_AND_ASSOCIATION.tsv",sep="\t",index=False)
    metrics.to_csv(DATA/"HELDOUT_ROLE_DISCRIMINATION.tsv",sep="\t",index=False)

    # Counts and coverage
    rc=out.refined_ring_topology.value_counts().rename_axis("refined_ring_topology").reset_index(name="compound_count")
    rc["fraction"]=rc.compound_count/len(out); rc.to_csv(DATA/"REFINED_RING_TOPOLOGY_COUNTS.tsv",sep="\t",index=False)
    fig,ax=plt.subplots(figsize=(6.8,4.2)); ax.barh(rc.refined_ring_topology[::-1],rc.compound_count[::-1],color=PALETTE[:len(rc)][::-1],alpha=.82)
    for i,v in enumerate(rc.compound_count[::-1]): ax.text(v,i,f" {v:,} ({v/len(out)*100:.1f}%)",va="center",fontsize=6)
    ax.set_xlabel("Structure-valid interaction-linked compounds"); ax.set_title("Refined ring topology",loc="left",fontweight="bold"); export(fig,FIG/"R1_refined_ring_topology")

    # Direct method comparison
    m=metrics.sort_values("macro_f1")
    fig,ax=plt.subplots(figsize=(6.8,4.2)); y0=np.arange(len(m)); ax.plot(m.macro_f1,y0,"o",color="#c85e62",label="Macro-F1"); ax.plot(m.balanced_accuracy,y0,"s",color="#7b95c6",label="Balanced accuracy")
    ax.set_yticks(y0,m.feature_set); ax.set_xlabel("Held-out dominant-role classification score"); ax.set_xlim(0,max(.35,m[["macro_f1","balanced_accuracy"]].to_numpy().max()*1.15)); ax.legend(); ax.set_title("Does finer chemistry improve target-role discrimination?",loc="left",fontweight="bold"); export(fig,FIG/"R2_heldout_role_discrimination")

    # Refined topology x role profiles
    tab=pd.crosstab(out.refined_ring_topology,d.dominant_membrane_role_descriptive); prof=tab.div(tab.sum(axis=1),axis=0); prof.to_csv(DATA/"REFINED_RING_ROLE_PROFILE.tsv",sep="\t")
    fig,ax=plt.subplots(figsize=(8.2,4.2)); sns.heatmap(prof,cmap="crest",linewidths=.35,linecolor="white",cbar_kws={"label":"Within-topology proportion"},ax=ax)
    ax.set_xlabel("Dominant formal membrane role"); ax.set_ylabel(""); ax.tick_params(axis="x",rotation=40); ax.set_title("Refined ring topology × membrane-target role",loc="left",fontweight="bold"); export(fig,FIG/"R3_refined_ring_role_profile")

    # Functional groups x role log2 enrichment
    fg_long=out[["compound_internal_id","functional_groups_multilabel"]].assign(functional_group=lambda x:x.functional_groups_multilabel.str.split(";")).explode("functional_group")
    fg_long["role"]=np.repeat(d.dominant_membrane_role_descriptive.to_numpy(),out.functional_groups_multilabel.str.count(";").add(1))
    ft=pd.crosstab(fg_long.functional_group,fg_long.role); expected=np.outer(ft.sum(axis=1),ft.sum(axis=0))/ft.values.sum(); log2or=np.log2((ft+0.5)/(expected+0.5)); support=ft.sum(axis=1); keep=support[support>=200].index; log2or=log2or.loc[keep]; log2or.to_csv(DATA/"FUNCTIONAL_GROUP_ROLE_LOG2_ENRICHMENT.tsv",sep="\t")
    fig,ax=plt.subplots(figsize=(8.4,5.4)); sns.heatmap(log2or.clip(-2,2),cmap="vlag",center=0,vmin=-2,vmax=2,linewidths=.3,linecolor="white",cbar_kws={"label":"log₂(observed / expected)"},ax=ax)
    ax.set_xlabel("Dominant formal membrane role"); ax.set_ylabel(""); ax.tick_params(axis="x",rotation=40); ax.set_title("Functional-group enrichment across membrane-target roles",loc="left",fontweight="bold"); export(fig,FIG/"R4_functional_group_role_enrichment")

    # Resolution/association comparison
    s=summary.sort_values("cramers_v_vs_dominant_role")
    fig,axs=plt.subplots(1,2,figsize=(8.2,4.0)); axs[0].barh(s.classification,s.effective_category_n,color="#a1d8e8"); axs[0].set_xlabel("Effective number of categories"); axs[0].set_title("Resolution",loc="left",fontweight="bold")
    axs[1].barh(s.classification,s.cramers_v_vs_dominant_role,color="#7b95c6"); axs[1].set_xlabel("Bias-corrected Cramér's V"); axs[1].set_title("Association with dominant target role",loc="left",fontweight="bold"); export(fig,FIG/"R5_resolution_vs_role_association")

    conclusions={"denominator":len(d),"old_ring":metrics[metrics.feature_set=="Old ring class"].iloc[0].to_dict(),"refined_ring":metrics[metrics.feature_set=="Refined ring topology"].iloc[0].to_dict(),"best_feature_set":metrics.sort_values("macro_f1",ascending=False).iloc[0].to_dict(),"interpretation":"A method is considered more useful only if held-out scores and role association improve, not merely category count."}
    (OUT/"REFINED_CLASSIFICATION_TRIAL_RESULTS.json").write_text(json.dumps(conclusions,indent=2,default=float),encoding="utf-8")
    print(json.dumps(conclusions,indent=2,default=float))

if __name__=="__main__": main()
