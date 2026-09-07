from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MultiLabelBinarizer, OneHotEncoder

ROOT=Path(r"D:\finale\compound_classification_exploration")
OUT=ROOT/"10_refined_classification_trial"; DATA=OUT/"data"
base=pd.read_csv(ROOT/"02_classification_tables/INTERACTION_LINKED_COMPOUND_ANALYSIS.tsv.gz",sep="\t",low_memory=False)
ref=pd.read_csv(DATA/"REFINED_COMPOUND_CLASSIFICATION.tsv.gz",sep="\t",low_memory=False)
d=base.merge(ref.drop(columns="dominant_membrane_role_descriptive"),on="compound_internal_id",validate="one_to_one")
counts=d.dominant_membrane_role_descriptive.value_counts()
roles=counts[(counts>=500)&(~counts.index.str.contains("unresolved",case=False))].index
x=d[d.dominant_membrane_role_descriptive.isin(roles)].reset_index(drop=True)
y=x.dominant_membrane_role_descriptive.to_numpy(); idx=np.arange(len(x)); tr,te=train_test_split(idx,test_size=.2,random_state=42,stratify=y)

def onehot(cols):
    enc=ColumnTransformer([("cat",OneHotEncoder(handle_unknown="ignore"),cols)])
    return enc.fit_transform(x)
def multihot(col):
    enc=MultiLabelBinarizer(sparse_output=True); return enc.fit_transform(x[col].str.split(";"))

features={
 "Old ring class":onehot(["old_ring_class"]),
 "Refined ring topology":onehot(["refined_ring_topology"]),
 "Local chemical regime":onehot(["chemical_regime_local"]),
 "Functional groups":multihot("functional_groups_multilabel"),
 "Named scaffold motifs":multihot("named_scaffold_motifs_multilabel"),
 "Multi-axis physicochemical":onehot(["refined_ring_topology","record_charge_class","mw_bin","logp_bin","tpsa_bin","rotatable_bond_bin","fraction_csp3_bin"]),
}
rows=[]
for name,X in features.items():
    model=LogisticRegression(max_iter=500,solver="lbfgs",class_weight="balanced",tol=1e-5)
    model.fit(X[tr],y[tr]); pred=model.predict(X[te])
    rows.append({"feature_set":name,"balanced_accuracy":balanced_accuracy_score(y[te],pred),"macro_f1":f1_score(y[te],pred,average="macro"),"feature_count":X.shape[1],"train_n":len(tr),"test_n":len(te),"model_iterations":int(np.max(model.n_iter_))})
    print(rows[-1],flush=True)
res=pd.DataFrame(rows); res.to_csv(DATA/"HELDOUT_DEFINED_ROLE_DISCRIMINATION.tsv",sep="\t",index=False)
meta={"source_n":len(d),"analysis_n":len(x),"excluded_n":len(d)-len(x),"inclusion_rule":"defined dominant membrane role with at least 500 compounds","included_roles":counts.loc[roles].to_dict(),"random_seed":42,"test_fraction":.2,"model":"class-weighted multinomial logistic regression; lbfgs; max_iter=500; tol=1e-5"}
(DATA/"HELDOUT_DEFINED_ROLE_DISCRIMINATION_METADATA.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
print(json.dumps(meta,indent=2))
