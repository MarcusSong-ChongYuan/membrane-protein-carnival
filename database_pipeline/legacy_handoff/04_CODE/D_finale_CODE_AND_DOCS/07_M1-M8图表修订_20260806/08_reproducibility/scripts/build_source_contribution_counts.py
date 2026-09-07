from __future__ import annotations
import os
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd

ROOT=Path(os.environ.get("MEMPRO_REVISION_ROOT",r"D:\finale\07_M1-M8图表修订_20260806"))
V62=Path(r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730")
def truth(s): return s.astype(str).str.lower().isin(["1","true","yes"])

def main():
    versions=pd.read_csv(ROOT/"01_source_audit/DATABASE_VERSION_AND_CONTRIBUTION_BASELINE.tsv",sep="\t")
    vmap=dict(zip(versions.Database,versions["Version / snapshot"]))
    p=pd.read_csv(V62/"human_membrane_protein_master_v6_2.tsv",sep="\t",usecols=["target_uniprot_id","hpa_mapping_status_v62","htp_present_v5","membranome_present_v5","opm_present_v5","pdbtm_present_v5","alphafolddb_ids"],low_memory=False)
    rows=[]
    def add(db,module,count,unit,proteins=None,compounds=None,sites=None,relations=None,note=""):
        rows.append({"module":module,"database":db,"version_snapshot":vmap.get(db,"see source registry"),"contributed_count":int(count),"primary_count_unit":unit,"unique_proteins":proteins,"unique_compounds":compounds,"binding_site_instances":sites,"relation_or_evidence_records":relations,"counting_note":note})
    add("UniProtKB","Protein",len(p),"canonical proteins",len(p),note="Canonical accession/sequence anchor for every protein row")
    add("Human Protein Atlas","Protein/expression",p.hpa_mapping_status_v62.eq("mapped").sum(),"HPA-mapped proteins",p.hpa_mapping_status_v62.eq("mapped").sum(),note="Expression matrices contain multiple tissue/cell records per protein")
    add("Human Transmembrane Proteome","Protein",truth(p.htp_present_v5).sum(),"proteins with HTP evidence",truth(p.htp_present_v5).sum())
    add("Membranome / BioMembHub","Protein",truth(p.membranome_present_v5).sum(),"proteins with Membranome evidence",truth(p.membranome_present_v5).sum())
    add("OPM","Structure",truth(p.opm_present_v5).sum(),"proteins with OPM structure",truth(p.opm_present_v5).sum())
    add("PDBTM / UniTmp","Structure",truth(p.pdbtm_present_v5).sum(),"proteins with PDBTM structure",truth(p.pdbtm_present_v5).sum())
    add("AlphaFoldDB","Structure",p.alphafolddb_ids.fillna("").ne("").sum(),"proteins with model ID",p.alphafolddb_ids.fillna("").ne("").sum(),note="Prediction coverage; not experimental validation")

    site=pd.read_csv(V62/"binding_site_instances_v6_2.tsv.gz",sep="\t",compression="gzip",usecols=["source_database","target_uniprot_id","compound_internal_id"],low_memory=False)
    for db,g in site.groupby("source_database"):
        add(db,"Binding site",len(g),"site instances",g.target_uniprot_id.nunique(),g.compound_internal_id.nunique(),len(g),note="One site instance may contain multiple contact residues")

    ev_count=Counter(); targets=defaultdict(set); compounds=defaultdict(set)
    for ch in pd.read_csv(V62/"binding_evidence_master_v6_2.tsv.gz",sep="\t",compression="gzip",usecols=["source_database","target_uniprot_id","compound_internal_id"],chunksize=250000,low_memory=False):
        for db,g in ch.groupby("source_database"):
            ev_count[db]+=len(g); targets[db].update(g.target_uniprot_id.dropna().astype(str)); compounds[db].update(g.compound_internal_id.dropna().astype(str))
    for db,n in ev_count.items():
        add(db,"Interaction",n,"evidence records",len(targets[db]),len(compounds[db]),relations=n,note="Evidence records, not unique protein–compound pairs")

    dis=pd.read_csv(V62/"protein_gene_disease_relations_v6_2.tsv",sep="\t",usecols=["source_database","target_uniprot_id","disease_relation_id"],low_memory=False)
    for db,g in dis.groupby("source_database"):
        add(db,"Disease",len(g),"source protein–disease relations",g.target_uniprot_id.nunique(),relations=len(g),note="V6.2 source relations before V6.3 exact-only disease canonicalization")
    out=pd.DataFrame(rows).sort_values(["module","database"])
    out.to_csv(ROOT/"01_source_audit/DATABASE_EXACT_CONTRIBUTION_COUNTS.tsv",sep="\t",index=False)
    print(out.to_string(index=False))
if __name__=="__main__": main()
