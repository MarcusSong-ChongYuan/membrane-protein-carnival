import csv,gzip,json
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813");OUT=ROOT/'04_negative_conflicts';OUT.mkdir(exist_ok=True)
AUTO=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data")
FULL=Path(r"D:\finale\02_V6.3_candidate_20260803\01_baseline_v6_2")
NEG=Path(r"D:\finale\negative_upgrade\upgraded_negative_evidence.tsv.gz")
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
proteins=set();positive_pairs=set()
with op(AUTO/'protein_master_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):proteins.add(r['canonical_uniprot_accession'])
with op(AUTO/'protein_compound_pair_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):positive_pairs.add((r['target_uniprot_id'],r['compound_internal_id']))

# Pass 1: strict eligibility independent of positive compound membership.
compound_ids=set();form_ids=set();eligible_ids=set();counts=Counter();negpairs=set()
with op(NEG) as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['target_uniprot_id'] not in proteins:continue
  counts['v7_protein_rows']+=1
  mapped=r.get('release_action_v62') in {'mapped_negative_release','mapped_negative_release_existing'}
  has_comp=bool(r.get('compound_internal_id'));has_form=bool(r.get('compound_form_id'))
  has_keys=bool(r.get('exact_inchikey_v62')) and bool(r.get('canonical_parent_inchikey_v62'))
  relation=r.get('canonical_relation_v62','')
  structurally_safe=(r.get('exact_inchikey_v62')==r.get('canonical_parent_inchikey_v62') or relation in {'stereochemically_defined_parent','canonical_parent','exact_form','baseline_form'})
  if mapped and has_comp and has_keys and structurally_safe:
   eligible_ids.add(r.get('negative_evidence_id_v62') or r['source_evidence_id']);compound_ids.add(r['compound_internal_id']);
   if has_form:form_ids.add(r['compound_form_id'])
   negpairs.add((r['target_uniprot_id'],r['compound_internal_id']));counts['eligible_public_negative']+=1
  else:
   counts['frozen_review']+=1

# Rebuild compound/form master as union of positive and negative canonical entities.
pos_comp=set()
with op(AUTO/'compound_master_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):pos_comp.add(r['compound_internal_id'])
union_comp=pos_comp|compound_ids
master_src=FULL/'small_molecule_master_v1_3.tsv';form_src=FULL/'compound_form_hierarchy_v1_3.tsv'
def filter_table(src,dst,pred):
 ni=no=0
 with op(src) as fi:
  rd=csv.DictReader(fi,delimiter='\t')
  with op(dst,'wt') as fo:
   w=csv.DictWriter(fo,fieldnames=rd.fieldnames,delimiter='\t');w.writeheader()
   for r in rd:
    ni+=1
    if pred(r):w.writerow(r);no+=1
 return ni,no
mi,mo=filter_table(master_src,OUT/'compound_master_positive_negative_union_v7.tsv.gz',lambda r:r['compound_internal_id'] in union_comp)
fi,fo=filter_table(form_src,OUT/'compound_form_positive_negative_union_v7.tsv.gz',lambda r:r['compound_internal_id'] in union_comp and (not form_ids or r['compound_form_id'] in form_ids or r['compound_internal_id'] in pos_comp))

# Pass 2 write public and frozen records.
pub=OUT/'negative_evidence_v7.tsv.gz';rev=OUT/'negative_evidence_frozen_v7.tsv.gz';np=nr=0
with op(NEG) as fi:
 rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['v7_disposition','v7_rule']
 with op(pub,'wt') as fp,op(rev,'wt') as fr:
  wp=csv.DictWriter(fp,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wp.writeheader();wr.writeheader()
  for r in rd:
   if r['target_uniprot_id'] not in proteins:continue
   nid=r.get('negative_evidence_id_v62') or r['source_evidence_id']
   if nid in eligible_ids:r['v7_disposition']='PUBLIC_NEGATIVE';r['v7_rule']='mapped release + canonical compound/form structural identity';wp.writerow(r);np+=1
   else:r['v7_disposition']='FROZEN_REVIEW';r['v7_rule']='identity or release rule not satisfied';wr.writerow(r);nr+=1

conf=sorted(negpairs & positive_pairs)
with op(OUT/'positive_negative_conflict_v7.tsv.gz','wt') as f:
 fields=['conflict_id','target_uniprot_id','compound_internal_id','conflict_type','automatic_interpretation','positive_retained','negative_retained'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for i,p in enumerate(conf,1):w.writerow({'conflict_id':f'PNC-V7-{i:08d}','target_uniprot_id':p[0],'compound_internal_id':p[1],'conflict_type':'positive_and_negative_across_assay_contexts','automatic_interpretation':'context_dependent; assay, concentration and endpoint must remain separate','positive_retained':'1','negative_retained':'1'})
report={'input_v7_protein_negative_rows':np+nr,'public_negative_rows':np,'frozen_negative_rows':nr,'unique_negative_pairs':len(negpairs),'positive_negative_conflict_pairs':len(conf),'positive_compounds':len(pos_comp),'negative_compounds':len(compound_ids),'union_compounds_requested':len(union_comp),'union_compound_master_rows':mo,'union_form_rows':fo,'policy':'Negative compounds need not occur in the positive master. Canonical identity is validated against the full V6.2 compound/form master.'}
(OUT/'T30_T31_NEGATIVE_REBUILD_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
