import csv,gzip,json
from pathlib import Path
from collections import Counter

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813");OUT=ROOT/'04_negative_conflicts';AUTO=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data");NEG=Path(r"D:\finale\negative_upgrade\upgraded_negative_evidence.tsv.gz");BASE=Path(r"D:\finale\negative_upgrade\backup\small_molecule_master_v1_3.tsv");APP=Path(r"D:\finale\negative_upgrade\new_compounds_v1_3_append.tsv");FORMBASE=Path(r"D:\finale\negative_upgrade\backup\compound_form_hierarchy_v1_3.tsv");FORMAPP=Path(r"D:\finale\negative_upgrade\new_forms_v1_3_append.tsv")
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
proteins=set();positive_pairs=set();positive_comp=set()
with op(AUTO/'protein_master_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):proteins.add(r['canonical_uniprot_accession'])
with op(AUTO/'protein_compound_pair_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):positive_pairs.add((r['target_uniprot_id'],r['compound_internal_id']));positive_comp.add(r['compound_internal_id'])

eligible=set();needed_comp=set();needed_form=set();negpairs=set();reason=Counter()
with op(NEG) as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['target_uniprot_id'] not in proteins:continue
  nid=r.get('negative_evidence_id_v62') or r['source_evidence_id'];mapped=r.get('release_action_v62') in {'mapped_negative_release','mapped_negative_release_existing'};keys=bool(r.get('exact_inchikey_v62')) and bool(r.get('canonical_parent_inchikey_v62'));rel=r.get('canonical_relation_v62','');safe=(r.get('exact_inchikey_v62')==r.get('canonical_parent_inchikey_v62') or rel in {'stereochemically_defined_parent','canonical_parent','exact_form','baseline_form'})
  if mapped and keys and safe and r.get('compound_internal_id'):
   eligible.add(nid);needed_comp.add(r['compound_internal_id']);negpairs.add((r['target_uniprot_id'],r['compound_internal_id']))
   if r.get('compound_form_id'):needed_form.add(r['compound_form_id'])
  else:reason['identity_or_release_rule']+=1

# Merge by internal primary key; duplicate IDs must be exact in identity_key/InChIKey.
comp_rows={};comp_conflict=[];fields=None
for p in [BASE,APP]:
 with op(p) as f:
  rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames
  for r in rd:
   cid=r['compound_internal_id']
   if cid not in needed_comp and cid not in positive_comp:continue
   if cid in comp_rows and (comp_rows[cid].get('identity_key'),comp_rows[cid].get('standard_inchikey'))!=(r.get('identity_key'),r.get('standard_inchikey')):comp_conflict.append(cid)
   else:comp_rows[cid]=r
form_rows={};form_conflict=[];fform=None
for p in [FORMBASE,FORMAPP]:
 with op(p) as f:
  rd=csv.DictReader(f,delimiter='\t');fform=rd.fieldnames
  for r in rd:
   fid=r['compound_form_id']
   if r['compound_internal_id'] not in comp_rows:continue
   if fid in form_rows and form_rows[fid].get('exact_inchikey')!=r.get('exact_inchikey'):form_conflict.append(fid)
   else:form_rows[fid]=r
valid_comp=set(comp_rows)-set(comp_conflict);valid_form=set(form_rows)-set(form_conflict)
with op(OUT/'compound_master_positive_negative_union_v7.tsv.gz','wt') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(comp_rows[x] for x in sorted(valid_comp))
with op(OUT/'compound_form_positive_negative_union_v7.tsv.gz','wt') as f:w=csv.DictWriter(f,fieldnames=fform,delimiter='\t');w.writeheader();w.writerows(form_rows[x] for x in sorted(valid_form))

pub=OUT/'negative_evidence_v7.tsv.gz';rev=OUT/'negative_evidence_frozen_v7.tsv.gz';np=nr=0;public_pairs=set();missing_comp=missing_form=0
with op(NEG) as fi:
 rd=csv.DictReader(fi,delimiter='\t');ff=list(rd.fieldnames)+['v7_disposition','v7_rule']
 with op(pub,'wt') as fp,op(rev,'wt') as fr:
  wp=csv.DictWriter(fp,fieldnames=ff,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=ff,delimiter='\t');wp.writeheader();wr.writeheader()
  for r in rd:
   if r['target_uniprot_id'] not in proteins:continue
   nid=r.get('negative_evidence_id_v62') or r['source_evidence_id'];cid=r.get('compound_internal_id');fid=r.get('compound_form_id')
   ok=nid in eligible and cid in valid_comp and (not fid or fid in valid_form)
   if ok:r['v7_disposition']='PUBLIC_NEGATIVE';r['v7_rule']='canonical compound and optional form FKs resolved from baseline+negative append';wp.writerow(r);np+=1;public_pairs.add((r['target_uniprot_id'],cid))
   else:
    if nid in eligible and cid not in valid_comp:missing_comp+=1
    elif nid in eligible and fid and fid not in valid_form:missing_form+=1
    r['v7_disposition']='FROZEN_REVIEW';r['v7_rule']='identity/release rule or compound/form FK unresolved';wr.writerow(r);nr+=1
conf=sorted(public_pairs & positive_pairs)
with op(OUT/'positive_negative_conflict_v7.tsv.gz','wt') as f:
 ff=['conflict_id','target_uniprot_id','compound_internal_id','conflict_type','automatic_interpretation','positive_retained','negative_retained'];w=csv.DictWriter(f,fieldnames=ff,delimiter='\t');w.writeheader()
 for i,p in enumerate(conf,1):w.writerow({'conflict_id':f'PNC-V7-{i:08d}','target_uniprot_id':p[0],'compound_internal_id':p[1],'conflict_type':'positive_and_negative_across_assay_contexts','automatic_interpretation':'context_dependent; preserve assay/concentration/endpoint','positive_retained':'1','negative_retained':'1'})
report={'public_negative_rows':np,'frozen_negative_rows':nr,'public_negative_pairs':len(public_pairs),'positive_negative_conflict_pairs':len(conf),'compound_master_union_rows':len(valid_comp),'compound_form_union_rows':len(valid_form),'compound_id_conflicts':len(set(comp_conflict)),'form_id_conflicts':len(set(form_conflict)),'eligible_but_missing_compound_fk':missing_comp,'eligible_but_missing_form_fk':missing_form,'source_compound_tables':['baseline small_molecule_master_v1_3','new_compounds_v1_3_append'],'source_form_tables':['baseline compound_form_hierarchy_v1_3','new_forms_v1_3_append']}
(OUT/'T30_T31_NEGATIVE_REBUILD_WITH_APPENDS_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
