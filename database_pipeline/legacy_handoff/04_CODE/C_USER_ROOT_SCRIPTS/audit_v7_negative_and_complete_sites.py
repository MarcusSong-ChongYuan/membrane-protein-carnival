import csv,gzip,json
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
REL=ROOT/'03_v7_relation_rebuild'; MOD=ROOT/'04_v7_remaining_modules'
NEG=Path(r"D:\finale\negative_upgrade\upgraded_negative_evidence.tsv.gz")
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
def vals(p,c):
 with op(p) as f:return {r[c] for r in csv.DictReader(f,delimiter='\t')}
proteins=vals(REL/'canonical_protein_entity_V7.tsv.gz','canonical_uniprot_accession')
compounds=vals(REL/'small_molecule_master_V7.tsv.gz','compound_internal_id')

stats=Counter(); combos=Counter(); eligible=[]
with op(NEG) as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['target_uniprot_id'] not in proteins:continue
  stats['protein_in_v7']+=1
  if r['compound_internal_id'] in compounds:stats['compound_in_v7']+=1
  else:continue
  stats['both_fk']+=1
  combos[(r.get('release_action_v62',''),r.get('identity_resolution_v62',''),r.get('canonical_relation_v62',''),r.get('compound_mapping_status',''))]+=1
  exact=bool(r.get('exact_inchikey_v62')) and bool(r.get('canonical_parent_inchikey_v62'))
  structure_consistent=(r.get('exact_inchikey_v62')==r.get('canonical_parent_inchikey_v62') or r.get('canonical_relation_v62') in {'canonical_parent','stereochemically_defined_parent','exact_form','baseline_form'})
  release_ok=r.get('release_action_v62')=='mapped_negative_release'
  if exact and structure_consistent and release_ok:
   eligible.append(r)
   stats['eligible_structure_unique_rule']+=1

top=[{'release_action':k[0],'identity_resolution':k[1],'canonical_relation':k[2],'compound_mapping_status':k[3],'count':v} for k,v in combos.most_common(30)]
audit={'counts':dict(stats),'top_field_combinations':top,'policy':'PUBLIC only when both FKs valid, exact and canonical parent InChIKeys present, structure relation consistent, and V6.2 release_action=mapped_negative_release'}
(MOD/'V7_NEGATIVE_ELIGIBILITY_AUDIT.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')

# Replace negative public/review outputs using audited deterministic rule.
pub=MOD/'negative_evidence_V7.tsv.gz'; rev=MOD/'negative_evidence_frozen_review_V7.tsv.gz'
with op(NEG) as fi:
 rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['v7_disposition','v7_disposition_reason']
 with op(pub,'wt') as fp,op(rev,'wt') as fr:
  wp=csv.DictWriter(fp,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wp.writeheader();wr.writeheader()
  np=nr=0;negpairs=set()
  for r in rd:
   if r['target_uniprot_id'] not in proteins:continue
   exact=bool(r.get('exact_inchikey_v62')) and bool(r.get('canonical_parent_inchikey_v62'))
   consistent=(r.get('exact_inchikey_v62')==r.get('canonical_parent_inchikey_v62') or r.get('canonical_relation_v62') in {'canonical_parent','stereochemically_defined_parent','exact_form','baseline_form'})
   ok=(r['compound_internal_id'] in compounds and r.get('release_action_v62')=='mapped_negative_release' and exact and consistent)
   if ok:r['v7_disposition']='PUBLIC_NEGATIVE';r['v7_disposition_reason']='both_FKs_valid;exact_structure_keys;canonical_relation_consistent;mapped_negative_release';wp.writerow(r);np+=1;negpairs.add((r['target_uniprot_id'],r['compound_internal_id']))
   else:r['v7_disposition']='FROZEN_REVIEW';r['v7_disposition_reason']='protein_valid_but_compound_or_structure_identity_not_uniquely_publishable';wr.writerow(r);nr+=1

# Complete membrane-side table with explicit unknown rows for every V7 site.
known={}
with op(MOD/'binding_site_membrane_side_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):known[r['binding_site_instance_id']]=r
with op(REL/'binding_site_instances_V7.tsv.gz') as fi:
 rd=csv.DictReader(fi,delimiter='\t')
 outfields=['binding_site_instance_id','target_uniprot_id','compound_internal_id','pdb_ids','pdb_chain_ids','source_database','residue_positions_uniprot','transmembrane_ranges_uniprot','residue_compartments','membrane_side','membrane_structure_support','coordinate_mapping_basis','orientation_basis','classification_scope']
 with op(MOD/'binding_site_membrane_side_complete_V7.tsv.gz','wt') as fo:
  w=csv.DictWriter(fo,fieldnames=outfields,delimiter='\t',extrasaction='ignore');w.writeheader();total=unknown=0
  for s in rd:
   total+=1
   if s['binding_site_instance_id'] in known:w.writerow(known[s['binding_site_instance_id']])
   else:
    unknown+=1;w.writerow({'binding_site_instance_id':s['binding_site_instance_id'],'target_uniprot_id':s['target_uniprot_id'],'compound_internal_id':s['compound_internal_id'],'pdb_ids':s['pdb_ids'],'pdb_chain_ids':s.get('pdb_chain_ids_v60',''),'source_database':s['source_database'],'residue_positions_uniprot':s.get('residue_or_site_description',''),'residue_compartments':'unknown','membrane_side':'unknown','membrane_structure_support':'none_or_not_mapped','coordinate_mapping_basis':'no_valid_orientation_mapping','orientation_basis':'not_inferred','classification_scope':'explicit_unknown_no_guess'})
audit['rebuilt_negative']={'public':np,'frozen_review':nr,'unique_public_negative_pairs':len(negpairs)}
audit['complete_binding_site_side']={'total':total,'classified':total-unknown,'unknown':unknown}
(MOD/'V7_NEGATIVE_AND_SITE_COMPLETION_REPORT.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
