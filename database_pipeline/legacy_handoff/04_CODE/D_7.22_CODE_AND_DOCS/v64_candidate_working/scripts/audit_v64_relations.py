from pathlib import Path
import csv,gzip,json,collections,re
V=Path(r'D:\finale\01_正式数据_V6.2');O=Path(r'D:\7.22\v64_candidate_working');R=Path(r'D:\finale\02_V6.3_candidate_20260803')
# load master keys
prot=set();cpid=set();forms=set();disease=set()
with open(V/'human_membrane_protein_master_v6_2.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):prot.add(r['target_uniprot_id'])
with open(V/'small_molecule_master_v1_3.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):cpid.add(r['compound_internal_id'])
with open(V/'compound_form_hierarchy_v1_3.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):forms.add(r['compound_form_id'])
with open(R/'04_disease'/'disease_entity_master_v0_1.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):disease.add(r['canonical_disease_id'])
report={'master_counts':{'proteins':len(prot),'compounds':len(cpid),'forms':len(forms),'diseases':len(disease)},'tables':{}}
def scan(p,checks,gz=False,pk=None):
 op=gzip.open if gz else open;c=collections.Counter();samp=collections.defaultdict(list);seen=set();rows=0
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   rows+=1
   if pk:
    v=r.get(pk,'');c['duplicate_primary_key']+=bool(v and v in seen);seen.add(v)
   for col,valid in checks:
    v=r.get(col,'').strip()
    if v and v not in valid:c[f'orphan:{col}']+=1;samp[col].append(v) if len(samp[col])<20 else None
 return {'rows':rows,'counts':dict(c),'samples':dict(samp)}
report['tables']['binding_evidence_v64']=scan(O/'candidate_tables'/'binding_evidence_master_v6_4_candidate.tsv.gz',[('target_uniprot_id',prot),('compound_internal_id',cpid),('compound_form_id',forms)],True,'evidence_id')
report['tables']['binding_sites_v64']=scan(O/'candidate_tables'/'binding_site_instances_v6_4_candidate.tsv.gz',[('target_uniprot_id',prot),('compound_internal_id',cpid),('compound_form_id',forms)],True,'binding_site_instance_id')
report['tables']['disease_v63']=scan(R/'04_disease'/'protein_disease_relation_v2_1.tsv',[('target_uniprot_id',prot),('canonical_disease_id',disease)],False,'disease_relation_id_v2')
for name in ['protein_tissue_expression_v2.tsv.gz','protein_cell_type_expression_v2.tsv.gz','protein_tissue_cell_ihc_expression_v2.tsv.gz','protein_subcellular_localization_v2.tsv.gz']:
 report['tables'][name]=scan(V/name,[('target_uniprot_id',prot)],True,'expression_record_id' if 'localization' not in name else 'localization_record_id')
# source registry coverage from evidence
src=collections.Counter()
with gzip.open(O/'candidate_tables'/'binding_evidence_master_v6_4_candidate.tsv.gz','rt',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):src[r['source_database']]+=1
with open(O/'stats'/'binding_evidence_source_counts_v64.tsv','w',encoding='utf-8',newline='') as f:w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['source_database','evidence_rows']);w.writerows(src.most_common())
report['binding_evidence_sources']=dict(src);(O/'qa'/'V64_RELATIONAL_INTEGRITY_AUDIT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
