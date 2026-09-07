from pathlib import Path
import csv,gzip,json,collections,hashlib,shutil
R=Path(r'D:\finale\02_V6.3_candidate_20260803');A=Path(r'D:\finale\04_V6.3.1_isoform_accession_audit_20260803');O=Path(r'D:\7.22\v64_candidate_working')
def readkeys(p,col,gz=False):
 op=gzip.open if gz else open;s=set();dup=0
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   v=r.get(col,'');dup+=v in s and bool(v);s.add(v) if v else None
 return s,dup
# core key sets
prot,dp=readkeys(R/'02_identity'/'canonical_protein_entity_v0_1.tsv.gz','canonical_protein_entity_id',True); up,_=readkeys(R/'02_identity'/'canonical_protein_entity_v0_1.tsv.gz','canonical_uniprot_accession',True); iso,di=readkeys(R/'02_identity'/'protein_isoform_v0_2.tsv.gz','protein_isoform_entity_id',True); gene,dg=readkeys(R/'02_identity'/'gene_entity_v0_1.tsv','gene_entity_id'); comp,dc=readkeys(R/'03_complex'/'complex_target_master_v0_2.tsv','complex_target_id'); disease,dd=readkeys(R/'04_disease'/'disease_entity_master_v0_1.tsv','canonical_disease_id')
report={'primary_key_duplicates':{'canonical_protein':dp,'protein_isoform':di,'gene':dg,'complex':dc,'disease':dd},'foreign_key_orphans':{}}
# generic FK checks
def check(p,checks,gz=False):
 op=gzip.open if gz else open;c=collections.Counter();samples=collections.defaultdict(list);rows=0
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   rows+=1
   for col,valid in checks:
    v=r.get(col,'').strip()
    if v and v not in valid:c[col]+=1;samples[col].append(v) if len(samples[col])<20 else None
 return {'rows':rows,'orphans':dict(c),'samples':dict(samples)}
report['foreign_key_orphans']['isoform']=check(R/'02_identity'/'protein_isoform_v0_2.tsv.gz',[('canonical_protein_entity_id',prot)],True)
report['foreign_key_orphans']['protein_gene_link']=check(R/'02_identity'/'canonical_protein_gene_link_v0_1.tsv',[('canonical_protein_entity_id',prot),('gene_entity_id',gene)])
report['foreign_key_orphans']['disease_relations']=check(R/'04_disease'/'protein_disease_relation_v2_1.tsv',[('target_uniprot_id',up),('canonical_disease_id',disease)])
report['foreign_key_orphans']['complex_components_v04']=check(R/'02_identity'/'complex_target_components_v0_4.tsv.gz',[('complex_target_id',comp)],True)
# apply v631 patch to component v04
patch={r['source_accession']:r for r in csv.DictReader(open(A/'complex_component_accession_patch_v631.tsv',encoding='utf-8-sig'),delimiter='\t')};src=R/'02_identity'/'complex_target_components_v0_4.tsv.gz';out=O/'candidate_tables'/'complex_target_components_v6_4_candidate.tsv.gz';counts=collections.Counter()
with gzip.open(src,'rt',encoding='utf-8-sig',newline='') as f,gzip.open(out,'wt',encoding='utf-8',newline='') as g:
 rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames+['source_accession_original_v64','accession_patch_status_v64','release_version_v64'];w=csv.DictWriter(g,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for r in rd:
  old=r.get('component_uniprot_isoform_id') or r.get('component_uniprot_id') or r.get('component_identifier');p=patch.get(old);r['source_accession_original_v64']='';r['accession_patch_status_v64']='unchanged'
  if p:
   r['source_accession_original_v64']=old;st=p['final_status'];r['accession_patch_status_v64']=st;counts[st]+=1
   if st=='resolved_high':
    if p['final_entity_type']=='protein_isoform':r['component_uniprot_isoform_id']=p['final_resolved_id'];r['component_uniprot_id']=p['final_canonical_parent_id'];r['component_entity_type']='protein_isoform'
    else:r['component_uniprot_id']=p['final_resolved_id'];r['component_uniprot_isoform_id']='';r['component_entity_type']='canonical_protein'
   elif st=='retain_ambiguous':r['component_identity_review_flag_v0_4']='1';r['component_reference_level_v0_4']='unresolved_gene_product_group'
   elif st=='exclude_pending_source_correction':r['component_identity_review_flag_v0_4']='1';r['component_reference_level_v0_4']='source_record_conflict'
  r['release_version_v64']='MemPro V6.4 candidate';w.writerow(r);counts['rows']+=1
report['v631_component_patch']=dict(counts);(O/'qa'/'V64_KEY_FK_ACCESSION_AUDIT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
