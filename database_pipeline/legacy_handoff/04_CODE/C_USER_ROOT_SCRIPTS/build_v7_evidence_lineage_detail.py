import csv,gzip,hashlib,json,re
from collections import Counter
from pathlib import Path
root=Path(r"D:\finale\09_V7_data_freeze_working_20260813");src=root/'03_v7_relation_rebuild'/'public_binding_evidence_master_V7.tsv.gz';outdir=root/'05_v7_automatic_release_candidate';out=outdir/'protein_compound_evidence_V7.tsv.gz'
def key(prefix,*parts):return prefix+'-'+hashlib.sha1('|'.join(str(x or '').strip().lower() for x in parts).encode()).hexdigest()[:24].upper()
def modality(r):
 if r.get('pdb_ids') or r.get('evidence_type','').lower().find('struct')>=0:return 'structure'
 if 'pubchem' in r.get('source_database','').lower():return 'bioassay'
 if r.get('standard_value_nM') or r.get('activity_type') in {'Kd','Ki'}:return 'quantitative_binding'
 if r.get('activity_type') in {'IC50','EC50'}:return 'functional_pharmacology'
 return 'functional_or_curated_relation'
counts=Counter();seen=set()
with gzip.open(src,'rt',encoding='utf-8',newline='') as fi:
 r=csv.DictReader(fi,delimiter='\t');fields=list(r.fieldnames)+['experiment_lineage_key_v7','structure_lineage_key_v7','evidence_modality_v7','lineage_resolution_v7']
 with gzip.open(out,'wt',encoding='utf-8',newline='') as fo:
  w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
  for x in r:
   if x['evidence_id'] in seen:counts['duplicate_evidence_id']+=1
   seen.add(x['evidence_id']);mod=modality(x)
   if x.get('pdb_ids'):
    sk=key('SLIN',x['target_uniprot_id'],x['pdb_ids'],x.get('ligand_het_id'),x['compound_internal_id']);ek=sk;res='structure_key'
   elif 'pubchem' in x.get('source_database','').lower():
    sk='';ek=key('PLIN',x.get('source_record_id'),x.get('relationship_context_id'),x['target_uniprot_id'],x['compound_internal_id'],x.get('assay_or_mechanism'));res='pubchem_assay_proxy_key'
   elif x.get('pubmed_ids') or x.get('doi'):
    sk='';ek=key('ELIN',x.get('pubmed_ids'),x.get('doi'),x['target_uniprot_id'],x['compound_internal_id'],x.get('activity_type'),x.get('activity_value'),x.get('activity_unit'),x.get('assay_or_mechanism'));res='literature_experiment_proxy_key'
   else:
    sk='';ek=key('DLIN',x.get('source_database'),x.get('source_record_id'),x['target_uniprot_id'],x['compound_internal_id'],x.get('activity_type'));res='database_record_lineage_not_independent_experiment'
   x.update({'experiment_lineage_key_v7':ek,'structure_lineage_key_v7':sk,'evidence_modality_v7':mod,'lineage_resolution_v7':res});counts[res]+=1;counts['rows']+=1;w.writerow(x)
report={'counts':dict(counts),'unique_evidence_ids':len(seen),'blocking_duplicate_evidence_ids':counts['duplicate_evidence_id'],'lineage_semantics':'Keys are deterministic experimental-lineage proxies. Database record lineage is explicitly not counted as an independent experiment.'}
(outdir/'V7_EVIDENCE_LINEAGE_DETAIL_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
