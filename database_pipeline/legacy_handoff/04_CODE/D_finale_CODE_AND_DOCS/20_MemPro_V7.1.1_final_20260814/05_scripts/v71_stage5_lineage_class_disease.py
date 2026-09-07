import csv,gzip,json,hashlib,shutil,sqlite3,itertools
from collections import Counter
from pathlib import Path

BASE=Path(r'D:\finale\16_MemPro_V7.0.2_final_20260814');ROOT=Path(r'D:\finale\17_MemPro_V7.1_candidate_20260814');OUT=ROOT/'01_release_tables';QA=ROOT/'04_QA';WORK=ROOT/'03_work'
DISEASE=Path(r'D:\finale\02_V6.3_candidate_20260803\04_disease')
def readp(p):
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f,delimiter='\t')
def read(rel):yield from readp(BASE/rel)
def write(path,fields,rows):
 path.parent.mkdir(parents=True,exist_ok=True);n=0
 with gzip.open(path,'wt',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
  for r in rows:w.writerow(r);n+=1
 return n
def hid(*xs):return hashlib.sha256('|'.join(str(x or '') for x in xs).encode()).hexdigest()[:24]

def lineage():
 WORK.mkdir(parents=True,exist_ok=True);db=WORK/'lineage_v71.sqlite'
 if db.exists():db.unlink()
 con=sqlite3.connect(db);con.execute('pragma journal_mode=off');con.execute('pragma synchronous=off');con.execute('create table l(e text primary key, x text, s text, d text)')
 batch=[];n=0
 for r in read('01_data/protein_compound_evidence_v7.tsv.gz'):
  batch.append((r['evidence_id'],r.get('experiment_lineage_key_v7',''),r.get('structure_lineage_key_v7',''),r.get('source_database','')));n+=1
  if len(batch)>=50000:con.executemany('insert into l values(?,?,?,?)',batch);batch=[]
 if batch:con.executemany('insert into l values(?,?,?,?)',batch)
 con.commit();con.execute('create index ix on l(x)');con.execute('create index isx on l(s)')
 dupx={k:(c,dc) for k,c,dc in con.execute("select x,count(*),count(distinct d) from l where x<>'' group by x having count(*)>1")}
 dups={k:(c,dc) for k,c,dc in con.execute("select s,count(*),count(distinct d) from l where s<>'' group by s having count(*)>1")}
 fields=['evidence_id','source_database','source_version','source_record_id','pubmed_ids','doi','experiment_lineage_key','structure_lineage_key','evidence_modality','experiment_key_record_count','experiment_key_database_count','structure_key_record_count','structure_key_database_count','cross_database_experiment_mirror_status','cross_database_structure_mirror_status','independence_semantics','release_version']
 def gen():
  for r in read('01_data/protein_compound_evidence_v7.tsv.gz'):
   x=r.get('experiment_lineage_key_v7','');s=r.get('structure_lineage_key_v7','');xc,xd=dupx.get(x,(1 if x else 0,1 if x else 0));sc,sd=dups.get(s,(1 if s else 0,1 if s else 0))
   yield {'evidence_id':r['evidence_id'],'source_database':r.get('source_database',''),'source_version':r.get('source_version',''),'source_record_id':r.get('source_record_id',''),'pubmed_ids':r.get('pubmed_ids',''),'doi':r.get('doi',''),'experiment_lineage_key':x,'structure_lineage_key':s,'evidence_modality':r.get('evidence_modality_v7',''),'experiment_key_record_count':xc,'experiment_key_database_count':xd,'structure_key_record_count':sc,'structure_key_database_count':sd,'cross_database_experiment_mirror_status':'CROSS_DATABASE_MIRROR' if xd>1 else ('WITHIN_DATABASE_DUPLICATE_KEY' if xc>1 else 'NOT_DETECTED'),'cross_database_structure_mirror_status':'CROSS_DATABASE_MIRROR' if sd>1 else ('WITHIN_DATABASE_DUPLICATE_KEY' if sc>1 else 'NOT_DETECTED'),'independence_semantics':'lineage_key_based_not_claimed_independent_when_key_missing','release_version':'V7.1-candidate'}
 out=write(OUT/'evidence_lineage_v71.tsv.gz',fields,gen());con.close()
 return {'rows':out,'duplicate_experiment_keys':len(dupx),'duplicate_structure_keys':len(dups)}

def classifications():
 fields=['classification_assertion_id','target_uniprot_id','classification_axis','classification_label','primary_flag','classification_source','source_record_or_annotation_ids','source_version','mapping_method','classification_level','is_inferred','release_version']
 def source(axis,support):
  t=support.lower()
  names=[]
  for key,label in [('interpro','InterPro'),('pfam','Pfam'),('reactome','Reactome'),('gpcr','GPCRdb'),('iuphar','IUPHAR/BPS Guide to PHARMACOLOGY'),('tcdb','TCDB'),('go:','Gene Ontology'),('ec:','EC')]:
   if key in t:names.append(label)
  if axis=='membrane_role':names.append('MemPro-derived')
  if not names:names.append('UniProtKB/GO/InterPro/Reactome integrated annotation')
  return ';'.join(dict.fromkeys(names))
 def gen():
  for r in read('01_data/protein_cross_classification_v7.tsv.gz'):
   axis=r['classification_axis'];sup=r['supporting_source_or_ids'];src=source(axis,sup)
   yield {'classification_assertion_id':'CLS-'+hid(r['target_uniprot_id'],axis,r['classification_label'],sup),'target_uniprot_id':r['target_uniprot_id'],'classification_axis':axis,'classification_label':r['classification_label'],'primary_flag':r['primary_flag'],'classification_source':src,'source_record_or_annotation_ids':sup,'source_version':r['source_release'],'mapping_method':r['classification_rule_version'],'classification_level':'PRIMARY' if r['primary_flag']=='1' else 'ADDITIONAL','is_inferred':'1' if 'MemPro-derived' in src else '0','release_version':'V7.1-candidate'}
 n=write(OUT/'protein_cross_classification_v71.tsv.gz',fields,gen())
 # Preserve the one-row-per-protein display summary.
 sg=read('01_data/protein_cross_classification_summary_v7.tsv.gz');first=next(sg);sf=list(first.keys())+['release_version_v71'];sn=write(OUT/'protein_cross_classification_summary_v71.tsv.gz',sf,({**r,'release_version_v71':'V7.1-candidate'} for r in itertools.chain([first],sg)))
 return {'assertions':n,'summary_rows':sn}

def copy_disease_table(srcname,outname,extra=None):
 p=DISEASE/srcname;g=readp(p);first=next(g,None)
 if first is None:return write(OUT/outname,[],[])
 fields=list(first.keys())+(['release_version_v71'] if 'release_version_v71' not in first else [])
 def gen():
  for r in itertools.chain([first],g):yield {**r,'release_version_v71':'V7.1-candidate'}
 return write(OUT/outname,fields,gen())
def diseases():
 counts={}
 # Current V7 canonical master and public exact protein relations are authoritative.
 for rel,out in [('01_data/disease_master_v7.tsv.gz','disease_entity_master_v71.tsv.gz'),('01_data/protein_disease_relation_v7.tsv.gz','protein_disease_relation_v71.tsv.gz')]:
  g=read(rel);first=next(g);fields=list(first.keys())+['release_version_v71'];counts[out]=write(OUT/out,fields,({**r,'release_version_v71':'V7.1-candidate'} for r in itertools.chain([first],g)))
 for src,out in [('disease_source_entity_v0_1.tsv','disease_source_entity_v71.tsv.gz'),('disease_xref_mapping_v0_1.tsv','disease_xref_mapping_v71.tsv.gz'),('disease_hierarchy_v0_1.tsv','disease_hierarchy_v71.tsv.gz'),('disease_anatomy_v0_1.tsv','disease_anatomy_v71.tsv.gz'),('disease_anatomical_system_v0_1.tsv','disease_anatomical_system_v71.tsv.gz'),('disease_therapeutic_area_v0_1.tsv','disease_therapeutic_area_v71.tsv.gz')]:counts[out]=copy_disease_table(src,out)
 return counts

def source_registry():
 rows=list(read('01_data/source_registry_v7.tsv'))
 extra=[
 {'module':'V7.1 rules','source_database':'MemPro-derived','source_versions':'V7.1-candidate','record_scope':'release semantics and derived classifications','registry_basis':'approved user rules'},
 {'module':'disease ontology','source_database':'MONDO/DO/Uberon','source_versions':'MONDO 2026-07-06; frozen DO/Uberon snapshots','record_scope':'exact identity, hierarchy and multi-label anatomy','registry_basis':'frozen V6.3 ontology module reused'},
 {'module':'expression','source_database':'Human Protein Atlas','source_versions':'inherited per-record HPA source_version','record_scope':'normal tissue/cell RNA, IHC and protein MS','registry_basis':'per-record provenance'},]
 fields=list(rows[0].keys())+['release_version_v71'];return write(OUT/'source_registry_v71.tsv.gz',fields,({**r,'release_version_v71':'V7.1-candidate'} for r in rows+extra))

def main():
 OUT.mkdir(parents=True,exist_ok=True);QA.mkdir(parents=True,exist_ok=True)
 lr=lineage();cl=classifications();ds=diseases();sr=source_registry()
 report={'lineage':lr,'classification':cl,'disease':ds,'source_registry_rows':sr}
 report['pass']=lr['rows']==942455 and cl['assertions']==59715 and cl['summary_rows']==7800 and ds['disease_entity_master_v71.tsv.gz']==3739 and ds['protein_disease_relation_v71.tsv.gz']==6753
 (QA/'V71_STAGE5_LINEAGE_CLASS_DISEASE_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
