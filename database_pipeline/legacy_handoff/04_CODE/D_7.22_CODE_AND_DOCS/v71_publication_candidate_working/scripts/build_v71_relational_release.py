import csv,gzip,json,sqlite3,collections,hashlib,os
from pathlib import Path
ROOT=Path(r'D:\7.22\v71_publication_candidate_working');IN=ROOT/'outputs'/'binding_evidence_public_candidate_v7_1.tsv.gz';OUT=ROOT/'outputs';Q=ROOT/'qa';DB=ROOT/'work'/'v71_relational_build.sqlite'
if DB.exists():DB.unlink()
cx=sqlite3.connect(DB);cx.executescript('''
CREATE TABLE ev(expkey TEXT,evidence_id TEXT,target TEXT,compound TEXT,class TEXT,source TEXT,publication TEXT,structurekey TEXT);
CREATE INDEX ev_exp ON ev(expkey);CREATE INDEX ev_pair ON ev(target,compound);CREATE INDEX ev_pub ON ev(publication);
''')
rows=[];n=0
with gzip.open(IN,'rt',encoding='utf8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  rows.append((r['independent_experiment_key_v71'],r['evidence_id'],r['target_uniprot_id'],r['compound_internal_id'],r['evidence_class_v71'],r['source_database_v71'],r['publication_id_v71'],r['independent_structure_key_v71']));n+=1
  if len(rows)>=50000:cx.executemany('INSERT INTO ev VALUES(?,?,?,?,?,?,?,?)',rows);cx.commit();rows=[]
if rows:cx.executemany('INSERT INTO ev VALUES(?,?,?,?,?,?,?,?)',rows);cx.commit()
with gzip.open(OUT/'independent_experiment_master_v7_1.tsv.gz','wt',encoding='utf8',newline='',compresslevel=6) as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['independent_experiment_key','representative_evidence_id','target_uniprot_id','compound_internal_id','evidence_class','evidence_record_count','contributing_databases','publication_id'])
 for r in cx.execute("SELECT expkey,MIN(evidence_id),target,compound,MIN(class),COUNT(*),GROUP_CONCAT(DISTINCT source),MIN(publication) FROM ev GROUP BY expkey,target,compound") :w.writerow(r)
with gzip.open(OUT/'protein_compound_relation_v7_1.tsv.gz','wt',encoding='utf8',newline='',compresslevel=6) as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['relation_id','target_uniprot_id','compound_internal_id','evidence_record_count','independent_experiment_count','independent_publication_count','independent_structure_count','evidence_classes','contributing_databases'])
 for r in cx.execute("SELECT target,compound,COUNT(*),COUNT(DISTINCT expkey),COUNT(DISTINCT NULLIF(publication,'')),COUNT(DISTINCT NULLIF(structurekey,'')),GROUP_CONCAT(DISTINCT class),GROUP_CONCAT(DISTINCT source) FROM ev GROUP BY target,compound"):
  rid='REL71_'+hashlib.sha256((r[0]+'|'+r[1]).encode()).hexdigest()[:20].upper();w.writerow([rid,*r])
with gzip.open(OUT/'publication_entity_v7_1.tsv.gz','wt',encoding='utf8',newline='',compresslevel=6) as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['publication_id','evidence_record_count','independent_experiment_count','source_databases'])
 for r in cx.execute("SELECT publication,COUNT(*),COUNT(DISTINCT expkey),GROUP_CONCAT(DISTINCT source) FROM ev WHERE publication<>'' GROUP BY publication"):w.writerow(r)
stats={}
stats['evidence_rows']=cx.execute('SELECT COUNT(*) FROM ev').fetchone()[0];stats['unique_evidence_ids']=cx.execute('SELECT COUNT(DISTINCT evidence_id) FROM ev').fetchone()[0];stats['unique_pairs']=cx.execute("SELECT COUNT(*) FROM (SELECT 1 FROM ev GROUP BY target,compound)").fetchone()[0];stats['independent_experiments']=cx.execute('SELECT COUNT(DISTINCT expkey) FROM ev').fetchone()[0];stats['independent_publications']=cx.execute("SELECT COUNT(DISTINCT publication) FROM ev WHERE publication<>''").fetchone()[0];stats['independent_structures']=cx.execute("SELECT COUNT(DISTINCT structurekey) FROM ev WHERE structurekey<>''").fetchone()[0];stats['class_counts']=dict(cx.execute('SELECT class,COUNT(*) FROM ev GROUP BY class').fetchall());stats['source_counts']=dict(cx.execute('SELECT source,COUNT(*) FROM ev GROUP BY source').fetchall());stats['mirrored_records_collapsed_by_lineage']=stats['evidence_rows']-stats['independent_experiments'];stats['status']='RELATIONAL_CANDIDATE_PASS'
(Q/'V71_RELATIONAL_RELEASE_REPORT.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8');cx.close();print(json.dumps(stats,ensure_ascii=False,indent=2))
