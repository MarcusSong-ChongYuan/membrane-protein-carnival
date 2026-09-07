import csv,gzip,hashlib,json,shutil
from datetime import datetime,timezone
from pathlib import Path
V7=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813");W=Path(r"D:\finale\12_V7_final_completion_20260813");DEST=Path(r"D:\finale\15_MemPro_V7.0.1_rescue_final_20260814");DATA=DEST/'01_data';REV=DEST/'02_frozen_review';QA=DEST/'03_QA'
for d in (DATA,REV,QA):d.mkdir(parents=True,exist_ok=True)
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
def merge_filtered(base,delta,out,idcol,base_status,delta_status):
 rows=[];fields=[];seen=set()
 for path,status in ((base,base_status),(delta,delta_status)):
  with op(path) as f:
   rd=csv.DictReader(f,delimiter='\t')
   for x in rd.fieldnames:
    if x not in fields:fields.append(x)
   for r in rd:
    if status and not r.get(status,'').startswith('PUBLIC'):continue
    if r[idcol] in seen:continue
    seen.add(r[idcol]);rows.append(r)
 with op(out,'wt') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
 return len(rows)
def copy_remaining(src,out):
 shutil.copy2(src,out)
 with op(src) as f:return sum(1 for _ in csv.DictReader(f,delimiter='\t'))
skip_data={'binding_site_coordinate_verified_v7.tsv.gz','protein_isoform_v7.tsv.gz','protein_processed_form_v7.tsv.gz'}
for p in (V7/'01_data').iterdir():
 if p.is_file() and p.name not in skip_data:shutil.copy2(p,DATA/p.name)
skip_rev={'binding_site_coordinate_frozen_v7.tsv.gz','isoform_frozen_review_v7.tsv.gz','processed_form_frozen_review_v7.tsv.gz'}
for p in (V7/'02_frozen_review').iterdir():
 if p.is_file() and p.name not in skip_rev:shutil.copy2(p,REV/p.name)
site_n=merge_filtered(V7/'01_data'/'binding_site_coordinate_verified_v7.tsv.gz',W/'03_binding_sites'/'binding_site_coordinate_rescued_v702.tsv.gz',DATA/'binding_site_coordinate_verified_v701.tsv.gz','binding_site_instance_id',None,'rescue_status_v702');site_rem=copy_remaining(W/'03_binding_sites'/'binding_site_coordinate_remaining_v702.tsv.gz',REV/'binding_site_coordinate_remaining_v701.tsv.gz')
iso_n=merge_filtered(V7/'01_data'/'protein_isoform_v7.tsv.gz',W/'06_remaining_modules'/'protein_isoform_v701.tsv.gz',DATA/'protein_isoform_v701.tsv.gz','protein_isoform_entity_id','isoform_disposition_v7','isoform_disposition_v701');iso_rem=copy_remaining(W/'06_remaining_modules'/'protein_isoform_v701.tsv.gz',REV/'isoform_remaining_v701.tsv.gz')
pro_n=merge_filtered(V7/'01_data'/'protein_processed_form_v7.tsv.gz',W/'06_remaining_modules'/'protein_processed_form_v701.tsv.gz',DATA/'protein_processed_form_v701.tsv.gz','processed_form_id','assignment_status','assignment_status_v701');pro_rem=copy_remaining(W/'06_remaining_modules'/'protein_processed_form_v701.tsv.gz',REV/'processed_form_remaining_v701.tsv.gz')
def ids(p,k):
 with op(p) as f:return [r[k] for r in csv.DictReader(f,delimiter='\t')]
proteins=set(ids(DATA/'protein_master_v7.tsv.gz','canonical_uniprot_accession'));evidence=set(ids(DATA/'protein_compound_evidence_v7.tsv.gz','evidence_id'));checks={}
for name,path,key,fk,parent in [('site',DATA/'binding_site_coordinate_verified_v701.tsv.gz','binding_site_instance_id','evidence_id',evidence),('isoform',DATA/'protein_isoform_v701.tsv.gz','protein_isoform_entity_id','canonical_uniprot_accession',proteins),('processed',DATA/'protein_processed_form_v701.tsv.gz','processed_form_id','target_uniprot_id',proteins)]:
 vals=[];bad=0
 with op(path) as f:
  for r in csv.DictReader(f,delimiter='\t'):vals.append(r[key]);bad+=r[fk] not in parent
 checks[name]={'rows':len(vals),'duplicate_primary_key':len(vals)-len(set(vals)),'bad_foreign_key':bad}
blocking=sum(x['duplicate_primary_key']+x['bad_foreign_key'] for x in checks.values());counts={'site_public_total':site_n,'site_remaining':site_rem,'isoform_public_total':iso_n,'isoform_remaining':iso_rem,'processed_public_total':pro_n,'processed_remaining':pro_rem};report={'status':'PASS' if blocking==0 else 'FAIL','blocking_errors':blocking,'checks':checks,'counts':counts,'site_rescue_breakdown':json.loads((W/'03_binding_sites'/'T26_T28_FROZEN_RESCUE_V702_REPORT.json').read_text(encoding='utf-8'))['status_counts'],'policy':'base PUBLIC rows only plus deterministic rescues; V7.0.0 immutable'};(QA/'V701_RESCUE_FINAL_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
manifest=[]
for p in sorted(DEST.rglob('*')):
 if not p.is_file() or p.name in {'MANIFEST.tsv','FREEZE_V701_FINAL.json'}:continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 manifest.append({'relative_path':str(p.relative_to(DEST)),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
with (DEST/'MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['relative_path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(manifest)
freeze={'release':'MemPro V7.0.1 deterministic rescue final','created_at_utc':datetime.now(timezone.utc).isoformat(),'status':'FROZEN_PASS' if blocking==0 else 'NOT_FROZEN_FAIL','blocking_errors':blocking,'manifest_entries':len(manifest),'counts':counts};(DEST/'FREEZE_V701_FINAL.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({'qa':report,'freeze':freeze},ensure_ascii=False,indent=2))
