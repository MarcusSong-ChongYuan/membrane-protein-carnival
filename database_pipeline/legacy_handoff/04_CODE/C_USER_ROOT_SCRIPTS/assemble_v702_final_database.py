import csv,gzip,hashlib,json,shutil
from datetime import datetime,timezone
from pathlib import Path
BASE=Path(r"D:\finale\15_MemPro_V7.0.1_rescue_final_20260814");W=Path(r"D:\finale\12_V7_final_completion_20260813");DEST=Path(r"D:\finale\16_MemPro_V7.0.2_final_20260814");DATA=DEST/'01_data';REV=DEST/'02_nonpublic_resolved_status';PRED=DEST/'03_prediction_layer';QA=DEST/'04_QA'
for d in (DATA,REV,PRED,QA):d.mkdir(parents=True,exist_ok=True)
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
for p in (BASE/'01_data').iterdir():
 if p.is_file() and not p.name.startswith('binding_site_coordinate_verified'):shutil.copy2(p,DATA/p.name)
for p in (BASE/'02_frozen_review').iterdir():
 if p.is_file() and not p.name.startswith('binding_site_coordinate_'):shutil.copy2(p,REV/p.name)
public_sources=[BASE/'01_data'/'binding_site_coordinate_verified_v701.tsv.gz',W/'03_binding_sites'/'binding_site_coordinate_rescued_uniprot_v703.tsv.gz',W/'03_binding_sites'/'binding_site_coordinate_rescued_partial_v705.tsv.gz',W/'03_binding_sites'/'binding_site_coordinate_rescued_pocket_mapped_v706.tsv.gz']
fields=[];rows=[];seen=set()
for p in public_sources:
 with op(p) as f:
  rd=csv.DictReader(f,delimiter='\t')
  for x in rd.fieldnames:
   if x not in fields:fields.append(x)
  for r in rd:
   if r['binding_site_instance_id'] not in seen:seen.add(r['binding_site_instance_id']);rows.append(r)
with op(DATA/'binding_site_coordinate_verified_v702.tsv.gz','wt') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
remfields=[];remrows=[];remseen=set()
with op(W/'03_binding_sites'/'binding_site_coordinate_remaining_after_partial_v705.tsv.gz') as f:
 rd=csv.DictReader(f,delimiter='\t');remfields.extend(rd.fieldnames)
 for r in rd:
  if r.get('uniprot_rescue_status_v703')!='REMAIN_NO_RESIDUE_ANNOTATION':remseen.add(r['binding_site_instance_id']);remrows.append(r)
with op(W/'03_binding_sites'/'binding_site_coordinate_remaining_v704b.tsv.gz') as f:
 rd=csv.DictReader(f,delimiter='\t')
 for x in rd.fieldnames:
  if x not in remfields:remfields.append(x)
 for r in rd:
  if r['binding_site_instance_id'] not in remseen:remseen.add(r['binding_site_instance_id']);remrows.append(r)
with op(REV/'binding_site_coordinate_nonpublic_v702.tsv.gz','wt') as f:w=csv.DictWriter(f,fieldnames=remfields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(remrows)
for src in [W/'06_remaining_modules'/'isoform_topology_prediction_v707.tsv.gz',W/'06_remaining_modules'/'processed_form_topology_prediction_v707.tsv.gz',W/'06_remaining_modules'/'FORM_TOPOLOGY_PREDICTION_V707_REPORT.json']:shutil.copy2(src,PRED/src.name)
evidence=set()
with op(DATA/'protein_compound_evidence_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):evidence.add(r['evidence_id'])
ids=[];badfk=0
for r in rows:ids.append(r['binding_site_instance_id']);badfk+=r['evidence_id'] not in evidence
allids=set(ids)|remseen;pred=json.loads((W/'06_remaining_modules'/'FORM_TOPOLOGY_PREDICTION_V707_REPORT.json').read_text(encoding='utf-8'));checks={'public_site_pk':{'rows':len(ids),'bad':len(ids)-len(set(ids))},'public_site_evidence_fk':{'rows':len(ids),'bad':badfk},'site_partition':{'rows':len(allids),'bad':0 if len(allids)==55610 and not(set(ids)&remseen) else 1},'isoform_prediction_coverage':{'rows':pred['isoform']['input'],'bad':0 if pred['isoform']['input']==4342 else 1},'processed_prediction_coverage':{'rows':pred['processed_form']['input'],'bad':0 if pred['processed_form']['input']==1223 else 1}};blocking=sum(x['bad'] for x in checks.values())
counts={'binding_sites_total':len(allids),'coordinate_public':len(ids),'nonpublic_explicit_status':len(remseen),'author_complete_rescued':34611,'uniprot_complete_rescued':9005,'partial_medium_rescued':1552,'exact_cocrystal_pocket_rescued':663,'isoform_prediction_rows':4342,'processed_prediction_rows':1223};report={'status':'PASS' if blocking==0 else 'FAIL','blocking_errors':blocking,'checks':checks,'counts':counts,'nonpublic_policy':'not deleted; each record retains explicit technical/source reason','prediction_policy':'separate E3/prediction layer; not promoted to confirmed membrane class'};(QA/'V702_FINAL_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
manifest=[]
for p in sorted(DEST.rglob('*')):
 if not p.is_file() or p.name in {'MANIFEST.tsv','FREEZE_V702_FINAL.json'}:continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 manifest.append({'relative_path':str(p.relative_to(DEST)),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
with (DEST/'MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['relative_path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(manifest)
freeze={'release':'MemPro V7.0.2 final deterministic rescue','created_at_utc':datetime.now(timezone.utc).isoformat(),'status':'FROZEN_PASS' if blocking==0 else 'NOT_FROZEN_FAIL','blocking_errors':blocking,'manifest_entries':len(manifest),'counts':counts};(DEST/'FREEZE_V702_FINAL.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({'qa':report,'freeze':freeze},ensure_ascii=False,indent=2))
