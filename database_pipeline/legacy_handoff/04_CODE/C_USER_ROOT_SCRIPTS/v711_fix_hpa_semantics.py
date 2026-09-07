import csv,gzip,json,hashlib,shutil,os,stat
from collections import Counter
from pathlib import Path

SRC=Path(r'D:\finale\18_MemPro_V7.1_final_20260814')
TMP=Path(r'D:\finale\19_MemPro_V7.1.1_candidate_20260814')
FINAL=Path(r'D:\finale\20_MemPro_V7.1.1_final_20260814')

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def writable(root):
 for p in root.rglob('*'):
  if p.is_file():p.chmod(stat.S_IREAD|stat.S_IWRITE)
def expected_layer(source):
 return {'HPA_tissue_RNA':'tissue_RNA','HPA_consensus_tissue_RNA':'tissue_RNA','HPA_single_cell_type_RNA':'cell_type_RNA','HPA_normal_tissue_IHC':'normal_tissue_IHC','HPA_tissue_MS':'tissue_protein_MS'}.get(source,'UNMAPPED_SOURCE_DATASET')
def main():
 if TMP.exists() or FINAL.exists():raise SystemExit('V7.1.1 candidate/final directory already exists')
 shutil.copytree(SRC,TMP);writable(TMP)
 p=TMP/'01_release_tables'/'expression_measurement_v71.tsv.gz';t=p.with_name('expression_measurement_v711.tsv.gz')
 n=0;changed=0;errors=Counter();layers=Counter();missing=Counter()
 with gzip.open(p,'rt',encoding='utf-8-sig',newline='') as fi:
  rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['legacy_v7_expression_layer_original','expression_layer_v711','measurement_status_v711','detection_status_v711','semantic_qc_v711','release_version_v711']
  with gzip.open(t,'wt',encoding='utf-8',newline='') as fo:
   w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
   for r in rd:
    n+=1;src=r['source_dataset'];layer=expected_layer(src);legacy=r.get('v7_expression_layer','')
    if legacy!=layer:changed+=1
    rawdet=r.get('detection_status','')
    if rawdet=='detected':ms='MEASURED';ds='DETECTED'
    elif rawdet=='not_detected':ms='MEASURED';ds='NOT_DETECTED'
    elif rawdet=='below_threshold':ms='MEASURED';ds='BELOW_THRESHOLD'
    else:ms='MISSING_OR_NOT_MEASURED_SOURCE_UNRESOLVED';ds='NOT_APPLICABLE'
    if src in {'HPA_tissue_RNA','HPA_consensus_tissue_RNA','HPA_single_cell_type_RNA'} and r.get('measurement_layer')!='RNA':errors['RNA_SOURCE_NOT_RNA_LAYER']+=1
    if src=='HPA_normal_tissue_IHC' and (r.get('measurement_layer')!='protein' or r.get('measurement_type')!='IHC_level'):errors['IHC_SEMANTIC_MISMATCH']+=1
    if src=='HPA_tissue_MS' and (r.get('measurement_layer')!='protein' or r.get('measurement_type')!='MS_intensity'):errors['MS_SEMANTIC_MISMATCH']+=1
    if src=='HPA_tissue_MS' and not r.get('value','').strip():missing['MS_EMPTY_VALUE']+=1
    if src=='HPA_tissue_MS' and not r.get('value','').strip() and (ms=='MEASURED' or ds in {'DETECTED','NOT_DETECTED'}):errors['EMPTY_MS_INTERPRETED_AS_MEASURED']+=1
    r.update({'legacy_v7_expression_layer_original':legacy,'v7_expression_layer':layer,'expression_layer_v711':layer,'measurement_status_v711':ms,'detection_status_v711':ds,'semantic_qc_v711':'PASS','release_version_v711':'V7.1.1'})
    w.writerow(r);layers[layer]+=1
 os.replace(t,p)
 report={'release':'V7.1.1','input_rows':n,'output_rows':n,'corrected_legacy_layer_rows':changed,'expression_layers':dict(layers),'missing_states':dict(missing),'semantic_errors':dict(errors),'blocking_errors':sum(errors.values()),'status':'PASS' if not errors and n==3046789 else 'FAIL'}
 q=TMP/'04_QA';q.mkdir(exist_ok=True);(q/'V711_HPA_SEMANTIC_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
 (TMP/'00_docs'/'V711_HPA_CORRECTION.md').write_text('''# V7.1.1 HPA semantic correction\n\nV7.1.1 corrects the inherited display-layer field for HPA tissue mass-spectrometry records. `HPA_tissue_MS` is now `tissue_protein_MS`, never `tissue_RNA`. The original value is preserved in `legacy_v7_expression_layer_original`. Empty MS intensity remains `MISSING_OR_NOT_MEASURED_SOURCE_UNRESOLVED` and is excluded from mapped+measured denominators; it is not interpreted as zero or not detected. RNA, IHC and MS cross-field consistency is now blocking QA.\n''',encoding='utf-8')
 if report['status']!='PASS':print(json.dumps(report,ensure_ascii=False));raise SystemExit(2)
 # Update final QA metadata without changing validated core counts/FKs.
 old=json.loads((q/'V71_FINAL_QA.json').read_text(encoding='utf-8'));old['release']='V7.1.1';old['status']='FROZEN_PASS';old['hpa_semantic_qa']='PASS';old['hpa_semantic_correction_rows']=changed
 (q/'V711_FINAL_QA.json').write_text(json.dumps(old,ensure_ascii=False,indent=2),encoding='utf-8')
 # Replace old release-level metadata and manifest.
 for name in ['MANIFEST.tsv','MANIFEST_VERIFICATION.json','FREEZE_V71_FINAL.json']:
  x=TMP/name
  if x.exists():x.unlink()
 rec=[]
 for x in sorted(TMP.rglob('*')):
  if x.is_file() and x.name not in {'MANIFEST.tsv','MANIFEST_VERIFICATION.json','FREEZE_V711_FINAL.json'}:rec.append({'relative_path':x.relative_to(TMP).as_posix(),'bytes':x.stat().st_size,'sha256':sha(x)})
 with (TMP/'MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['relative_path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(rec)
 ver={'checked':0,'missing':0,'mismatched':0,'status':'PASS'}
 for r in rec:
  x=TMP/r['relative_path'];ver['checked']+=1
  if not x.exists():ver['missing']+=1
  elif x.stat().st_size!=r['bytes'] or sha(x)!=r['sha256']:ver['mismatched']+=1
 if ver['missing'] or ver['mismatched']:ver['status']='FAIL'
 (TMP/'MANIFEST_VERIFICATION.json').write_text(json.dumps(ver,indent=2),encoding='utf-8')
 freeze={'release':'MemPro V7.1.1','status':'FROZEN_PASS','baseline':'MemPro V7.1','correction':'HPA expression-layer semantic consistency','corrected_rows':changed,'blocking_errors':0,'manual_sampling':'WAIVED_BY_USER','manifest_entries':len(rec),'manifest_verification':ver['status']}
 (TMP/'FREEZE_V711_FINAL.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2),encoding='utf-8')
 TMP.rename(FINAL)
 for x in FINAL.rglob('*'):
  if x.is_file():x.chmod(stat.S_IREAD)
 print(json.dumps({'final':str(FINAL),'report':report,'manifest':ver},ensure_ascii=False))
if __name__=='__main__':main()
