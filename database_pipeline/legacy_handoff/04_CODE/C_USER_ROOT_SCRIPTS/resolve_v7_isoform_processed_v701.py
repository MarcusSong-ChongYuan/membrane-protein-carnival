import csv, gzip, json
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813")
WORK=Path(r"D:\finale\12_V7_final_completion_20260813\06_remaining_modules")
JSONDIR=Path(r"D:\finale\09_V7_data_freeze_working_20260813\02_protein_audit\uniprot_current_jsonl")

def op(p,m='rt'):
 return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
def pos(loc,key):
 x=loc.get(key,{}) or {};return x.get('value'),x.get('modifier','')

proteins={}
with op(ROOT/'01_data'/'protein_master_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):proteins[r['canonical_uniprot_accession']]=r
records={}
for p in JSONDIR.glob('batch_*.json'):
 d=json.loads(p.read_text(encoding='utf-8'))
 for r in d.get('results',[]):records[r.get('primaryAccession','')]=r

# Two-pass feature extraction: all membrane mechanisms are collected before any form is classified.
mechanisms={}
for acc,rec in records.items():
 seq=rec.get('sequence',{}).get('value',''); tm=[]; lipid=[]
 for ft in rec.get('features',[]):
  typ=ft.get('type','');s,sm=pos(ft.get('location',{}),'start');e,em=pos(ft.get('location',{}),'end')
  if not (s and e and sm=='EXACT' and em=='EXACT'):continue
  if typ in {'Transmembrane','Intramembrane'}:tm.append((typ,int(s),int(e),seq[int(s)-1:int(e)]))
  if typ in {'Lipidation','Glycosylphosphatidylinositol anchor'}:lipid.append((typ,int(s),int(e),seq[int(s)-1:int(e)],ft.get('description','')))
 mechanisms[acc]={'tm':tm,'lipid':lipid,'length':len(seq),'sequence':seq}

# Reclassify processed forms from their exact range and complete feature set.
processed=[];pc=Counter()
with op(ROOT/'02_frozen_review'/'processed_form_frozen_review_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  acc=r['target_uniprot_id'];m=mechanisms.get(acc,{});exact=r.get('sequence_range_status')=='exact'
  try:s=int(r['start']);e=int(r['end'])
  except:s=e=0;exact=False
  tm=[x for x in m.get('tm',[]) if exact and s<=x[1] and e>=x[2]]
  lip=[x for x in m.get('lipid',[]) if exact and s<=x[1] and e>=x[2]]
  parent=proteins.get(acc,{});full=exact and s==1 and e==m.get('length')
  if tm:
   cls='A';status='PUBLIC_EXACT_TM_OR_INTRAMEMBRANE_RETAINED';basis=';'.join(f'{x[0]}:{x[1]}-{x[2]}' for x in tm);confidence='DIRECT_UNIPROT_FEATURE_COORDINATE'
  elif lip:
   cls='B';status='PUBLIC_EXACT_LIPID_ANCHOR_RETAINED';basis=';'.join(f'{x[0]}:{x[1]}-{x[2]}:{x[4]}' for x in lip);confidence='DIRECT_UNIPROT_FEATURE_COORDINATE'
  elif full and parent:
   cls=parent['membrane_class_v7'];status='PUBLIC_FULL_LENGTH_FORM_INHERITS_CANONICAL';basis='processed range spans full canonical sequence';confidence='EXACT_SEQUENCE_RANGE'
  else:
   cls='UNRESOLVED_PROCESSED_FORM';status='FROZEN_SOURCE_INSUFFICIENT_FORM_SPECIFIC_MECHANISM';basis='exact range does not retain a coordinate-annotated TM/intramembrane/lipid-anchor feature';confidence='SOURCE_INSUFFICIENT'
  x=dict(r);x.update({'membrane_class_assignment_v701':cls,'assignment_status_v701':status,'assignment_basis_v701':basis,'assignment_confidence_v701':confidence});processed.append(x);pc[status]+=1

with op(WORK/'protein_processed_form_v701.tsv.gz','wt') as f:
 w=csv.DictWriter(f,fieldnames=list(processed[0]),delimiter='\t');w.writeheader();w.writerows(processed)

# Reclassify frozen isoforms using exact retention of any complete canonical membrane feature.
iso=[];ic=Counter()
with op(ROOT/'02_frozen_review'/'isoform_frozen_review_v7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  acc=r['canonical_uniprot_accession'];seq=r.get('isoform_sequence','');m=mechanisms.get(acc,{});parent=proteins.get(acc,{})
  tm=[];lip=[]
  for typ,s,e,frag in m.get('tm',[]):
   if frag and seq.count(frag)==1:tm.append((typ,s,e))
  for typ,s,e,frag,desc in m.get('lipid',[]):
   # Lipid anchors are single-residue features; require the local 15-aa context to be uniquely retained.
   cseq=m.get('sequence','');lo=max(0,s-8);hi=min(len(cseq),e+7);context=cseq[lo:hi]
   if context and len(context)>=8 and seq.count(context)==1:lip.append((typ,s,e,desc))
  if tm:
   cls='A';status='PUBLIC_EXACT_UNIQUE_TM_SEGMENT_TRANSFER';basis=';'.join(f'{x[0]}:{x[1]}-{x[2]}' for x in tm);confidence='EXACT_SEQUENCE_FEATURE_TRANSFER'
  elif lip:
   cls='B';status='PUBLIC_EXACT_LIPID_ANCHOR_CONTEXT_TRANSFER';basis=';'.join(f'{x[0]}:{x[1]}-{x[2]}:{x[3]}' for x in lip);confidence='EXACT_LOCAL_SEQUENCE_CONTEXT'
  elif seq and seq==m.get('sequence') and parent:
   cls=parent['membrane_class_v7'];status='PUBLIC_IDENTICAL_SEQUENCE_INHERITS_CANONICAL';basis='isoform sequence identical to canonical';confidence='EXACT_SEQUENCE_IDENTITY'
  else:
   cls='UNRESOLVED_ISOFORM';status='FROZEN_SOURCE_INSUFFICIENT_ISOFORM_SPECIFIC_MECHANISM';basis='isoform lacks an exactly transferable coordinate-annotated membrane feature';confidence='SOURCE_INSUFFICIENT'
  x=dict(r);x.update({'isoform_membrane_class_v701':cls,'isoform_disposition_v701':status,'isoform_membrane_basis_v701':basis,'isoform_membrane_confidence_v701':confidence});iso.append(x);ic[status]+=1

with op(WORK/'protein_isoform_v701.tsv.gz','wt') as f:
 w=csv.DictWriter(f,fieldnames=list(iso[0]),delimiter='\t');w.writeheader();w.writerows(iso)
report={'uniprot_records':len(records),'processed_input':len(processed),'processed_status':dict(pc),'processed_rescued':sum(v for k,v in pc.items() if k.startswith('PUBLIC')),'isoform_input':len(iso),'isoform_status':dict(ic),'isoform_rescued':sum(v for k,v in ic.items() if k.startswith('PUBLIC')),'policy':'direct coordinate or exact sequence-context transfer only; no hydropathy prediction promoted to confirmed class'}
(WORK/'ISOFORM_PROCESSED_V701_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
