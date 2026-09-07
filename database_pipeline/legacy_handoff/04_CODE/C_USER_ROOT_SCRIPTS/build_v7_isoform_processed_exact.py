import csv,gzip,hashlib,json
from pathlib import Path

ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
REL=ROOT/'03_v7_relation_rebuild'; MOD=ROOT/'04_v7_remaining_modules'
JSONDIR=ROOT/'02_protein_audit'/'uniprot_current_jsonl'
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')

protein={}
with op(ROOT/'02_protein_audit'/'V7_PROTEIN_FREEZE_CANDIDATE'/'human_membrane_protein_master_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):protein[r['target_uniprot']] = r
records={}
for p in sorted(JSONDIR.glob('batch_*.json')):
 d=json.loads(p.read_text(encoding='utf-8'))
 for r in d.get('results',[]):records[r['primaryAccession']]=r

def pos(loc,key):
 x=loc.get(key,{})
 return x.get('value'),x.get('modifier','')
def evidences(feat):
 return ';'.join(sorted({('|'.join(filter(None,[e.get('evidenceCode',''),e.get('source',''),e.get('id','')]))) for e in feat.get('evidences',[])}))

processed=[]; local_segments={}
for acc,p in protein.items():
 rec=records.get(acc,{});seq=rec.get('sequence',{}).get('value','');segments=[]
 for ft in rec.get('features',[]):
  typ=ft.get('type','');loc=ft.get('location',{});start,sm=pos(loc,'start');end,em=pos(loc,'end')
  if typ in {'Transmembrane','Intramembrane'} and start and end and sm=='EXACT' and em=='EXACT':segments.append((typ,start,end,seq[start-1:end]))
  if typ not in {'Chain','Peptide','Propeptide','Signal','Transit peptide'}:continue
  exact=bool(start and end and sm=='EXACT' and em=='EXACT')
  contains=[s for s in segments if exact and start<=s[1] and end>=s[2]]
  if contains:memclass='A';status='PUBLIC_EXACT_MECHANISM_WITHIN_PROCESSED_RANGE'
  elif typ in {'Signal','Transit peptide','Propeptide'}:memclass='not_applicable';status='PROCESSING_FEATURE_NOT_INDEPENDENT_PROTEIN_ENTITY'
  else:memclass='UNRESOLVED_PROCESSED_FORM';status='FROZEN_REVIEW_NO_FORM_SPECIFIC_MEMBRANE_MECHANISM'
  raw=f"{acc}|{typ}|{start}|{end}|{ft.get('description','')}"
  processed.append({'processed_form_id':'PFORM-'+hashlib.sha1(raw.encode()).hexdigest()[:16].upper(),'target_uniprot_id':acc,'feature_id':ft.get('featureId',''),'processed_form_type':typ,'description':ft.get('description',''),'start':start or '','end':end or '','start_modifier':sm,'end_modifier':em,'sequence':seq[start-1:end] if exact else '','sequence_range_status':'exact' if exact else 'non_exact_or_missing','membrane_class_assignment':memclass,'assignment_status':status,'evidence_codes_and_sources':evidences(ft),'source_database':'UniProtKB current audit'})
 local_segments[acc]=segments
fields=list(processed[0])
with op(MOD/'protein_processed_form_exact_V7.tsv.gz','wt') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(processed)

# Conservative isoform status based on identical canonical sequence or exact,
# unique occurrence of a canonical TM/intramembrane segment in the isoform.
iso=[]
with op(REL/'protein_isoform_V7.tsv.gz') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  acc=r['canonical_uniprot_accession'];p=protein[acc];iseq=r.get('isoform_sequence','');cseq=records.get(acc,{}).get('sequence',{}).get('value','')
  same=r.get('same_sequence_as_canonical','').lower() in {'1','true','yes'} or (iseq and cseq and iseq==cseq)
  matched=[]
  if not same and iseq:
   for typ,s,e,frag in local_segments.get(acc,[]):
    if frag and iseq.count(frag)==1:matched.append(f'{typ}:{s}-{e}')
  if same:
   cls=p['final_class_v7_3'];disp='PUBLIC_INHERITED_IDENTICAL_SEQUENCE';basis='isoform sequence identical to canonical'
  elif matched:
   cls='A';disp='PUBLIC_EXACT_UNIQUE_TM_SEGMENT_TRANSFER';basis=';'.join(matched)
  else:
   cls='UNRESOLVED_ISOFORM';disp='FROZEN_REVIEW_NO_UNIQUE_FORM_SPECIFIC_MECHANISM';basis='no exact unique canonical membrane segment transfer; no inference from gene name'
  x=dict(r);x.update({'isoform_membrane_class_v7':cls,'isoform_disposition_v7':disp,'isoform_membrane_basis_v7':basis});iso.append(x)
if iso:
 fields=list(iso[0])
 with op(MOD/'protein_isoform_membrane_state_V7.tsv.gz','wt') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(iso)
report={'uniprot_records_loaded':len(records),'public_proteins':len(protein),'processed_feature_records':len(processed),'processed_public_A':sum(x['membrane_class_assignment']=='A' for x in processed),'processed_frozen_review':sum(x['assignment_status'].startswith('FROZEN') for x in processed),'isoform_records':len(iso),'isoform_public_identical':sum(x['isoform_disposition_v7']=='PUBLIC_INHERITED_IDENTICAL_SEQUENCE' for x in iso),'isoform_public_exact_tm_transfer':sum(x['isoform_disposition_v7']=='PUBLIC_EXACT_UNIQUE_TM_SEGMENT_TRANSFER' for x in iso),'isoform_frozen_review':sum(x['isoform_disposition_v7'].startswith('FROZEN') for x in iso)}
(MOD/'V7_ISOFORM_PROCESSED_EXACT_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
