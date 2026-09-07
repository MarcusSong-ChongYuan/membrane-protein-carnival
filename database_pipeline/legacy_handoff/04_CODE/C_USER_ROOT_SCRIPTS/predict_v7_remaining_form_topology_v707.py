import csv,gzip,json
from collections import Counter
from pathlib import Path
ROOT=Path(r"D:\finale\15_MemPro_V7.0.1_rescue_final_20260814");OUT=Path(r"D:\finale\12_V7_final_completion_20260813\06_remaining_modules");KD={'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5};HYD=set('IVLFCMAGWY');CHG=set('DEKR')
def segments(seq,w=19):
 hits=[]
 for i in range(max(0,len(seq)-w+1)):
  s=seq[i:i+w]
  if sum(KD.get(a,0) for a in s)/w>=1.6 and sum(a in HYD for a in s)>=13 and sum(a in CHG for a in s)<=2:hits.append((i+1,i+w))
 merged=[]
 for a,b in hits:
  if merged and a<=merged[-1][1]+2:merged[-1]=(merged[-1][0],max(merged[-1][1],b))
  else:merged.append((a,b))
 return merged
proteins={}
with gzip.open(ROOT/'01_data'/'protein_master_v7.tsv.gz','rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):proteins[r['canonical_uniprot_accession']]=r['membrane_class_v7']
def run(inp,out,kind,idfield,accfield,seqfield):
 rows=[];cnt=Counter()
 with gzip.open(inp,'rt',encoding='utf-8',newline='') as f:
  rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames
  for r in rd:
   seq=r.get(seqfield,'');seg=segments(seq);internal=[x for x in seg if not (x[0]<=10 and x[1]<=40)];parent=proteins.get(r[accfield],'')
   if internal:status='PREDICTED_TM_FORM_E3';pred='A';basis='internal hydrophobic segment(s)'
   elif seg:status='N_TERMINAL_SIGNAL_OR_TM_AMBIGUOUS';pred='UNRESOLVED';basis='only N-terminal hydrophobic segment; signal peptide cannot be excluded'
   elif parent=='A':status='PREDICTED_TM_LOSS_FORM';pred='NON_TM_PREDICTED';basis='canonical A but no high-confidence hydrophobic segment in form sequence'
   else:status='NO_TM_PREDICTED_BC_MECHANISM_UNRESOLVED';pred=parent or 'UNRESOLVED';basis='absence of TM does not distinguish lipid-anchor/peripheral mechanisms'
   x=dict(r);x.update({'topology_prediction_status_v707':status,'predicted_membrane_class_v707':pred,'predicted_tm_segments_v707':';'.join(f'{a}-{b}' for a,b in seg),'topology_prediction_method_v707':'Kyte-Doolittle 19-aa window; mean>=1.6; hydrophobic>=13; charged<=2','topology_prediction_basis_v707':basis,'prediction_release_scope_v707':'PREDICTION_ONLY_NOT_CONFIRMED_PUBLIC_CLASS'});rows.append(x);cnt[status]+=1
 with gzip.open(out,'wt',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=base+['topology_prediction_status_v707','predicted_membrane_class_v707','predicted_tm_segments_v707','topology_prediction_method_v707','topology_prediction_basis_v707','prediction_release_scope_v707'],delimiter='\t');w.writeheader();w.writerows(rows)
 return {'input':len(rows),'status':dict(cnt)}
iso=run(ROOT/'02_frozen_review'/'isoform_remaining_v701.tsv.gz',OUT/'isoform_topology_prediction_v707.tsv.gz','isoform','protein_isoform_entity_id','canonical_uniprot_accession','isoform_sequence');pro=run(ROOT/'02_frozen_review'/'processed_form_remaining_v701.tsv.gz',OUT/'processed_form_topology_prediction_v707.tsv.gz','processed','processed_form_id','target_uniprot_id','sequence');report={'isoform':iso,'processed_form':pro,'policy':'prediction layer only; never promoted to E1/E2 or confirmed A/B/C without source evidence'};(OUT/'FORM_TOPOLOGY_PREDICTION_V707_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
