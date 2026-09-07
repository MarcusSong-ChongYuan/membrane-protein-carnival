import csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813");OUT=ROOT/'05_hpa_expression';OUT.mkdir(exist_ok=True)
BASE=Path(r"D:\finale\02_V6.3_candidate_20260803\01_baseline_v6_2")
MASTER=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data\protein_master_v7.tsv.gz")
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
proteins=set()
with op(MASTER) as f:
 for r in csv.DictReader(f,delimiter='\t'):proteins.add(r['canonical_uniprot_accession'])
layers=[('tissue_RNA',BASE/'protein_tissue_expression_v2.tsv.gz'),('cell_type_RNA',BASE/'protein_cell_type_expression_v2.tsv.gz'),('IHC',BASE/'protein_tissue_cell_ihc_expression_v2.tsv.gz')]
counts={};covered=defaultdict(set);mapping=Counter();seen=set();dups=0
out=OUT/'expression_measurement_v7_complete.tsv.gz'
with op(out,'wt') as fo:
 writer=None
 for layer,p in layers:
  ni=no=0
  with op(p) as fi:
   rd=csv.DictReader(fi,delimiter='\t')
   if writer is None:
    fields=list(rd.fieldnames)+['v7_expression_layer','v7_mapping_state'];writer=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');writer.writeheader()
   for r in rd:
    ni+=1
    if r['target_uniprot_id'] not in proteins:continue
    eid=r['expression_record_id'];
    if eid in seen:dups+=1
    seen.add(eid);no+=1;covered[r['target_uniprot_id']].add(layer)
    ps=r.get('protein_mapping_status','')
    if ps.startswith('E1_') or ps.startswith('mapped'):state='mapped'
    elif ps:state='unmapped_or_ambiguous'
    else:state='missing_mapping_status'
    mapping[(layer,state)]+=1;r['v7_expression_layer']=layer;r['v7_mapping_state']=state;writer.writerow(r)
  counts[layer]={'input':ni,'v7_output':no}

# Localization remains a separate measurement modality but uses the same target FK.
locout=OUT/'subcellular_localization_v7_complete.tsv.gz';locn=0;loccovered=set();locmap=Counter()
with op(BASE/'protein_subcellular_localization_v2.tsv.gz') as fi:
 rd=csv.DictReader(fi,delimiter='\t');fields=list(rd.fieldnames)+['v7_mapping_state']
 with op(locout,'wt') as fo:
  w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
  for r in rd:
   if r['target_uniprot_id'] not in proteins:continue
   ps=r.get('protein_mapping_status','');state='mapped' if ps.startswith('E1_') or ps.startswith('mapped') else ('unmapped_or_ambiguous' if ps else 'missing_mapping_status')
   r['v7_mapping_state']=state;w.writerow(r);locn+=1;loccovered.add(r['target_uniprot_id']);locmap[state]+=1

# One row per protein distinguishes mapped data from absent records; absent is not
# interpreted as biological non-expression.
covout=OUT/'hpa_protein_coverage_state_v7.tsv';covcounts=Counter()
with covout.open('w',encoding='utf-8',newline='') as fo:
 fields=['target_uniprot_id','tissue_RNA_state','cell_type_RNA_state','IHC_state','subcellular_localization_state','interpretation'];w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t');w.writeheader()
 for p in sorted(proteins):
  states={x:('mapped_records_available' if x in covered[p] else 'missing_or_unmapped_no_record') for x,_ in layers};loc='mapped_records_available' if p in loccovered else 'missing_or_unmapped_no_record'
  row={'target_uniprot_id':p,**{x+'_state':states[x] for x,_ in layers},'subcellular_localization_state':loc,'interpretation':'Missing/unmapped is not biological zero expression'};w.writerow(row)
  for k,v in states.items():covcounts[(k,v)]+=1
  covcounts[('subcellular_localization',loc)]+=1

report={'proteins':len(proteins),'layer_counts':counts,'localization_records':locn,'duplicate_expression_record_ids':dups,'mapping_state_counts':{f'{a}|{b}':n for (a,b),n in mapping.items()},'localization_mapping_counts':dict(locmap),'protein_coverage_counts':{f'{a}|{b}':n for (a,b),n in covcounts.items()},'source_version':'HPA 25.1 inherited frozen raw standardized modules; re-filtered and re-QA against V7 master','manual_review':'waived; ambiguous/missing explicitly labeled'}
(OUT/'T36_T38_HPA_COMPLETE_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
