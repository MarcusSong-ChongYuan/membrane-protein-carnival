from pathlib import Path
import csv,gzip,json,collections,numpy as np
V=Path(r'D:\finale\01_正式数据_V6.2');R=Path(r'D:\finale\02_V6.3_candidate_20260803');O=Path(r'D:\7.22\v64_candidate_working');S=O/'stats'
def write(n,h,rows):
 with open(S/n,'w',encoding='utf-8',newline='') as f:w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(h);w.writerows(rows)
# protein breadth from V63 master
T=[];C=[];I=[];loc=collections.Counter();srcd=collections.Counter()
with gzip.open(R/'05_protein_annotation'/'human_membrane_protein_master_v6_3_candidate.tsv.gz','rt',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  for arr,k in [(T,'hpa_tissue_detected_count_v62'),(C,'hpa_cell_type_detected_count_v62'),(I,'hpa_ihc_detected_tissue_count_v62')]:
   try:arr.append(float(r.get(k) or 0))
   except:arr.append(0)
  for k in ['hpa_main_locations_v62','hpa_additional_locations_v62','hpa_extracellular_locations_v62']:
   for x in (r.get(k) or '').split(';'):
    if x.strip():loc[x.strip()]+=1
write('expression_breadth_summary_v64.tsv',['layer','n','median','p25','p75','zero_count'],[(n,len(a),np.median(a),np.percentile(a,25),np.percentile(a,75),sum(x==0 for x in a)) for n,a in [('tissue_RNA_detected_count',T),('cell_type_RNA_detected_count',C),('IHC_detected_tissue_count',I)]])
write('subcellular_location_counts_v64.tsv',['location','protein_count'],loc.most_common())
# IHC candidate ordinal/detection/reliability
det=collections.Counter();val=collections.Counter();rel=collections.Counter();tissues=collections.Counter();cells=collections.Counter()
with gzip.open(O/'candidate_tables'/'protein_tissue_cell_ihc_expression_v6_4_candidate.tsv.gz','rt',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):det[r['detection_status']]+=1;val[r['value']]+=1;rel[r['reliability']]+=1;tissues[r['tissue']]+=1;cells[r['cell_type']]+=1
write('ihc_detection_counts_v64.tsv',['detection_status','records'],det.most_common());write('ihc_reliability_counts_v64.tsv',['reliability','records'],rel.most_common())
# site types/sources/readiness
st=collections.Counter();ss=collections.Counter();spec=collections.Counter();pdbq=collections.Counter()
with gzip.open(O/'candidate_tables'/'binding_site_instances_v6_4_candidate.tsv.gz','rt',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):st[r['site_type'] or 'missing']+=1;ss[r['source_database'] or 'missing']+=1;spec[r['site_compound_specificity'] or 'missing']+=1;pdbq[r['pdb_identity_status_v64']]+=1
write('binding_site_type_counts_v64.tsv',['site_type','site_rows'],st.most_common());write('binding_site_source_counts_v64.tsv',['source_database','site_rows'],ss.most_common());write('binding_site_specificity_counts_v64.tsv',['specificity','site_rows'],spec.most_common())
# disease sources
src=collections.Counter()
with open(R/'04_disease'/'protein_disease_relation_v2_1.tsv',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  for x in (r.get('source_databases') or '').split(';'):
   if x.strip():src[x.strip()]+=1
write('disease_source_counts_v64.tsv',['source_database','protein_disease_pairs'],src.most_common())
print('supplementary stats complete')
