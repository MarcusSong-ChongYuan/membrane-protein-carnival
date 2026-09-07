import csv,gzip,collections,json,hashlib
p=r'D:\finale\01_正式数据_V6.2\protein_tissue_cell_ihc_expression_v2.tsv.gz';out=r'D:\7.22\v64_candidate_working\review_queues\ihc_duplicate_record_id_review_v64.tsv';seen={};dups=[];classes=collections.Counter()
with gzip.open(p,'rt',encoding='utf-8-sig',newline='') as f:
 rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames
 for r in rd:
  k=r['expression_record_id'];sig=tuple(r[x] for x in fields if x!='expression_record_id')
  if k in seen:
   same=sig==seen[k][0];classes['exact_duplicate_rows' if same else 'same_id_different_content']+=1
   if len(dups)<100000:dups.append([k,'exact_duplicate' if same else 'content_conflict',seen[k][1],r])
  else:seen[k]=(sig,r)
with open(out,'w',encoding='utf-8',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['expression_record_id','duplicate_class','first_record','duplicate_record'])
 for k,c,a,b in dups:w.writerow([k,c,json.dumps(a,ensure_ascii=False),json.dumps(b,ensure_ascii=False)])
print(dict(classes),'unique_ids',len(seen),'rows',len(seen)+sum(classes.values()))
