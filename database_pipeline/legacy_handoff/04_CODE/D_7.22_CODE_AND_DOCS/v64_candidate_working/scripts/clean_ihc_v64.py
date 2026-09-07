import csv,gzip,collections,json,hashlib
p=r'D:\finale\01_正式数据_V6.2\protein_tissue_cell_ihc_expression_v2.tsv.gz';out=r'D:\7.22\v64_candidate_working\candidate_tables\protein_tissue_cell_ihc_expression_v6_4_candidate.tsv.gz';seen=set();old_sig={};counts=collections.Counter()
with gzip.open(p,'rt',encoding='utf-8-sig',newline='') as f,gzip.open(out,'wt',encoding='utf-8',newline='') as g:
 rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;fields=base+['source_expression_record_id_v62','ihc_identity_status_v64','release_version_v64'];w=csv.DictWriter(g,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for r in rd:
  old=r['expression_record_id'];content='\x1f'.join(r[x] for x in base if x!='expression_record_id');key=hashlib.sha256(content.encode()).hexdigest()
  if key in seen:counts['exact_duplicate_excluded']+=1;continue
  seen.add(key);status='original_unique';new=old
  if old in old_sig and old_sig[old]!=key:
   new='IHC64_'+key[:20].upper();status='same_legacy_id_distinct_content_rekeyed';counts['content_distinct_rekeyed']+=1
  else:old_sig[old]=key
  r['source_expression_record_id_v62']=old;r['expression_record_id']=new;r['ihc_identity_status_v64']=status;r['release_version_v64']='MemPro V6.4 candidate';w.writerow(r);counts['rows_output']+=1
print(dict(counts))
