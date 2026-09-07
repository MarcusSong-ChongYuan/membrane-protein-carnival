import csv,gzip,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path

BASE=Path(r'D:\finale\16_MemPro_V7.0.2_final_20260814');ROOT=Path(r'D:\finale\17_MemPro_V7.1_candidate_20260814');OUT=ROOT/'01_release_tables';ARC=ROOT/'02_archive';QA=ROOT/'04_QA'
def read(rel):
 p=BASE/rel;op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f,delimiter='\t')
def write(path,fields,rows):
 path.parent.mkdir(parents=True,exist_ok=True);n=0
 with gzip.open(path,'wt',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
  for r in rows:w.writerow(r);n+=1
 return n
def hid(*xs):return hashlib.sha256('|'.join(str(x or '') for x in xs).encode()).hexdigest()[:24]

def expression_rows(stats):
 for r in read('01_data/expression_measurement_v7.tsv.gz'):
  rawmap=r.get('v7_mapping_state',''); pmap=r.get('protein_mapping_status','')
  mapping='MAPPED' if rawmap=='mapped' else ('AMBIGUOUS' if 'ambig' in rawmap.lower() or pmap.startswith('E3') else 'UNMAPPED')
  rawdet=r.get('detection_status','')
  if rawdet=='detected':det='DETECTED';meas='MEASURED'
  elif rawdet=='not_detected':det='NOT_DETECTED';meas='MEASURED'
  elif rawdet=='below_threshold':det='BELOW_THRESHOLD';meas='MEASURED'
  else:det='NOT_APPLICABLE';meas='MISSING_OR_NOT_MEASURED_SOURCE_UNRESOLVED'
  layer=r.get('measurement_layer','');mtype=r.get('measurement_type','')
  assay='RNA_ABUNDANCE' if layer=='RNA' else ('IHC_PROTEIN_STAINING' if mtype=='IHC_level' else 'PROTEIN_MASS_SPECTROMETRY')
  context='NORMAL_TISSUE_OR_CELL'
  denom='1' if mapping=='MAPPED' and meas=='MEASURED' else '0'
  x={**r,'mapping_status_v71':mapping,'measurement_status_v71':meas,'detection_status_v71':det,'assay_semantics_v71':assay,'biological_context_v71':context,'mapped_measured_denominator_eligible_v71':denom,'missingness_limitation_v71':'source_combines_missing_and_not_measured' if meas.startswith('MISSING') else '','release_version_v71':'V7.1-candidate'}
  stats['mapping'][mapping]+=1;stats['measurement'][meas]+=1;stats['detection'][det]+=1;stats['assay'][assay]+=1;stats['denom'][denom]+=1
  yield x

def localization_rows(stats):
 for r in read('01_data/subcellular_localization_v7.tsv.gz'):
  mapping='MAPPED' if r.get('v7_mapping_state')=='mapped' else ('AMBIGUOUS' if r.get('protein_mapping_status','').startswith('E3') else 'UNMAPPED')
  direct='DIRECT_IF_OBSERVATION' if r.get('direct_observation_flag')=='1' else 'HPA_ANNOTATION_NOT_DIRECT'
  x={**r,'mapping_status_v71':mapping,'observation_semantics_v71':direct,'biological_context_v71':'NORMAL_CELL_LINE','isoform_specificity_v71':'NOT_ISOFORM_SPECIFIC','release_version_v71':'V7.1-candidate'}
  stats[mapping]+=1;yield x

def main():
 OUT.mkdir(parents=True,exist_ok=True);ARC.mkdir(parents=True,exist_ok=True);QA.mkdir(parents=True,exist_ok=True)
 # Positive core remains unchanged in row semantics; add explicit V7.1 release role.
 pos_fields=None
 def positives():
  nonlocal pos_fields
  for r in read('01_data/protein_compound_evidence_v7.tsv.gz'):
   if pos_fields is None:pos_fields=list(r.keys())+['release_role_v71','site_required_for_inclusion_v71','release_version_v71']
   yield {**r,'release_role_v71':'POSITIVE_INTERACTION_CORE','site_required_for_inclusion_v71':'0','release_version_v71':'V7.1-candidate'}
 # initialize generator once so fields are known
 pg=positives();first=next(pg)
 def pgen():yield first;yield from pg
 pc=write(OUT/'positive_interaction_evidence_v71.tsv.gz',pos_fields,pgen())
 pairs={(r['target_uniprot_id'],r['compound_internal_id']) for r in read('01_data/protein_compound_pair_v7.tsv.gz')}
 neg_fields=None;ctx=[]; # contextual rows are expected to be manageable
 conflict_pairs=defaultdict(lambda:{'count':0,'sources':set(),'outcomes':set()})
 archive_tmp=ARC/'archived_isolated_negative_evidence_v71.tsv.gz';archive_count=0;context_count=0;all_neg=0
 streams=[read('01_data/negative_evidence_v7.tsv.gz'),read('02_nonpublic_resolved_status/negative_evidence_frozen_v7.tsv.gz')]
 aw=None;af=None
 with gzip.open(archive_tmp,'wt',encoding='utf-8',newline='') as fo:
  for stream in streams:
   for r in stream:
    all_neg+=1
    if neg_fields is None:
     neg_fields=list(r.keys())+['negative_evidence_id_v71','negative_release_layer_v71','positive_pair_match_v71','conflict_interpretation_v71','release_version_v71'];af=neg_fields;aw=csv.DictWriter(fo,fieldnames=af,delimiter='\t',extrasaction='ignore');aw.writeheader()
    nid=r.get('negative_evidence_id_v62') or 'NEG-'+hid(r.get('source_evidence_id'),r.get('source_database'),r.get('source_record_id'))
    key=(r.get('target_uniprot_id',''),r.get('compound_internal_id',''));matched=bool(key[0] and key[1] and key in pairs)
    x={**r,'negative_evidence_id_v71':nid,'negative_release_layer_v71':'CONTEXTUAL_PUBLIC' if matched else 'ISOLATED_ARCHIVE','positive_pair_match_v71':'1' if matched else '0','conflict_interpretation_v71':'SAME_CANONICAL_PAIR_CONTEXT_REQUIRES_ASSAY_COMPARISON' if matched else 'NO_POSITIVE_CANONICAL_PAIR_IN_V71','release_version_v71':'V7.1-candidate'}
    if matched:
     ctx.append(x);context_count+=1;g=conflict_pairs[key];g['count']+=1;g['sources'].add(r.get('source_database',''));g['outcomes'].add(r.get('activity_outcome',''))
    else:aw.writerow(x);archive_count+=1
 if context_count:
  write(OUT/'contextual_negative_evidence_v71.tsv.gz',neg_fields,ctx)
 else:
  write(OUT/'contextual_negative_evidence_v71.tsv.gz',neg_fields or [],[])
 def conflicts():
  for i,(k,g) in enumerate(sorted(conflict_pairs.items()),1):
   yield {'conflict_id':f'PNC-V71-{i:08d}','target_uniprot_id':k[0],'compound_internal_id':k[1],'negative_evidence_count':g['count'],'negative_sources':';'.join(sorted(g['sources'])),'conflict_type':'POSITIVE_AND_INACTIVE_EVIDENCE_SAME_CANONICAL_PAIR','automatic_interpretation':'CONTEXT_DEPENDENT_UNTIL_ASSAY_CONDITIONS_MATCHED','positive_retained':'1','negative_retained_contextually':'1','release_version':'V7.1-candidate'}
 cc=write(OUT/'evidence_conflict_v71.tsv.gz',['conflict_id','target_uniprot_id','compound_internal_id','negative_evidence_count','negative_sources','conflict_type','automatic_interpretation','positive_retained','negative_retained_contextually','release_version'],conflicts())
 estats={'mapping':Counter(),'measurement':Counter(),'detection':Counter(),'assay':Counter(),'denom':Counter()}
 eg=expression_rows(estats);firste=next(eg);ef=list(firste.keys());ec=write(OUT/'expression_measurement_v71.tsv.gz',ef,(x for x in [firste]));
 # append remaining rows without recompressing through a second file
 # Re-run once to preserve simple deterministic streaming and overwrite complete table.
 estats={'mapping':Counter(),'measurement':Counter(),'detection':Counter(),'assay':Counter(),'denom':Counter()};eg=expression_rows(estats);firste=next(eg);ef=list(firste.keys());ec=write(OUT/'expression_measurement_v71.tsv.gz',ef,(x for x in [firste]+list(eg)))
 lstats=Counter();lg=localization_rows(lstats);firstl=next(lg);lf=list(firstl.keys());lc=write(OUT/'subcellular_localization_v71.tsv.gz',lf,(x for x in [firstl]+list(lg)))
 (ROOT/'00_docs'/'HPA_V71_INTERPRETATION.md').write_text('# HPA V7.1 interpretation\n\nDetection denominators use only MAPPED + MEASURED records. RNA abundance, IHC protein staining and mass-spectrometry measurements are separate. The current frozen HPA inputs represent normal tissue/cell contexts; no disease-tissue expression matrix is inferred. Source status `missing_or_not_measured` cannot be safely split and remains explicit.\n',encoding='utf-8')
 report={'positive_evidence':pc,'all_negative_evidence':all_neg,'contextual_negative_evidence':context_count,'isolated_negative_archive':archive_count,'conflict_pairs':cc,'negative_partition_pass':context_count+archive_count==all_neg,'expression_rows':ec,'expression':{k:dict(v) for k,v in estats.items()},'localization_rows':lc,'localization_mapping':dict(lstats),'pass':pc==942455 and all_neg==1936309 and context_count+archive_count==all_neg and ec==3046789 and lc==19091}
 (QA/'V71_STAGE4_EXPRESSION_NEGATIVE_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
