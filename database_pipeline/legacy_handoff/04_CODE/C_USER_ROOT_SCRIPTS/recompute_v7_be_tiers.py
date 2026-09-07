import csv,gzip,json,math,re
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813")
SRC=Path(r"D:\finale\11_MemPro_V7_automatic_data_release_20260813\01_data\protein_compound_evidence_v7.tsv.gz")
ASSIGN=ROOT/'01_entity_assignment'/'evidence_target_entity_assignment_v7_final.tsv.gz'
OUT=ROOT/'02_be_recompute';OUT.mkdir(exist_ok=True)
def op(p,m='rt'):return gzip.open(p,m,encoding='utf-8',newline='') if str(p).endswith('.gz') else p.open(m,encoding='utf-8',newline='')
assign={}
with op(ASSIGN) as f:
 for r in csv.DictReader(f,delimiter='\t'):assign[r['evidence_id']]=r

struct_types={'experimental_structure_residue_contact','compound_specific_structure_context','experimental_complex_with_affinity','compound_specific_structure'}
direct_types={'direct_quantitative_binding','pubchem_direct_quantitative_binding'}
direct_measures={'kd','ki','kb','ka','affinity constant','kd(app)','kd,app','kdapp','kd/ki'}
functional_measures={'ic50','ec50','ac50','ed50','gi50','potency','activity','inhibition','ic90','dc50','km','kinact','ke','inh','fIC50','fEC50','log ic50','log ki'}
valid_target={'single_human_uniprot','single_uniprot_mapping','single_protein_direct'}

def positive_number(x):
 try:return math.isfinite(float(x)) and float(x)>0
 except:return False
def has_anchor(r):return bool(r.get('source_record_id') or r.get('pubmed_ids') or r.get('doi') or r.get('pdb_ids'))

counts=Counter();transitions=Counter();reasons=Counter()
pub=OUT/'protein_compound_evidence_BE_recomputed_v7.tsv.gz';review=OUT/'evidence_BE_frozen_review_v7.tsv.gz'
with op(SRC) as fi:
 rd=csv.DictReader(fi,delimiter='\t');extra=['evidence_tier_recomputed_v7','be_rule_v7','be_disposition_v7','be_limitations_v7'];fields=list(rd.fieldnames)+extra
 with op(pub,'wt') as fp,op(review,'wt') as fr:
  wp=csv.DictWriter(fp,fieldnames=fields,delimiter='\t');wr=csv.DictWriter(fr,fieldnames=fields,delimiter='\t');wp.writeheader();wr.writeheader()
  for r in rd:
   a=assign.get(r['evidence_id'],{});etype=r.get('evidence_type','');direct=r.get('evidence_directness','');act=r.get('activity_type','').strip().lower();tier='';rule='';limit=[]
   safe_target=(a.get('assignment_disposition_v7','').startswith('PUBLIC') and r.get('target_assignment_status') in valid_target)
   if not safe_target:limit.append('target entity not uniquely public')
   if not has_anchor(r):limit.append('no source/PubMed/DOI/PDB anchor')
   structural=(etype in struct_types or 'structur' in direct) and bool(r.get('pdb_ids'))
   direct_quant=(etype in direct_types or direct in {'direct_quantitative','direct'}) and act in direct_measures and positive_number(r.get('standard_value_nM'))
   functional=(etype not in struct_types and (act in functional_measures or direct in {'functional_or_pharmacological','indirect_or_functional','curated_or_functional','curated_target_relation','binding_assay_activity'}))
   if safe_target and has_anchor(r) and structural:
    tier='BE1';rule='BE1_STRUCTURAL_DIRECT_PDB_CONTEXT'
   elif safe_target and has_anchor(r) and direct_quant:
    tier='BE2';rule='BE2_DIRECT_QUANTITATIVE_KD_KI_OR_AFFINITY_WITH_STANDARD_NM'
   elif safe_target and has_anchor(r) and functional:
    tier='BE3';rule='BE3_TARGET_SPECIFIC_FUNCTIONAL_PHARMACOLOGY_OR_CURATED_RELATION'
   elif safe_target and has_anchor(r) and etype in {'pubchem_target_specific_quantitative_activity','pubchem_target_specific_activity','functional_or_bioactivity_measurement','binding_assay_activity','curated_pharmacological_interaction','curated_pharmacological_target_assertion','target_specific_functional_pharmacology'}:
    tier='BE3';rule='BE3_TARGET_SPECIFIC_ACTIVITY_FALLBACK_NO_DIRECT_BINDING_CLAIM';limit.append('activity type not a direct Kd/Ki binding measure')
   else:
    if not limit:limit.append('evidence type/value combination does not satisfy BE1-BE3 publication rule')
   if tier:
    r.update({'evidence_tier_recomputed_v7':tier,'be_rule_v7':rule,'be_disposition_v7':'PUBLIC','be_limitations_v7':';'.join(limit)});wp.writerow(r);counts['public']+=1;counts[tier]+=1;transitions[(r['evidence_tier'],tier)]+=1
   else:
    r.update({'evidence_tier_recomputed_v7':'UNRESOLVED_BE','be_rule_v7':'NO_VALID_BE_RULE','be_disposition_v7':'FROZEN_REVIEW','be_limitations_v7':';'.join(limit)});wr.writerow(r);counts['frozen_review']+=1
    for x in limit:reasons[x]+=1

report={'input_rows':sum(counts[k] for k in ['public','frozen_review']),'counts':dict(counts),'tier_transitions':{f'{a}->{b}':n for (a,b),n in transitions.items()},'frozen_reasons':dict(reasons),'rules':{'BE1':'direct experimental structural context with PDB','BE2':'unique human target plus direct Kd/Ki/affinity-type quantitative value standardized to positive nM','BE3':'target-specific functional/pharmacological/curated relation; IC50/EC50/AC50 etc do not claim direct equilibrium affinity'},'manual_review':'waived; unresolved rows frozen'}
(OUT/'T24_BE_RECOMPUTE_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
