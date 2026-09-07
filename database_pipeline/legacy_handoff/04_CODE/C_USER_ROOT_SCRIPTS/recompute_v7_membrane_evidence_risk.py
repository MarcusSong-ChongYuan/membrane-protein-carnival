import csv,gzip,json,re
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813");OUT=ROOT/'06_remaining_modules';SRC=Path(r"D:\finale\09_V7_data_freeze_working_20260813\02_protein_audit\final_adjudication_v7_3_1\protein_identity_final_adjudication_all_v7_3_1.tsv")
def iv(x):
 try:return int(float(x or 0))
 except:return 0
rows=[];levels=Counter();risks=Counter();blocking=0
with SRC.open(encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  acc=r['target_uniprot'];cls=r['final_class_v7_3'];included=r['public_inclusion_v7_3']=='1';mode=r['final_mode_v7_3'];text=' | '.join([r.get('uniprot_transmembrane_features_current',''),r.get('uniprot_intramembrane_features_current',''),r.get('uniprot_lipidation_features_current',''),r.get('uniprot_subcellular_locations_current',''),r.get('final_evidence_v7_3','')])
  experimental=bool(re.search(r'ECO:0000269|ECO:0007744|PubMed:',text)) or 'VALIDATED' in r['final_decision_v7_3']
  structural=bool(r.get('opm_integral')=='1' or r.get('pdbtm_present')=='1' or 'structure_supported' in mode or 'beta_barrel' in mode)
  curated=iv(r.get('uniprot_transmembrane_count_current'))>0 or iv(r.get('uniprot_intramembrane_count_current'))>0 or iv(r.get('uniprot_lipidation_count_current'))>0 or any(x in mode for x in ['uniprot_peripheral','explicit_uniprot','covalent_lipid_anchor','immunoglobulin_heavy_chain_membrane_isoform'])
  if not included:level='E0';rule='E0_NOT_IN_PUBLIC_MEMBRANE_LAYER'
  elif experimental or structural:level='E1';rule='E1_DIRECT_EXPERIMENT_OR_STRUCTURE_SUPPORTED_MEMBRANE_MECHANISM'
  elif curated:level='E2';rule='E2_CURATED_UNIPROT_TOPOLOGY_ANCHOR_OR_PERIPHERAL_ANNOTATION'
  else:level='E3';rule='E3_INFERRED_OR_EXTERNAL_NONEXPERIMENTAL_MEMBRANE_MECHANISM'
  flags=[]
  if included and cls=='A' and not any(x in mode for x in ['transmembrane','intramembrane','integral','beta_barrel','membrane_isoform']):flags.append('A_WITHOUT_INTEGRAL_MODE')
  if included and cls=='B' and not any(x in mode for x in ['anchor','lipid']):flags.append('B_WITHOUT_ANCHOR_MODE')
  if included and cls=='C' and not any(x in mode for x in ['peripheral','complex_mediated']):flags.append('C_WITHOUT_PERIPHERAL_MODE')
  if r.get('cross_source_conflict')=='1':flags.append('CROSS_SOURCE_CONFLICT')
  if r.get('form_specificity_v7_2') and r.get('form_specificity_v7_2')!='not_applicable':flags.append('FORM_SPECIFIC_STATE')
  if not included:priority='P0' if r['final_decision_v7_3']=='EXCLUDED_NON_MEMBRANE' else 'P4';disp='EXCLUDED_AUDIT'
  elif any(x.endswith('WITHOUT_INTEGRAL_MODE') or x.endswith('WITHOUT_ANCHOR_MODE') or x.endswith('WITHOUT_PERIPHERAL_MODE') for x in flags):priority='P0';disp='FROZEN_RULE_FAILURE';blocking+=1
  elif 'CROSS_SOURCE_CONFLICT' in flags:priority='P1';disp='PUBLIC_WITH_CONFLICT_FLAG'
  elif level=='E3':priority='P2';disp='PUBLIC_E3_INFERRED'
  elif cls in {'B','C'} or 'FORM_SPECIFIC_STATE' in flags:priority='P3';disp='PUBLIC_MECHANISM_SPECIFIC'
  else:priority='P4';disp='PUBLIC_LOW_RISK'
  levels[level]+=1;risks[priority]+=1
  rows.append({'target_uniprot_id':acc,'membrane_class_v7':cls,'public_inclusion_v7':int(included),'membrane_evidence_level_v7_recomputed':level,'evidence_rule_triggered_v7':rule,'evidence_source_ids_v7':';'.join(filter(None,[r.get('external_integral_supports_v7',''),r.get('external_localization_supports_v7','')])), 'confidence_limitations_v7':';'.join(flags),'risk_priority_v7':priority,'automatic_disposition_v7':disp,'manual_review_status':'WAIVED_BY_USER'})
with gzip.open(OUT/'protein_membrane_evidence_and_risk_v7.tsv.gz','wt',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
report={'rows':len(rows),'evidence_levels':dict(levels),'risk_priorities':dict(risks),'blocking_rule_failures':blocking,'manual_review_status':'WAIVED_BY_USER','status':'PASS' if blocking==0 else 'FAIL'}
(OUT/'T07_T08_MEMBRANE_EVIDENCE_RISK_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
