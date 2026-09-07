from __future__ import annotations

import csv, gzip, json
from collections import Counter
from pathlib import Path

BASE = Path(r"D:\finale\16_MemPro_V7.0.2_final_20260814")
OUT = Path(r"D:\finale\17_MemPro_V7.1_candidate_20260814\00_baseline\V702_MODULE_PROFILE.json")

def rows(rel):
    p=BASE/rel
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt',encoding='utf-8-sig',newline='') as f:
        yield from csv.DictReader(f,delimiter='\t')

def count_values(rel, fields):
    cs={f:Counter() for f in fields}; n=0
    for r in rows(rel):
        n+=1
        for f in fields: cs[f][r.get(f,'')]+=1
    return {'rows':n,'fields':{f:dict(c.most_common()) for f,c in cs.items()}}

def main():
    specs={
      'sites_public':('01_data/binding_site_coordinate_verified_v702.tsv.gz',['site_type','coordinate_qc_stratum_v7','coordinate_mapping_status_v7','coordinate_mapping_confidence_v7','actual_chain_used_v7','pocket_rescue_status_v704b','partial_rescue_status_v705']),
      'sites_nonpublic':('02_nonpublic_resolved_status/binding_site_coordinate_nonpublic_v702.tsv.gz',['site_type','coordinate_qc_stratum_v7','coordinate_mapping_status_v7','coordinate_mapping_confidence_v7','coordinate_mapping_reason_v7','actual_chain_used_v7']),
      'assignment':('01_data/evidence_target_entity_assignment_v7.tsv.gz',['target_entity_type_v7','assignment_resolution_v7','assignment_disposition_v7','historical_isoform_specificity','identity_review_flag_source']),
      'negative_public':('01_data/negative_evidence_v7.tsv.gz',['positive_negative_conflict_flag_v62','positive_negative_conflict_type_v62','v7_disposition','v7_rule','activity_outcome']),
      'negative_frozen':('02_nonpublic_resolved_status/negative_evidence_frozen_v7.tsv.gz',['positive_negative_conflict_flag_v62','positive_negative_conflict_type_v62','v7_disposition','v7_rule','activity_outcome']),
      'expression':('01_data/expression_measurement_v7.tsv.gz',['source_dataset','measurement_layer','measurement_type','detection_status','protein_mapping_status','v7_expression_layer','v7_mapping_state','tissue_mapping_status','cell_type_mapping_status']),
      'localization':('01_data/subcellular_localization_v7.tsv.gz',['source_dataset','hpa_reliability','direct_observation_flag','protein_mapping_status','v7_mapping_state']),
      'complexes':('01_data/protein_complex_v7.tsv.gz',['complex_class','membrane_relevance','stoichiometry_status','required_optional_subunit_status','complex_disposition_v7_final']),
      'components':('01_data/complex_component_v7.tsv.gz',['component_entity_type','stoichiometry_status','component_relation','required_subunit_status','optional_subunit_status','component_identity_review_flag_v0_4']),
      'classification':('01_data/protein_cross_classification_v7.tsv.gz',['classification_axis','classification_rule_version','source_release']),
      'disease':('01_data/protein_disease_relation_v7.tsv.gz',['best_evidence_level','canonical_mapping_statuses','mapping_review_flag','release_disposition_v63','disease_disposition_v7_final']),
    }
    out={k:count_values(*v) for k,v in specs.items()}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v['rows'] for k,v in out.items()},ensure_ascii=False))
if __name__=='__main__': main()
