from pathlib import Path
import csv, json, shutil

ROOT=Path(r"D:\finale\17_MemPro_V7.1_candidate_20260814")
DOC=ROOT/'00_docs'; SCHEMA=ROOT/'00_schema'

RULES='''# MemPro V7.1 approved release rules

V7.0.2 is immutable. V7.1 is an incremental data-model and publication-semantics release.

## Interaction and site independence

A reliable human membrane-protein–small-molecule interaction does not require residue or coordinate data. Site status is represented independently as RESIDUE_COMPLETE, RESIDUE_PARTIAL, RESIDUE_SOURCE_ONLY, NO_RESIDUE_REPORTED, PREDICTED_POCKET, or NO_SITE_INFORMATION. Missing site coordinates never remove an otherwise valid interaction.

## Protein entity granularity

Evidence targets exactly one explicit level when supported: gene product unspecified, canonical protein, protein isoform, processed protein form, protein complex, or unresolved target. No isoform or subunit is selected from a name alone.

## Membrane classes

A = transmembrane or strongly integral; B = directly embedded or lipid/GPI anchored without a conventional transmembrane span; C = peripheral/state- or complex-mediated membrane association. Class is stored per canonical protein, isoform, processed form, and state. The most direct mechanism is primary; other mechanisms remain secondary.

## Structure and chain

Source-reported residues are immutable. Mapped and unmapped residues are separate. SIFTS residue mapping is preferred, followed by sequence-confirmed and segment-safe DBREF mapping. Multiple chains remain a candidate set unless the source or structure context selects one uniquely.

## Negative evidence

Positive evidence is the public core. Negative evidence is public-contextual only when it conflicts with or informs selectivity for a positive pair. Isolated inactive records are archived outside website, API, figure, and manuscript default statistics.

## Expression

Detection denominators use mapped and measured records only. NOT_DETECTED, UNMAPPED, NOT_MEASURED, and MISSING are distinct. RNA abundance and IHC protein staining are separate; normal and disease contexts are separate.

## Provenance

Distinct contributing databases, independent experiments, independent structures, PubChem assays, and evidence modalities are separate quantities. Database mirrors do not increase independent-experiment counts.

## Classification and disease

Five independent classification axes are retained with source/version/method: structural family, molecular function, biological process, membrane role, and specialist classification. Disease identity uses MONDO with official exact/equivalent mappings only; broad/narrow/related links do not merge identities. Anatomy is ontology-derived and multi-label.

## Manual validation

Manual sampling is WAIVED_BY_USER for this release. The release may report automated QA, but may not claim that all records were manually reviewed or quote a manual accuracy estimate.
'''

TABLES=[
('interaction_site_v71','site_record_id','One row per source/derived site assertion; interaction inclusion is independent of coordinates.'),
('gene_membrane_summary_v71','gene_entity_id','Gene-level summary of membrane-related protein forms.'),
('canonical_protein_membrane_class_v71','canonical_uniprot_accession','Canonical protein A/B/C classification.'),
('isoform_membrane_class_v71','protein_isoform_entity_id','Isoform-specific classification; predictions remain explicitly predicted.'),
('processed_form_membrane_class_v71','processed_form_id','Processed-form classification and membrane-retention state.'),
('state_dependent_membrane_association_v71','membrane_state_id','State/form-specific membrane association.'),
('protein_complex_v71','complex_target_id','Complex entity.'),
('complex_component_v71','component_assertion_id','Complex component assertion.'),
('complex_stoichiometry_v71','stoichiometry_assertion_id','Per-component stoichiometry.'),
('complex_assembly_v71','complex_assembly_id','Biological assembly assertion or explicit not-reported state.'),
('complex_state_v71','complex_state_id','State-dependent complex assembly.'),
('complex_binding_evidence_v71','complex_binding_candidate_id','Ligand evidence attached to complex or resolved component.'),
('expression_measurement_v71','expression_record_id','HPA RNA/protein measurement with mapping and missingness semantics.'),
('subcellular_localization_v71','localization_record_id','HPA localization observation.'),
('positive_interaction_evidence_v71','evidence_id','Public positive core.'),
('contextual_negative_evidence_v71','negative_evidence_id','Negative evidence relevant to a positive pair/selectivity/conflict.'),
('archived_isolated_negative_evidence_v71','negative_evidence_id','Isolated inactive evidence excluded from default publication.'),
('evidence_lineage_v71','evidence_id','Database, experiment, structure and modality provenance.'),
('protein_cross_classification_v71','classification_assertion_id','Five-axis classification with explicit provenance.'),
('disease_source_entity_v71','source_disease_entity_id','Original disease identity.'),
('disease_entity_master_v71','canonical_disease_id','Canonical MONDO disease identity.'),
('disease_hierarchy_v71','disease_hierarchy_edge_id','Disease parent-child edge.'),
('disease_anatomy_v71','disease_anatomy_assertion_id','Ontology-derived multi-label disease anatomy.'),
('protein_disease_relation_v71','disease_relation_id_v2','Protein to canonical disease relation retaining source identities.'),
]

ENUMS={
'site_residue_status':['RESIDUE_COMPLETE','RESIDUE_PARTIAL','RESIDUE_SOURCE_ONLY','NO_RESIDUE_REPORTED','PREDICTED_POCKET','NO_SITE_INFORMATION'],
'chain_assignment_status':['UNIQUE_CHAIN','MULTIPLE_EQUIVALENT_CHAINS','MULTIPLE_NON_EQUIVALENT_CHAINS','CHAIN_NOT_REPORTED','CHAIN_UNRESOLVED'],
'structure_ligand_relationship':['EXACT_COCRYSTAL_LIGAND','SAME_PARENT_DIFFERENT_FORM','ANALOG_LIGAND','OTHER_BIOLOGICAL_LIGAND','STRUCTURE_ONLY_NO_QUERY_LIGAND','ADDITIVE_OR_ARTIFACT','LIGAND_IDENTITY_UNRESOLVED','NOT_APPLICABLE'],
'target_entity_type':['gene_product_unspecified','canonical_protein','protein_isoform','processed_protein_form','protein_complex','unresolved_target'],
'membrane_class':['A','B','C','NON_MEMBRANE','UNRESOLVED'],
'measurement_status':['MEASURED','NOT_MEASURED','MISSING'],
'mapping_status':['MAPPED','UNMAPPED','AMBIGUOUS','NOT_APPLICABLE'],
'detection_status':['DETECTED','NOT_DETECTED','BELOW_THRESHOLD','NOT_APPLICABLE'],
'negative_release_layer':['CONTEXTUAL_PUBLIC','ISOLATED_ARCHIVE'],
}

def main():
    DOC.mkdir(parents=True,exist_ok=True); SCHEMA.mkdir(parents=True,exist_ok=True)
    (DOC/'APPROVED_RELEASE_RULES_V71.md').write_text(RULES,encoding='utf-8')
    with (SCHEMA/'TABLE_REGISTRY_V71.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['table_name','primary_key','description']); w.writerows(TABLES)
    with (SCHEMA/'CONTROLLED_VOCABULARY_V71.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['vocabulary','term']);
        for k,vals in ENUMS.items():
            for v in vals:w.writerow([k,v])
    old=Path(r'C:\Users\Administrator\v71_data_dictionary.tsv')
    if old.exists(): shutil.copy2(old,SCHEMA/'LEGACY_PUBLICATION_FIELD_DICTIONARY.tsv')
    (DOC/'V71_SCOPE_AND_LIMITATIONS.md').write_text('''# Scope and limitations\n\nV7.1 restructures existing frozen evidence and approved ontology snapshots. It does not invent residues, chains, isoforms, subunits, disease equivalences, or experimental independence when the source lacks them. Manual sampling was waived. Prediction-layer topology and pockets remain separate from confirmed evidence.\n''',encoding='utf-8')
    print(json.dumps({'tables':len(TABLES),'vocabularies':len(ENUMS)},ensure_ascii=False))
if __name__=='__main__':main()
