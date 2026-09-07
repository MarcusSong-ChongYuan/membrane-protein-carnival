import csv,gzip,json
from pathlib import Path
root=Path(r"D:\finale\12_V7_final_completion_20260813");out=root/'00_baseline';out.mkdir(exist_ok=True)
rules=[
('membrane_class','A','transmembrane/intramembrane feature OR validated integral structure','A overrides B/C; secondary mechanisms retained','PUBLIC if target identity valid'),
('membrane_class','B','covalent GPI/lipid anchor with no A mechanism','not inferred from localization alone','PUBLIC if target identity valid'),
('membrane_class','C','direct/state-dependent peripheral binding OR validated complex-mediated attachment','GO/HPA location alone prohibited','PUBLIC if mechanism explicit'),
('membrane_evidence','E1','direct experimental ECO/PubMed or structure-supported membrane mechanism','strongest satisfied rule wins','PUBLIC'),
('membrane_evidence','E2','curated topology, anchor, or peripheral annotation without direct experiment','manual annotation accepted','PUBLIC'),
('membrane_evidence','E3','inferred/external nonexperimental mechanism','limitations retained','PUBLIC_E3'),
('membrane_evidence','E0','not in final A/B/C layer','not a membrane confidence grade','EXCLUDED_AUDIT'),
('binding_evidence','BE1','direct experimental PDB structural context','PDB required','PUBLIC'),
('binding_evidence','BE2','unique human target + direct Kd/Ki/affinity + positive standardized nM + source anchor','IC50/EC50 prohibited from BE2','PUBLIC'),
('binding_evidence','BE3','target-specific functional/pharmacological/curated relationship','does not claim direct equilibrium affinity','PUBLIC'),
('target_entity','protein_isoform','explicit recoverable isoform foreign key','gene-name inference prohibited','PUBLIC_EXACT_ISOFORM'),
('target_entity','canonical_protein','explicit/current single human UniProt with no safe isoform specificity','historical isoform loss stated','PUBLIC_CANONICAL'),
('target_entity','complex','source-asserted complex + unique compound + direct relation','otherwise frozen','PUBLIC or FROZEN_REVIEW'),
('negative_evidence','PUBLIC_NEGATIVE','V7 protein + canonical compound + optional form FK + mapped negative release + safe structure identity','negative compound need not occur in positive table','PUBLIC'),
('binding_site','HIGH','current RCSB PDB + explicit chain + unique SIFTS mapping + all residues mapped/observed','no chain guessing','PUBLIC_COORDINATE_VERIFIED'),
('binding_site','MEDIUM','unique best SIFTS PDB-chain with partial/unobserved mapping','limitations explicit','PUBLIC_COORDINATE_PARTIAL'),
('binding_site','FROZEN','missing chain, zero match, tied best mapping, invalid PDB','cannot generate docking box','FROZEN_REVIEW'),
('general','missing_value','empty means missing/unavailable','never convert to numeric zero','preserve'),
('general','manual_review','WAIVED_BY_USER','unresolved states become frozen final data states','not claimed completed'),
]
with (out/'V7_CLASSIFICATION_AND_RELEASE_RULES.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.writer(f,delimiter='\t');w.writerow(['rule_domain','label','trigger','prohibition_or_limitation','disposition']);w.writerows(rules)

# Generate a field-level dictionary from every current module output available.
folders=[root/'01_entity_assignment',root/'02_be_recompute',root/'03_binding_sites',root/'04_negative_conflicts',root/'05_hpa_expression',root/'06_remaining_modules'];rows=[]
for folder in folders:
 for p in sorted(folder.glob('*.tsv*')):
  try:
   op=gzip.open if p.suffix=='.gz' else open
   kw={'mode':'rt','encoding':'utf-8','newline':''} if p.suffix=='.gz' else {'mode':'r','encoding':'utf-8','newline':''}
   with op(p,**kw) as f:header=next(csv.reader(f,delimiter='\t'))
  except Exception:continue
  for i,c in enumerate(header,1):rows.append([folder.name,p.name,i,c,'empty means missing/not available; see source module and rules','MemPro V7 final completion'])
with (out/'V7_DATA_DICTIONARY_DRAFT.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.writer(f,delimiter='\t');w.writerow(['module','table','column_order','column_name','null_semantics','release']);w.writerows(rows)
report={'rule_rows':len(rules),'dictionary_rows':len(rows),'manual_review':'WAIVED_BY_USER'};(out/'T02_RULES_DICTIONARY_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
