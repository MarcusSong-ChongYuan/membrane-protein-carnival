#!/usr/bin/env python3
import csv,json
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(r"C:\Users\Administrator\MemPro_Docking_D1_AE1_20260813")
SRC=ROOT/'hetatm_rcsb_enrichment'/'hetatm_entities_rcsb_enriched.tsv'
OUT=ROOT/'hetatm_rcsb_enrichment'
with SRC.open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
by=defaultdict(list)
for r in rows:by[r['task_id']].append(r)
fields=['task_id','pdb_id','author_chain','target_uniprot','compound_internal_id','compound_name','hetero_entity_count','rcsb_codes','rcsb_names','has_metal_or_ion','has_cofactor','has_glycan','has_additive_candidate','has_ligand_or_unknown','task_hetatm_decision','decision_reason','automatic_deletion_allowed']
out=[]
for tid,xs in sorted(by.items()):
 cats={x['category_candidate'] for x in xs}; codes=sorted({x['het_resname'] for x in xs});names=sorted({x['rcsb_name'] for x in xs if x['rcsb_name']})
 metal='common_ion_or_metal' in cats;cof='cofactor_candidate' in cats;gly='glycan_or_sugar' in cats;add='buffer_or_crystallization_additive_candidate' in cats;lig='ligand_or_unknown_heteroentity' in cats
 if metal or cof:dec='HOLD_SPECIAL_CHEMISTRY';reason='Pocket contains metal/ion or cofactor; requires retain/remove decision and compatible receptor parameterization.'
 elif gly:dec='HOLD_GLYCAN_CONTEXT';reason='Pocket contains glycan/sugar; biological or crystallographic role must be reviewed.'
 elif lig:dec='HOLD_COGNATE_LIGAND_IDENTIFICATION';reason='Pocket contains ligand/unknown heteroentity; must identify cognate crystallographic ligand versus essential component.'
 elif add:dec='REMOVAL_CANDIDATE_AFTER_STRUCTURE_CONFIRMATION';reason='Only known additive candidates detected; verify structure context before receptor regeneration.'
 else:dec='HOLD_UNCLASSIFIED';reason='No safe automatic rule.'
 z=xs[0];out.append({'task_id':tid,'pdb_id':z['pdb_id'],'author_chain':z['author_chain'],'target_uniprot':z['target_uniprot'],'compound_internal_id':z['compound_internal_id'],'compound_name':z['compound_name'],'hetero_entity_count':len(xs),'rcsb_codes':';'.join(codes),'rcsb_names':';'.join(names),'has_metal_or_ion':int(metal),'has_cofactor':int(cof),'has_glycan':int(gly),'has_additive_candidate':int(add),'has_ligand_or_unknown':int(lig),'task_hetatm_decision':dec,'decision_reason':reason,'automatic_deletion_allowed':0})
with (OUT/'HETATM_TASK_DECISIONS_V4.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
summary={'tasks':len(out),'decisions':dict(Counter(x['task_hetatm_decision'] for x in out)),'rule':'No task permits automatic HETATM deletion. Removal candidates require structure-context confirmation.'}
(OUT/'HETATM_TASK_DECISIONS_V4_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
