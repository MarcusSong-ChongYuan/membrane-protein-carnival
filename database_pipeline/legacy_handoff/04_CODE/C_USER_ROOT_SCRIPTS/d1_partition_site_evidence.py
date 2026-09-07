#!/usr/bin/env python3
import csv,json
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
SRC=ROOT/'docking_manifest_D1_sifts_qc.tsv'
with SRC.open(encoding='utf-8',newline='') as f: rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames;rows=list(rd)
for r in rows:
    has=(r.get('has_binding_residues')=='True' and int(r.get('requested_residue_count') or 0)>0)
    if has and r.get('pdb_source')=='direct': tier='G1_EXPERIMENTAL_RESIDUES'
    elif not has and r.get('pdb_source')=='direct': tier='G2_DIRECT_PDB_LIGAND_RECOVERY'
    elif r.get('pdb_source')=='fallback': tier='G3_FALLBACK_STRUCTURE_POCKET_TRANSFER_OR_PREDICTION'
    else: tier='G4_NO_STRUCTURE_POCKET_PREDICTION'
    r['site_evidence_route_v31']=tier
for tier in ['G1_EXPERIMENTAL_RESIDUES','G2_DIRECT_PDB_LIGAND_RECOVERY','G3_FALLBACK_STRUCTURE_POCKET_TRANSFER_OR_PREDICTION','G4_NO_STRUCTURE_POCKET_PREDICTION']:
    path=ROOT/f'docking_manifest_D1_{tier}.tsv'
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields+['site_evidence_route_v31'],delimiter='\t');w.writeheader();w.writerows(r for r in rows if r['site_evidence_route_v31']==tier)
summary={'tasks':len(rows),'routes':dict(Counter(r['site_evidence_route_v31'] for r in rows))}
(ROOT/'metadata'/'D1_SITE_EVIDENCE_PARTITION.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
