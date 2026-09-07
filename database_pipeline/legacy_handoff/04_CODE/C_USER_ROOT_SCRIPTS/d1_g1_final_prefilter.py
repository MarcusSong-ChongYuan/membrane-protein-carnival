#!/usr/bin/env python3
import csv,json,re
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
with (ROOT/'D1_identity_whitelist_v7_2_1.tsv').open(encoding='utf-8',newline='') as f: allowed={r['target_uniprot'] for r in csv.DictReader(f,delimiter='\t')}
src=ROOT/'docking_manifest_D1_G1_coordinate_pass.tsv'
with src.open(encoding='utf-8',newline='') as f: rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames;rows=list(rd)

def ligand_flags(r):
    s=(r.get('smiles') or '').strip(); name=(r.get('preferred_name') or '').upper(); flags=[]
    # Conservative standard-Vina suitability gates; no Lipinski filtering.
    if not s: flags.append('MISSING_SMILES')
    atoms=re.findall(r'Br|Cl|Si|Se|Na|Ca|Mg|Zn|Fe|Mn|Cu|Co|Ni|[A-Z][a-z]?',s)
    organic=sum(a in {'C','N','O','S','P','F','Cl','Br','I','B','Si','Se'} for a in atoms)
    carbons=sum(a=='C' for a in atoms)
    if atoms and carbons==0: flags.append('NO_CARBON_INORGANIC_OR_ION')
    if len(atoms)<=2: flags.append('VERY_SMALL_FRAGMENT_OR_ION')
    if len(atoms)>120: flags.append('TOO_MANY_ATOMS_FOR_STANDARD_VINA')
    if s.count('.')>0: flags.append('MULTICOMPONENT_FORM_REQUIRES_PARENT_SELECTION')
    if any(x in name for x in ['SULFATE','BROMIDE','CHLORIDE','NITRATE','PHOSPHATE']) and carbons==0: flags.append('LIKELY_BUFFER_OR_INORGANIC')
    return sorted(set(flags))

extra=['identity_whitelist_pass','ligand_standard_vina_status','ligand_standard_vina_flags','d1_prefilter_status','d1_prefilter_reason']
out=[]
for r in rows:
    flags=ligand_flags(r); reasons=[]
    ident=r['target_uniprot'] in allowed
    if not ident: reasons.append('MEMBRANE_IDENTITY_NOT_WHITELISTED')
    if flags: reasons.extend(flags)
    x=dict(r);x.update({'identity_whitelist_pass':'1' if ident else '0','ligand_standard_vina_status':'PASS' if not flags else 'HOLD','ligand_standard_vina_flags':';'.join(flags),'d1_prefilter_status':'PASS_TO_RECEPTOR_QC' if not reasons else 'HOLD','d1_prefilter_reason':';'.join(reasons)});out.append(x)
for status in ['PASS_TO_RECEPTOR_QC','HOLD']:
    fn='docking_manifest_D1_G1_prefilter_pass.tsv' if status.startswith('PASS') else 'docking_manifest_D1_G1_prefilter_hold.tsv'
    with (ROOT/fn).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields+extra,delimiter='\t');w.writeheader();w.writerows(r for r in out if r['d1_prefilter_status']==status)
summary={'coordinate_pass_input':len(out),'pass_to_receptor_qc':sum(r['d1_prefilter_status']=='PASS_TO_RECEPTOR_QC' for r in out),'hold':sum(r['d1_prefilter_status']=='HOLD' for r in out),'hold_reasons':dict(Counter(x for r in out for x in r['d1_prefilter_reason'].split(';') if x))}
(ROOT/'metadata'/'D1_G1_PREFILTER_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
