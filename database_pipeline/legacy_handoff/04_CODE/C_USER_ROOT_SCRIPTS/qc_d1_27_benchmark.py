#!/usr/bin/env python3
import csv,json,re
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parent
MAN=ROOT/'docking_manifest_D1_G1_chain_receptor_qc_v2.tsv'
OUT=ROOT/'benchmark_qc_chain_v2';OUT.mkdir(exist_ok=True)
with MAN.open(encoding='utf-8',newline='') as f:rows=[r for r in csv.DictReader(f,delimiter='\t') if r['config_status']=='READY_FOR_PARAMETER_BENCHMARK']
def atom_stats(p):
 total=blank=0
 if not p.exists() or not p.stat().st_size:return total,blank
 for line in p.read_text(encoding='utf-8',errors='ignore').splitlines():
  if line.startswith(('ATOM  ','HETATM')):
   total+=1
   if len(line)<79 or not line[77:79].strip():blank+=1
 return total,blank
fields=['task_id','pdb_id','actual_chain_used','target_uniprot','compound_internal_id','preferred_name','config_exists','receptor_exists','receptor_atom_count','receptor_blank_type_count','ligand_exists','ligand_atom_count','ligand_blank_type_count','box_x','box_y','box_z','box_max_axis_A','matched_residue_fraction','identity_whitelist_pass','qc_status','qc_reasons']
out=[]
for r in rows:
 conf=ROOT/'vina_configs_chain_v2'/(r['task_id']+'.conf');rec=ROOT/'receptors_pdbqt_clean_chain_v2'/(r['chain_receptor_key']+'.pdbqt');lig=ROOT/'ligands_pdbqt_chain_v2'/(r['compound_internal_id']+'.pdbqt')
 ra,rb=atom_stats(rec);la,lb=atom_stats(lig);sizes=[float(r['grid_size_'+a]) for a in 'xyz']; reasons=[]
 if not conf.exists():reasons.append('CONFIG_MISSING')
 if not ra:reasons.append('RECEPTOR_MISSING_OR_EMPTY')
 if rb:reasons.append('RECEPTOR_BLANK_AUTODOCK_TYPE')
 if not la:reasons.append('LIGAND_MISSING_OR_EMPTY')
 if lb:reasons.append('LIGAND_BLANK_AUTODOCK_TYPE')
 if max(sizes)>30:reasons.append('BOX_AXIS_GT_30A')
 if float(r.get('matched_residue_fraction_coord') or 0)<0.8:reasons.append('RESIDUE_COVERAGE_LT_80PCT')
 if r.get('identity_whitelist_pass')!='1':reasons.append('IDENTITY_NOT_WHITELISTED')
 out.append({'task_id':r['task_id'],'pdb_id':r['pdb_id'],'actual_chain_used':r['actual_chain_used'],'target_uniprot':r['target_uniprot'],'compound_internal_id':r['compound_internal_id'],'preferred_name':r['preferred_name'],'config_exists':int(conf.exists()),'receptor_exists':int(rec.exists()),'receptor_atom_count':ra,'receptor_blank_type_count':rb,'ligand_exists':int(lig.exists()),'ligand_atom_count':la,'ligand_blank_type_count':lb,'box_x':sizes[0],'box_y':sizes[1],'box_z':sizes[2],'box_max_axis_A':max(sizes),'matched_residue_fraction':r.get('matched_residue_fraction_coord',''),'identity_whitelist_pass':r.get('identity_whitelist_pass',''),'qc_status':'PASS_STATIC_QC' if not reasons else 'HOLD','qc_reasons':';'.join(reasons)})
with (OUT/'D1_27_BENCHMARK_STATIC_QC.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
summary={'input':len(out),'status':dict(Counter(x['qc_status'] for x in out)),'reasons':dict(Counter(y for x in out for y in x['qc_reasons'].split(';') if y)),'note':'Static QC only. PASS does not launch docking and does not validate protonation, stereochemistry, missing pocket residues, or redocking RMSD.'}
(OUT/'D1_27_BENCHMARK_STATIC_QC_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
