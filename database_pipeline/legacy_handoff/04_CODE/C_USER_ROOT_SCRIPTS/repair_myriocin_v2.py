#!/usr/bin/env python3
import csv,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent;MAN=ROOT/'docking_manifest_D1_G1_chain_receptor_qc_v2.tsv';OUT=ROOT/'ligands_pdbqt_repaired_v2';OUT.mkdir(exist_ok=True)
with MAN.open(encoding='utf-8',newline='') as f:row=next(r for r in csv.DictReader(f,delimiter='\t') if r['task_id']=='DOCK-002555')
cid=row['compound_internal_id'];smi=row['smiles'];sdf=OUT/f'{cid}.sdf';pdbqt=OUT/f'{cid}.pdbqt';log=OUT/f'{cid}.prepare.log'
commands=[['/home/csong/miniconda3/envs/docking/bin/obabel',f'-:{smi}','--gen3d','-p','7.4','-O',str(sdf)],['/home/csong/miniconda3/envs/docking/bin/obabel',str(sdf),'-O',str(pdbqt)]];messages=[]
for cmd in commands:
 r=subprocess.run(cmd,text=True,capture_output=True,timeout=600);messages.append({'cmd':cmd,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
 if r.returncode:break
log.write_text(json.dumps(messages,indent=2),encoding='utf-8');print(json.dumps({'compound':cid,'pdbqt':str(pdbqt),'exists':pdbqt.exists(),'bytes':pdbqt.stat().st_size if pdbqt.exists() else 0,'steps':[x['returncode'] for x in messages]},indent=2))
