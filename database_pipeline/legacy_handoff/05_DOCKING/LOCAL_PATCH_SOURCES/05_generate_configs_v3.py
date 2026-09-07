#!/usr/bin/env python3
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'vina_configs'; META=ROOT/'metadata'
OUT.mkdir(exist_ok=True); META.mkdir(exist_ok=True)
for p in OUT.glob('*.conf'): p.unlink()
with (META/'task_prep_qc.tsv').open(encoding='utf-8') as f:
    prep={r['task_id']:r for r in csv.DictReader(f,delimiter='\t')}
ok=[]; blocked=[]
with (ROOT/'docking_manifest_with_grid.tsv').open(encoding='utf-8') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        tid=r['task_id']; status=prep.get(tid,{}).get('prep_qc_status','MISSING_PREP_QC')
        grid=(r.get('grid_confidence','').upper() in {'HIGH','MEDIUM'} and all((r.get(k) or '').strip() for k in ['grid_center_x','grid_center_y','grid_center_z','grid_size_x','grid_size_y','grid_size_z']))
        pdb=r['pdb_id'].lower(); cmpd=r['compound_internal_id']
        receptor=ROOT/'receptors_pdbqt_clean'/f'{pdb}.pdbqt'; ligand=ROOT/'ligands_pdbqt'/f'{cmpd}.pdbqt'
        ligand_atom_count=0
        if ligand.exists():
            with ligand.open(encoding='utf-8',errors='ignore') as lf:
                ligand_atom_count=sum(1 for line in lf if line.startswith(('ATOM  ','HETATM')))
        if status!='PASS' or not grid:
            blocked.append((tid,status,r.get('grid_qc_status',''))); continue
        if not receptor.exists() or receptor.stat().st_size == 0:
            blocked.append((tid,'MISSING_CLEAN_RECEPTOR',r.get('grid_qc_status',''))); continue
        if not ligand.exists() or ligand.stat().st_size == 0:
            blocked.append((tid,'MISSING_LIGAND_PDBQT',r.get('grid_qc_status',''))); continue
        if ligand_atom_count < 5:
            blocked.append((tid,'LIGAND_TOO_SMALL_FOR_POSE_DOCKING',r.get('grid_qc_status',''))); continue
        seed=20260000 + int(tid.rsplit('-',1)[-1])
        text=(f'receptor = receptors_pdbqt/{pdb}.pdbqt\nligand = ligands_pdbqt/{cmpd}.pdbqt\n'
              f'center_x = {r["grid_center_x"]}\ncenter_y = {r["grid_center_y"]}\ncenter_z = {r["grid_center_z"]}\n'
              f'size_x = {r["grid_size_x"]}\nsize_y = {r["grid_size_y"]}\nsize_z = {r["grid_size_z"]}\n'
              f'seed = {seed}\nthread = 8000\nnum_modes = 9\nenergy_range = 3\n'
              f'out = results/{tid}_out.pdbqt\n')
        (OUT/f'{tid}.conf').write_text(text,encoding='ascii'); ok.append(tid)
with (META/'runnable_after_production_qc.tsv').open('w',encoding='utf-8') as f: f.write('task_id\n'+'\n'.join(ok)+'\n')
with (META/'blocked_after_production_qc.tsv').open('w',encoding='utf-8') as f:
    f.write('task_id\tprep_qc_status\tgrid_qc_status\n'); f.writelines('\t'.join(x)+'\n' for x in blocked)
print(f'configs={len(ok)} blocked={len(blocked)} thread=8000')
