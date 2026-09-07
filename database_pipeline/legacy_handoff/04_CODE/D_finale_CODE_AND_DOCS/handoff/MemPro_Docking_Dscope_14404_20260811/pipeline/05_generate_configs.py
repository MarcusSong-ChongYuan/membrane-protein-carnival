#!/usr/bin/env python3
from pathlib import Path
import csv
r=Path(__file__).resolve().parents[1];o=r/'vina_configs';o.mkdir(exist_ok=True);ok=[];bad=[]
with open(r/'docking_manifest_with_grid.tsv',encoding='utf-8') as f:
 for a in csv.DictReader(f,delimiter='\t'):
  s=(a.get('grid_qc_status') or '').strip();fs=['grid_center_x','grid_center_y','grid_center_z','grid_size_x','grid_size_y','grid_size_z'];ready=s.upper() in {'HIGH','MEDIUM','PASS'} and all((a.get(x) or '').strip() for x in fs)
  if not ready:bad.append((a.get('task_id',''),s,a.get('grid_method','')));continue
  t=a['task_id'];p=a['pdb_id'].lower();c=a['compound_internal_id'];txt=f"receptor = receptors_pdbqt/{p}.pdbqt\nligand = ligands_pdbqt/{c}.pdbqt\ncenter_x = {a['grid_center_x']}\ncenter_y = {a['grid_center_y']}\ncenter_z = {a['grid_center_z']}\nsize_x = {a['grid_size_x']}\nsize_y = {a['grid_size_y']}\nsize_z = {a['grid_size_z']}\nnum_modes = 9\nenergy_range = 3\nout = results/{t}_out.pdbqt\n";(o/f'{t}.conf').write_text(txt,encoding='ascii');ok.append(t)
with open(r/'metadata'/'runnable_after_grid.tsv','w',encoding='utf-8') as f:f.write('task_id\n'+'\n'.join(ok)+'\n')
with open(r/'metadata'/'blocked_after_grid.tsv','w',encoding='utf-8') as f:f.write('task_id\tgrid_qc_status\tgrid_method\n');f.writelines(f'{a}\t{b}\t{c}\n' for a,b,c in bad)
print(f'configs={len(ok)} blocked={len(bad)}');raise SystemExit(0 if ok else 2)
