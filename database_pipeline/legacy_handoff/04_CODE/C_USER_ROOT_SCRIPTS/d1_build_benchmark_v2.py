#!/usr/bin/env python3
import csv,json,math,re
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent
MAN=ROOT/'docking_manifest_D1_G1_chain_receptor_qc_v2.tsv'
OUT=ROOT/'benchmark_v2_adjusted';CONF=OUT/'configs';OUT.mkdir(exist_ok=True);CONF.mkdir(exist_ok=True)
def ligand_qc(path):
 atoms={};types=[];charges=[];branches=[];issues=[]
 for line in path.read_text(encoding='utf-8',errors='ignore').splitlines():
  if line.startswith(('ATOM  ','HETATM')):
   try:serial=int(line[6:11]);xyz=tuple(float(line[a:b]) for a,b in ((30,38),(38,46),(46,54)))
   except:issues.append('UNPARSEABLE_ATOM');continue
   atoms[serial]=xyz;typ=line[77:79].strip() if len(line)>=79 else '';types.append(typ)
   try:charges.append(float(line[69:76]))
   except:charges.append(0.0)
  elif line.startswith('BRANCH'):
   x=line.split()
   if len(x)>=3:
    try:branches.append((int(x[1]),int(x[2])))
    except:issues.append('UNPARSEABLE_BRANCH')
 if not atoms:issues.append('NO_ATOMS')
 if any(not x for x in types):issues.append('BLANK_AUTODOCK_TYPE')
 coords=list(atoms.values());dup=len(coords)-len({tuple(round(v,4) for v in x) for x in coords})
 if dup:issues.append('DUPLICATE_ATOM_COORDINATES')
 deg=0;missing=0
 for a,b in branches:
  if a not in atoms or b not in atoms:missing+=1;continue
  d=math.dist(atoms[a],atoms[b])
  if d<0.1:deg+=1
 if missing:issues.append('BRANCH_REFERENCES_MISSING_ATOM')
 if deg:issues.append('DEGENERATE_TORSION_AXIS')
 spans=[max(x[i] for x in coords)-min(x[i] for x in coords) for i in range(3)] if coords else [0,0,0]
 return {'atoms':len(atoms),'rotors':len(branches),'duplicate_coords':dup,'degenerate_axes':deg,'net_charge':round(sum(charges),3),'span_x':spans[0],'span_y':spans[1],'span_z':spans[2],'max_span':max(spans),'issues':issues}

with MAN.open(encoding='utf-8',newline='') as f:rd=csv.DictReader(f,delimiter='\t');base=rd.fieldnames;rows=[r for r in rd if r['config_status']=='READY_FOR_PARAMETER_BENCHMARK']
extra=['ligand_atom_count_v2','ligand_rotor_count_v2','ligand_net_charge_v2','ligand_max_span_A_v2','ligand_qc_flags_v2','method_class_v2','adaptive_box_x_v2','adaptive_box_y_v2','adaptive_box_z_v2','thread_v2','benchmark_v2_disposition','benchmark_v2_reason']
out=[]
for r in rows:
 lig=ROOT/'ligands_pdbqt_chain_v2'/f"{r['compound_internal_id']}.pdbqt";q=ligand_qc(lig);reason=[]
 if q['issues']:reason.extend(q['issues'])
 if q['atoms']>80 or q['rotors']>20 or abs(q['net_charge'])>6 or q['max_span']>26:method='SPECIAL_OUT_OF_DOMAIN';reason.append('ROUTINE_VINA_DOMAIN_EXCEEDED')
 elif q['rotors']>12 or q['atoms']>55 or abs(q['net_charge'])>3:method='DIFFICULT'
 elif q['rotors']>6 or q['atoms']>35:method='MEDIUM'
 else:method='STANDARD'
 sizes=[]
 for axis in 'xyz':sizes.append(max(22.0,float(r['grid_size_'+axis]),q['span_'+axis]+10.0))
 if max(sizes)>32:reason.append('ADAPTIVE_BOX_AXIS_GT_32A')
 volume=sizes[0]*sizes[1]*sizes[2]
 thread=1000
 if method=='MEDIUM' or volume>20000:thread=3000
 if method=='DIFFICULT' or volume>26000:thread=5000
 if method=='SPECIAL_OUT_OF_DOMAIN':thread=0
 blocking={'NO_ATOMS','BLANK_AUTODOCK_TYPE','DUPLICATE_ATOM_COORDINATES','BRANCH_REFERENCES_MISSING_ATOM','DEGENERATE_TORSION_AXIS','UNPARSEABLE_ATOM','UNPARSEABLE_BRANCH'}
 if blocking & set(reason):disp='LIGAND_PREP_FAILED'
 elif method=='SPECIAL_OUT_OF_DOMAIN' or 'ADAPTIVE_BOX_AXIS_GT_32A' in reason:disp='METHOD_OUT_OF_DOMAIN'
 else:disp='READY_FOR_ADJUSTED_BENCHMARK'
 x=dict(r);x.update({'ligand_atom_count_v2':q['atoms'],'ligand_rotor_count_v2':q['rotors'],'ligand_net_charge_v2':q['net_charge'],'ligand_max_span_A_v2':round(q['max_span'],3),'ligand_qc_flags_v2':';'.join(q['issues']),'method_class_v2':method,'adaptive_box_x_v2':round(sizes[0],2),'adaptive_box_y_v2':round(sizes[1],2),'adaptive_box_z_v2':round(sizes[2],2),'thread_v2':thread,'benchmark_v2_disposition':disp,'benchmark_v2_reason':';'.join(dict.fromkeys(reason))});out.append(x)
 if disp=='READY_FOR_ADJUSTED_BENCHMARK':
  src=ROOT/'vina_configs_chain_v2'/f"{r['task_id']}.conf";lines=[]
  for line in src.read_text().splitlines():
   key=line.split('=',1)[0].strip() if '=' in line else ''
   if key in {'size_x','size_y','size_z','thread','num_modes','energy_range','out'}:continue
   lines.append(line)
  lines += [f'size_x = {sizes[0]:.2f}',f'size_y = {sizes[1]:.2f}',f'size_z = {sizes[2]:.2f}',f'thread = {thread}','num_modes = 20','energy_range = 4',f'out = {OUT}/results/{r["task_id"]}_out.pdbqt']
  (CONF/f"{r['task_id']}.conf").write_text('\n'.join(lines)+'\n',encoding='ascii')
with (OUT/'D1_ADJUSTED_BENCHMARK_TASKS.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=base+extra,delimiter='\t');w.writeheader();w.writerows(out)
summary={'input':len(out),'disposition':dict(Counter(x['benchmark_v2_disposition'] for x in out)),'method_class':dict(Counter(x['method_class_v2'] for x in out)),'reasons':dict(Counter(y for x in out for y in x['benchmark_v2_reason'].split(';') if y)),'policy':'thread by ligand/box complexity; adaptive box; malformed torsion and routine-Vina out-of-domain tasks blocked'}
(OUT/'D1_ADJUSTED_BENCHMARK_BUILD_REPORT.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
