#!/usr/bin/env python3
"""Create paired thread=5000/8000 configs for representative runnable tasks."""
import csv
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]; SRC=ROOT/'vina_configs'; OUT=ROOT/'benchmark_configs'; META=ROOT/'metadata'
OUT.mkdir(exist_ok=True); META.mkdir(exist_ok=True)
for p in OUT.glob('*.conf'): p.unlink()
rows=[]
with (ROOT/'docking_manifest_with_grid.tsv').open(encoding='utf-8') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
by_id={r['task_id']:r for r in rows}
runnable=[p.stem for p in sorted(SRC.glob('*.conf'))]
def ligand_complexity(r):
    smi=r.get('smiles',''); return (sum(smi.count(x) for x in ['=','#','(',')','@']),len(smi))
ordered=sorted(runnable,key=lambda x:ligand_complexity(by_id[x]))
if not ordered: raise SystemExit('No runnable configs')
n=min(12,len(ordered)); idx=sorted(set(round(i*(len(ordered)-1)/(n-1)) if n>1 else 0 for i in range(n)))
selected=[ordered[i] for i in idx]
with (META/'thread_benchmark_tasks.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['task_id','smiles_length','box_x','box_y','box_z'])
    for tid in selected:
        r=by_id[tid]; w.writerow([tid,len(r.get('smiles','')),r['grid_size_x'],r['grid_size_y'],r['grid_size_z']])
        base=(SRC/f'{tid}.conf').read_text(encoding='ascii')
        for thread in (5000,8000):
            text='\n'.join((f'thread = {thread}' if line.startswith('thread = ') else f'out = results/{tid}_t{thread}_out.pdbqt' if line.startswith('out = ') else line) for line in base.splitlines())+'\n'
            (OUT/f'{tid}_t{thread}.conf').write_text(text,encoding='ascii')
print(f'benchmark_tasks={len(selected)} configs={len(selected)*2}')
