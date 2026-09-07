#!/usr/bin/env python3
import csv,glob,json,statistics
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'benchmark_v2_adjusted';rows=[]
for p in glob.glob(str(OUT/'status'/'*.tsv')):
 with open(p,encoding='utf-8') as f:rows.extend(csv.DictReader(f,delimiter='\t'))
by=defaultdict(list)
for r in rows:by[r['task_id']].append(r)
task=[]
for t,rs in by.items():
 vals=[float(r['best_affinity_kcal_mol']) for r in rs if r['exit_code']=='0' and r['best_affinity_kcal_mol']]
 task.append({'task_id':t,'successful_seeds':len(vals),'affinity_range':max(vals)-min(vals) if len(vals)>1 else None,'positive_affinity':any(v>0 for v in vals)})
success=sum(r['exit_code']=='0' and bool(r['best_affinity_kcal_mol']) for r in rows);complete3=sum(x['successful_seeds']==3 for x in task);stable=sum(x['affinity_range'] is not None and x['affinity_range']<=1 for x in task);positive=sum(x['positive_affinity'] for x in task)
technical=(success/len(rows)>=.99 if rows else False) and complete3==len(task) and stable/max(1,len(task))>=.9 and positive==0
report={'runs':len(rows),'successful_runs':success,'tasks':len(task),'tasks_all_three_seeds':complete3,'tasks_affinity_range_le_1':stable,'tasks_positive_affinity':positive,'technical_gate':'PASS' if technical else 'FAIL','scientific_gate':'PENDING_REDOCKING_RMSD_AND_CONTACT_RECALL','formal_submission_ready':False,'task_results':sorted(task,key=lambda x:x['task_id'])}
(OUT/'D1_ADJUSTED_BENCHMARK_EVALUATION.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='task_results'},indent=2))
