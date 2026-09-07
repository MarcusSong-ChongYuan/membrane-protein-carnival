import csv, glob, json, math, os, statistics
from collections import Counter, defaultdict

root = "/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813/production_batch_v1"
rows=[]
for path in glob.glob(root+"/status/*.status.tsv"):
    with open(path,encoding="utf-8") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            r["exit_code"]=int(r["exit_code"] or -999)
            r["thread"]=int(r["thread"] or 0)
            r["elapsed_seconds"]=float(r["elapsed_seconds"] or 0)
            r["affinity"]=float(r["best_affinity_kcal_mol"]) if r["best_affinity_kcal_mol"] else None
            r["status_file"]=os.path.basename(path); rows.append(r)
success=[r for r in rows if r["exit_code"]==0 and r["affinity"] is not None]
failed=[r for r in rows if r not in success]
def quantile(xs,q):
    xs=sorted(xs);p=(len(xs)-1)*q;lo=int(math.floor(p));hi=int(math.ceil(p))
    return xs[lo] if lo==hi else xs[lo]+(xs[hi]-xs[lo])*(p-lo)
aff=[r["affinity"] for r in success]
bytask=defaultdict(list)
for r in rows:bytask[r["task_id"]].append(r)
taskstats=[]
for task,rs in sorted(bytask.items()):
    vals=[r["affinity"] for r in rs if r["exit_code"]==0 and r["affinity"] is not None]
    taskstats.append({"task_id":task,"success_seeds":len(vals),"failed_seeds":len(rs)-len(vals),"min_affinity":min(vals) if vals else None,"max_affinity":max(vals) if vals else None,"range":max(vals)-min(vals) if len(vals)>1 else None,"values":vals,"threads":sorted(set(r['thread'] for r in rs))})
ranges=[x["range"] for x in taskstats if x["range"] is not None]
threadstats={}
for t in sorted(set(r['thread'] for r in rows)):
    rs=[r for r in rows if r['thread']==t];ok=[r for r in rs if r in success]
    threadstats[str(t)]={"runs":len(rs),"success":len(ok),"failed":len(rs)-len(ok),"median_seconds":statistics.median([r['elapsed_seconds'] for r in ok]) if ok else None,"median_affinity":statistics.median([r['affinity'] for r in ok]) if ok else None}
report={
 "runs":{"planned":len(rows),"success":len(success),"failed":len(failed),"success_rate":len(success)/len(rows) if rows else 0},
 "affinity":{"min":min(aff),"q25":quantile(aff,.25),"median":statistics.median(aff),"q75":quantile(aff,.75),"max":max(aff),"mean":statistics.mean(aff)},
 "seed_reproducibility":{"tasks":len(taskstats),"tasks_all_3_success":sum(x['success_seeds']==3 for x in taskstats),"tasks_with_failure":sum(x['failed_seeds']>0 for x in taskstats),"median_affinity_range":statistics.median(ranges),"range_le_0_5":sum(x<=.5 for x in ranges),"range_le_1_0":sum(x<=1 for x in ranges),"range_gt_1_0":sum(x>1 for x in ranges)},
 "thread_stats":threadstats,
 "failed_runs":[{k:r[k] for k in ('task_id','seed','gpu','thread','exit_code','elapsed_seconds','status_file')} for r in failed],
 "least_reproducible_tasks":sorted([x for x in taskstats if x['range'] is not None],key=lambda x:x['range'],reverse=True)[:10],
 "strongest_tasks":sorted([x for x in taskstats if x['min_affinity'] is not None],key=lambda x:x['min_affinity'])[:10],
 "weakest_tasks":sorted([x for x in taskstats if x['max_affinity'] is not None],key=lambda x:x['max_affinity'],reverse=True)[:10],
}
print(json.dumps(report,ensure_ascii=False,indent=2))
