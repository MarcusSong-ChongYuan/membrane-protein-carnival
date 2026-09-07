#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813
VDIR=/home/csong/Vina-GPU
BIN=$VDIR/Vina-GPU
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}
OUT=$ROOT/multiseed_benchmark_v5
mkdir -p "$OUT"/{configs,results,logs}
cd "$ROOT"
cp thread_benchmark_v4/SELECTED_TASKS.tsv "$OUT/SELECTED_TASKS.tsv"
printf 'task_id\treceptor_atoms\tthread\tseed\texit_code\telapsed_seconds\tbest_affinity_kcal_mol\n' > "$OUT/MULTISEED_RUNS.tsv"
while IFS=$'\t' read -r task pdb chain uniprot cmpd name cfgex recex atoms rb ligex la lb bx by bz bmax frac white status reasons; do
  [[ "$task" == task_id ]] && continue
  src="$ROOT/vina_configs_chain_v2/${task}.conf"
  for thread in 1000 5000 8000; do
    for seed in 20260813 20260814 20260815; do
      cfg="$OUT/configs/${task}_t${thread}_s${seed}.conf"
      awk -v root="$ROOT" -v out="$OUT/results/${task}_t${thread}_s${seed}_out.pdbqt" -v th="$thread" '
        /^receptor[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
        /^ligand[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
        /^thread[[:space:]]*=/ {$0="thread = " th}
        /^search_depth[[:space:]]*=/ {next}
        /^out[[:space:]]*=/ {$0="out = " out}
        {print}
      ' "$src" > "$cfg"
      start=$(date +%s); code=0
      cd "$VDIR"
      timeout 1800 "$BIN" --config "$cfg" --seed "$seed" --log "$OUT/logs/${task}_t${thread}_s${seed}.vina.log" > "$OUT/logs/${task}_t${thread}_s${seed}.stdout.log" 2>&1 || code=$?
      cd "$ROOT"; end=$(date +%s)
      affinity=$(awk '/^[[:space:]]*1[[:space:]]+[-0-9.]+/{print $2;exit}' "$OUT/logs/${task}_t${thread}_s${seed}.stdout.log")
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$atoms" "$thread" "$seed" "$code" "$((end-start))" "${affinity:-}" >> "$OUT/MULTISEED_RUNS.tsv"
    done
  done
done < "$OUT/SELECTED_TASKS.tsv"
python - <<'PY'
import csv, json, statistics
from collections import defaultdict
p='multiseed_benchmark_v5/MULTISEED_RUNS.tsv'
with open(p,encoding='utf-8') as f: rows=list(csv.DictReader(f,delimiter='\t'))
g=defaultdict(list)
for r in rows:
    if r['exit_code']=='0' and r['best_affinity_kcal_mol']:
        g[(r['task_id'],r['thread'])].append(float(r['best_affinity_kcal_mol']))
summary=[]
for (task,thread),vals in sorted(g.items()):
    summary.append({'task_id':task,'thread':int(thread),'n':len(vals),'mean_affinity':round(statistics.mean(vals),3),'sd_affinity':round(statistics.stdev(vals),3) if len(vals)>1 else 0.0,'range_affinity':round(max(vals)-min(vals),3)})
with open('multiseed_benchmark_v5/MULTISEED_SUMMARY.json','w',encoding='utf-8') as f:json.dump(summary,f,indent=2)
print(json.dumps(summary,indent=2))
PY
