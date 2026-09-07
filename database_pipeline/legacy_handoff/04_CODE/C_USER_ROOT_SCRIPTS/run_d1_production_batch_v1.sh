#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813
VDIR=/home/csong/Vina-GPU
BIN=$VDIR/Vina-GPU
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}
OUT=$ROOT/production_batch_v1
mkdir -p "$OUT"/{configs,results,logs,status}
cd "$ROOT"
QC=benchmark_qc_chain_v2/D1_27_BENCHMARK_STATIC_QC.tsv
printf 'task_id\treceptor_atoms\tthread\tseed\tgpu\tconfig\n' > "$OUT/PRODUCTION_TASKS.tsv"

make_config() {
  local task=$1 atoms=$2 seed=$3 gpu=$4
  local thread=1000
  # Roughly >=300 aa receptors generally exceed ~2400 heavy atoms. These receive
  # at least 5000 threads as requested; smaller receptors retain the benchmarked
  # 1000-thread baseline. Three fixed seeds quantify stochastic variation.
  if (( atoms >= 2400 )); then thread=5000; fi
  local src="$ROOT/vina_configs_chain_v2/${task}.conf"
  local cfg="$OUT/configs/${task}_s${seed}.conf"
  awk -v root="$ROOT" -v out="$OUT/results/${task}_s${seed}_out.pdbqt" -v th="$thread" '
    /^receptor[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
    /^ligand[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
    /^thread[[:space:]]*=/ {$0="thread = " th}
    /^search_depth[[:space:]]*=/ {next}
    /^num_modes[[:space:]]*=/ {$0="num_modes = 20"}
    /^energy_range[[:space:]]*=/ {$0="energy_range = 4"}
    /^out[[:space:]]*=/ {$0="out = " out}
    {print}
  ' "$src" > "$cfg"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$atoms" "$thread" "$seed" "$gpu" "$cfg" >> "$OUT/PRODUCTION_TASKS.tsv"
}

i=0
while IFS=$'\t' read -r task pdb chain uniprot cmpd name cfgex recex atoms rb ligex la lb bx by bz bmax frac white status reasons; do
  [[ "$task" == task_id ]] && continue
  [[ "$status" == PASS_STATIC_QC ]] || continue
  for seed in 20260813 20260814 20260815; do
    gpu=$((i % 2)); make_config "$task" "$atoms" "$seed" "$gpu"; i=$((i+1))
  done
done < "$QC"

worker() {
  local gpu=$1
  awk -F '\t' -v g="$gpu" 'NR>1 && $5==g {print}' "$OUT/PRODUCTION_TASKS.tsv" | while IFS=$'\t' read -r task atoms thread seed gp cfg; do
    statusfile="$OUT/status/${task}_s${seed}.status.tsv"
    if [[ -s "$OUT/results/${task}_s${seed}_out.pdbqt" ]]; then continue; fi
    start=$(date +%s); code=0
    cd "$VDIR"
    CUDA_VISIBLE_DEVICES="$gpu" timeout 3600 "$BIN" --config "$cfg" --seed "$seed" --log "$OUT/logs/${task}_s${seed}.vina.log" > "$OUT/logs/${task}_s${seed}.stdout.log" 2>&1 || code=$?
    cd "$ROOT"; end=$(date +%s)
    affinity=$(awk '/^[[:space:]]*1[[:space:]]+[-0-9.]+/{print $2;exit}' "$OUT/logs/${task}_s${seed}.stdout.log")
    printf 'task_id\tseed\tgpu\tthread\texit_code\telapsed_seconds\tbest_affinity_kcal_mol\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$seed" "$gpu" "$thread" "$code" "$((end-start))" "${affinity:-}" > "$statusfile"
  done
}
worker 0 & p0=$!
worker 1 & p1=$!
wait "$p0" "$p1"
python - <<'PY'
import csv, glob, json, os
rows=[]
for p in glob.glob('production_batch_v1/status/*.status.tsv'):
    with open(p,encoding='utf-8') as f: rows.extend(csv.DictReader(f,delimiter='\t'))
summary={'planned_runs':81,'completed_status_rows':len(rows),'successful_runs':sum(r['exit_code']=='0' for r in rows),'failed_runs':sum(r['exit_code']!='0' for r in rows),'result_files':len(glob.glob('production_batch_v1/results/*_out.pdbqt'))}
with open('production_batch_v1/PRODUCTION_SUMMARY.json','w') as f:json.dump(summary,f,indent=2)
print(json.dumps(summary,indent=2))
PY
