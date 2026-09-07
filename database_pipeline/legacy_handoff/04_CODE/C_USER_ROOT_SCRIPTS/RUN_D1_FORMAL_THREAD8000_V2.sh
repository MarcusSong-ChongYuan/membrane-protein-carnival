#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813
VDIR=/home/csong/Vina-GPU
BIN=$VDIR/Vina-GPU
QC=$ROOT/benchmark_qc_chain_v2/D1_27_BENCHMARK_STATIC_QC.tsv
OUT=$ROOT/formal_thread8000_v2
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}
mkdir -p "$OUT"/{configs,results,logs,status}
printf 'task_id\treceptor_atoms\tthread\tseed\tgpu\tconfig\n' > "$OUT/FORMAL_TASKS.tsv"

make_config() {
  local task=$1 atoms=$2 seed=$3 gpu=$4
  local src="$ROOT/vina_configs_chain_v2/${task}.conf"
  local cfg="$OUT/configs/${task}_s${seed}.conf"
  awk -v root="$ROOT" -v out="$OUT/results/${task}_s${seed}_out.pdbqt" '
    /^receptor[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
    /^ligand[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
    /^thread[[:space:]]*=/ {$0="thread = 8000"}
    /^search_depth[[:space:]]*=/ {next}
    /^num_modes[[:space:]]*=/ {$0="num_modes = 20"}
    /^energy_range[[:space:]]*=/ {$0="energy_range = 4"}
    /^out[[:space:]]*=/ {$0="out = " out}
    {print}
  ' "$src" > "$cfg"
  printf '%s\t%s\t8000\t%s\t%s\t%s\n' "$task" "$atoms" "$seed" "$gpu" "$cfg" >> "$OUT/FORMAL_TASKS.tsv"
}

i=0
while IFS=$'\t' read -r task pdb chain uniprot cmpd name cfgex recex atoms rb ligex la lb bx by bz bmax frac white status reasons; do
  [[ "$task" == task_id ]] && continue
  [[ "$status" == PASS_STATIC_QC ]] || continue
  for seed in 20260813 20260814 20260815; do
    gpu=$((i % 2))
    make_config "$task" "$atoms" "$seed" "$gpu"
    i=$((i+1))
  done
done < "$QC"

planned=$(( $(wc -l < "$OUT/FORMAL_TASKS.tsv") - 1 ))
[[ "$planned" -eq 81 ]] || { echo "Expected 81 runs, got $planned" >&2; exit 8; }

worker() {
  local gpu=$1
  awk -F '\t' -v g="$gpu" 'NR>1 && $5==g {print}' "$OUT/FORMAL_TASKS.tsv" |
  while IFS=$'\t' read -r task atoms thread seed gp cfg; do
    result="$OUT/results/${task}_s${seed}_out.pdbqt"
    statusfile="$OUT/status/${task}_s${seed}.status.tsv"
    [[ -s "$result" ]] && continue
    start=$(date +%s); code=0
    cd "$VDIR"
    CUDA_VISIBLE_DEVICES="$gpu" timeout 3600 "$BIN" --config "$cfg" --seed "$seed" \
      --log "$OUT/logs/${task}_s${seed}.vina.log" \
      > "$OUT/logs/${task}_s${seed}.stdout.log" 2>&1 || code=$?
    cd "$ROOT"; end=$(date +%s)
    affinity=$(awk '/^[[:space:]]*1[[:space:]]+[-0-9.]+/{print $2;exit}' "$OUT/logs/${task}_s${seed}.stdout.log")
    printf 'task_id\tseed\tgpu\tthread\texit_code\telapsed_seconds\tbest_affinity_kcal_mol\n%s\t%s\t%s\t8000\t%s\t%s\t%s\n' \
      "$task" "$seed" "$gpu" "$code" "$((end-start))" "${affinity:-}" > "$statusfile"
  done
}

printf 'RUNNING pid=%s planned=%s start=%s\n' "$$" "$planned" "$(date -Iseconds)" > "$OUT/STATUS"
worker 0 & p0=$!
worker 1 & p1=$!
wait "$p0" "$p1"

python3 - <<'PY'
import csv, glob, json
rows=[]
for path in glob.glob('formal_thread8000_v2/status/*.status.tsv'):
    with open(path,encoding='utf-8') as f: rows.extend(csv.DictReader(f,delimiter='\t'))
summary={
    'planned_runs':81,
    'completed_status_rows':len(rows),
    'successful_runs':sum(r['exit_code']=='0' for r in rows),
    'failed_runs':sum(r['exit_code']!='0' for r in rows),
    'result_files':len(glob.glob('formal_thread8000_v2/results/*_out.pdbqt')),
    'thread':8000,
    'seeds':[20260813,20260814,20260815],
}
with open('formal_thread8000_v2/FORMAL_SUMMARY.json','w') as f: json.dump(summary,f,indent=2)
print(json.dumps(summary,indent=2))
PY
printf 'COMPLETE time=%s\n' "$(date -Iseconds)" > "$OUT/STATUS"
