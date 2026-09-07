#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813
BIN=/home/csong/Vina-GPU/Vina-GPU
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}
cd "$ROOT"
mkdir -p benchmark_smoke_v2/{results,logs,configs}
cfg=$(find vina_configs_chain_v2 -maxdepth 1 -name '*.conf' | sort | head -1)
id=$(basename "$cfg" .conf)
awk -v out="benchmark_smoke_v2/results/${id}_t1000_out.pdbqt" '
  /^thread[[:space:]]*=/ {$0="thread = 1000"}
  /^search_depth[[:space:]]*=/ {next}
  /^out[[:space:]]*=/ {$0="out = " out}
  {print}
' "$cfg" > "benchmark_smoke_v2/configs/${id}_t1000.conf"
printf '%s\n' "$id" > benchmark_smoke_v2/SMOKE_TASK_ID.txt
start=$(date +%s); code=0
timeout 1800 "$BIN" --config "benchmark_smoke_v2/configs/${id}_t1000.conf" --seed 20260813 --log "benchmark_smoke_v2/logs/${id}_t1000.vina.log" > "benchmark_smoke_v2/logs/${id}_t1000.stdout.log" 2>&1 || code=$?
end=$(date +%s)
printf 'task_id\tthread\tseed\texit_code\telapsed_seconds\n%s\t1000\t20260813\t%s\t%s\n' "$id" "$code" "$((end-start))" > benchmark_smoke_v2/SMOKE_RUN_STATUS.tsv
cat benchmark_smoke_v2/SMOKE_RUN_STATUS.tsv
exit "$code"
