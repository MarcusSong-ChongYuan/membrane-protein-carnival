#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813
VDIR=/home/csong/Vina-GPU
BIN=$VDIR/Vina-GPU
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}
mkdir -p "$ROOT"/benchmark_smoke_v3/{results,logs,configs}
cfg=$(find "$ROOT/vina_configs_chain_v2" -maxdepth 1 -name '*.conf' | sort | head -1)
id=$(basename "$cfg" .conf)
# Config paths must resolve from Vina-GPU's installation directory.
awk -v root="$ROOT" -v out="$ROOT/benchmark_smoke_v3/results/${id}_t1000_out.pdbqt" '
  /^receptor[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
  /^ligand[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
  /^thread[[:space:]]*=/ {$0="thread = 1000"}
  /^search_depth[[:space:]]*=/ {next}
  /^out[[:space:]]*=/ {$0="out = " out}
  {print}
' "$cfg" > "$ROOT/benchmark_smoke_v3/configs/${id}_t1000.conf"
printf '%s\n' "$id" > "$ROOT/benchmark_smoke_v3/SMOKE_TASK_ID.txt"
cd "$VDIR"
start=$(date +%s); code=0
timeout 1800 "$BIN" --config "$ROOT/benchmark_smoke_v3/configs/${id}_t1000.conf" --seed 20260813 --log "$ROOT/benchmark_smoke_v3/logs/${id}_t1000.vina.log" > "$ROOT/benchmark_smoke_v3/logs/${id}_t1000.stdout.log" 2>&1 || code=$?
end=$(date +%s)
printf 'task_id\tthread\tseed\texit_code\telapsed_seconds\n%s\t1000\t20260813\t%s\t%s\n' "$id" "$code" "$((end-start))" > "$ROOT/benchmark_smoke_v3/SMOKE_RUN_STATUS.tsv"
exit "$code"
