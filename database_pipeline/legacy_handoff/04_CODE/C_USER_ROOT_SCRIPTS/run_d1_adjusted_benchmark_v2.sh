#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813
OUT=$ROOT/benchmark_v2_adjusted
BIN=/home/csong/Vina-GPU/Vina-GPU
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}
mkdir -p "$OUT"/{results,logs,status}
if nvidia-smi --query-compute-apps=used_memory --format=csv,noheader,nounits | awk '$1>10000{bad=1} END{exit bad?0:1}'; then
  echo 'GPU_BUSY_OVER_10GB; benchmark not launched' >&2; exit 75
fi
python3 "$ROOT/d1_build_benchmark_v2.py"
printf 'task_id\tseed\tgpu\tthread\tconfig\n' > "$OUT/RUNS.tsv"
i=0
while IFS=$'\t' read -r task thread disp; do
 [[ "$task" == task_id ]] && continue
 [[ "$disp" == READY_FOR_ADJUSTED_BENCHMARK ]] || continue
 for seed in 20260813 20260814 20260815; do gpu=$((i%2)); printf '%s\t%s\t%s\t%s\t%s\n' "$task" "$seed" "$gpu" "$thread" "$OUT/configs/$task.conf" >> "$OUT/RUNS.tsv";i=$((i+1));done
done < <(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="task_id")a=i;if($i=="thread_v2")b=i;if($i=="benchmark_v2_disposition")c=i};print "task_id\tthread\tdisp";next}{print $a"\t"$b"\t"$c}' "$OUT/D1_ADJUSTED_BENCHMARK_TASKS.tsv")
worker(){
 local gpu=$1
 awk -F '\t' -v g="$gpu" 'NR>1&&$3==g{print}' "$OUT/RUNS.tsv" | while IFS=$'\t' read -r task seed gp thread cfg; do
  status="$OUT/status/${task}_s${seed}.tsv";code=1;attempt=0;start=$(date +%s)
  for attempt in 1 2 3; do
   log="$OUT/logs/${task}_s${seed}_a${attempt}.stdout.log";vlog="$OUT/logs/${task}_s${seed}_a${attempt}.vina.log"
   code=0;CUDA_VISIBLE_DEVICES="$gpu" timeout 7200 "$BIN" --config "$cfg" --seed "$seed" --log "$vlog" >"$log" 2>&1 || code=$?
   [[ $code -eq 0 ]] && break
   grep -qE 'CL_BUILD_PROGRAM_FAILURE|Err-11' "$log" || break
   sleep 2
  done
  end=$(date +%s);log="$OUT/logs/${task}_s${seed}_a${attempt}.stdout.log";aff=$(awk '/^[[:space:]]*1[[:space:]]+[-0-9.]+/{print $2;exit}' "$log")
  printf 'task_id\tseed\tgpu\tthread\tattempts\texit_code\telapsed_seconds\tbest_affinity_kcal_mol\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$seed" "$gpu" "$thread" "$attempt" "$code" "$((end-start))" "${aff:-}" > "$status"
 done
}
worker 0 & p0=$!;worker 1 & p1=$!;wait "$p0" "$p1"
python3 "$ROOT/evaluate_d1_adjusted_benchmark_v2.py"
