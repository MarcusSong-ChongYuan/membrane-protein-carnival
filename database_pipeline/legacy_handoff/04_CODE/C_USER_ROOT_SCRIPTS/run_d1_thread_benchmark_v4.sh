#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813
VDIR=/home/csong/Vina-GPU
BIN=$VDIR/Vina-GPU
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}
OUT=$ROOT/thread_benchmark_v4
mkdir -p "$OUT"/{configs,results,logs}
cd "$ROOT"
# Deterministically select small/median/large receptor atom-count tasks among 27 static-QC passes.
python - <<'PY'
import csv
p='benchmark_qc_chain_v2/D1_27_BENCHMARK_STATIC_QC.tsv'
with open(p,encoding='utf-8') as f:r=sorted(csv.DictReader(f,delimiter='\t'),key=lambda x:int(x['receptor_atom_count']))
sel=[r[0],r[len(r)//2],r[-1]]
with open('thread_benchmark_v4/SELECTED_TASKS.tsv','w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=r[0].keys(),delimiter='\t');w.writeheader();w.writerows(sel)
PY
printf 'task_id\treceptor_atoms\tthread\tseed\texit_code\telapsed_seconds\tbest_affinity_kcal_mol\n' > "$OUT/BENCHMARK_RUNS.tsv"
while IFS=$'\t' read -r task pdb chain uniprot cmpd name cfgex recex atoms rb ligex la lb bx by bz bmax frac white status reasons; do
  [[ "$task" == task_id ]] && continue
  src="$ROOT/vina_configs_chain_v2/${task}.conf"
  for thread in 1000 5000 8000; do
    cfg="$OUT/configs/${task}_t${thread}.conf"
    awk -v root="$ROOT" -v out="$OUT/results/${task}_t${thread}_out.pdbqt" -v th="$thread" '
      /^receptor[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
      /^ligand[[:space:]]*=/ {sub(/=[[:space:]]*/,"= " root "/",$0)}
      /^thread[[:space:]]*=/ {$0="thread = " th}
      /^search_depth[[:space:]]*=/ {next}
      /^out[[:space:]]*=/ {$0="out = " out}
      {print}
    ' "$src" > "$cfg"
    start=$(date +%s);code=0
    cd "$VDIR"
    timeout 1800 "$BIN" --config "$cfg" --seed 20260813 --log "$OUT/logs/${task}_t${thread}.vina.log" > "$OUT/logs/${task}_t${thread}.stdout.log" 2>&1 || code=$?
    cd "$ROOT";end=$(date +%s)
    affinity=$(awk '/^[[:space:]]*1[[:space:]]+[-0-9.]+/{print $2;exit}' "$OUT/logs/${task}_t${thread}.stdout.log")
    printf '%s\t%s\t%s\t20260813\t%s\t%s\t%s\n' "$task" "$atoms" "$thread" "$code" "$((end-start))" "${affinity:-}" >> "$OUT/BENCHMARK_RUNS.tsv"
  done
done < "$OUT/SELECTED_TASKS.tsv"
cat "$OUT/BENCHMARK_RUNS.tsv"
