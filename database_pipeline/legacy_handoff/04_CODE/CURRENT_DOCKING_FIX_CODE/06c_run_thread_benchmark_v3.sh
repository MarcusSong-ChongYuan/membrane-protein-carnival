#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; VDIR="${VINA_DIR:-$HOME/Vina-GPU}"; BIN="${VINA_BIN:-$VDIR/Vina-GPU}"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}"
[ -x "$BIN" ] || { echo "Missing Vina executable: $BIN" >&2; exit 3; }
mkdir -p "$ROOT/benchmark_results" "$ROOT/benchmark_logs"
ln -sfn "$ROOT/receptors_pdbqt_clean" "$VDIR/receptors_pdbqt"; ln -sfn "$ROOT/ligands_pdbqt" "$VDIR/ligands_pdbqt"; ln -sfn "$ROOT/benchmark_results" "$VDIR/results"
printf 'config\texit_code\telapsed_seconds\n' > "$ROOT/metadata/thread_benchmark_run.tsv"
for cfg in "$ROOT"/benchmark_configs/*.conf; do
  id="$(basename "$cfg" .conf)"; start=$(date +%s); code=0
  (cd "$VDIR" && timeout 1800 "$BIN" --config "$cfg") >"$ROOT/benchmark_logs/${id}.log" 2>&1 || code=$?
  end=$(date +%s); printf '%s\t%s\t%s\n' "$id" "$code" "$((end-start))" >> "$ROOT/metadata/thread_benchmark_run.tsv"
done
echo 'Benchmark complete. Compare paired scores/poses before creating PRODUCTION_READY.'
