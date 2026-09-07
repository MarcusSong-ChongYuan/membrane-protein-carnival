#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; VDIR="${VINA_DIR:-$HOME/Vina-GPU}"; BIN="${VINA_BIN:-$VDIR/Vina-GPU}"; LIMIT="${TASK_TIMEOUT:-1800}"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}"
[ -f "$ROOT/metadata/PRODUCTION_READY" ] || { echo 'Missing reviewed metadata/PRODUCTION_READY gate' >&2; exit 4; }
[ -x "$BIN" ] || { echo "Missing Vina executable: $BIN" >&2; exit 3; }
mkdir -p "$ROOT/results" "$ROOT/logs"; : > "$ROOT/logs/failed.current.txt"
ln -sfn "$ROOT/receptors_pdbqt_clean" "$VDIR/receptors_pdbqt"; ln -sfn "$ROOT/ligands_pdbqt" "$VDIR/ligands_pdbqt"; ln -sfn "$ROOT/results" "$VDIR/results"
mapfile -t cfgs < <(find "$ROOT/vina_configs" -type f -name '*.conf' | sort)
done_n=0; fail_n=0; skip_n=0
for cfg in "${cfgs[@]}"; do
  id="$(basename "$cfg" .conf)"; [ -s "$ROOT/results/${id}_out.pdbqt" ] && { skip_n=$((skip_n+1)); continue; }
  if (cd "$VDIR" && timeout "$LIMIT" "$BIN" --config "$cfg") >"$ROOT/logs/${id}.log" 2>&1; then done_n=$((done_n+1)); else fail_n=$((fail_n+1)); echo "$id" >> "$ROOT/logs/failed.current.txt"; fi
done
echo "done=$done_n fail=$fail_n skip=$skip_n total=${#cfgs[@]}"
