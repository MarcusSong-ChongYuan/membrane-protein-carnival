#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2
PREP_PID="${1:?usage: A100_AUTOCHAIN_SAFE_V4.sh PREP_PID}"
LOG="$ROOT/logs/autochain_v4.log"
STATUS="$ROOT/metadata/AUTOCHAIN_V4_STATUS"
BIN="$HOME/Vina-GPU/Vina-GPU"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}"
cd "$ROOT"

printf 'WAITING_FOR_PREP pid=%s time=%s\n' "$PREP_PID" "$(date -Iseconds)" > "$STATUS"
printf '[%s] waiting for preparation pid=%s\n' "$(date -Iseconds)" "$PREP_PID" >> "$LOG"
while kill -0 "$PREP_PID" 2>/dev/null; do
  sleep 30
done

printf '[%s] preparation process ended; validating artifacts\n' "$(date -Iseconds)" >> "$LOG"
configs=$(find vina_configs -maxdepth 1 -type f -name '*.conf' -size +0c | wc -l)
clean_receptors=$(find receptors_pdbqt_clean -maxdepth 1 -type f -name '*.pdbqt' -size +0c | wc -l)
ligands=$(find ligands_pdbqt -maxdepth 1 -type f -name '*.pdbqt' -size +0c | wc -l)
if [ "$configs" -lt 5 ] || [ "$clean_receptors" -lt 1 ] || [ "$ligands" -lt 1 ]; then
  printf 'BLOCKED_PREP configs=%s clean_receptors=%s ligands=%s time=%s\n' "$configs" "$clean_receptors" "$ligands" "$(date -Iseconds)" > "$STATUS"
  exit 10
fi

printf 'SMOKE_RUNNING configs=%s time=%s\n' "$configs" "$(date -Iseconds)" > "$STATUS"
bash pipeline/06_smoke_test_v3.sh >> "$LOG" 2>&1

# The original benchmark configs reused the production output filename. Rewrite
# each derived benchmark config so the 5000- and 8000-thread poses remain distinct.
python3 - <<'PY'
from pathlib import Path
root=Path('/home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2')
for path in (root/'benchmark_configs').glob('*.conf'):
    lines=[]
    for line in path.read_text(encoding='ascii').splitlines():
        if line.startswith('out = '):
            line=f'out = results/{path.stem}_out.pdbqt'
        lines.append(line)
    path.write_text('\n'.join(lines)+'\n',encoding='ascii')
PY

printf 'BENCHMARK_RUNNING time=%s\n' "$(date -Iseconds)" > "$STATUS"
bash pipeline/06c_run_thread_benchmark_v3.sh >> "$LOG" 2>&1
python3 pipeline/a100_validate_benchmark_v4.py >> "$LOG" 2>&1

binary_hash=$(sha256sum "$BIN" | awk '{print $1}')
cat > metadata/PRODUCTION_READY <<EOF
review=automated_paired_thread_validation_v4
date=$(date -Iseconds)
vina_binary=$BIN
vina_binary_sha256=$binary_hash
accepted_profile=Vina-GPU thread=8000 num_modes=9 energy_range=3 max_box_axis=30A
benchmark_report=metadata/thread_benchmark_validation_v4.tsv
benchmark_summary=metadata/thread_benchmark_validation_v4.txt
EOF

printf 'DOCKING_STARTING configs=%s time=%s\n' "$configs" "$(date -Iseconds)" > "$STATUS"
nohup bash pipeline/07_run_docking_v3.sh > logs/production_v4_docking.log 2>&1 < /dev/null &
dock_pid=$!
printf 'DOCKING_RUNNING pid=%s configs=%s time=%s\n' "$dock_pid" "$configs" "$(date -Iseconds)" > "$STATUS"
printf '[%s] production docking launched pid=%s configs=%s\n' "$(date -Iseconds)" "$dock_pid" "$configs" >> "$LOG"
