#!/usr/bin/env bash
set -u
cd /home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813 || exit 1
echo TASK_METADATA
grep -R -n -E 'DOCK-001624|DOCK-002555|DOCK-000166' --include='*.tsv' --include='*.csv' --include='*.json' . | head -60 || true
for task in DOCK-001624 DOCK-002555 DOCK-000166; do
  echo "LOGS_${task}"
  for f in production_batch_v1/logs/${task}*; do
    [ -f "$f" ] || continue
    echo "==== $f"
    tail -40 "$f"
  done
done
