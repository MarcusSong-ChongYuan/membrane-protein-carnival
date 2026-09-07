#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2
cd "$ROOT"
prep_pid="${1:?prep pid required}"
echo "WAITING_FOR_RECOVERY prep_pid=$prep_pid time=$(date -Iseconds)" > metadata/AUTOCHAIN_V4_STATUS
while [ ! -f metadata/RECOVERY_20260817_STATUS ]; do
  sleep 30
done
# The recovery checkpoint is authoritative. Passing the now-finished preparation
# PID makes the existing validated autochain continue immediately.
bash A100_AUTOCHAIN_SAFE_V4.sh "$prep_pid"
