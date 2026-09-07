#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2
cd "$ROOT"
stamp=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "benchmark_archive/$stamp"
find benchmark_results -maxdepth 1 -type f -exec mv -t "benchmark_archive/$stamp" {} + 2>/dev/null || true
find benchmark_logs -maxdepth 1 -type f -exec mv -t "benchmark_archive/$stamp" {} + 2>/dev/null || true
bash pipeline/06c_run_thread_benchmark_v3.sh
if python3 pipeline/a100_validate_benchmark_v4.py; then
  mkdir -p metadata production_archive
  {
    echo "status=PRODUCTION_READY"
    echo "created_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "benchmark=thread_5000_vs_8000_same_seed"
    echo "config_count=$(find vina_configs -maxdepth 1 -type f -name '*.conf' | wc -l)"
  } > metadata/PRODUCTION_READY
  prodstamp=$(date -u +%Y%m%dT%H%M%SZ)
  mkdir -p "production_archive/$prodstamp"
  find results -maxdepth 1 -type f -exec mv -t "production_archive/$prodstamp" {} + 2>/dev/null || true
  CUDA_VISIBLE_DEVICES=0 bash pipeline/07_run_docking_v3.sh
else
  echo "Benchmark gate failed; production not started."
  exit 5
fi
