set -Eeuo pipefail
cd /home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2
stamp=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "benchmark_archive/$stamp"
find benchmark_results -maxdepth 1 -type f -exec mv -t "benchmark_archive/$stamp" {} + 2>/dev/null || true
find benchmark_logs -maxdepth 1 -type f -exec mv -t "benchmark_archive/$stamp" {} + 2>/dev/null || true
python3 pipeline/06b_benchmark_threads_v3.py
head -12 "$(find benchmark_configs -type f | sort | head -1)"
bash pipeline/06c_run_thread_benchmark_v3.sh
python3 pipeline/a100_validate_benchmark_v4.py
