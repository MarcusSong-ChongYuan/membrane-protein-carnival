#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2
cd "$ROOT"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
bash pipeline/02_prepare_proteins_v3.sh
bash pipeline/03_prepare_ligands.sh
python3 pipeline/04_compute_grid_boxes_v3.py
python3 pipeline/04b_finalize_receptor_qc.py
python3 pipeline/05_generate_configs_v3.py
python3 pipeline/06b_benchmark_threads_v3.py
