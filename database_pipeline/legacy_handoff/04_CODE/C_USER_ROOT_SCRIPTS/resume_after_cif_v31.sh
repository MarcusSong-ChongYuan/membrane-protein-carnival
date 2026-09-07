#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"; cd "$ROOT"
source /home/csong/miniconda3/etc/profile.d/conda.sh
conda activate docking
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:${LD_LIBRARY_PATH:-}
bash pipeline/02_prepare_proteins_v3.sh
bash pipeline/03_prepare_ligands.sh
python3 pipeline/04_compute_grid_boxes_v3.py
python3 pipeline/04b_finalize_receptor_qc_v31.py
python3 pipeline/05_generate_configs_v3.py
python3 pipeline/06b_benchmark_threads_v3.py
