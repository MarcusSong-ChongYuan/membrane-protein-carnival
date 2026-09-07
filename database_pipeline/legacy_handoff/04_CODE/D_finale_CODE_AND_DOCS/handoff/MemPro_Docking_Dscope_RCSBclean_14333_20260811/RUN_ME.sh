#!/usr/bin/env bash
set -Eeuo pipefail
R="$(cd "$(dirname "$0")"&&pwd)";M="${1:---prepare}";cd "$R";set +u;source "${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}";conda activate "${CONDA_ENV:-docking}";set -u;python3 pipeline/00_validate_and_stage.py
if [ "$M" != '--dock' ];then bash pipeline/01_download_pdbs.sh;bash pipeline/03_prepare_ligands.sh;bash pipeline/02_prepare_proteins.sh;python3 pipeline/download_sifts_mappings.py||echo 'SIFTS incomplete; resumable';python3 pipeline/04_compute_grid_boxes.py;python3 pipeline/05_generate_configs.py;bash pipeline/06_smoke_test.sh;fi
if [ "$M" = '--dock' ]||[ "$M" = '--all' ];then bash pipeline/07_run_docking.sh;else echo 'Prepared. Review metadata/blocked_after_grid.tsv; run bash RUN_ME.sh --dock';fi
