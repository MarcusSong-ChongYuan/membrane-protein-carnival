#!/usr/bin/env bash
set -Eeuo pipefail
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_DIR"
PY="$PROJECT_DIR/.venv/bin/python"
LOG="$PROJECT_DIR/pipeline_a100.log"
exec >>"$LOG" 2>&1

echo "$(date -Is) installing analysis packages"
"$PY" -m pip install fair-esm hdbscan umap-learn pandas scikit-learn matplotlib seaborn plotly pyyaml
echo "$(date -Is) package installation passed"

# GPU 1 is intentionally selected because GPU 0 belongs to an unrelated job.
export CUDA_VISIBLE_DEVICES=1
export TORCH_HOME="$PROJECT_DIR/model_cache"
echo "$(date -Is) starting restart-safe pipeline on CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
exec "$PY" run_pipeline.py
