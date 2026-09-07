#!/bin/bash
# ============================================================
# Compute Vina grid boxes from binding site residues
# ============================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
conda activate docking
python3 compute_grid_boxes.py
