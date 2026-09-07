#!/usr/bin/env bash
#SBATCH --job-name=mempro_v631
#SBATCH --array=0-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
set -euo pipefail

: "${MEMPRO_ROOT:?Set MEMPRO_ROOT to the transferred project root}"
: "${VINA_BIN:?Set VINA_BIN to the AutoDock Vina executable}"

MAP="${MEMPRO_ROOT}/04_docking_hpc/docking_array_map_v631.tsv"
ROW=$(awk -F '\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==i {print; exit}' "$MAP")
test -n "$ROW"
IFS=$'\t' read -r ARRAY_INDEX TARGET_ID COMPOUND_ID BATCH_ID OUTPUT_STEM <<< "$ROW"

RECEPTOR="${MEMPRO_ROOT}/04_docking_hpc/receptors/${TARGET_ID}.pdbqt"
LIGAND="${MEMPRO_ROOT}/04_docking_hpc/ligands/${COMPOUND_ID}.pdbqt"
CONFIG="${MEMPRO_ROOT}/04_docking_hpc/configs/${TARGET_ID}.txt"
OUTDIR="${MEMPRO_ROOT}/04_docking_hpc/results/batch_${BATCH_ID}"
mkdir -p "$OUTDIR"

test -s "$RECEPTOR" && test -s "$LIGAND" && test -s "$CONFIG"
"$VINA_BIN" --receptor "$RECEPTOR" --ligand "$LIGAND" --config "$CONFIG" --out "$OUTDIR/${OUTPUT_STEM}.pdbqt" --log "$OUTDIR/${OUTPUT_STEM}.log"
