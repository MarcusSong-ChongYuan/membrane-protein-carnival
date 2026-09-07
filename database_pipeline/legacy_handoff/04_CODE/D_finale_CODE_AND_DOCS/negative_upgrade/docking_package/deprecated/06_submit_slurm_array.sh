#!/bin/bash
#SBATCH --job-name=vina_dock
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --array=1-390264%1000
#SBATCH --output=logs/vina_%A_%a.out
#SBATCH --error=logs/vina_%A_%a.err

set -e
cd "D:\finale\negative_upgrade\docking_package"
mkdir -p logs docking_results

CONF_FILE=$(printf "vina_configs/DOCK-%06d.conf" $SLURM_ARRAY_TASK_ID)
if [ ! -f "$CONF_FILE" ]; then
    echo "SKIP: $CONF_FILE not found"
    exit 0
fi
echo "=== $(date) Task $SLURM_ARRAY_TASK_ID ==="
head -1 "$CONF_FILE"

vina_gpu --config "$CONF_FILE" --gpu_id 0
echo "Done: $(date)"
