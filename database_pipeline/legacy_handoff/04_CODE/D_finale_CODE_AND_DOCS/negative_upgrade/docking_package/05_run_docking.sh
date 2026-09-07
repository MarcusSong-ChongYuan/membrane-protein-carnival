#!/bin/bash
# Run Vina-GPU docking tasks (resume-safe)
# Skips tasks whose receptor/ligand are not yet prepared
# IMPORTANT: Vina-GPU must be run from its source directory (OpenCL kernels)
# IMPORTANT: Uses receptors_pdbqt_clean (UNK/UNX atoms removed)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate docking

VINA_DIR="$HOME/Vina-GPU"
VINA_BIN="$VINA_DIR/Vina-GPU"
PROJ_DIR="$SCRIPT_DIR"
mkdir -p results logs

# Setup symlinks in Vina-GPU dir so relative paths in configs resolve
cd "$VINA_DIR"
ln -sf "$PROJ_DIR/receptors_pdbqt_clean" receptors_pdbqt
ln -sf "$PROJ_DIR/ligands_pdbqt" ligands_pdbqt
ln -sf "$PROJ_DIR/results" results

TOTAL=$(ls "$PROJ_DIR/vina_configs"/DOCK-*.conf 2>/dev/null | wc -l)
echo "Docking up to $TOTAL tasks..."
DONE=0; SKIP=0; FAIL=0

for conf in $(ls "$PROJ_DIR/vina_configs"/DOCK-*.conf | sort); do
    task_id=$(basename "$conf" .conf)
    out_file="$PROJ_DIR/results/${task_id}_out.pdbqt"
    if [ -f "$out_file" ] && [ -s "$out_file" ]; then
        SKIP=$((SKIP+1)); continue
    fi

    # Skip if receptor not ready
    rec=$(grep "^receptor =" "$conf" | awk '{print $3}')
    [ ! -f "$VINA_DIR/$rec" ] && SKIP=$((SKIP+1)) && continue

    # Skip if ligand not ready
    lig=$(grep "^ligand =" "$conf" | awk '{print $3}')
    [ ! -f "$VINA_DIR/$lig" ] && SKIP=$((SKIP+1)) && continue

    PROGRESS=$((DONE+SKIP+FAIL+1))
    if [ $((PROGRESS % 200)) -eq 0 ]; then
        echo "  [$PROGRESS/$TOTAL] done=$DONE fail=$FAIL skip=$SKIP"
    fi

    # Must run from Vina-GPU directory for kernel files
    cd "$VINA_DIR"
    timeout 60 "$VINA_BIN" --config "$conf" >> "$PROJ_DIR/logs/${task_id}.log" 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then
        DONE=$((DONE+1))
    else
        echo "  FAILED: $task_id (exit=$rc)"
        echo "$task_id" >> "$PROJ_DIR/logs/failed.txt"
        FAIL=$((FAIL+1))
    fi
done

echo ""
echo "========================================="
echo " DOCKING COMPLETE"
echo " Done: $DONE  Failed: $FAIL  Skipped: $SKIP"
echo " Results in: results/"
echo "========================================="
