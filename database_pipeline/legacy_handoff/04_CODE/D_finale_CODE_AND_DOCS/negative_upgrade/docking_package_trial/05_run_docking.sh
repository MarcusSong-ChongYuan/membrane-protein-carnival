#!/bin/bash
# ============================================================
# Run 7,227 Vina-GPU docking tasks
# Resume-safe: re-run to skip completed, retry failed
# ============================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
conda activate docking
mkdir -p results logs

TOTAL=$(ls vina_configs/DOCK-*.conf 2>/dev/null | wc -l)
DONE=0
SKIP=0
FAIL=0

echo "Docking $TOTAL tasks on GPU 0..."
echo ""

for conf in $(ls vina_configs/DOCK-*.conf | sort); do
    task_id=$(basename "$conf" .conf)
    out_file="results/${task_id}_out.pdbqt"

    if [ -f "$out_file" ] && [ -s "$out_file" ]; then
        SKIP=$((SKIP+1))
        continue
    fi

    PROGRESS=$((DONE+SKIP+FAIL+1))
    if [ $((PROGRESS % 100)) -eq 0 ]; then
        echo "  [$PROGRESS/$TOTAL] done=$DONE fail=$FAIL skip=$SKIP"
    fi

    if vina_gpu --config "$conf" --gpu_id 0 \
        >> "logs/${task_id}.log" 2>&1; then
        DONE=$((DONE+1))
    else
        echo "  FAILED: $task_id"
        echo "$task_id" >> logs/failed.txt
        FAIL=$((FAIL+1))
    fi
done

echo ""
echo "========================================="
echo " DOCKING COMPLETE"
echo " Done: $DONE  Failed: $FAIL  Skipped: $SKIP"
echo " Results in: results/"
echo "========================================="
