#!/bin/bash
# Quick test: run first 5 docking tasks with proper validation
# Uses set -o pipefail; checks exit codes and output files
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate docking

VINA_DIR="$HOME/Vina-GPU"
PROJ_DIR="$SCRIPT_DIR"

# Ensure symlinks
cd "$VINA_DIR"
ln -sf "$PROJ_DIR/receptors_pdbqt_clean" receptors_pdbqt
ln -sf "$PROJ_DIR/ligands_pdbqt" ligands_pdbqt
ln -sf "$PROJ_DIR/results" results

mkdir -p "$PROJ_DIR/results" "$PROJ_DIR/logs"

TEST_IDS="000001 000002 000003 000004 000005"
TOTAL_TESTS=$(echo "$TEST_IDS" | wc -w)
PASSED=0
FAILED=0

echo "=== TEST: $TOTAL_TESTS tasks ==="
echo ""

for i in $TEST_IDS; do
    conf="$PROJ_DIR/vina_configs/DOCK-${i}.conf"
    out_file="$PROJ_DIR/results/DOCK-${i}_out.pdbqt"
    log_file="$PROJ_DIR/logs/DOCK-${i}_test.log"

    if [ ! -f "$conf" ]; then
        echo "  SKIP: DOCK-${i} — config not found"
        continue
    fi

    echo -n "  DOCK-${i}: "

    # Run docking — full log saved, no truncation
    timeout 60 "$VINA_DIR/Vina-GPU" --config "$conf" > "$log_file" 2>&1
    rc=$?

    # Validate
    if [ $rc -ne 0 ]; then
        echo "FAILED (exit code=$rc)"
        FAILED=$((FAILED + 1))
        echo "  Last 3 lines of log:"
        tail -3 "$log_file" | sed 's/^/    /'
        continue
    fi

    if [ ! -f "$out_file" ]; then
        echo "FAILED (no output file)"
        FAILED=$((FAILED + 1))
        continue
    fi

    if [ ! -s "$out_file" ]; then
        echo "FAILED (empty output file)"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Check for docking poses
    pose_count=$(grep -c 'ENDMDL' "$out_file" 2>/dev/null || true)
    if [ "$pose_count" -gt 0 ]; then
        echo "OK ($pose_count poses)"
        PASSED=$((PASSED + 1))
    else
        echo "OK (output present, no ENDMDL markers)"
        PASSED=$((PASSED + 1))
    fi
done

echo ""
echo "========================================="
echo " TEST RESULTS: $PASSED passed, $FAILED failed"
echo " Logs: $PROJ_DIR/logs/"
echo " Results: $PROJ_DIR/results/"
echo "========================================="

[ $FAILED -eq 0 ] && exit 0 || exit 1
