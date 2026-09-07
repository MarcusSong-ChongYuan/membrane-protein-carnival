#!/bin/bash
# ============================================================
# Prepare 3,226 receptor PDBQT files using ADFR suite
# Resume-safe: skips already-prepared files
# ============================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
conda activate docking
mkdir -p receptors_pdbqt

COUNT=0
TOTAL=$(ls pdb_structures/*.pdb 2>/dev/null | wc -l)
echo "Preparing $TOTAL receptors..."

for pdb_file in pdb_structures/*.pdb; do
    pdb_id=$(basename "$pdb_file" .pdb)
    out="receptors_pdbqt/${pdb_id}.pdbqt"
    if [ -f "$out" ] && [ -s "$out" ]; then
        continue
    fi
    COUNT=$((COUNT+1))
    if [ $((COUNT % 200)) -eq 0 ]; then
        echo "  [$COUNT] $pdb_id"
    fi
    prepare_receptor -r "$pdb_file" -o "$out" -A hydrogens -U nphs 2>/dev/null
done
echo "Prepared $COUNT receptors"
