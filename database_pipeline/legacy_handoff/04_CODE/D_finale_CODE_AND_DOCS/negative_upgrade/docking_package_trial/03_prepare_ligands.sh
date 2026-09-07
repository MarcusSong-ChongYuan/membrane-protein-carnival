#!/bin/bash
# ============================================================
# Prepare 547 ligand PDBQT files using OpenBabel
# Resume-safe: skips already-prepared files
# ============================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
conda activate docking
mkdir -p ligands_pdbqt

TOTAL=$(wc -l < ligands.smi)
echo "Converting $TOTAL compounds from SMILES to PDBQT..."
COUNT=0
FAIL=0

while IFS=$'\t' read -r smiles cpd_id name; do
    out="ligands_pdbqt/${cpd_id}.pdbqt"
    if [ -f "$out" ] && [ -s "$out" ]; then
        continue
    fi
    COUNT=$((COUNT+1))
    if [ $((COUNT % 100)) -eq 0 ]; then
        echo "  [$COUNT/$TOTAL]"
    fi
    if ! echo "$smiles" | obabel -ismi -opdbqt -O "$out" \
        --gen3d --conformers --nconf 1 2>/dev/null; then
        FAIL=$((FAIL+1))
        rm -f "$out"
    fi
done < ligands.smi
echo "Prepared $COUNT ligands ($FAIL failed)"
