#!/bin/bash
# Prepare ligand PDBQT files from SMILES using OpenBabel
# Resume-safe: skips already-prepared files
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate docking
mkdir -p ligands_pdbqt logs

TOTAL=$(wc -l < ligands.smi)
echo "Converting $TOTAL compounds from SMILES to PDBQT..."
COUNT=0; SKIP=0; FAIL=0

while IFS=$'\t' read -r smiles cpd_id name; do
    # Strip CR
    smiles=$(echo "$smiles" | tr -d '\r')
    cpd_id=$(echo "$cpd_id" | tr -d '\r' | xargs)
    [ -z "$cpd_id" ] && continue

    out="ligands_pdbqt/${cpd_id}.pdbqt"
    if [ -f "$out" ] && [ -s "$out" ]; then
        SKIP=$((SKIP+1)); continue
    fi
    COUNT=$((COUNT+1))
    if [ $((COUNT % 200)) -eq 0 ]; then
        echo "  [$(ls ligands_pdbqt/*.pdbqt 2>/dev/null | wc -l)/$TOTAL] $cpd_id..."
    fi

    timeout 60 obabel -:"$smiles" -opdbqt -O "$out" --gen3d --fastest 2>&1
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "  FAILED: $cpd_id (exit=$rc)"
        rm -f "$out"
        echo "$cpd_id" >> logs/ligand_failed.txt
        FAIL=$((FAIL+1))
    fi
done < ligands.smi
echo "Done: $(ls ligands_pdbqt/*.pdbqt 2>/dev/null | wc -l) ligands ($FAIL failed)"
