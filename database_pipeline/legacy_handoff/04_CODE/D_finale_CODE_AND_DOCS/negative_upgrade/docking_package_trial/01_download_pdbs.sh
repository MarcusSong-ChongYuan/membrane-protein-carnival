#!/bin/bash
# ============================================================
# Download 3,226 PDB structures from PDBe/RCSB
# Resume-safe: skips already-downloaded files
# ============================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p pdb_structures

TOTAL=$(wc -l < pdb_download_list.txt)
echo "Downloading $TOTAL PDB files..."
DONE=0
FAIL=0

while read -r pdb; do
    if [ -f "pdb_structures/${pdb}.pdb" ] && [ -s "pdb_structures/${pdb}.pdb" ]; then
        DONE=$((DONE+1))
        continue
    fi
    DONE=$((DONE+1))
    if [ $((DONE % 100)) -eq 0 ]; then
        echo "  [$DONE/$TOTAL] $pdb"
    fi
    # Try PDBe first, fallback to RCSB
    curl -sf --connect-timeout 10 --max-time 60 \
        -o "pdb_structures/${pdb}.pdb" \
        "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${pdb}.ent" || \
    curl -sf --connect-timeout 10 --max-time 60 \
        -o "pdb_structures/${pdb}.pdb" \
        "https://files.rcsb.org/download/${pdb}.pdb" || {
        echo "  FAILED: $pdb"
        FAIL=$((FAIL+1))
        rm -f "pdb_structures/${pdb}.pdb"
    }
done < pdb_download_list.txt

DOWNLOADED=$(ls pdb_structures/*.pdb 2>/dev/null | wc -l)
echo ""
echo "Download complete: $DOWNLOADED PDBs ($FAIL failed)"
