#!/bin/bash
# Download PDB structures from RCSB (PDB + CIF fallback)
# Resume-safe: skips already-downloaded files
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate docking
mkdir -p pdb_structures logs

TOTAL=$(wc -l < pdb_download_list.txt)
echo "Downloading $TOTAL PDB files..."
DONE=0; FAIL=0

while read -r pdb; do
    # Strip CR and whitespace
    pdb=$(echo "$pdb" | tr -d '\r' | xargs)
    [ -z "$pdb" ] && continue

    out="pdb_structures/${pdb}.pdb"
    if [ -f "$out" ] && [ -s "$out" ]; then
        DONE=$((DONE+1)); continue
    fi

    DONE=$((DONE+1))
    if [ $((DONE % 200)) -eq 0 ]; then
        echo "  [$DONE/$TOTAL] $pdb (fail=$FAIL)"
    fi

    # Try PDB format first (RCSB)
    if curl -sf --connect-timeout 10 --max-time 60 \
        -o "$out" "https://files.rcsb.org/download/${pdb}.pdb" 2>/dev/null; then
        continue
    fi

    # Try CIF format (newer structures - RCSB)
    tmp_cif="/tmp/${pdb}.cif"
    if curl -sf --connect-timeout 10 --max-time 60 \
        -o "$tmp_cif" "https://files.rcsb.org/download/${pdb}.cif" 2>/dev/null; then
        if timeout 120 obabel "$tmp_cif" -O "$out" 2>/dev/null && [ -s "$out" ]; then
            rm -f "$tmp_cif"
            continue
        fi
        rm -f "$tmp_cif"
    fi

    # Try PDBe (PDB format)
    if curl -sf --connect-timeout 10 --max-time 60 \
        -o "$out" "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${pdb}.ent" 2>/dev/null; then
        continue
    fi

    # Try PDBe (CIF format)
    if curl -sf --connect-timeout 10 --max-time 60 \
        -o "$tmp_cif" "https://www.ebi.ac.uk/pdbe/entry-files/download/${pdb}.cif" 2>/dev/null; then
        if timeout 120 obabel "$tmp_cif" -O "$out" 2>/dev/null && [ -s "$out" ]; then
            rm -f "$tmp_cif"
            continue
        fi
        rm -f "$tmp_cif"
    fi

    echo "  FAILED: $pdb"
    echo "$pdb" >> logs/failed_ids.txt
    FAIL=$((FAIL+1))
    rm -f "$out" "$tmp_cif"
done < pdb_download_list.txt

echo ""
echo "Download complete: $(ls pdb_structures/*.pdb 2>/dev/null | wc -l) PDBs ($FAIL failed)"
