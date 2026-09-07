#!/bin/bash
# Prepare receptor PDBQT files from PDB structures
# Uses obabel for conversion, strips torsion tree tags,
# and removes UNK/UNX atoms with blank atom types (Vina-GPU compatible)
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCRIPT_DIR"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate docking

mkdir -p receptors_pdbqt receptors_pdbqt_clean logs
COUNT=0; SKIP=0; FAIL=0
TOTAL=$(ls pdb_structures/*.pdb 2>/dev/null | wc -l)
echo "Preparing $TOTAL receptors..."

for pdb_file in pdb_structures/*.pdb; do
    pdb_id=$(basename "$pdb_file" .pdb)
    out="receptors_pdbqt/${pdb_id}.pdbqt"
    [ -f "$out" ] && [ -s "$out" ] && SKIP=$((SKIP+1)) && continue
    COUNT=$((COUNT+1))
    if [ $((COUNT % 200)) -eq 0 ]; then
        echo "  [$(ls receptors_pdbqt/*.pdbqt 2>/dev/null | wc -l)/$TOTAL] $pdb_id..."
    fi

    tmp="${out}.tmp"
    if timeout 120 obabel "$pdb_file" -O "$tmp" -h 2>/dev/null; then
        # Strip torsion tree tags (ROOT/ENDBRANCH/BRANCH/TORSDOF) for receptor
        grep -P '^(ATOM|HETATM|TER)' "$tmp" > "$out"
        rm -f "$tmp"
    else
        echo "  FAILED: $pdb_id"
        rm -f "$out" "$tmp"
        echo "$pdb_id" >> logs/receptor_failed.txt
        FAIL=$((FAIL+1))
    fi
done

# Step 2: Clean receptors — remove UNK/UNX atoms with blank AutoDock atom types
# See: parse_pdbqt.cpp:71 assertion crash when atom type column (78-79) is blank
echo ""
echo "=== Cleaning receptors (removing UNK/UNX with blank atom types) ==="
mkdir -p receptors_pdbqt_clean
printf "pdbqt_file\tremoved_atom_count\n" > logs/receptor_cleanup.tsv
cleaned=0; total_bad=0

for f in receptors_pdbqt/*.pdbqt; do
    name=$(basename "$f")
    clean="receptors_pdbqt_clean/$name"

    count=$(awk '/^(ATOM  |HETATM)/ && substr($0,78,2) ~ /^[[:space:]]*$/' "$f" | wc -l)

    if [ "$count" -gt 0 ]; then
        total_bad=$((total_bad + count))
        awk '/^(ATOM  |HETATM)/ && substr($0,78,2) ~ /^[[:space:]]*$/ {next} {print}' "$f" > "$clean"
        cleaned=$((cleaned+1))
    else
        cp "$f" "$clean"
    fi
    printf "%s\t%s\n" "$name" "$count" >> logs/receptor_cleanup.tsv
done

echo "Cleaned $cleaned receptors ($total_bad atoms removed)"

# Step 3: Audit UNK/UNX removal (distance to pocket, coordinate logging)
if command -v python3 &>/dev/null && [ -f audit_unk_removal.py ]; then
    echo ""
    echo "=== Auditing UNK/UNX removal ==="
    python3 audit_unk_removal.py
fi

echo "Done: $(ls receptors_pdbqt_clean/*.pdbqt 2>/dev/null | wc -l) clean receptors, $FAIL failed"
