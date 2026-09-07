#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source "${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
conda activate "${DOCKING_ENV:-docking}"
mkdir -p receptors_pdbqt receptors_pdbqt_clean logs receptor_cleanup_logs
find receptors_pdbqt_clean -maxdepth 1 -type f -name '*.pdbqt' -delete
: > logs/receptor_failed_v3.txt
printf 'pdb_id\tstatus\tblank_atom_count\n' > logs/receptor_prepare_v3.tsv
for pdb_file in pdb_structures/*.pdb; do
  pdb_id="$(basename "$pdb_file" .pdb)"
  out="receptors_pdbqt/${pdb_id}.pdbqt"
  tmp="${out}.tmp"
  if timeout 180 obabel "$pdb_file" -opdbqt -O "$tmp" -h >/dev/null 2>&1; then
    awk '/^(ATOM|HETATM|TER)/{print}' "$tmp" > "$out"
    rm -f "$tmp"
    blank=$(awk '/^(ATOM  |HETATM)/ && substr($0,78,2) ~ /^[[:space:]]*$/{n++} END{print n+0}' "$out")
    printf '%s\tPREPARED_UNMODIFIED\t%s\n' "$pdb_id" "$blank" >> logs/receptor_prepare_v3.tsv
  else
    rm -f "$out" "$tmp"
    printf '%s\tCONVERSION_FAILED\t\n' "$pdb_id" >> logs/receptor_prepare_v3.tsv
    printf '%s\n' "$pdb_id" >> logs/receptor_failed_v3.txt
  fi
done
echo 'Receptors were prepared without deleting unsupported atoms.'
echo 'Run 04b_finalize_receptor_qc.py after grid construction.'
