#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2
cd "$ROOT"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate docking
mkdir -p receptors_pdbqt receptors_pdbqt_clean logs metadata
required=metadata/required_receptors_from_configs_20260817.txt
grep -h '^receptor[[:space:]]*=' vina_configs/*.conf | cut -d= -f2- | sed 's/^[[:space:]]*//' | sort -u > "$required"
printf 'receptor\tstatus\n' > logs/receptor_recovery_20260817.tsv

prep_one() {
  rel="$1"; pdb_id="$(basename "$rel" .pdbqt)"; src="$ROOT/pdb_structures/${pdb_id}.pdb"
  out="$ROOT/receptors_pdbqt/${pdb_id}.pdbqt"; tmp="${out}.tmp"
  if [ ! -s "$src" ]; then printf '%s\tMISSING_PDB\n' "$pdb_id" >> "$ROOT/logs/receptor_recovery_20260817.tsv"; return; fi
  if timeout 240 obabel -ipdb "$src" -opdbqt -O "$tmp" -h >/dev/null 2>&1; then
    awk '/^(ATOM|HETATM|TER)/{print}' "$tmp" > "$out"; rm -f "$tmp"
    if [ -s "$out" ]; then printf '%s\tPREPARED_UNMODIFIED\n' "$pdb_id" >> "$ROOT/logs/receptor_recovery_20260817.tsv"
    else rm -f "$out"; printf '%s\tEMPTY_OUTPUT\n' "$pdb_id" >> "$ROOT/logs/receptor_recovery_20260817.tsv"; fi
  else rm -f "$tmp" "$out"; printf '%s\tCONVERSION_FAILED\n' "$pdb_id" >> "$ROOT/logs/receptor_recovery_20260817.tsv"; fi
}
export -f prep_one
export ROOT
xargs -r -P "${PREP_JOBS:-8}" -n 1 bash -c 'prep_one "$1"' _ < "$required"

python3 pipeline/04b_finalize_receptor_qc_v31.py >> logs/receptor_recovery_20260817.log 2>&1
python3 - <<'PY'
import csv, re, shutil
from pathlib import Path
root=Path('/home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2')
with (root/'metadata/task_prep_qc.tsv').open(encoding='utf-8',newline='') as f:
    qc={r['task_id']:r['prep_qc_status'] for r in csv.DictReader(f,delimiter='\t')}
blocked=root/'vina_configs_blocked_recovery_20260817'; blocked.mkdir(exist_ok=True)
kept=moved=0
for cfg in sorted((root/'vina_configs').glob('*.conf')):
    text=cfg.read_text(encoding='ascii',errors='ignore'); rec=re.search(r'^receptor\s*=\s*(.+)$',text,re.M); lig=re.search(r'^ligand\s*=\s*(.+)$',text,re.M)
    ok=(qc.get(cfg.stem)=='PASS' and rec and lig and (root/'receptors_pdbqt_clean'/Path(rec.group(1).strip()).name).is_file() and (root/'ligands_pdbqt'/Path(lig.group(1).strip()).name).is_file())
    if ok: kept+=1
    else: shutil.move(str(cfg),str(blocked/cfg.name)); moved+=1
(root/'metadata/recovery_config_filter_20260817.txt').write_text(f'kept={kept}\nblocked_or_missing={moved}\n',encoding='utf-8')
print(f'kept={kept} blocked_or_missing={moved}')
PY
find vina_configs -maxdepth 1 -type f -name '*.conf' -size +0c | wc -l > metadata/recovered_config_count_20260817.txt
find receptors_pdbqt_clean -maxdepth 1 -type f -name '*.pdbqt' -size +0c | wc -l > metadata/recovered_clean_receptor_count_20260817.txt
echo "RECOVERY_COMPLETE $(date -Iseconds)" > metadata/RECOVERY_20260817_STATUS
