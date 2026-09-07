#!/usr/bin/env bash
set -Eeuo pipefail
R="$(cd "$(dirname "$0")/.."&&pwd)";D="${VINA_DIR:-$HOME/Vina-GPU}";B="${VINA_BIN:-$D/Vina-GPU}";mapfile -t C < <(find "$R/vina_configs" -name '*.conf' -type f|sort|head -n 5);[ ${#C[@]} -gt 0 ]||exit 2;ln -sfn "$R/receptors_pdbqt_clean" "$D/receptors_pdbqt";ln -sfn "$R/ligands_pdbqt" "$D/ligands_pdbqt";ln -sfn "$R/results" "$D/results";cd "$D";for c in "${C[@]}";do timeout 180 "$B" --config "$c";done
