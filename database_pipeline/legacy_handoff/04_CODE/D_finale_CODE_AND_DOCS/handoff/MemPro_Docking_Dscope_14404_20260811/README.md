# MemPro D-scope Docking pipeline ? 14,404 pairs

Frozen scope SHA256: `b2c8e07dbe38eb658d3955d8dc0eae7080af5b518dfa5edc69db6b7205b5e8a0`.

## One-command use
```bash
tar -xzf MemPro_Docking_Dscope_14404_20260811.tar.gz
cd MemPro_Docking_Dscope_14404_20260811
bash RUN_ME.sh --prepare
# Review smoke test and blocked queue
bash RUN_ME.sh --dock
# Or preparation + docking: bash RUN_ME.sh --all
```
Override paths with `VINA_DIR`, `VINA_BIN`, `CONDA_SH`, `CONDA_ENV`, `TASK_TIMEOUT`.

## Pocket routing
- G1 experimental residues: 2,871; automatically mapped and QC gated.
- G2 structural context without residues: 847; requires co-crystal ligand pocket extraction/validation.
- G3 PDB but no pocket: 9,429; requires P2Rank or equivalent prediction.
- G4 no experimental PDB: 1,257; requires frozen predicted structure plus pocket prediction.

No default whole-protein box is generated for G2?G4. Unsafe tasks stay in `metadata/blocked_after_grid.tsv`. Therefore this package processes all 14,404 pairs reproducibly, while only QC-approved tasks enter Vina. Exact pair/PDB/chain data are frozen in `docking_manifest.tsv`.

Verify before use: `sha256sum -c SHA256SUMS.txt`.
