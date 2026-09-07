# Vina-GPU Docking Package (A100) — Fixed for Deployment

Generated: 2026-08-08 | Fixed: 2026-08-10
Source: MemPro V6.2 — Biological Status Compounds × Membrane Protein Binding Evidence

## Scope

| Metric | Count |
|--------|-------|
| Docking tasks (with PDB) | 124,332 |
| Unique compounds | 5,563 |
| Unique PDB structures | 3,961 |
| Membrane protein targets | 3,385 |
| AlphaFold-needed (no PDB) | 27,771 tasks / 1,106 proteins |
| Total pairs in manifest | 152,103 |

## Server Requirements

- **OS**: Ubuntu Linux
- **GPU**: NVIDIA A100-SXM4-80GB (2×)
- **CUDA**: 12.4
- **Driver**: 550+
- **No sudo required**

## Software Stack (pre-installed via conda)

| Tool | Purpose | Source |
|------|---------|--------|
| Miniconda3 | Package management | `~/miniconda3/` |
| Python 3.11 | Grid box computation | conda env `docking` |
| OpenBabel (obabel) | Protein + ligand PDBQT preparation | conda env `docking` |
| Vina-GPU | GPU-accelerated docking | compiled in `~/Vina-GPU/` |
| Boost 1.85 | Vina-GPU dependency | conda env `docking` |
| meeko | Alternative PDBQT prep (fallback) | pip in conda env |

## Key Environment Variable

Automatically set by conda env activate:
```bash
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:$LD_LIBRARY_PATH
```

## Files

```
docking_package/
├── docking_manifest.tsv          # 152,103 pairs with all metadata
├── ligands.smi                   # 5,563 SMILES (tab-separated: SMILES\tID\tName)
├── pdb_download_list.txt         # 3,961 PDB IDs to download
├── vina_configs/                 # 124,332 .conf files (grid box pre-computed)
├── 01_download_pdbs.sh           # Download PDBs (PDB + CIF fallback, ~4h)
├── 02_prepare_proteins.sh        # PDB → receptor PDBQT (obabel, ~2-4h)
├── 03_prepare_ligands.sh         # SMILES → ligand PDBQT (obabel, ~30min)
├── 04_compute_grid_boxes.py      # Residue → grid box coordinates
├── 05_run_docking.sh             # Sequential GPU docking (10-20h)
├── 06_test_run.sh                # Quick 5-task smoke test
├── fix_all_configs.py           # Config file patcher (run locally before transfer)
└── README_DOCKING.md             # This file
```

## Workflow (A100 Server)

All commands from `~/docking_package/` directory.

```bash
# 0. Transfer folder to A100 first!
#    scp -r docking_package/ csong@star:~/docking_package/

cd ~/docking_package

# 1. Download PDB structures (~4 hours, 3,961 PDBs, ~12 GB)
nohup bash 01_download_pdbs.sh > logs/download.log 2>&1 &

# 2. Prepare proteins (~2-4 hours, CPU)
nohup bash 02_prepare_proteins.sh > logs/receptors.log 2>&1 &

# 3. Prepare ligands (~30 min, CPU)
bash 03_prepare_ligands.sh

# 4. Compute grid boxes (~10 min, CPU)
python3 04_compute_grid_boxes.py

# 5. Smoke test (check GPU availability first with nvidia-smi)
bash 06_test_run.sh

# 6. Full docking (~10-40 hours, GPU)
nohup bash 05_run_docking.sh > logs/docking.log 2>&1 &
```

## Changes from Original (2026-08-06 version)

| Issue | Original | Fixed |
|-------|----------|-------|
| Protein prep | ADFR `prepare_receptor` (unavailable) | OpenBabel `obabel -h` + strip torsions |
| Ligand prep | `echo $smiles \| obabel -ismi` (broken) | `obabel -:"$smiles" -opdbqt` |
| conda activate | `conda activate docking` (fails in scripts) | `source ~/miniconda3/etc/profile.d/conda.sh` + `conda activate` |
| Vina binary | `vina_gpu` (not in PATH) | `$HOME/Vina-GPU/Vina-GPU` (absolute) |
| GPU selection | `--gpu_id 0` (not supported) | Auto-detected via OpenCL |
| Config: exhaustiveness | `exhaustiveness = 32` | Removed (Vina-GPU uses `--thread` + `--search_depth`) |
| Config: log | `log = ...` | Removed (Vina-GPU uses `--log` argument) |
| Config: paths | `../receptors_pdbqt/` etc. | Relative to project root |
| Config: exhaustiveness | CRLF line endings | LF |
| PDB download | `.pdb` only (404 for new structures) | `.cif` fallback via obabel conversion |

## Vina-GPU Compilation Notes

Located at `~/Vina-GPU/Vina-GPU`. Compiled with Boost 1.85 + patches:
1. `boost/filesystem/convenience.hpp` shim (header removed in 1.85)
2. `-DBOOST_TIMER_ENABLE_DEPRECATED` flag (deprecated→error in 1.85)
3. Include path: `-I$CONDA_PREFIX/include`

See `A100_docking_setup_guide.md` for full compilation details.

## Docking Parameters

| Parameter | Value |
|-----------|-------|
| num_modes | 9 |
| energy_range | 3 kcal/mol |
| Grid box | Residue-based with 10 Å margin, min 22 Å |
| Default grid | 20×20×20 Å (for tasks without residue data) |
| Scoring | Vina empirical scoring function (gauss1, gauss2, repulsion, hydrophobic, hydrogen) |
