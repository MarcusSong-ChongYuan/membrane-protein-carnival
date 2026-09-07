# Deprecated Scripts

These scripts were part of an earlier Windows-based pipeline and contain
hardcoded paths, unsupported flags, and incompatible assumptions.

## Why Removed

| Script | Issues |
|--------|--------|
| `06_submit_slurm_array.sh` | Windows path `D:\finale\...`; array range 390,264 mismatched; `%1000` unrealistic; uses `--gpu_id 0` (unsupported) |
| `06_run_batch_vina_gpu.py` | Windows paths; wrong working directory; `vina_gpu` command not installed; uses `--gpu_id 0` |

## Replacement

- For single-task testing: use `06_test_run.sh`
- For batch docking: use `05_run_docking.sh`
- For Slurm submission: rewrite from scratch against A100 cluster config
