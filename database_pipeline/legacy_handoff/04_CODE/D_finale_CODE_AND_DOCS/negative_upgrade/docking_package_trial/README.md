# MemPro Docking - Trial Package

## What is this?

This package contains everything needed to run **7,227 molecular docking calculations** on an A100 GPU server. Each calculation predicts how well a small molecule (drug compound) binds to a membrane protein target.

This is the **simplest, highest-quality subset** of the full MemPro database — every single task has:
- **BE1 evidence** (PDB co-crystal structure confirmation)
- **Binding site residue data** (we know exactly which amino acids form the pocket)
- **One PDB structure per task** (no cross-docking complexity)

## Package Stats

| Item | Count |
|------|-------|
| Docking tasks | 7,227 |
| Unique compounds (ligands) | 547 |
| Unique PDB structures | 3,226 |
| Unique protein targets | 707 |
| Estimated disk needed | ~21 GB |
| Package size (without PDBs) | ~5 MB |

## Files

```
docking_package_trial/
  docking_manifest.tsv          Master task list with all metadata
  ligands.smi                   SMILES strings for 547 compounds
  pdb_download_list.txt         3,226 PDB IDs to download

  vina_configs/                 7,227 .conf files (grid boxes filled in step 4)

  compute_grid_boxes.py         Parses residues -> computes grid boxes
  00_setup_env.sh               Install conda + Vina-GPU (run ONCE)
  01_download_pdbs.sh           Download PDB structures (~8 GB)
  02_prepare_proteins.sh        Convert PDB -> PDBQT receptors
  03_prepare_ligands.sh         Convert SMILES -> PDBQT ligands
  04_compute_grids.sh           Compute binding site grid boxes
  05_run_docking.sh             Run all docking on GPU

  README.md                     This file
```

---

## Beginner's Guide: Step by Step

You are on a Windows machine, connected to an A100 Linux server via VSCode Remote-SSH.

### Step 0: Transfer this folder to the A100 server

In VSCode, you can simply **drag and drop** the `docking_package_trial` folder from your local file explorer into the VSCode file panel (which shows the remote server's files).

Or using the terminal (run this on your **local** Git Bash):
```bash
scp -r "D:/finale/negative_upgrade/docking_package_trial" \
    YOUR_USERNAME@YOUR_SERVER_IP:~/docking_trial
```

Then on the **A100 server** (in your VSCode SSH terminal):
```bash
cd ~/docking_trial   # or wherever you put the folder
ls                    # you should see all the files listed above
```

### Step 1: Install software (ONE TIME, ~10 minutes)

```bash
bash 00_setup_env.sh
```

This installs:
- **Miniconda** (Python package manager) with Python 3.11
- **OpenBabel** (converts SMILES chemical strings to 3D structures)
- **ADFR Suite** (prepares protein structures for docking)
- **Vina-GPU** (the actual docking software, compiled for A100 GPU)

After this finishes, close and reopen your terminal, then:
```bash
conda activate docking
```

### Step 2: Download PDB structures (~2-4 hours)

```bash
conda activate docking
bash 01_download_pdbs.sh
```

This downloads 3,226 protein structures from the Protein Data Bank (PDBe/RCSB).
Each PDB is ~2.5 MB, so total is ~8 GB.

**You can interrupt this anytime** (Ctrl+C) — when you re-run, it skips already-downloaded files.

### Step 3: Prepare protein structures (~30 minutes)

```bash
conda activate docking
bash 02_prepare_proteins.sh
```

This converts each raw PDB file into PDBQT format:
- Adds hydrogen atoms
- Merges non-polar hydrogens (saves computation)
- Assigns AutoDock atom types

### Step 4: Prepare ligand structures (~20 minutes)

```bash
conda activate docking
bash 03_prepare_ligands.sh
```

This converts each compound's SMILES string into a 3D PDBQT structure:
- Generates 3D coordinates from the SMILES
- Finds the lowest-energy conformer
- Assigns AutoDock atom types and rotatable bonds

### Step 5: Compute grid boxes (~5 minutes)

```bash
conda activate docking
bash 04_compute_grids.sh
```

This is the key step! For each task, it:
1. Reads the binding site residue list from the manifest
2. Finds those residues in the PDB structure (by chain + name + number)
3. Gets their CA (alpha carbon) coordinates
4. Computes the bounding box center + size with 8 Angstrom padding
5. Writes the grid box into the Vina config file

### Step 6: Run docking (~5-10 seconds per task on A100)

```bash
conda activate docking
bash 05_run_docking.sh
```

This runs all 7,227 docking calculations sequentially on GPU 0. Each task takes ~5-10 seconds, so total time is ~10-20 hours.

**Resume-safe**: If interrupted, re-run the script — it skips tasks that already have results.

---

## After Docking: Understanding Results

Results are in the `results/` folder. Each `DOCK-XXXXX_out.pdbqt` file contains up to 9 predicted binding poses with their Vina affinity scores (kcal/mol).

To extract the best (lowest/most negative) score for each task:
```bash
grep "REMARK VINA RESULT:" results/*_out.pdbqt | \
    awk '{print $4, FILENAME}' | \
    sort -n > best_scores.txt
head -20 best_scores.txt
```

More negative = better predicted binding. Typical ranges:
- **-12 to -9**: Excellent (likely true binder)
- **-9 to -7**: Good
- **-7 to -5**: Moderate
- **>-5**: Weak (likely non-binder)

---

## FAQ

**Q: Why prepare files locally but run docking on A100?**
A: The "preparation" we did locally was generating the task list, config files, and scripts from the MemPro database. The computationally heavy steps (PDB download, protein/ligand preparation, docking) all run on the A100 where the GPU is.

**Q: What if some PDB downloads fail?**
A: Some PDB IDs may be obsolete. The docking script skips tasks whose receptor PDBQT is missing. Check the download log and `logs/failed.txt` afterwards.

**Q: Can I run this on multiple GPUs?**
A: Yes! Split the config files across GPUs:
```bash
# GPU 0: first half
ls vina_configs/DOCK-0*.conf | while read conf; do
    vina_gpu --config "$conf" --gpu_id 0
done

# GPU 1: second half (in another terminal)
ls vina_configs/DOCK-5*.conf | while read conf; do
    vina_gpu --config "$conf" --gpu_id 1
done
```

**Q: How do I check progress?**
A: Count completed results:
```bash
ls results/*_out.pdbqt | wc -l
```
