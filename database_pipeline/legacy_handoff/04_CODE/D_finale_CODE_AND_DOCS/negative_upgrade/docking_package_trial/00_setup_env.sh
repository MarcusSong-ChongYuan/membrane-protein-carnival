#!/bin/bash
# ============================================================
# Environment Setup for MemPro Docking
# Run ONCE on the A100 server before anything else.
# ============================================================
set -e

echo "=== Step 1: Create conda environment ==="
conda create -n docking -c conda-forge \
    python=3.11 openbabel adfr-suite curl -y

echo ""
echo "=== Step 2: Build Vina-GPU ==="
cd $HOME
if [ ! -d "Vina-GPU" ]; then
    git clone https://github.com/DeltaGroupNJUPT/Vina-GPU.git
fi
cd Vina-GPU
make clean 2>/dev/null || true
make
echo 'export PATH=$HOME/Vina-GPU:$PATH' >> ~/.bashrc

echo ""
echo "=== Done! ==="
echo "Now run: bash 01_download_pdbs.sh"
