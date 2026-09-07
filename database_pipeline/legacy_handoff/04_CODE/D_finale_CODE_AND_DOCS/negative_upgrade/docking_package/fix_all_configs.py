#!/usr/bin/env python3
"""
Fix all Vina config files for A100 deployment:
1. Remove exhaustiveness line (Vina-GPU doesn't support it)
2. Remove log line (Vina-GPU writes to --log arg, not config)
3. Fix relative paths: ../receptors_pdbqt/ -> receptors_pdbqt/
                         ../ligands_pdbqt/ -> ligands_pdbqt/
                         ../docking_results/ -> results/
4. Convert CRLF -> LF
"""
import os
import glob

CONFIG_DIR = r"D:\finale\negative_upgrade\docking_package\vina_configs"

configs = sorted(glob.glob(os.path.join(CONFIG_DIR, "DOCK-*.conf")))
print(f"Found {len(configs)} config files")

fixed = 0
for i, path in enumerate(configs):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.read()

    original = lines

    # 1. Remove exhaustiveness line
    lines = '\n'.join(
        l for l in lines.split('\n')
        if not l.strip().startswith('exhaustiveness')
    )

    # 2. Remove log = line
    lines = '\n'.join(
        l for l in lines.split('\n')
        if not l.strip().startswith('log =')
    )

    # 3. Fix paths
    lines = lines.replace('../receptors_pdbqt/', 'receptors_pdbqt/')
    lines = lines.replace('../ligands_pdbqt/', 'ligands_pdbqt/')
    lines = lines.replace('../docking_results/', 'results/')

    # 4. Strip trailing whitespace
    lines = '\n'.join(l.rstrip() for l in lines.split('\n'))

    if '\r' in lines:
        lines = lines.replace('\r', '')
        fixed += 1

    if i == 0:
        print(f"\n--- BEFORE (first file) ---")
        print(original[:300])
        print(f"\n--- AFTER (first file) ---")
        print(lines[:300])
        print()

    if lines != original:
        with open(path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(lines)
        fixed += 1

    if (i + 1) % 50000 == 0:
        print(f"  Processed {i+1}/{len(configs)}... fixed={fixed}")

print(f"\nDone: {len(configs)} files processed, {fixed} files modified")
