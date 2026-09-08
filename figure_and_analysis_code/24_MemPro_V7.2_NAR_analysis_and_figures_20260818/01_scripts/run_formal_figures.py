"""Run the archival V7.2 NAR figure workflow against a portable FORMAL release.

This runner deliberately writes to a new output directory.  It never edits the
frozen release or historical figure package in place.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPTS = (
    "00_build_analysis_layer.py",
    "01_compute_statistics.py",
    "01b_compute_evidence_flow.py",
    "02_compute_chemical_space.py",
    "03_draw_main_figures.py",
    "04_draw_graphics.py",
    "08_run_final_qa.py",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path, help="Extracted 01_database_FORMAL directory")
    parser.add_argument("--output-root", required=True, type=Path, help="New directory for derived analysis and figures")
    parser.add_argument("--asset-root", type=Path, help="Directory containing human_anatomy_base_v1.png from the migration package")
    parser.add_argument("--overwrite", action="store_true", help="Allow writing into an existing empty/derived output directory")
    args = parser.parse_args()

    required = args.data_root / "01_core_v72" / "01_release_tables" / "protein_master_v72.tsv.gz"
    if not required.exists():
        raise SystemExit(f"FORMAL release table not found: {required}")
    if args.output_root.exists() and any(args.output_root.iterdir()) and not args.overwrite:
        raise SystemExit("Output directory already contains files. Choose a new directory or pass --overwrite.")
    args.output_root.mkdir(parents=True, exist_ok=True)

    if args.asset_root:
        source = args.asset_root / "human_anatomy_base_v1.png"
        destination = args.output_root / "05_graphics" / source.name
        if not source.exists():
            raise SystemExit(f"Required anatomy asset not found: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    environment = os.environ.copy()
    environment["MEMPRO_DATA_ROOT"] = str(args.data_root.resolve())
    environment["MEMPRO_FIGURE_OUTPUT_ROOT"] = str(args.output_root.resolve())
    script_root = Path(__file__).resolve().parent
    for script in SCRIPTS:
        subprocess.run([sys.executable, str(script_root / script)], check=True, cwd=script_root, env=environment)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

