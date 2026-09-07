from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
from pathlib import Path


def configure_module(module, root: Path) -> None:
    if hasattr(module, "ROOT"):
        module.ROOT = root
    mappings = {
        "V62": root / "inputs/v62",
        "V63": root / "inputs/v63",
        "OLD": root / "inputs/old_figures",
        "MM": root / "inputs/membrane_master",
        "SOURCE": root / "inputs/v62/small_molecule_master_v1_3.tsv",
        "AUDIT": root / "01_source_audit",
        "DATA": root / "02_data",
        "FIG": root / "03_figures",
        "QA": root / "09_qa",
        "OUT": root / "04_pagtn",
        "DATA_OUT": root / "02_data",
    }
    for name, value in mappings.items():
        if hasattr(module, name):
            setattr(module, name, value)


def ensure_docking_fallback(root: Path) -> None:
    source = root / "inputs/docking/docking_pilot_manifest_v631.tsv"
    target = root / "inputs/v63/07_docking/docking_priority_v6_3_candidate.tsv.gz"
    if target.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as src, gzip.open(target, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst)


def run_figures(root: Path) -> None:
    import build_revision_figures as figures
    configure_module(figures, root)
    for folder in [root / "01_source_audit", root / "02_data", root / "03_figures/png", root / "03_figures/svg", root / "03_figures/pdf", root / "09_qa"]:
        folder.mkdir(parents=True, exist_ok=True)
    figures.AUDIT = root / "01_source_audit"
    figures.DATA = root / "02_data"
    figures.FIG = root / "03_figures"
    figures.QA = root / "09_qa"
    figures.main()

    import build_m2_fixed as m2
    configure_module(m2, root)
    m2.ROOT = root
    m2.OUT = root / "03_figures"
    m2.main()

    import build_source_contribution_counts as source_counts
    configure_module(source_counts, root)
    source_counts.main()


def run_full_augmentations(root: Path) -> None:
    import classify_binding_site_side as side
    configure_module(side, root)
    side.OUT = root / "06_binding_site_side"
    side.main()

    import fetch_alphafold_confidence as alphafold
    configure_module(alphafold, root)
    alphafold.OUT = root / "05_alphafold"
    alphafold.main()

    import run_pagtn_sample_experiment_v2 as pagtn_wrapper
    exp = pagtn_wrapper.experiment
    configure_module(exp, root)
    exp.OUT = root / "04_pagtn"
    exp.DATA_OUT = root / "02_data"
    exp.main()


def run_final_panels(root: Path) -> None:
    import build_m5_pagtn_final as m5
    configure_module(m5, root)
    m5.main()

    import build_revision_figures as figures
    configure_module(figures, root)
    figures.DATA = root / "02_data"
    figures.FIG = root / "03_figures"
    figures.make_m6()


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce the MemPro M1-M8 figure revision from bundled inputs.")
    parser.add_argument("--full", action="store_true", help="Rerun membrane-side, AlphaFold and PAGTN augmentations; AlphaFold needs internet.")
    args = parser.parse_args()
    root = Path(os.environ.get("MEMPRO_BUNDLE_ROOT", Path(__file__).resolve().parents[2])).resolve()
    scripts = root / "08_reproducibility/scripts"
    sys.path.insert(0, str(scripts))
    os.environ["MEMPRO_BUNDLE_ROOT"] = str(root)
    os.environ["MEMPRO_REVISION_ROOT"] = str(root)
    ensure_docking_fallback(root)
    run_figures(root)
    if args.full:
        run_full_augmentations(root)
    run_final_panels(root)
    import validate_agent_reproduction as validation
    validation.validate(root)


if __name__ == "__main__":
    main()
