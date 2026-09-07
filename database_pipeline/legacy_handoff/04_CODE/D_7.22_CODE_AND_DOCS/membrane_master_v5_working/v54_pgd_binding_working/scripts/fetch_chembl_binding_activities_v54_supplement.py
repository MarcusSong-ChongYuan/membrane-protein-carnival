#!/usr/bin/env python3
"""Run the ChEMBL Ki/Kd downloader for v5.3 accessions absent from v4."""

from __future__ import annotations

import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve().parent
module_path = HERE / "fetch_chembl_binding_activities_v54.py"
spec = importlib.util.spec_from_file_location("chembl_fetch_v54", module_path)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

module.TARGET_INDEX = (
    module.WORK
    / "intermediate"
    / "external_chembl_target_index_v54_supplement.tsv"
)
module.RAW = module.WORK / "raw" / "chembl_37_direct_binding_v54_supplement_batches"
module.OUT = (
    module.WORK
    / "intermediate"
    / "chembl_37_binding_activities_v54_supplement.tsv"
)
module.REPORT = (
    module.WORK
    / "reports"
    / "CHEMBL_BINDING_ACTIVITY_V54_SUPPLEMENT_REPORT.json"
)
module.MAX_WORKERS = 8

module.main()
