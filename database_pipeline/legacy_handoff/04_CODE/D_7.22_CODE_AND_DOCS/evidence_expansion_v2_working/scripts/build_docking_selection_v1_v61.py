#!/usr/bin/env python3
"""Refresh docking selection V1 from the formally frozen MemPro V6.1 release."""

from __future__ import annotations

import json
from pathlib import Path

import build_docking_selection_v1 as base


base.RELEASE_DIR = Path(
    r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_1_20260729"
)
base.RUN_DIR = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\docking_selection_v1_v61_20260729"
)
base.STAGING_DIR = base.RUN_DIR / "staging"
base.QA_DIR = base.RUN_DIR / "qa"
base.PAIR_FILE = base.RELEASE_DIR / "protein_compound_summary_v6_1.tsv.gz"
base.SITE_FILE = base.RELEASE_DIR / "binding_site_instances_v6_1.tsv.gz"
base.PROTEIN_FILE = (
    base.RELEASE_DIR / "human_membrane_protein_master_v6_1.tsv"
)
base.COMPOUND_FILE = base.RELEASE_DIR / "small_molecule_master_v1_2.tsv"
base.ALL_OUTPUT = (
    base.STAGING_DIR / "docking_pair_priority_v1_v61.tsv.gz"
)
base.RECOMMENDED_OUTPUT = (
    base.STAGING_DIR / "docking_recommended_pairs_v1_v61.tsv.gz"
)
base.TARGET_OUTPUT = (
    base.STAGING_DIR / "docking_target_summary_v1_v61.tsv"
)
base.VALIDATION_OUTPUT = (
    base.QA_DIR / "DOCKING_SELECTION_V1_V61_VALIDATION.json"
)
base.DICTIONARY_OUTPUT = (
    base.RUN_DIR / "DOCKING_SELECTION_V1_V61_DATA_DICTIONARY.md"
)
base.EXPECTED_TIER_COUNTS = None
base.RELEASE_BASIS = "the formally frozen MemPro V6.1 release"
base.SELECTION_BASIS_TAG = "MemPro_V6.1_formal_frozen_release"
base.DICTIONARY_TITLE = "Docking selection V1 (V6.1 formal refresh)"


if __name__ == "__main__":
    report = base.build_outputs()
    base.write_dictionary(report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
