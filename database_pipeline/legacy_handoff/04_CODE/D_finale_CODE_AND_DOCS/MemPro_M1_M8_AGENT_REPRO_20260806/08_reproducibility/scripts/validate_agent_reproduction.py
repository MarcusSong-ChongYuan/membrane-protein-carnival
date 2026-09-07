from __future__ import annotations

import json
from pathlib import Path


STEMS = [
    "M1_database_architecture_quality_revised",
    "M2_membrane_proteome_landscape_revised",
    "M3_expression_localization_atlas_revised",
    "M4_disease_association_landscape_revised",
    "M5_chemical_identity_space_revised",
    "M6_binding_evidence_sites_revised",
    "M7_negative_conflicts_quality_revised",
    "M8_docking_prioritization_revised",
]


def validate(root: Path) -> None:
    checks = []
    for stem in STEMS:
        for ext, folder in [("png", "png"), ("svg", "svg"), ("pdf", "pdf")]:
            path = root / f"03_figures/{folder}/{stem}.{ext}"
            checks.append({"file": path.relative_to(root).as_posix(), "pass": path.exists() and path.stat().st_size > 5000})
    for rel in [
        "02_data/M2D_spearman_main_and_bootstrap.tsv",
        "02_data/M5E_pagtn_cluster_summary.tsv",
        "05_alphafold/ALPHAFOLD_CONFIDENCE_QA.json",
        "06_binding_site_side/BINDING_SITE_MEMBRANE_SIDE_QA.json",
        "07_ppt/MemPro_V6.3.1_M1-M8科学图表修订版_20260806.pptx",
    ]:
        path = root / rel
        checks.append({"file": rel, "pass": path.exists() and path.stat().st_size > 0})
    status = "PASS" if all(row["pass"] for row in checks) else "FAIL"
    output = {"status": status, "checks": checks}
    out_path = root / "09_qa/AGENT_REPRODUCTION_VALIDATION.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": status, "checks": len(checks)}))
    if status != "PASS":
        raise RuntimeError("Agent reproduction validation failed")


if __name__ == "__main__":
    validate(Path(__file__).resolve().parents[2])
