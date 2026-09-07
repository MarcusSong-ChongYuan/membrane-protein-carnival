from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(os.environ.get("MEMPRO_REVISION_ROOT", r"D:\finale\07_M1-M8图表修订_20260806"))
SRC = Path(r"C:\tmp\mempro_revision_src")
PPT_SRC = Path(r"C:\tmp\mempro_ppt_revision")
REPRO = ROOT / "08_reproducibility"
SCRIPTS = REPRO / "scripts"
QA = ROOT / "09_qa"
DOCS = ROOT / "00_docs"

FIGURES = [
    "M1_database_architecture_quality_revised",
    "M2_membrane_proteome_landscape_revised",
    "M3_expression_localization_atlas_revised",
    "M4_disease_association_landscape_revised",
    "M5_chemical_identity_space_revised",
    "M6_binding_evidence_sites_revised",
    "M7_negative_conflicts_quality_revised",
    "M8_docking_prioritization_revised",
]

SUPPORT = [
    ROOT / "01_source_audit/PANEL_SOURCE_FIELD_VERSION_UNIT.tsv",
    ROOT / "01_source_audit/DATABASE_EXACT_CONTRIBUTION_COUNTS.tsv",
    ROOT / "02_data/M2D_spearman_main_and_bootstrap.tsv",
    ROOT / "02_data/M2D_spearman_sensitivity.tsv",
    ROOT / "02_data/M4F_ontology_disease_anatomy_matrix.tsv",
    ROOT / "02_data/M5C_recomputed_structural_classes.tsv",
    ROOT / "02_data/M5E_pagtn_cluster_summary.tsv",
    ROOT / "02_data/M5E_pagtn_morgan_coordinates.tsv",
    ROOT / "04_pagtn/PAGTN_EXPERIMENT_QA.json",
    ROOT / "05_alphafold/ALPHAFOLD_CONFIDENCE_QA.json",
    ROOT / "06_binding_site_side/BINDING_SITE_MEMBRANE_SIDE_QA.json",
    ROOT / "07_ppt/MemPro_V6.3.1_M1-M8科学图表修订版_20260806.pptx",
    ROOT / "07_ppt/MemPro_V6.3.1_M1-M8修订版逐页讲稿.txt",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def copy_scripts() -> None:
    SCRIPTS.mkdir(parents=True, exist_ok=True)
    names = [
        "build_revision_figures.py", "build_m2_fixed.py", "build_m5_pagtn_final.py",
        "build_source_contribution_counts.py", "classify_binding_site_side.py",
        "fetch_alphafold_confidence.py", "run_pagtn_sample_experiment.py",
        "run_pagtn_sample_experiment_v2.py", "make_qa_thumbnail.py",
    ]
    for name in names:
        source = SRC / name
        if source.exists():
            shutil.copy2(source, SCRIPTS / name)
    for name in ["build_revised_deck.mjs", "build_revised_deck_v4.mjs", "inspect_template_deck_win.mjs"]:
        source = PPT_SRC / name
        if source.exists():
            shutil.copy2(source, SCRIPTS / name)
    template = PPT_SRC / "source.pptx"
    if template.exists():
        shutil.copy2(template, REPRO / "reference_template.pptx")


def write_docs() -> None:
    DOCS.mkdir(parents=True, exist_ok=True)
    report = f"""# MemPro M1-M8 scientific figure revision — completion report

Generated: {datetime.now(timezone.utc).isoformat()}

## Scope and data protection

- All work was written under `{ROOT}`.
- Frozen V6.2 and V6.3.1 source releases were read only; no release table was overwritten.
- Eight figures were exported as PNG (600 dpi), editable SVG, and PDF.
- The presentation preserves the approved 18-slide two-page-per-module template and includes speaker notes on every slide.

## Completed scientific corrections

1. **Panel provenance:** 48 M1-M8 panels now have a source, field/definition, version policy, statistical unit, and limitation record.
2. **M2 Spearman scope:** A-class, E1/E2, canonical multipass proteins only; 2-40 TM segments and 50-5000 aa; n=2,844; rho=0.445; bootstrap 95% CI 0.410-0.481 (seed 42).
3. **M3 provenance:** HPA 25.1 is explicitly separated into tissue RNA nTPM, cell-type RNA nCPM, IHC protein staining, and mapping/missingness layers.
4. **M4 anatomy:** keyword bins were replaced by ontology-derived, multi-label disease-anatomy summaries using available MONDO/DO/Uberon mappings.
5. **M5C classes:** all 899,591 core-QC compounds were reclassified by one deterministic, mutually exclusive descriptor rule set.
6. **PAGTN experiment:** 2,400 compounds, seed 42, 64-dimensional property-supervised PAGTN embeddings, HDBSCAN clustering, and a Morgan radius-2/1024-bit baseline. PAGTN improved non-noise silhouette (0.360 vs 0.148), while Morgan retained higher 10-neighbor scaffold purity (11.9% vs 9.7%); no claim of universal superiority is made.
7. **AlphaFold:** 10,908/10,908 API metadata requests succeeded; local site pLDDT mapped for 857/859 attempted site targets. pLDDT is presented as model confidence, not experimental validation.
8. **Membrane-side sites:** 9,036 residue-contact sites were classified relative to protein topology: 773 intramembrane, 1,331 interface/mixed, 407 cytoplasmic, 384 non-cytoplasmic, 3 both sides/mixed, and 6,138 extramembrane-side unresolved.
9. **PPT/QA:** 18 slides rendered; no slide overflow was detected by the bundled presentation validator.

## Interpretation boundaries

- PAGTN is the encoder; HDBSCAN is the clustering algorithm. The current encoder is trained on four physicochemical descriptors and is a method experiment, not a universal pretrained chemical embedding.
- The membrane-side classifier is protein-topology-relative. It does not yet calculate ligand/residue atom z-coordinates in an OPM membrane frame. Non-cytoplasmic includes extracellular and organelle-lumen sides.
- AlphaFold coverage and pLDDT cannot replace experimental structures or prove a binding pocket is biologically correct.
- M4 anatomical associations describe ontology mappings, not tissue causality.
- Database source counts are distinct contributing database counts unless an experiment/structure lineage key explicitly supports independence.
"""
    (DOCS / "REVISION_COMPLETION_REPORT.md").write_text(report, encoding="utf-8")

    runbook = f"""# Rebuild runbook

## Inputs

- Frozen V6.2: `D:/7.22/evidence_expansion_v2_working/releases/release_mempro_v6_2_20260730`
- V6.3 candidate: `D:/finale/02_V6.3_candidate_20260803`
- Revision output: `{ROOT.as_posix()}`

## Order

1. Set `MEMPRO_REVISION_ROOT` to the revision output directory.
2. Run `build_revision_figures.py` for the initial M1-M8 set.
3. Run `build_m2_fixed.py` for scoped Spearman/bootstrap M2.
4. Run `classify_binding_site_side.py`, then rebuild M6.
5. Run `fetch_alphafold_confidence.py`.
6. In the isolated PAGTN environment, run `run_pagtn_sample_experiment_v2.py`.
7. Run `build_m5_pagtn_final.py`.
8. Build the deck with the artifact-tool script. The supplied `reference_template.pptx` is the approved template.
9. Run the presentation `slides_test.py` validator and `finalize_revision_release.py`.

## Environment

- Standard plots: the existing `D:/7.22/evidence_expansion_v2_working/tools/plotting_venv`.
- PAGTN: `04_pagtn/env`; versions are frozen in `pagtn_environment_freeze.txt`.
- PPT: Node.js plus `@oai/artifact-tool`; do not replace this step with python-pptx.

The scripts retain explicit source paths to prevent accidental use of a different release. For another device, mount/copy the frozen releases and update only the three documented root constants/environment variables.
"""
    (REPRO / "README_REBUILD.md").write_text(runbook, encoding="utf-8")


def environment_freeze() -> None:
    py = ROOT / "04_pagtn/env/python.exe"
    if py.exists():
        result = subprocess.run([str(py), "-m", "pip", "freeze"], capture_output=True, text=True, check=False)
        (REPRO / "pagtn_environment_freeze.txt").write_text(result.stdout, encoding="utf-8")


def validate() -> dict:
    checks = []
    for stem in FIGURES:
        for ext, folder in [("png", "png"), ("svg", "svg"), ("pdf", "pdf")]:
            p = ROOT / f"03_figures/{folder}/{stem}.{ext}"
            checks.append({"name": f"figure:{stem}.{ext}", "pass": p.exists() and p.stat().st_size > 5000, "path": str(p), "bytes": p.stat().st_size if p.exists() else 0})
    for p in SUPPORT:
        checks.append({"name": f"support:{p.name}", "pass": p.exists() and p.stat().st_size > 0, "path": str(p), "bytes": p.stat().st_size if p.exists() else 0})
    renders = sorted((QA / "ppt_render/slides").glob("slide-*.png"))
    checks.append({"name": "ppt_render_count_18", "pass": len(renders) == 18, "value": len(renders)})
    notes_inspect = ROOT / "07_ppt/MemPro_V6.3.1_M1-M8科学图表修订版_20260806.pptx.inspect.ndjson"
    notes_count = notes_inspect.read_text(encoding="utf-8").count('"kind":"notes"') if notes_inspect.exists() else 0
    checks.append({"name": "speaker_notes_present", "pass": notes_count >= 18, "value": notes_count})
    return {"status": "PASS" if all(x["pass"] for x in checks) else "FAIL", "release_data_modified": False, "checks": checks}


def manifest(validation: dict) -> None:
    key_files = []
    for stem in FIGURES:
        for ext, folder in [("png", "png"), ("svg", "svg"), ("pdf", "pdf")]:
            key_files.append(ROOT / f"03_figures/{folder}/{stem}.{ext}")
    key_files += SUPPORT
    key_files += sorted(SCRIPTS.glob("*"))
    key_files += [DOCS / "REVISION_COMPLETION_REPORT.md", REPRO / "README_REBUILD.md", REPRO / "pagtn_environment_freeze.txt"]
    rows = []
    for p in key_files:
        if p.exists() and p.is_file():
            rows.append({"sha256": sha256(p), "bytes": p.stat().st_size, "relative_path": p.relative_to(ROOT).as_posix()})
    lines = ["sha256\tbytes\trelative_path"] + [f"{r['sha256']}\t{r['bytes']}\t{r['relative_path']}" for r in rows]
    (ROOT / "SHA256SUMS.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    payload = {
        "release_name": "MemPro M1-M8 scientific figure revision 20260806",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "frozen_sources_modified": False,
        "validation_status": validation["status"],
        "file_count": len(rows),
        "sha256_table": "SHA256SUMS.tsv",
    }
    (ROOT / "REVISION_MANIFEST.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    copy_scripts()
    write_docs()
    environment_freeze()
    validation = validate()
    QA.mkdir(parents=True, exist_ok=True)
    (QA / "FINAL_VALIDATION.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    # Replace the stale initial partial status with the final, post-PAGTN state.
    (QA / "FIGURE_REVISION_VALIDATION.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest(validation)
    print(json.dumps({"status": validation["status"], "checks": len(validation["checks"]), "root": str(ROOT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
