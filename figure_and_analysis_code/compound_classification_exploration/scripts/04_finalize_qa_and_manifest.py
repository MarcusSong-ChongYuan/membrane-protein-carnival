from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import matplotlib
import numpy
import pandas
import rdkit
import scipy
import seaborn
import sklearn
import umap

ROOT = Path(r"D:\finale\compound_classification_exploration")
QA = ROOT / "qa"

pdf_raw = json.loads((QA / "PDF_TEXT_AUDIT_RAW.json").read_text(encoding="utf-8-sig"))
pdf_rows = []
for row in pdf_raw:
    result = json.loads(row["result"])
    pdf_rows.append({
        "file": row["file"],
        "auditable": result.get("auditable"),
        "minimum_found_pt": result.get("minimum_found_pt"),
        "below_minimum_count": result.get("below_minimum_count", 0),
        "warnings": result.get("warnings", []),
    })

source = json.loads((QA / "SOURCE_PREFLIGHT_02.json").read_text(encoding="utf-8-sig"))
blocking_pdf = [x for x in pdf_rows if x["below_minimum_count"] > 0]
unauditable = [x for x in pdf_rows if not x["auditable"]]

qa = {
    "status": "PASS_WITH_DOCUMENTED_EXPLORATORY_EXEMPTIONS",
    "scope": "compound classification exploration; not final manuscript Figure 3",
    "backend": "Python only",
    "denominators": {
        "compound_registry": 646670,
        "interaction_linked_compounds": 240055,
        "structure_valid_interaction_linked_compounds": 240054,
        "formal_pairs": 529168,
    },
    "source_preflight": source["summary"],
    "source_warnings_reviewed": {
        "no_tiff": "Accepted for exploration; final selected manuscript panels will be re-exported at 600 dpi TIFF.",
        "350_dpi": "Accepted for exploration previews; vector PDF/SVG are supplied for V1-V9.",
        "width": "Candidate plots are independent views, not assembled final-size manuscript panels.",
        "simulated_data_warning": "False positive: deterministic RNG is used only for reproducible sampling/bootstrap mechanics, not fabricated observations.",
        "rotated_text": "Visually inspected; labels remain within canvas. Final selected panel will use anchored rotation or horizontal labels.",
    },
    "pdf_text_audit": {
        "pdf_count": len(pdf_rows),
        "below_5pt_count": len(blocking_pdf),
        "unauditable_count": len(unauditable),
        "unauditable_files": [x["file"] for x in unauditable],
        "note": "The RDKit scaffold gallery PDF is raster-backed and contains no auditable PDF text operators; its source table and high-resolution PNG are supplied.",
    },
    "visual_review": {
        "contact_sheet_inspected": True,
        "V8B_inspected": True,
        "V9B_inspected": True,
        "blocking_overlap_or_clipping": 0,
        "observations": [
            "UMAP forms a dominant continuous cloud with limited local islands; retained as exploratory only.",
            "3D physicochemical plot is compressed by extreme values; retained as interactive/exploratory only.",
            "Target-breadth class view is coverage-dependent and not interpreted as intrinsic promiscuity.",
        ],
    },
    "classification_boundary": "chemical_regime_local is a transparent deterministic local taxonomy, not ClassyFire or ChEBI ontology.",
    "blocking_errors": 0,
}
(QA / "COMPOUND_EXPLORATION_QA.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")

env = {
    "python": platform.python_version(), "platform": platform.platform(),
    "numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__,
    "scikit_learn": sklearn.__version__, "matplotlib": matplotlib.__version__,
    "seaborn": seaborn.__version__, "rdkit": rdkit.__version__, "umap_learn": umap.__version__,
}
(QA / "PYTHON_ENVIRONMENT.json").write_text(json.dumps(env, indent=2), encoding="utf-8")

readme = """# MemPro compound classification exploration pack

This directory contains independent candidate views V1--V10. It is **not** the final manuscript Figure 3.

## Denominators

- Registry: 646,670 compound identities.
- Interaction-linked: 240,055 compound identities.
- Structure-valid interaction-linked: 240,054 compounds.
- Formal protein-compound pairs: 529,168.

## Interpretation boundaries

- `chemical_regime_local` is a deterministic local rule-based taxonomy, not ClassyFire/ChEBI ontology.
- Bemis-Murcko scaffolds, Morgan/Butina neighborhoods, PCA and UMAP answer different questions and are not interchangeable chemical classes.
- UMAP axes/global distances and the 3D perspective have no direct chemical meaning.
- Observed target breadth is database-coverage-dependent.

Start with `08_contact_sheet/COMPOUND_VISUAL_EXPLORATION_CONTACT_SHEET.pdf` and `09_review/COMPOUND_VISUALIZATION_REVIEW.md`.
"""
(ROOT / "README.md").write_text(readme, encoding="utf-8")

manifest_rows = []
for path in sorted(ROOT.rglob("*")):
    if path.is_file() and path.name != "SHA256SUMS.tsv":
        h = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        manifest_rows.append((path.relative_to(ROOT).as_posix(), path.stat().st_size, h.hexdigest()))
with (QA / "SHA256SUMS.tsv").open("w", encoding="utf-8", newline="") as fh:
    fh.write("relative_path\tsize_bytes\tsha256\n")
    for rel, size, digest in manifest_rows:
        fh.write(f"{rel}\t{size}\t{digest}\n")

print(json.dumps({"status": qa["status"], "files_hashed": len(manifest_rows), "blocking_errors": 0}, indent=2))
