from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from figure_style_v62 import CAPTION_DIR, DATA_DIR, PDF_DIR, PNG_DIR, QA_DIR, ROOT, SVG_DIR


EXPECTED = [
    "M1_database_architecture_quality",
    "M2_membrane_proteome_landscape",
    "M3_expression_localization_atlas",
    "M4_disease_association_landscape",
    "M5_chemical_identity_space",
    "M6_binding_evidence_sites",
    "M7_negative_conflicts_quality",
    "M8_docking_prioritization",
    "S1_source_integration",
    "S2_topology_classification_assembly",
    "S3_expression_localization_detail",
    "S4_disease_ontology_evidence",
    "S5_compound_identity_qc",
    "S6_docking_robustness",
]


def make_contact_sheet(stems: list[str], output: Path, columns: int = 2, width: int = 2400) -> None:
    thumbs = []
    cell_width = width // columns
    for stem in stems:
        image = Image.open(PNG_DIR / f"{stem}.png").convert("RGB")
        scale = (cell_width - 24) / image.width
        image = image.resize((int(image.width * scale), int(image.height * scale)), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (cell_width, image.height + 46), "white")
        canvas.paste(image, ((cell_width - image.width) // 2, 40))
        draw = ImageDraw.Draw(canvas)
        draw.text((12, 12), stem, fill="#17365D", font=ImageFont.load_default())
        thumbs.append(canvas)
    rows = (len(thumbs) + columns - 1) // columns
    row_height = max(item.height for item in thumbs)
    sheet = Image.new("RGB", (width, rows * row_height), "white")
    for index, item in enumerate(thumbs):
        x = (index % columns) * cell_width
        y = (index // columns) * row_height
        sheet.paste(item, (x, y))
    sheet.save(output, dpi=(150, 150), quality=94)


def main() -> None:
    records = []
    failures = []
    for stem in EXPECTED:
        row = {"figure": stem}
        for extension, directory in [("svg", SVG_DIR), ("pdf", PDF_DIR), ("png", PNG_DIR)]:
            path = directory / f"{stem}.{extension}"
            row[f"{extension}_exists"] = path.exists()
            row[f"{extension}_bytes"] = path.stat().st_size if path.exists() else 0
            if not path.exists() or path.stat().st_size == 0:
                failures.append(f"Missing/empty {path}")
        png_path = PNG_DIR / f"{stem}.png"
        if png_path.exists():
            with Image.open(png_path) as image:
                row["png_width_px"] = image.width
                row["png_height_px"] = image.height
                row["png_dpi_x"] = round(float(image.info.get("dpi", (0, 0))[0]), 1)
                row["png_dpi_y"] = round(float(image.info.get("dpi", (0, 0))[1]), 1)
                if image.width < 7000 or image.height < 4000:
                    failures.append(f"Low pixel dimensions {png_path}: {image.size}")
        svg_path = SVG_DIR / f"{stem}.svg"
        if svg_path.exists():
            svg_text = svg_path.read_text(encoding="utf-8", errors="ignore")
            row["svg_editable_text"] = "<text" in svg_text
            if not row["svg_editable_text"]:
                failures.append(f"No editable SVG text in {svg_path}")
        records.append(row)

    captions = [
        CAPTION_DIR / "CAPTIONS_M1_M4.md",
        CAPTION_DIR / "CAPTIONS_M5_M8.md",
        CAPTION_DIR / "CAPTIONS_SUPPLEMENTARY_FIGURES.md",
    ]
    for path in captions:
        if not path.exists() or path.stat().st_size == 0:
            failures.append(f"Missing caption file {path}")

    combined = "# MemPro V6.2 figure captions\n\n"
    for path in captions:
        if path.exists():
            combined += path.read_text(encoding="utf-8").strip() + "\n\n"
    (CAPTION_DIR / "CAPTIONS_ALL_FIGURES.md").write_text(combined, encoding="utf-8")

    make_contact_sheet(EXPECTED[:8], QA_DIR / "MAIN_FIGURES_CONTACT_SHEET.png")
    make_contact_sheet(EXPECTED[8:], QA_DIR / "SUPPLEMENTARY_FIGURES_CONTACT_SHEET.png")
    make_contact_sheet(EXPECTED, QA_DIR / "ALL_FIGURES_CONTACT_SHEET.png")

    report = {
        "status": "PASS" if not failures else "FAIL",
        "expected_figure_count": len(EXPECTED),
        "formats": ["editable SVG", "vector PDF", "600 dpi PNG"],
        "figure_records": records,
        "data_table_count": len(list(DATA_DIR.glob("*.tsv"))),
        "caption_files": [str(path) for path in captions],
        "failures": failures,
    }
    (QA_DIR / "FIGURE_PACKAGE_VALIDATION.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    readme = f"""# MemPro V6.2 scientific figure package

This package contains {len(EXPECTED)} concentrated 2×3 composite figures: eight main figures (M1–M8) and six supplementary figures (S1–S6).

## Formats

- `svg/`: editable vector artwork with live text.
- `pdf/`: publication-ready vector output.
- `png/`: 600 dpi raster output for review and submission systems.
- `data/`: panel-level plotted statistics and deterministic embedding coordinates.
- `captions/`: detailed captions, including interpretation limits.
- `qa/`: input profiles, validation reports and contact sheets.

## Visual semantics

Membrane classes A/B/C use blue/green/amber; evidence tiers E1/BE1, E2/BE2 and E3/BE3 use dark blue, green and orange; positive, negative, conflict and review states use green, grey, vermilion and magenta. Docking tiers R0/D1/D2/D3 use purple, green, ochre and slate blue. The palette is color-vision-aware and remains interpretable through labels in grayscale.

## Reproducibility

All figures were generated from the frozen V6.2 release. UMAP and rendered compound-space samples use random seed 42. Chemical-space coordinates are saved in `data/`. The full categorical totals use complete source tables; large scatter/hexbin panels use deterministic sampling solely for rendering.

Package QA status: **{report["status"]}**.
"""
    (ROOT / "README_FIGURES.md").write_text(readme, encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
