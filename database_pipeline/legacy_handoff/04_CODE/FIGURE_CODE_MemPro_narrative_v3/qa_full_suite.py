from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pdfplumber
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(r"C:\Users\Administrator\MemPro_narrative_v3")
FIG = ROOT / "figures"
QA = ROOT / "qa"
QA.mkdir(parents=True, exist_ok=True)

STEMS = {
    "M1": "M1_source_overlap_and_marginal_gain",
    "M2": "M2_multi_axis_annotation_space_cleanfinal",
    "M3": "M3_expression_space_and_RNA_IHC",
    "M4": "M4_disease_expression_concordance_cleanfinal",
    "M5": "M5_chemical_space_embedding_comparison_cleanfinal",
    "M6": "M6_polypharmacology_and_membrane_side_cleanfinal",
    "M7": "M7_evidence_profiles_and_lineage_cleanfinal",
    "M8": "M8_value_readiness_and_knowledge_gaps",
}


def overlap_area(a: dict, b: dict) -> float:
    x = max(0.0, min(a["x1"], b["x1"]) - max(a["x0"], b["x0"]))
    y = max(0.0, min(a["bottom"], b["bottom"]) - max(a["top"], b["top"]))
    return x * y


def pdf_checks(path: Path) -> dict:
    with pdfplumber.open(path) as pdf:
        assert len(pdf.pages) == 1
        page = pdf.pages[0]
        words = page.extract_words(use_text_flow=True, keep_blank_chars=False)
        outside = [w for w in words if w["x0"] < -0.5 or w["top"] < -0.5 or w["x1"] > page.width + 0.5 or w["bottom"] > page.height + 0.5]
        overlaps = []
        for i, a in enumerate(words):
            for b in words[i + 1:]:
                # Text on different baselines is allowed to be close, but not to occupy the same area.
                area = overlap_area(a, b)
                if area > 4.0:
                    overlaps.append((a["text"], b["text"], round(area, 2)))
        text = " ".join(w["text"] for w in words)
        bad = [c for c in ["�", "□", "锟", "顖"] if c in text]
        return {
            "pdf_pages": 1,
            "pdf_width_pt": round(page.width, 2),
            "pdf_height_pt": round(page.height, 2),
            "outside_word_count": len(outside),
            "overlap_count": len(overlaps),
            "overlap_examples": overlaps[:8],
            "bad_encoding_tokens": bad,
        }


def png_checks(path: Path) -> dict:
    image = Image.open(path).convert("RGB")
    arr = np.asarray(image)
    gray = arr.mean(axis=2)
    nonwhite = float(np.mean(gray < 248))
    dark = float(np.mean(gray < 90))
    return {
        "png_width_px": image.width,
        "png_height_px": image.height,
        "aspect_ratio": round(image.width / image.height, 3),
        "nonwhite_fraction": round(nonwhite, 4),
        "dark_fraction": round(dark, 4),
    }


def contact_sheet() -> Path:
    thumb_w, thumb_h = 920, 650
    margin, label_h = 30, 45
    canvas = Image.new("RGB", (2 * thumb_w + 3 * margin, 4 * (thumb_h + label_h) + 5 * margin), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("arialbd.ttf", 27)
    except OSError:
        font = ImageFont.load_default()
    for idx, (code, stem) in enumerate(STEMS.items()):
        source = Image.open(FIG / f"{stem}.png").convert("RGB")
        source.thumbnail((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        row, col = divmod(idx, 2)
        x = margin + col * (thumb_w + margin)
        y = margin + row * (thumb_h + label_h + margin)
        px = x + (thumb_w - source.width) // 2
        py = y + label_h + (thumb_h - source.height) // 2
        canvas.paste(source, (px, py))
        draw.text((x, y + 5), code, fill="#27313d", font=font)
    out = QA / "M1-M8_contact_sheet.png"
    canvas.save(out, quality=95)
    return out


def main() -> None:
    rows = []
    details = {}
    for code, stem in STEMS.items():
        paths = {ext: FIG / f"{stem}.{ext}" for ext in ("png", "svg", "pdf")}
        missing = [ext for ext, path in paths.items() if not path.exists() or path.stat().st_size == 0]
        result = {"figure": code, "stem": stem, "missing_formats": ";".join(missing)}
        if not missing:
            result.update(png_checks(paths["png"]))
            result.update(pdf_checks(paths["pdf"]))
        passed = (
            not missing
            and result.get("png_width_px", 0) >= 3000
            and result.get("png_height_px", 0) >= 1800
            and 0.9 <= result.get("aspect_ratio", 0) <= 2.0
            and result.get("nonwhite_fraction", 0) >= .05
            and result.get("outside_word_count", 1) == 0
            and result.get("overlap_count", 1) == 0
            and not result.get("bad_encoding_tokens", ["bad"])
        )
        result["status"] = "PASS" if passed else "FAIL"
        rows.append({k: v for k, v in result.items() if k != "overlap_examples"})
        details[code] = result
    sheet = contact_sheet()
    table = pd.DataFrame(rows)
    table.to_csv(QA / "M1-M8_QA_report.tsv", sep="\t", index=False)
    summary = {
        "all_pass": bool((table.status == "PASS").all()),
        "pass_count": int((table.status == "PASS").sum()),
        "figure_count": len(table),
        "contact_sheet": str(sheet),
        "figures": details,
    }
    (QA / "M1-M8_QA_report.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(table[["figure", "status", "png_width_px", "png_height_px", "overlap_count", "outside_word_count"]].to_string(index=False))
    print(json.dumps({"all_pass": summary["all_pass"], "pass_count": summary["pass_count"]}))


if __name__ == "__main__":
    main()
