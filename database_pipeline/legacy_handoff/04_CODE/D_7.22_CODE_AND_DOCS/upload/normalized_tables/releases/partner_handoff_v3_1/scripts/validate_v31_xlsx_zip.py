import hashlib
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from openpyxl import load_workbook

BASE = Path(r"C:\github-repos\upload\normalized_tables\outputs\v3_1_completion")
XLSX = BASE / "normalized_tables_v3_1.xlsx"
EXPECTED = [
    ("1_Binary_Relationships", 101802, 38),
    ("2_Protein_Database", 1773, 42),
    ("3_Small_Molecules", 87611, 37),
    ("4_Pocket_Instances", 6815, 16),
]
NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pr": "http://schemas.openxmlformats.org/package/2006/relationships",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def scan_sheet(z, entry):
    rows = 0
    formulas = 0
    has_pane = False
    has_filter = False
    first = b""
    carry = b""
    with z.open(entry) as f:
        while True:
            chunk = f.read(4 * 1024 * 1024)
            if not chunk:
                break
            if len(first) < 2 * 1024 * 1024:
                first += chunk[: 2 * 1024 * 1024 - len(first)]
            data = carry + chunk
            rows += data.count(b"<row ")
            formulas += len(re.findall(br"<f(?:\s|>)", data))
            has_pane = has_pane or b"<pane " in data
            has_filter = has_filter or b"<autoFilter " in data
            carry = data[-16:]
    row1 = re.search(br"<row[^>]*\br=\"1\"[^>]*>(.*?)</row>", first, re.S)
    if not row1:
        raise RuntimeError(f"No header row in {entry}")
    cells = re.findall(br"<c\s+([^>]*)>", row1.group(1))
    styles = []
    for attrs in cells:
        m = re.search(br"\bs=\"(\d+)\"", attrs)
        styles.append(int(m.group(1)) if m else 0)
    pane = re.search(br"<pane\s+([^>]*)/>", first)
    return {
        "total_rows": rows,
        "data_rows": rows - 1,
        "header_cols": len(cells),
        "formula_tags": formulas,
        "freeze_pane_xml": pane.group(1).decode("utf-8", "replace") if pane else "",
        "has_freeze_pane": has_pane,
        "has_auto_filter": has_filter,
        "header_style_nonzero": all(x > 0 for x in styles),
    }


def main():
    if not XLSX.exists():
        raise FileNotFoundError(XLSX)
    with zipfile.ZipFile(XLSX) as z:
        bad_member = z.testzip()
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        targets = {x.get("Id"): x.get("Target") for x in rels.findall("pr:Relationship", NS)}
        actual_sheets = []
        for s in wb.find("m:sheets", NS):
            name = s.get("name")
            rid = s.get(f"{{{NS['r']}}}id")
            target = targets[rid].lstrip("/")
            entry = target if target.startswith("xl/") else "xl/" + target
            actual_sheets.append((name, entry))
        scans = []
        for (expected_name, expected_rows, expected_cols), (actual_name, entry) in zip(EXPECTED, actual_sheets):
            scan = scan_sheet(z, entry)
            scan.update({
                "name": actual_name,
                "expected_name": expected_name,
                "expected_data_rows": expected_rows,
                "expected_cols": expected_cols,
            })
            scans.append(scan)

    read = load_workbook(XLSX, read_only=True, data_only=False)
    for scan in scans:
        cell = read[scan["name"]]["A1"]
        scan["header_fill_rgb"] = cell.fill.fgColor.rgb
        scan["header_font_bold"] = bool(cell.font.bold)
        scan["header_font_rgb"] = cell.font.color.rgb if cell.font.color and cell.font.color.type == "rgb" else ""
    read.close()

    passed = (
        bad_member is None
        and len(scans) == len(EXPECTED)
        and all(x["name"] == x["expected_name"] for x in scans)
        and all(x["data_rows"] == x["expected_data_rows"] for x in scans)
        and all(x["header_cols"] == x["expected_cols"] for x in scans)
        and all(x["formula_tags"] == 0 for x in scans)
        and all(x["has_freeze_pane"] and x["has_auto_filter"] for x in scans)
        and all(x["header_style_nonzero"] and x["header_font_bold"] for x in scans)
    )
    report = {
        "output": str(XLSX),
        "bytes": XLSX.stat().st_size,
        "sha256": sha256(XLSX),
        "zip_integrity_error": bad_member,
        "sheets": scans,
        "validation_passed": passed,
    }
    (BASE / "V3_1_XLSX_VALIDATION_REPORT.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
