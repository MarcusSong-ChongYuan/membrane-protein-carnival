import json
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "p": "http://schemas.openxmlformats.org/package/2006/relationships",
}

def col_num(ref):
    letters = re.match(r"[A-Z]+", ref).group()
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n

def analyze(path):
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", NS):
                shared.append("".join(t.text or "" for t in si.iterfind(".//m:t", NS)))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {r.attrib["Id"]: r.attrib["Target"] for r in rels}
        result = {"file": str(path), "sheets": []}
        for s in wb.find("m:sheets", NS):
            name = s.attrib["name"]
            target = relmap[s.attrib[f"{{{NS['r']}}}id"]].lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            root = ET.fromstring(z.read(target))
            dim = root.find("m:dimension", NS)
            rows = []
            sheet_data = root.find("m:sheetData", NS)
            for row in list(sheet_data)[:8]:
                vals = {}
                for c in row.findall("m:c", NS):
                    ref = c.attrib["r"]
                    typ = c.attrib.get("t")
                    v = c.find("m:v", NS)
                    inline = c.find("m:is", NS)
                    value = ""
                    if typ == "s" and v is not None:
                        value = shared[int(v.text)]
                    elif typ == "inlineStr" and inline is not None:
                        value = "".join(t.text or "" for t in inline.iterfind(".//m:t", NS))
                    elif v is not None:
                        value = v.text or ""
                    vals[col_num(ref)] = value
                if vals:
                    rows.append([vals.get(i, "") for i in range(1, min(max(vals), 60) + 1)])
            result["sheets"].append({"name": name, "dimension": dim.attrib.get("ref") if dim is not None else "", "sample": rows})
        return result

for f in [
    Path(r"C:\Users\Administrator\Desktop\semifinal\normalized_tables_v3_1.xlsx"),
    Path(r"C:\Users\Administrator\Desktop\semifinal\protein_gene_disease_rule_b_summary_v1.xlsx"),
]:
    print(json.dumps(analyze(f), ensure_ascii=False))
