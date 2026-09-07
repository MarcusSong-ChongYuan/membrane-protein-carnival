from __future__ import annotations

import io
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image

SOURCE = Path(r"C:\Users\Administrator\Downloads\不冤.pptx")
OUTPUT = Path(r"D:\finale\figures_nar_final\09_assets\human_anatomy_highres_v2.png")
CANVAS_W = 3000

NS = {
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "p": "http://schemas.openxmlformats.org/presentationml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


def main() -> None:
    with zipfile.ZipFile(SOURCE) as zf:
        pres = ET.fromstring(zf.read("ppt/presentation.xml"))
        size = pres.find("p:sldSz", NS)
        slide_w = int(size.attrib["cx"])
        slide_h = int(size.attrib["cy"])
        canvas_h = round(CANVAS_W * slide_h / slide_w)

        rel_root = ET.fromstring(zf.read("ppt/slides/_rels/slide1.xml.rels"))
        rels = {x.attrib["Id"]: x.attrib["Target"].split("/")[-1] for x in rel_root}
        slide = ET.fromstring(zf.read("ppt/slides/slide1.xml"))

        canvas = Image.new("RGBA", (CANVAS_W, canvas_h), (255, 255, 255, 0))
        for pic in slide.findall(".//p:pic", NS):
            blip = pic.find("p:blipFill/a:blip", NS)
            rid = blip.attrib[f"{{{NS['r']}}}embed"]
            media = rels[rid]
            off = pic.find("p:spPr/a:xfrm/a:off", NS)
            ext = pic.find("p:spPr/a:xfrm/a:ext", NS)
            x = round(int(off.attrib["x"]) / slide_w * CANVAS_W)
            y = round(int(off.attrib["y"]) / slide_h * canvas_h)
            w = round(int(ext.attrib["cx"]) / slide_w * CANVAS_W)
            h = round(int(ext.attrib["cy"]) / slide_h * canvas_h)
            image = Image.open(io.BytesIO(zf.read(f"ppt/media/{media}"))).convert("RGBA")
            image = image.resize((w, h), Image.Resampling.LANCZOS)
            canvas.alpha_composite(image, (x, y))

    # Crop the empty slide margins while preserving the complete hands and body.
    alpha_bbox = canvas.getchannel("A").getbbox()
    pad = 28
    crop = (
        max(0, alpha_bbox[0] - pad),
        max(0, alpha_bbox[1] - pad),
        min(canvas.width, alpha_bbox[2] + pad),
        min(canvas.height, alpha_bbox[3] + pad),
    )
    result = canvas.crop(crop)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.save(OUTPUT, optimize=True)
    print({"output": str(OUTPUT), "size_px": result.size, "source_layers": 16})


if __name__ == "__main__":
    main()
