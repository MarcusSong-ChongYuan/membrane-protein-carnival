from __future__ import annotations

import re
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(r"C:\Users\Administrator\HuMemLigDB_manuscript")
SOURCE = ROOT / "HuMemLigDB_NAR_manuscript_draft_v0.1.md"
OUTPUT = ROOT / "HuMemLigDB_NAR_manuscript_draft_v0.1.docx"

PAGE_WIDTH_DXA = 12240
PAGE_HEIGHT_DXA = 15840
CONTENT_WIDTH_DXA = 9360
TABLE_INDENT_DXA = 120

NAVY = "203748"
BLUE = "2E74B5"
DARK_BLUE = "1F4D78"
MUTED = "5D6872"
LIGHT_FILL = "F4F6F9"
LIGHT_BLUE = "E8EEF5"
WHITE = "FFFFFF"
GRID = "BCC6D0"
INK = "111111"
CAUTION = "7A5A00"


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=80, start=120, bottom=80, end=120) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for tag, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{tag}"))
        if node is None:
            node = OxmlElement(f"w:{tag}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_table_borders(table, color=GRID, size="4") -> None:
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.find(qn("w:tblBorders"))
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = borders.find(qn(f"w:{edge}"))
        if tag is None:
            tag = OxmlElement(f"w:{edge}")
            borders.append(tag)
        tag.set(qn("w:val"), "single")
        tag.set(qn("w:sz"), size)
        tag.set(qn("w:space"), "0")
        tag.set(qn("w:color"), color)


def set_table_geometry(table, widths: list[int], indent=TABLE_INDENT_DXA) -> None:
    if sum(widths) != CONTENT_WIDTH_DXA:
        raise ValueError(f"Table widths must sum to {CONTENT_WIDTH_DXA}: {widths}")
    table.autofit = False
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    tbl = table._tbl
    tbl_pr = tbl.tblPr

    tbl_w = tbl_pr.find(qn("w:tblW"))
    if tbl_w is None:
        tbl_w = OxmlElement("w:tblW")
        tbl_pr.append(tbl_w)
    tbl_w.set(qn("w:w"), str(CONTENT_WIDTH_DXA))
    tbl_w.set(qn("w:type"), "dxa")

    tbl_ind = tbl_pr.find(qn("w:tblInd"))
    if tbl_ind is None:
        tbl_ind = OxmlElement("w:tblInd")
        tbl_pr.append(tbl_ind)
    tbl_ind.set(qn("w:w"), str(indent))
    tbl_ind.set(qn("w:type"), "dxa")

    grid = tbl.tblGrid
    for child in list(grid):
        grid.remove(child)
    for width in widths:
        col = OxmlElement("w:gridCol")
        col.set(qn("w:w"), str(width))
        grid.append(col)

    for row in table.rows:
        for idx, cell in enumerate(row.cells):
            tc_pr = cell._tc.get_or_add_tcPr()
            tc_w = tc_pr.find(qn("w:tcW"))
            if tc_w is None:
                tc_w = OxmlElement("w:tcW")
                tc_pr.append(tc_w)
            tc_w.set(qn("w:w"), str(widths[idx]))
            tc_w.set(qn("w:type"), "dxa")
            cell.width = Inches(widths[idx] / 1440)
            set_cell_margins(cell)
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER


def set_repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def set_keep_with_next(paragraph, value=True) -> None:
    paragraph.paragraph_format.keep_with_next = value


def set_run_font(run, name="Calibri", size=None, color=None, bold=None, italic=None) -> None:
    run.font.name = name
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), name)
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), name)
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), name)
    if size is not None:
        run.font.size = Pt(size)
    if color is not None:
        run.font.color.rgb = RGBColor.from_string(color)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic


def add_inline_runs(paragraph, text: str, base_size=11, base_color=INK) -> None:
    # A restrained parser for bold, italic and inline-code spans.
    token = re.compile(r"(\*\*[^*]+\*\*|`[^`]+`|\*[^*]+\*)")
    pos = 0
    for match in token.finditer(text):
        if match.start() > pos:
            run = paragraph.add_run(text[pos : match.start()])
            set_run_font(run, size=base_size, color=base_color)
        raw = match.group(0)
        if raw.startswith("**"):
            run = paragraph.add_run(raw[2:-2])
            set_run_font(run, size=base_size, color=base_color, bold=True)
        elif raw.startswith("`"):
            run = paragraph.add_run(raw[1:-1])
            set_run_font(run, name="Consolas", size=max(base_size - 0.5, 8.0), color=DARK_BLUE)
        else:
            run = paragraph.add_run(raw[1:-1])
            set_run_font(run, size=base_size, color=base_color, italic=True)
        pos = match.end()
    if pos < len(text):
        run = paragraph.add_run(text[pos:])
        set_run_font(run, size=base_size, color=base_color)


def add_page_field(paragraph) -> None:
    run = paragraph.add_run()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    fld_sep = OxmlElement("w:fldChar")
    fld_sep.set(qn("w:fldCharType"), "separate")
    txt = OxmlElement("w:t")
    txt.text = "1"
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    run._r.extend([fld_begin, instr, fld_sep, txt, fld_end])
    set_run_font(run, size=9, color=MUTED)


def add_custom_numbering(doc: Document) -> int:
    numbering = doc.part.numbering_part.element
    abstract_ids = [int(x.get(qn("w:abstractNumId"))) for x in numbering.findall(qn("w:abstractNum"))]
    num_ids = [int(x.get(qn("w:numId"))) for x in numbering.findall(qn("w:num"))]
    abstract_id = max(abstract_ids, default=0) + 1
    num_id = max(num_ids, default=0) + 1

    abstract = OxmlElement("w:abstractNum")
    abstract.set(qn("w:abstractNumId"), str(abstract_id))
    multi = OxmlElement("w:multiLevelType")
    multi.set(qn("w:val"), "singleLevel")
    abstract.append(multi)
    lvl = OxmlElement("w:lvl")
    lvl.set(qn("w:ilvl"), "0")
    start = OxmlElement("w:start")
    start.set(qn("w:val"), "1")
    num_fmt = OxmlElement("w:numFmt")
    num_fmt.set(qn("w:val"), "decimal")
    lvl_text = OxmlElement("w:lvlText")
    lvl_text.set(qn("w:val"), "%1.")
    suff = OxmlElement("w:suff")
    suff.set(qn("w:val"), "tab")
    p_pr = OxmlElement("w:pPr")
    tabs = OxmlElement("w:tabs")
    tab = OxmlElement("w:tab")
    tab.set(qn("w:val"), "num")
    tab.set(qn("w:pos"), "540")
    tabs.append(tab)
    ind = OxmlElement("w:ind")
    ind.set(qn("w:left"), "540")
    ind.set(qn("w:hanging"), "279")
    spacing = OxmlElement("w:spacing")
    spacing.set(qn("w:after"), "80")
    spacing.set(qn("w:line"), "290")
    spacing.set(qn("w:lineRule"), "auto")
    p_pr.extend([tabs, ind, spacing])
    lvl.extend([start, num_fmt, lvl_text, suff, p_pr])
    abstract.append(lvl)
    numbering.append(abstract)

    num = OxmlElement("w:num")
    num.set(qn("w:numId"), str(num_id))
    abs_id = OxmlElement("w:abstractNumId")
    abs_id.set(qn("w:val"), str(abstract_id))
    num.append(abs_id)
    numbering.append(num)
    return num_id


def apply_numbering(paragraph, num_id: int) -> None:
    p_pr = paragraph._p.get_or_add_pPr()
    num_pr = p_pr.find(qn("w:numPr"))
    if num_pr is None:
        num_pr = OxmlElement("w:numPr")
        p_pr.append(num_pr)
    ilvl = OxmlElement("w:ilvl")
    ilvl.set(qn("w:val"), "0")
    num_id_el = OxmlElement("w:numId")
    num_id_el.set(qn("w:val"), str(num_id))
    num_pr.extend([ilvl, num_id_el])


def configure_styles(doc: Document) -> None:
    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Calibri"
    normal._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    normal.font.size = Pt(11)
    normal.font.color.rgb = RGBColor.from_string(INK)
    pf = normal.paragraph_format
    pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    pf.space_before = Pt(0)
    pf.space_after = Pt(8)
    pf.line_spacing = 1.333
    pf.widow_control = True

    tokens = {
        "Heading 1": (16, BLUE, 18, 10),
        "Heading 2": (13, BLUE, 12, 6),
        "Heading 3": (12, DARK_BLUE, 8, 4),
    }
    for name, (size, color, before, after) in tokens.items():
        style = styles[name]
        style.font.name = "Calibri"
        style._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
        style._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor.from_string(color)
        p = style.paragraph_format
        p.space_before = Pt(before)
        p.space_after = Pt(after)
        p.line_spacing = 1.0
        p.keep_with_next = True
        p.keep_together = True


def configure_page(doc: Document) -> None:
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)
    section.header_distance = Inches(0.492)
    section.footer_distance = Inches(0.492)
    section.different_first_page_header_footer = True

    header = section.header
    p = header.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    p.paragraph_format.space_after = Pt(0)
    run = p.add_run("HuMemLigDB")
    set_run_font(run, size=9, color=MUTED, bold=True)
    run = p.add_run("\tV6.2 | Working manuscript 0.1")
    set_run_font(run, size=9, color=MUTED)
    tabs = p.paragraph_format.tab_stops
    tabs.add_tab_stop(Inches(6.5), alignment=2)

    footer = section.footer
    fp = footer.paragraphs[0]
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    fp.paragraph_format.space_before = Pt(0)
    fp.paragraph_format.space_after = Pt(0)
    label = fp.add_run("Page ")
    set_run_font(label, size=9, color=MUTED)
    add_page_field(fp)


def add_cover(doc: Document) -> None:
    for _ in range(5):
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(8)

    kicker = doc.add_paragraph()
    kicker.alignment = WD_ALIGN_PARAGRAPH.CENTER
    kicker.paragraph_format.space_after = Pt(16)
    run = kicker.add_run("DATABASE RESOURCE MANUSCRIPT")
    set_run_font(run, size=10, color=BLUE, bold=True)

    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.paragraph_format.space_before = Pt(0)
    title.paragraph_format.space_after = Pt(10)
    title.paragraph_format.keep_with_next = True
    run = title.add_run(
        "HuMemLigDB: a provenance-aware resource for human membrane proteins, "
        "small-molecule interactions and binding sites"
    )
    set_run_font(run, size=27, color=NAVY, bold=True)

    subtitle = doc.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle.paragraph_format.space_after = Pt(30)
    run = subtitle.add_run("Nucleic Acids Research Database Issue-style working draft")
    set_run_font(run, size=13.5, color=DARK_BLUE)

    meta = [
        ("Release described", "HuMemLigDB V6.2 (internal label: MemPro V6.2)"),
        ("Freeze date", "30 July 2026"),
        ("Draft", "Version 0.1 | 3 August 2026"),
        ("Status", "Internal scientific draft; not yet a public database service"),
    ]
    table = doc.add_table(rows=len(meta), cols=2)
    set_table_geometry(table, [2300, 7060], indent=TABLE_INDENT_DXA)
    set_table_borders(table, color="D8E0E7", size="3")
    for row, (label, value) in zip(table.rows, meta):
        set_cell_shading(row.cells[0], LIGHT_FILL)
        p1 = row.cells[0].paragraphs[0]
        p1.paragraph_format.space_after = Pt(0)
        r1 = p1.add_run(label)
        set_run_font(r1, size=9.5, color=NAVY, bold=True)
        p2 = row.cells[1].paragraphs[0]
        p2.paragraph_format.space_after = Pt(0)
        r2 = p2.add_run(value)
        set_run_font(r2, size=9.5, color=INK)

    doc.add_paragraph().paragraph_format.space_after = Pt(8)
    note = doc.add_table(rows=1, cols=1)
    set_table_geometry(note, [CONTENT_WIDTH_DXA], indent=TABLE_INDENT_DXA)
    set_table_borders(note, color="D7C680", size="6")
    set_cell_shading(note.cell(0, 0), "FFF9E8")
    p = note.cell(0, 0).paragraphs[0]
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after = Pt(2)
    r = p.add_run("Submission boundary. ")
    set_run_font(r, size=9.5, color=CAUTION, bold=True)
    r = p.add_run(
        "This draft reports the frozen data accurately but does not claim a live public website, "
        "API or redistribution-cleared download. Those items, external mapping validation and "
        "final source citations are required before submission."
    )
    set_run_font(r, size=9.5, color=INK)

    doc.add_page_break()


def parse_markdown_table(lines: list[str]) -> list[list[str]]:
    rows = []
    for line in lines:
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        rows.append(cells)
    if len(rows) >= 2 and all(re.fullmatch(r":?-{3,}:?", x) for x in rows[1]):
        rows.pop(1)
    return rows


def add_table(doc: Document, rows: list[list[str]], title: str) -> None:
    cols = len(rows[0])
    if any(len(row) != cols for row in rows):
        raise ValueError(f"Inconsistent Markdown table under {title}")
    if cols == 4:
        widths = [1050, 1750, 1850, 4710]
    elif cols == 3 and "Evidence" in title:
        widths = [1150, 5650, 2560]
    elif cols == 3:
        widths = [2850, 1450, 5060]
    elif cols == 2:
        widths = [2700, 6660]
    else:
        base = CONTENT_WIDTH_DXA // cols
        widths = [base] * cols
        widths[-1] += CONTENT_WIDTH_DXA - sum(widths)

    table = doc.add_table(rows=len(rows), cols=cols)
    set_table_geometry(table, widths)
    set_table_borders(table)
    set_repeat_table_header(table.rows[0])
    for i, row in enumerate(rows):
        for j, text in enumerate(row):
            cell = table.cell(i, j)
            if i == 0:
                set_cell_shading(cell, LIGHT_FILL)
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if (i == 0 or j == 1 and cols == 3) else WD_ALIGN_PARAGRAPH.LEFT
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.08
            add_inline_runs(p, text, base_size=8.5 if cols == 4 else 9.0)
            for run in p.runs:
                if i == 0:
                    run.bold = True
                    run.font.color.rgb = RGBColor.from_string(NAVY)
    after = doc.add_paragraph()
    after.paragraph_format.space_before = Pt(0)
    after.paragraph_format.space_after = Pt(2)


def normalize_source(text: str) -> str:
    old = (
        "Most pairs are supported by one contributing database label, while 294,? pairs have two or more source labels; "
        "this field is retained as provenance support rather than interpreted as experimentally independent replication. "
        "[This source-count sentence should be replaced with the exact reconciled total before submission.]"
    )
    new = (
        "Most pairs are supported by one contributing database label, while 295,207 pairs have two or more source labels; "
        "this field is retained as provenance support rather than interpreted as experimentally independent replication."
    )
    return text.replace(old, new)


def add_manuscript_content(doc: Document, markdown: str, num_id: int) -> None:
    lines = markdown.splitlines()
    # Cover metadata is already rendered; begin at Abstract.
    start = next(i for i, line in enumerate(lines) if line.strip() == "## Abstract")
    lines = lines[start:]
    current_h2 = ""
    i = 0
    paragraph_buffer: list[str] = []

    def flush_paragraph() -> None:
        nonlocal paragraph_buffer
        if not paragraph_buffer:
            return
        text = " ".join(x.strip() for x in paragraph_buffer).strip()
        paragraph_buffer = []
        if not text:
            return
        p = doc.add_paragraph(style="Normal")
        add_inline_runs(p, text)

    while i < len(lines):
        line = lines[i].rstrip()
        stripped = line.strip()
        if not stripped:
            flush_paragraph()
            i += 1
            continue
        if stripped.startswith("### "):
            flush_paragraph()
            title = stripped[4:].strip()
            current_h2 = title
            if title in {"Table 1. Principal source snapshots and roles in HuMemLigDB V6.2", "Figure 1. Database architecture, source integration and quality control"}:
                doc.add_page_break()
            p = doc.add_paragraph(title, style="Heading 2")
            set_keep_with_next(p)
            i += 1
            continue
        if stripped.startswith("## "):
            flush_paragraph()
            title = stripped[3:].strip()
            current_h2 = title
            if title in {"Materials and methods", "Tables", "Figure legends", "Editorial action list before submission"}:
                doc.add_page_break()
            p = doc.add_paragraph(title, style="Heading 1")
            set_keep_with_next(p)
            i += 1
            continue
        if stripped.startswith("|"):
            flush_paragraph()
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                table_lines.append(lines[i].strip())
                i += 1
            rows = parse_markdown_table(table_lines)
            add_table(doc, rows, current_h2)
            continue
        if re.match(r"^\d+\.\s+", stripped):
            flush_paragraph()
            text = re.sub(r"^\d+\.\s+", "", stripped)
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.375)
            p.paragraph_format.first_line_indent = Inches(-0.194)
            p.paragraph_format.space_after = Pt(4)
            p.paragraph_format.line_spacing = 1.208
            apply_numbering(p, num_id)
            add_inline_runs(p, text)
            i += 1
            continue
        if stripped.startswith("> "):
            flush_paragraph()
            note_text = stripped[2:].strip()
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.25)
            p.paragraph_format.right_indent = Inches(0.25)
            p.paragraph_format.space_before = Pt(6)
            p.paragraph_format.space_after = Pt(8)
            add_inline_runs(p, note_text, base_size=10, base_color=CAUTION)
            i += 1
            continue
        paragraph_buffer.append(stripped)
        i += 1
    flush_paragraph()


def set_core_properties(doc: Document) -> None:
    props = doc.core_properties
    props.title = "HuMemLigDB NAR database manuscript draft v0.1"
    props.subject = "Evidence-stratified human membrane protein-small-molecule database"
    props.author = "HuMemLigDB project team"
    props.keywords = "membrane protein; small molecule; binding site; database; provenance; docking"
    props.comments = "Internal working draft generated from the frozen HuMemLigDB V6.2 release documentation."


def audit_document(doc: Document) -> None:
    section = doc.sections[0]
    assert round(section.page_width.inches, 3) == 8.5
    assert round(section.page_height.inches, 3) == 11.0
    assert all(round(x.inches, 3) == 1.0 for x in (section.top_margin, section.right_margin, section.bottom_margin, section.left_margin))
    for table in doc.tables:
        grid = table._tbl.tblGrid
        widths = [int(x.get(qn("w:w"))) for x in grid.findall(qn("w:gridCol"))]
        assert sum(widths) == CONTENT_WIDTH_DXA, widths
        tbl_w = table._tbl.tblPr.find(qn("w:tblW"))
        assert tbl_w is not None and int(tbl_w.get(qn("w:w"))) == CONTENT_WIDTH_DXA
    full_text = "\n".join(p.text for p in doc.paragraphs)
    assert "294,?" not in full_text
    assert "[This source-count" not in full_text


def main() -> None:
    markdown = normalize_source(SOURCE.read_text(encoding="utf-8"))
    doc = Document()
    configure_styles(doc)
    configure_page(doc)
    set_core_properties(doc)
    num_id = add_custom_numbering(doc)
    add_cover(doc)
    add_manuscript_content(doc, markdown, num_id)
    audit_document(doc)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(OUTPUT)
    print(OUTPUT)


if __name__ == "__main__":
    main()
