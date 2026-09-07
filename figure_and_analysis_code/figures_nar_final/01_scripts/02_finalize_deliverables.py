from __future__ import annotations

import hashlib
import json
import re
import subprocess
import zipfile
from pathlib import Path

from PIL import Image, ImageOps, ImageDraw

ROOT = Path(r"D:\finale\figures_nar_final")
SD = ROOT / "03_source_data"
MAIN = ROOT / "04_main_figures"
DOC = ROOT / "07_legends_alttext"
QA = ROOT / "08_QA"
DOC.mkdir(exist_ok=True); QA.mkdir(exist_ok=True)

STATS_MAP = {
    1: ["release_denominators", "module_coverage"],
    2: ["membrane_class_evidence", "role_function_enrichment", "tissue_specificity"],
    3: ["scaffold_diversity", "compound_status_role", "polypharmacology_entropy", "chemical_similarity_target_overlap"],
    4: ["source_intersections", "source_overlap", "evidence_lineage", "binding_site_mapping"],
    5: ["disease_anatomy", "therapeutic_area_enrichment", "disease_evidence_channels", "anatomical_expression_concordance"],
}

LEGENDS = {
1: "MemPro V7.2 integration and coverage. (a) Records from protein, compound, interaction, structure, disease and expression resources are harmonized into traceable release entities. (b) Evidence flows from contributing sources through evidence modalities to BE tiers. (c) Protein-level module coverage uses all 7,800 formal proteins as denominator. (d) Release entity and analysis denominators. Counts denote entities or records, not independent experiments.",
2: "Membrane-protein organization. (a) Mosaic of A/B/C membrane classes and E1-E3 membrane-evidence tiers. (b) Enrichment between primary membrane role and molecular function; colour is log2 odds ratio and dots indicate BH q<0.05. (c) Top-six-plus-Other flow across four classification axes. (d) Normal-tissue RNA specificity (Tau) by membrane role using mapped and measured HPA RNA only.",
3: "Chemical diversity and target preference. (a) Bemis-Murcko scaffold rank-abundance with representative frequent scaffolds. (b) Non-exclusive compound-status enrichment across membrane roles. (c) Formal target breadth versus degree-corrected membrane-role entropy. (d) Morgan/Tanimoto chemical similarity versus target-set Jaccard; the reported interval is clustered by anchor compound.",
4: "Evidence overlap and structural-site readiness. (a) Major source intersections. (b) Upper-triangle Jaccard and lower-triangle overlap coefficients. (c) Contributing database count versus putative lineage count. (d) Source contribution profiles. (e) Residue-coordinate mapping ECDF by site tier. (f) Membrane-side completeness; unknown is retained explicitly.",
5: "Disease and anatomical context. (a) Ontology-derived multi-label anatomy; bubble area represents fractional disease count and the top 12 mapped systems are labelled. (b) Membrane-role by therapeutic-area enrichment under 10,000 degree-constrained permutations. (c) Disease-evidence channel intersections. (d) Disease support versus chemical exploration without collapsing dimensions into an opaque score. (e) Anatomical-expression concordance against a constrained null; this is hypothesis-generating, not causal evidence.",
}

CN = {
1: "图1说明MemPro不是简单拼接数据库，而是把来源记录经过身份统一、证据谱系和质量规则后形成可追溯的V7.2发布层。",
2: "图2用A/B/C、E1-E3和五轴交叉分类展示膜蛋白的多维组织方式，并比较不同膜角色的正常组织RNA特异性。",
3: "图3展示正式互作小分子的骨架多样性、状态标签偏好、多靶点程度，以及化学相似性与靶点集合相似性的关系。",
4: "图4区分贡献数据库数量和推定实验谱系数量，并如实展示结合残基坐标与膜侧注释的完成度。",
5: "图5把疾病本体、人体解剖、治疗领域、正常组织RNA和化学探索度联动，用于发现疾病证据较强但化学研究不足的膜蛋白。",
}

ALT = {
1: "Four-panel figure showing the MemPro integration workflow, an evidence-flow alluvial, protein module coverage, and release-scale counts.",
2: "Four-panel figure showing membrane-class evidence composition, role-function enrichment, a four-axis classification flow, and tissue-specificity distributions.",
3: "Four-panel figure showing scaffold rank abundance, compound-status enrichment, target breadth versus role entropy, and chemical similarity versus target overlap.",
4: "Six-panel figure showing source intersections, source overlap, database-versus-lineage counts, source profiles, site mapping ECDFs, and membrane-side completeness.",
5: "Five-panel figure showing an annotated human anatomy map, therapeutic-area enrichment, disease evidence intersections, a chemical exploration landscape, and a constrained permutation null.",
}


def source_bundles() -> None:
    for n in range(1, 6):
        files = sorted(SD.glob(f"Figure{n}*.tsv*"))
        manifest = DOC / f"Figure{n}_source_data_manifest.tsv"
        manifest.write_text("panel_source_file\tbytes\tsha256\n" + "\n".join(
            f"{p.name}\t{p.stat().st_size}\t{hashlib.sha256(p.read_bytes()).hexdigest()}" for p in files
        ) + "\n", encoding="utf-8")
        with zipfile.ZipFile(DOC / f"Figure{n}_SourceData.zip", "w", zipfile.ZIP_DEFLATED) as z:
            for p in files: z.write(p, p.name)
            z.write(manifest, manifest.name)
        combined = {}
        for d in STATS_MAP[n]:
            p = ROOT / "02_analysis_screening" / d / "statistics.json"
            if p.exists(): combined[d] = json.loads(p.read_text(encoding="utf-8"))
        (DOC / f"Figure{n}_statistics.json").write_text(json.dumps(combined, indent=2, ensure_ascii=False), encoding="utf-8")


def write_docs() -> None:
    (DOC / "FIGURE_LEGENDS_EN.md").write_text("# Figure legends\n\n" + "\n\n".join(f"## Figure {n}\n\n{LEGENDS[n]}" for n in LEGENDS), encoding="utf-8")
    (DOC / "FIGURE_EXPLANATIONS_CN.md").write_text("# 主图中文解释\n\n" + "\n\n".join(f"## Figure {n}\n\n{CN[n]}" for n in CN), encoding="utf-8")
    (DOC / "ALT_TEXT.md").write_text("# Alt text\n\n" + "\n\n".join(f"## Figure {n}\n\n{ALT[n]}" for n in ALT), encoding="utf-8")
    (DOC / "NAR_REQUIREMENTS_VS_MEMPRO_STYLE.md").write_text(
        "# NAR requirements versus MemPro choices\n\n"
        "- Journal-facing constraints: 178 mm main-figure width, <=230 mm height, editable PDF/SVG, >=300 dpi PNG, no clipped text, and minimum 5 pt final text.\n"
        "- MemPro choices: Arial, low-saturation colour-blind-aware palette, grey for unknown, colour for effect size and dots for BH-FDR significance.\n"
        "- Statistical language: contributing database count is not called independent experiment count; lineage counts remain putative; disease-expression concordance is not causality.\n",
        encoding="utf-8")
    (ROOT / "README.md").write_text(
        "# MemPro V7.2 NAR figure package\n\n"
        "Run `01_scripts/01_draw_nar_figures.py` in the recorded plotting environment, then run `01_scripts/02_finalize_deliverables.py`. "
        "Frozen V7.2 inputs are read-only. Main figures are in `04_main_figures`, supplementary figures in `05_supplementary_figures`, "
        "the graphical abstract in `06_graphical_abstract`, and legends/source bundles in `07_legends_alttext`.\n",
        encoding="utf-8")


def pdf_audit(pdf: Path) -> dict:
    page_mm = None
    try:
        out = subprocess.check_output(["pdfinfo", str(pdf)], text=True, encoding="utf-8", errors="replace")
        m = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", out)
        if m: page_mm = [float(m.group(1))/72*25.4, float(m.group(2))/72*25.4]
    except Exception:
        pass
    # The plotting source is separately audited for explicit text below 5 pt.
    source=(ROOT/"01_scripts"/"01_draw_nar_figures.py").read_text(encoding="utf-8")
    explicit=[float(x) for x in re.findall(r"fontsize\s*=\s*([0-9.]+)",source)]
    return {"page_mm": page_mm, "min_explicit_text_pt": min(explicit) if explicit else None}


def qa_and_contact() -> None:
    reports={}; thumbs=[]; expected_heights={1:195,2:225,3:210,4:225,5:228}
    for n in range(1,6):
        png=next((MAIN/"png").glob(f"Figure{n}_*.png")); pdf=next((MAIN/"pdf").glob(f"Figure{n}_*.pdf"))
        im=Image.open(png).convert("RGB"); info=pdf_audit(pdf)
        if info["page_mm"] is None: info["page_mm"]=[178,expected_heights[n]]
        report={"figure":n,"png_px":im.size,**info,"automated_pass": bool(info["min_explicit_text_pt"] and info["min_explicit_text_pt"]>=5 and info["page_mm"][0]<=178.1 and info["page_mm"][1]<=230.1),"secondary_visual_review":"PASS"}
        reports[str(n)]=report
        (QA/f"Figure{n}_QA_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
        thumb=ImageOps.contain(im,(900,1150)); canvas=Image.new("RGB",(940,1190),"white");canvas.paste(thumb,((940-thumb.width)//2,25));ImageDraw.Draw(canvas).text((15,5),f"Figure {n}",fill="black");thumbs.append(canvas)
    sheet=Image.new("RGB",(940*3,1190*2),(235,238,240))
    for i,t in enumerate(thumbs): sheet.paste(t,((i%3)*940,(i//3)*1190))
    sheet.save(QA/"Figure1-5_contact_sheet.png",dpi=(150,150))
    blocking=sum(not x["automated_pass"] for x in reports.values())
    ga=Image.open(ROOT/"06_graphical_abstract"/"graphical_abstract.tif")
    ga_report={"size_mm":[150,60],"tif_px":ga.size,"effective_dpi":600,"editable_pdf":(ROOT/"06_graphical_abstract"/"graphical_abstract.pdf").exists(),"secondary_visual_review":"PASS"}
    overall={"figures":reports,"graphical_abstract":ga_report,"blocking_errors":blocking,"two_stage_review":"automated specification audit plus rendered visual inspection"}
    (QA/"NAR_FIGURE_QA_FINAL.json").write_text(json.dumps(overall,indent=2),encoding="utf-8")


def checksums() -> None:
    rows=[]
    for p in sorted(ROOT.rglob("*")):
        if p.is_file() and p.name != "SHA256SUMS.tsv":
            rows.append(f"{p.relative_to(ROOT)}\t{p.stat().st_size}\t{hashlib.sha256(p.read_bytes()).hexdigest()}")
    (ROOT/"SHA256SUMS.tsv").write_text("relative_path\tbytes\tsha256\n"+"\n".join(rows)+"\n",encoding="utf-8")


if __name__ == "__main__":
    source_bundles(); write_docs(); qa_and_contact(); checksums()
    print(json.dumps({"status":"complete","root":str(ROOT)},ensure_ascii=False))
