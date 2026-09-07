import fs from "node:fs/promises";
import path from "node:path";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const TMP = "D:/7.22/evidence_expansion_v2_working/deliverables/mempro_ppt_20260803";
const OUT = "D:/finale/MemPro_V6.2_数据库建立流程与科学意义.pptx";
const FIG = "D:/finale/02_展示图表_V6.2/png";

const W = 1280;
const H = 720;
const C = {
  navy: "#12324A",
  navy2: "#1B4B66",
  blue: "#4D87C7",
  teal: "#2AA79B",
  green: "#68B984",
  coral: "#E97862",
  gold: "#D5A600",
  purple: "#7B61A8",
  ink: "#19313D",
  muted: "#58707C",
  pale: "#EAF3F4",
  paleBlue: "#E9F1FA",
  paleCoral: "#FBEDE9",
  line: "#D4E1E5",
  bg: "#F8FBFA",
  white: "#FFFFFF",
};
const FONT = "Microsoft YaHei";

async function bytes(file) {
  const b = await fs.readFile(file);
  return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);
}

function addBox(slide, name, pos, fill = "none", lineFill = "none", radius = 0) {
  return slide.shapes.add({
    geometry: radius ? "roundRect" : "rect",
    name,
    position: pos,
    fill,
    line: { style: "solid", fill: lineFill, width: lineFill === "none" ? 0 : 1 },
    ...(radius ? { borderRadius: radius } : {}),
  });
}

function addText(slide, name, text, pos, size = 24, color = C.ink, bold = false, align = "left", valign = "top") {
  const s = slide.shapes.add({
    geometry: "textbox",
    name,
    position: pos,
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  s.text = text;
  s.text.style = {
    fontSize: size,
    color,
    bold,
    alignment: align,
    verticalAlignment: valign,
    typeface: FONT,
    autoFit: "shrinkText",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  };
  return s;
}

function addTitle(slide, title, kicker = "") {
  if (kicker) addText(slide, "section-kicker", kicker, { left: 56, top: 28, width: 320, height: 24 }, 15, C.teal, true);
  addText(slide, "slide-title", title, { left: 56, top: 52, width: 1168, height: 60 }, 48, C.navy, true);
  addBox(slide, "title-accent", { left: 56, top: 113, width: 88, height: 5 }, C.teal, "none", 2);
}

function addFooter(slide, n) {
  addText(slide, "footer-left", "MemPro V6.2 · 人类膜蛋白—小分子互作数据库", { left: 56, top: 690, width: 520, height: 18 }, 12, "#78909A");
  addText(slide, "footer-page", String(n).padStart(2, "0"), { left: 1175, top: 687, width: 48, height: 20 }, 13, "#78909A", true, "right");
}

function baseSlide(p, n, title, kicker = "") {
  const slide = p.slides.add();
  slide.background.fill = C.bg;
  addTitle(slide, title, kicker);
  addFooter(slide, n);
  return slide;
}

function addNotes(slide, body, sources) {
  slide.speakerNotes.textFrame.setText(`${body}\n\n[Sources]\n${sources.map((x) => `- ${x}`).join("\n")}`);
  slide.speakerNotes.setVisible(true);
}

function addRule(slide, x, y, w, color = C.line, h = 2) {
  addBox(slide, `rule-${x}-${y}`, { left: x, top: y, width: w, height: h }, color, "none");
}

function addMetric(slide, x, y, w, number, label, color, detail = "") {
  addText(slide, `metric-number-${x}-${y}`, number, { left: x, top: y, width: w, height: 60 }, 52, color, true);
  addText(slide, `metric-label-${x}-${y}`, label, { left: x, top: y + 62, width: w, height: 34 }, 22, C.navy, true);
  if (detail) addText(slide, `metric-detail-${x}-${y}`, detail, { left: x, top: y + 99, width: w, height: 56 }, 17, C.muted);
}

async function addFigure(slide, file, alt, pos, crop = undefined) {
  return slide.images.add({
    blob: await bytes(file),
    contentType: "image/png",
    alt,
    fit: "contain",
    position: pos,
    ...(crop ? { crop } : {}),
  });
}

function addFlatColumn(slide, x, y, w, heading, body, color, index) {
  addBox(slide, `col-accent-${index}`, { left: x, top: y, width: 8, height: 54 }, color, "none", 3);
  addText(slide, `col-head-${index}`, heading, { left: x + 22, top: y, width: w - 22, height: 54 }, 30, C.navy, true, "left", "middle");
  addText(slide, `col-body-${index}`, body, { left: x, top: y + 78, width: w, height: 250 }, 22, C.ink);
}

async function main() {
  await fs.mkdir(path.dirname(OUT), { recursive: true });
  await fs.mkdir(path.join(TMP, "previews"), { recursive: true });
  const p = Presentation.create({ slideSize: { width: W, height: H } });

  // 1. Title
  {
    const s = p.slides.add();
    s.background.fill = C.navy;
    addBox(s, "title-membrane-band", { left: 0, top: 0, width: 24, height: H }, C.teal, "none");
    addBox(s, "title-accent-coral", { left: 24, top: 0, width: 8, height: H }, C.coral, "none");
    addText(s, "title-kicker", "DATABASE CONSTRUCTION · EVIDENCE INTEGRATION · TRANSLATIONAL USE", { left: 80, top: 84, width: 1030, height: 28 }, 18, "#77D4C9", true);
    addText(s, "deck-title", "MemPro V6.2", { left: 80, top: 158, width: 820, height: 90 }, 72, C.white, true);
    addText(s, "deck-subtitle", "人类膜蛋白—小分子互作数据库\n建立流程、科学意义与下一步", { left: 80, top: 260, width: 970, height: 150 }, 48, C.white, true);
    addText(s, "deck-meta", "冻结版本：2026-07-30  ·  汇报稿：2026-08-03", { left: 82, top: 594, width: 670, height: 28 }, 18, "#BDD1DB");
    addText(s, "deck-tagline", "从分散证据到可追溯、可计算、可验证的数据底座", { left: 708, top: 556, width: 480, height: 72 }, 24, "#E7F4F2", true, "right");
    addNotes(s, "开场强调：这不是单一来源的条目集合，而是围绕身份、证据与可追溯性构建的研究型数据库。", [
      "Internal: D:\\finale\\01_正式数据_V6.2\\RELEASE_INFO_v6_2.json",
      "Internal: D:\\finale\\01_正式数据_V6.2\\V62_VALIDATION_REPORT.json",
    ]);
  }

  // 2. Background
  {
    const s = baseSlide(p, 2, "膜蛋白连接细胞环境，也连接疾病机制与药物作用", "01 · 研究背景");
    addText(s, "background-lead", "人类蛋白质组约有 2 万个 canonical 蛋白条目；其中膜蛋白承担信号感知、物质转运、细胞黏附和能量转换等关键任务。", { left: 56, top: 150, width: 1120, height: 76 }, 29, C.ink, true);
    addRule(s, 56, 248, 1168);
    addFlatColumn(s, 70, 286, 338, "生物学界面", "膜蛋白决定细胞如何感知配体、离子、营养与机械信号，是器官功能与疾病表型之间的重要桥梁。", C.teal, 1);
    addFlatColumn(s, 470, 286, 338, "药理学界面", "受体、离子通道、转运体与膜酶是小分子调控的核心类别；构象、膜环境和复合体状态直接影响结合。", C.blue, 2);
    addFlatColumn(s, 870, 286, 320, "计算挑战", "结构异质、身份不统一、证据粒度不同，使跨数据库比较和大规模 docking 难以直接开展。", C.coral, 3);
    addNotes(s, "先建立必要性：膜蛋白是环境—细胞—药物的共同界面，但其结构状态和膜环境让数据整合与计算设计比可溶性蛋白更复杂。", [
      "UniProt human proteome overview: https://www.uniprot.org/release-notes/homo_sapiens",
      "Li H et al. Nat Biotechnol 2024. https://doi.org/10.1038/s41587-023-01987-2",
    ]);
  }

  // 3. Gap
  {
    const s = baseSlide(p, 3, "瓶颈不是缺少数据，而是证据分散且不能直接相加", "01 · 研究背景");
    addText(s, "gap-lead", "同一个蛋白、同一个化合物、同一次实验，在不同数据库中可能使用不同 ID、粒度和结论。", { left: 56, top: 148, width: 1120, height: 60 }, 30, C.navy, true);
    const xs = [60, 355, 650, 945];
    const heads = ["蛋白身份", "化学身份", "证据语义", "结构与位置"];
    const bodies = [
      "canonical/isoform\n基因名/复合物\n旧 ID/新 ID",
      "名称/CID/InChIKey\n母体/盐型\n立体异构体/电荷",
      "结合/活性/功能\n定量/定性\n阳性/阴性/冲突",
      "PDB/链/组装体\n配体 HET\n残基坐标与膜方向",
    ];
    const colors = [C.blue, C.teal, C.coral, C.purple];
    for (let i = 0; i < 4; i++) {
      addText(s, `gap-head-${i}`, heads[i], { left: xs[i], top: 258, width: 245, height: 46 }, 29, colors[i], true);
      addRule(s, xs[i], 316, 230, colors[i], 4);
      addText(s, `gap-body-${i}`, bodies[i], { left: xs[i], top: 346, width: 245, height: 160 }, 23, C.ink, false, "left");
    }
    addBox(s, "gap-conclusion", { left: 56, top: 548, width: 1168, height: 90 }, C.pale, "none", 16);
    addText(s, "gap-conclusion-text", "如果跳过身份统一与证据门控：重复记录会被误认为独立支持，盐型会被误认为新化合物，活性与直接结合会被混为一谈。", { left: 84, top: 568, width: 1112, height: 52 }, 24, C.navy, true, "center", "middle");
    addNotes(s, "这一页解释为什么不能简单下载后拼表。数据库价值主要来自 identity resolution、证据语义标准化和冲突保留。", [
      "Internal: D:\\finale\\01_正式数据_V6.2\\METHODS_v6_2.md",
      "PubChem BioAssay documentation: https://pubchem.ncbi.nlm.nih.gov/docs/bioassay",
      "ChEMBL: https://www.ebi.ac.uk/chembl/",
    ]);
  }

  // 4. Data model
  {
    const s = baseSlide(p, 4, "MemPro用四张主表回答四类相互关联的问题", "02 · 建库目标");
    const nodes = [
      { x: 70, y: 210, w: 250, h: 120, c: C.paleBlue, t: "膜蛋白主表", d: "它是谁？\n属于哪类膜蛋白？\n证据有多强？" },
      { x: 960, y: 210, w: 250, h: 120, c: C.paleCoral, t: "小分子主表", d: "它是什么结构？\n母体与 form 如何对应？" },
      { x: 70, y: 450, w: 250, h: 120, c: "#EEF5EA", t: "蛋白—基因—疾病", d: "关联什么疾病？\n证据来自哪里？" },
      { x: 960, y: 450, w: 250, h: 120, c: "#F1EDF7", t: "结合与位点证据", d: "是否互作？强度如何？\n在哪个残基附近结合？" },
    ];
    // Relationship arrows first, so they remain behind the entity nodes.
    addBox(s, "arrow-left-top", { left: 310, top: 255, width: 170, height: 30 }, "#BFD4DD", "none", 10);
    addBox(s, "arrow-right-top", { left: 800, top: 255, width: 170, height: 30 }, "#BFD4DD", "none", 10);
    addBox(s, "arrow-left-bottom", { left: 310, top: 495, width: 170, height: 30 }, "#BFD4DD", "none", 10);
    addBox(s, "arrow-right-bottom", { left: 800, top: 495, width: 170, height: 30 }, "#BFD4DD", "none", 10);
    addBox(s, "arrow-vertical", { left: 625, top: 350, width: 30, height: 90 }, "#BFD4DD", "none", 10);
    for (let i = 0; i < nodes.length; i++) {
      const n = nodes[i];
      addBox(s, `model-node-${i}`, { left: n.x, top: n.y, width: n.w, height: n.h }, n.c, C.line, 18);
      addText(s, `model-title-${i}`, n.t, { left: n.x + 20, top: n.y + 18, width: n.w - 40, height: 34 }, 27, C.navy, true, "center");
      addText(s, `model-body-${i}`, n.d, { left: n.x + 22, top: n.y + 58, width: n.w - 44, height: 54 }, 18, C.ink, false, "center");
    }
    addBox(s, "model-center", { left: 460, top: 286, width: 360, height: 170 }, C.navy, "none", 22);
    addText(s, "model-center-title", "统一身份层", { left: 505, top: 318, width: 270, height: 46 }, 34, C.white, true, "center");
    addText(s, "model-center-body", "UniProt accession\ncanonical compound / form\n来源、证据与版本可追溯", { left: 505, top: 372, width: 270, height: 70 }, 21, "#DCECEF", false, "center");
    addText(s, "model-bottom", "四张表通过稳定主键连接；复杂表达矩阵、负证据和生物学组装保留为标准化附表。", { left: 220, top: 614, width: 840, height: 38 }, 21, C.muted, true, "center");
    addNotes(s, "四张主表是项目主线。表达定位、负证据、亚基/生物学组装体等高维信息单独标准化，只把一行一蛋白摘要合并回主表。", [
      "Internal: D:\\finale\\01_正式数据_V6.2\\README_v6_2.md",
      "Internal: D:\\finale\\00_交接说明\\01_项目与数据交接说明.md",
    ]);
  }

  // 5. Sources
  {
    const s = baseSlide(p, 5, "V6.2整合14个来源，并固定版本与贡献边界", "03 · 数据来源");
    await addFigure(s, `${FIG}/M0_source_versions_contributions.png`, "MemPro source versions and contributions", { left: 52, top: 128, width: 820, height: 530 });
    addText(s, "source-side-head", "四类来源", { left: 914, top: 158, width: 270, height: 40 }, 30, C.navy, true);
    const groups = [
      ["蛋白", "UniProt · HPA · HTP · Membranome", C.blue],
      ["结构", "PDBe · OPM · PDBTM · PDBbind", C.green],
      ["互作", "ChEMBL · BindingDB · PubChem · BRENDA", C.coral],
      ["疾病", "Open Targets · UniProt疾病注释", C.purple],
    ];
    groups.forEach((g, i) => {
      addBox(s, `source-line-${i}`, { left: 916, top: 224 + i * 92, width: 8, height: 58 }, g[2], "none", 4);
      addText(s, `source-group-${i}`, g[0], { left: 940, top: 222 + i * 92, width: 220, height: 28 }, 22, g[2], true);
      addText(s, `source-db-${i}`, g[1], { left: 940, top: 254 + i * 92, width: 250, height: 52 }, 17, C.ink);
    });
    addNotes(s, "说明版本固定的重要性：数据库内容会持续变化，必须把来源版本、抓取日期和字段贡献写入 release，才能让统计与论文可复现。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M0_source_versions_contributions.png",
      "Internal source catalog: D:\\finale\\00_交接说明\\05_SOURCE_VERSIONS.tsv",
      "ChEMBL release portal: https://www.ebi.ac.uk/chembl/",
      "HPA downloads: https://www.proteinatlas.org/about/download",
      "PDBe: https://www.ebi.ac.uk/pdbe/",
    ]);
  }

  // 6. Integration pipeline
  {
    const s = baseSlide(p, 6, "每条记录经过身份、层级、去重、证据分级与发布门控", "04 · 数据处理流程");
    await addFigure(s, `${FIG}/M1_database_architecture_quality.png`, "Database architecture, integration and quality control", { left: 42, top: 128, width: 880, height: 530 });
    addText(s, "pipeline-head", "核心原则", { left: 952, top: 150, width: 250, height: 38 }, 30, C.navy, true);
    const principles = [
      ["ID mapping", "统一 UniProt、CID、ChEMBL ID 与完整 InChIKey"],
      ["Parent / form", "母体与盐型、立体化学、质子化形式分层"],
      ["Cross-source dedup", "同一实体与同一实验语境下跨来源去重"],
      ["Evidence tier", "膜证据 E1–E3；结合证据 BE1–BE3"],
      ["Release gate", "主键、外键、数量、冲突与哈希全部通过才冻结"],
    ];
    principles.forEach((v, i) => {
      addText(s, `p-title-${i}`, v[0], { left: 952, top: 212 + i * 82, width: 250, height: 24 }, 19, i === 4 ? C.coral : C.teal, true);
      addText(s, `p-body-${i}`, v[1], { left: 952, top: 240 + i * 82, width: 260, height: 48 }, 16, C.ink);
    });
    addNotes(s, "M1A/M1B 中的 integration 指的是完整处理链，而不是简单合并。名称只产生候选；结构身份唯一时才自动合并。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M1_database_architecture_quality.png",
      "Internal methods: D:\\finale\\01_正式数据_V6.2\\METHODS_v6_2.md",
    ]);
  }

  // 7. Release overview
  {
    const s = baseSlide(p, 7, "冻结版本已经形成可验证、可追溯的研究数据资产", "05 · V6.2成果");
    addMetric(s, 60, 170, 270, "10,997", "膜蛋白审计行", C.blue, "其中 7,904 条进入默认 E1/E2 高置信层");
    addMetric(s, 360, 170, 270, "2,016,064", "canonical 小分子", C.teal, "另有 2,020,898 条 form—parent 层级记录");
    addMetric(s, 660, 170, 270, "1,502,456", "阳性蛋白—小分子关系", C.coral, "由 3,003,306 条阳性证据汇总");
    addMetric(s, 960, 170, 260, "7,784,223", "正式负证据", C.purple, "汇总为 6,331,307 组负关系");
    addRule(s, 56, 358, 1168);
    addMetric(s, 110, 404, 300, "95,598", "结合位点实例", C.green, "结构、链、配体与残基级记录");
    addMetric(s, 490, 404, 300, "10,264", "蛋白—基因—疾病关系", C.gold, "保留来源 ID 与证据渠道");
    addText(s, "pass-big", "PASS", { left: 915, top: 418, width: 230, height: 76 }, 58, C.teal, true, "center");
    addText(s, "pass-small", "V6.2 validation\n阻断错误 = 0", { left: 915, top: 500, width: 230, height: 60 }, 22, C.navy, true, "center");
    addBox(s, "overview-bottom", { left: 56, top: 610, width: 1168, height: 48 }, C.pale, "none", 12);
    addText(s, "overview-bottom-text", "这些数量是冻结 release 总数；图中任何筛选后数字都不能替代正式版本统计。", { left: 80, top: 620, width: 1120, height: 28 }, 19, C.navy, true, "center");
    addNotes(s, "强调审计总行与默认高置信层的区别。10,997 包含 E3 候选与 E0 排除审计记录；网站/论文核心统计应明确筛选层。", [
      "Internal: D:\\finale\\01_正式数据_V6.2\\RELEASE_INFO_v6_2.json",
      "Internal: D:\\finale\\01_正式数据_V6.2\\V62_VALIDATION_REPORT.json",
      "Internal protein master: D:\\finale\\01_正式数据_V6.2\\human_membrane_protein_master_v6_2.tsv",
    ]);
  }

  // 8. Protein landscape
  {
    const s = baseSlide(p, 8, "蛋白库同时区分膜类型、证据强度、功能与结构可用性", "06 · 膜蛋白主表");
    await addFigure(s, `${FIG}/M2_membrane_proteome_landscape.png`, "Human membrane proteome landscape", { left: 48, top: 122, width: 1184, height: 550 });
    addNotes(s, "A/B/C 是膜结合方式；E1–E3 是确定性。A=5,608，B=442，C=4,506，另有 441 条 E0/unknown 审计记录。E1=3,451，E2=4,453，E3=2,652。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M2_membrane_proteome_landscape.png",
      "Internal protein master: D:\\finale\\01_正式数据_V6.2\\human_membrane_protein_master_v6_2.tsv",
    ]);
  }

  // 9. Expression
  {
    const s = baseSlide(p, 9, "表达与定位回答“在哪里、在哪类细胞中存在”", "07 · 表达与定位");
    await addFigure(s, `${FIG}/M3_expression_localization_atlas.png`, "Tissue, cell-type and subcellular localization atlas", { left: 48, top: 122, width: 1184, height: 550 });
    addNotes(s, "组织 nTPM、细胞类型 nCPM、IHC 染色和亚细胞定位保留原始量纲；缺失值与测得为零严格区分。主表只合并一行一蛋白摘要。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M3_expression_localization_atlas.png",
      "Internal: D:\\finale\\01_正式数据_V6.2\\protein_tissue_expression_v2.tsv.gz",
      "Internal: D:\\finale\\01_正式数据_V6.2\\protein_cell_type_expression_v2.tsv.gz",
      "HPA downloads: https://www.proteinatlas.org/about/download",
    ]);
  }

  // 10. Disease
  {
    const s = baseSlide(p, 10, "疾病关系保留来源、证据渠道与本体粒度", "08 · 蛋白—基因—疾病");
    await addFigure(s, `${FIG}/M4_disease_association_landscape.png`, "Disease association landscape", { left: 48, top: 122, width: 1184, height: 550 });
    addNotes(s, "疾病系统图是项目标准化映射，用于统计展示；原始来源 ID 仍保留。人类遗传、临床遗传、体细胞突变、功能实验、动物模型、表达和文献等渠道主要来自 Open Targets，UniProt 疾病注释单独保留。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M4_disease_association_landscape.png",
      "Internal: D:\\finale\\01_正式数据_V6.2\\protein_gene_disease_relations_v6_2.tsv",
      "Open Targets Platform: https://platform.opentargets.org/",
      "UniProt human disease annotation overview: https://www.uniprot.org/release-notes/homo_sapiens",
    ]);
  }

  // 11. Compound identity
  {
    const s = baseSlide(p, 11, "canonical母体与form双层身份避免盐型重复计数", "09 · 小分子主表");
    await addFigure(s, `${FIG}/M5_chemical_identity_space.png`, "Chemical identity and small-molecule space", { left: 48, top: 122, width: 1184, height: 550 });
    addNotes(s, "完整 InChIKey 是结构身份键：前14位反映连接关系，后续部分反映立体化学、同位素和质子化等层面。只有结构身份唯一才自动合并；名称相同不足以合并。ECDF 展示物化性质分布，Morgan 指纹嵌入展示结构相似性空间。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M5_chemical_identity_space.png",
      "Internal: D:\\finale\\01_正式数据_V6.2\\small_molecule_master_v1_3.tsv",
      "Internal: D:\\finale\\01_正式数据_V6.2\\compound_form_hierarchy_v1_3.tsv",
      "InChI documentation: https://www.inchi-trust.org/",
    ]);
  }

  // 12. Binding evidence
  {
    const s = baseSlide(p, 12, "结合证据同时保存来源、定量强度与残基级位点", "10 · 结合与位点证据");
    await addFigure(s, `${FIG}/M6_binding_evidence_sites.png`, "Binding evidence and binding-site landscape", { left: 48, top: 122, width: 1184, height: 550 });
    addNotes(s, "BE1通常代表结构或直接结合层面的强证据；BE2代表可量化直接作用；BE3包含药理或功能支持。独立来源数当前表示不同数据库标签数，不等同于完成原始实验谱系去重。位点覆盖只针对拥有结构/残基信息的关系，不覆盖所有蛋白—小分子对。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M6_binding_evidence_sites.png",
      "Internal: D:\\finale\\01_正式数据_V6.2\\binding_evidence_master_v6_2.tsv.gz",
      "Internal: D:\\finale\\01_正式数据_V6.2\\binding_site_instances_v6_2.tsv.gz",
      "PDBe: https://www.ebi.ac.uk/pdbe/",
    ]);
  }

  // 13. Negative and conflict
  {
    const s = baseSlide(p, 13, "负证据和正负冲突被显式保留，而不是被静默覆盖", "11 · 冲突与可靠性");
    await addFigure(s, `${FIG}/M7_negative_conflicts_quality.png`, "Negative evidence, conflicts and data reliability", { left: 48, top: 122, width: 1184, height: 550 });
    addNotes(s, "阳性与阴性并非简单互斥：不同浓度、体系、构象、终点和 assay 都可能造成状态依赖。正式负证据必须先有可靠的小分子外键；未解析记录继续留在 review 表。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M7_negative_conflicts_quality.png",
      "Internal: D:\\finale\\01_正式数据_V6.2\\negative_binding_evidence_v1_2.tsv.gz",
      "Internal: D:\\finale\\01_正式数据_V6.2\\negative_binding_evidence_unmapped_review_v1_2.tsv.gz",
      "Internal: D:\\finale\\01_正式数据_V6.2\\METHODS_v6_2.md",
    ]);
  }

  // 14. QA interpretation
  {
    const s = baseSlide(p, 14, "质量门控保证可计算性，同时保留不确定性", "12 · 质量控制");
    addText(s, "qa-lead", "V6.2的“PASS”意味着结构完整、规则可重放、审计层可追踪；它不意味着所有生物学问题已经解决。", { left: 56, top: 152, width: 1120, height: 70 }, 30, C.navy, true);
    const qa = [
      ["完整性", "主键唯一；外键有效；文件数量、行数与来源统计对账。", C.blue],
      ["可追溯性", "来源、版本、原始 ID、转换规则与稳定关系 ID 被保留。", C.teal],
      ["不确定性", "E3、冲突、未映射、复杂复合物和授权边界进入独立复核层。", C.coral],
    ];
    qa.forEach((v, i) => addFlatColumn(s, 70 + i * 400, 282, 330, v[0], v[1], v[2], i + 10));
    addBox(s, "qa-gate", { left: 258, top: 570, width: 764, height: 72 }, C.navy, "none", 18);
    addText(s, "qa-gate-text", "Preflight PASS  +  Blocking errors = 0  +  SHA256 manifest  →  Frozen release", { left: 290, top: 590, width: 700, height: 34 }, 24, C.white, true, "center");
    addNotes(s, "解释 QA 的边界：数据工程正确性和生物学真值不是同一件事。审核层的存在是质量设计，而不是失败记录。", [
      "Internal: D:\\finale\\01_正式数据_V6.2\\V62_PREFLIGHT_QA_GATE.json",
      "Internal: D:\\finale\\01_正式数据_V6.2\\V62_VALIDATION_REPORT.json",
      "Internal: D:\\finale\\01_正式数据_V6.2\\SHA256SUMS_v6_2.txt",
    ]);
  }

  // 15. Significance
  {
    const s = baseSlide(p, 15, "MemPro把分散证据转化为可检索、可比较的研究底座", "13 · 科学意义");
    const sig = [
      ["发现", "从蛋白、疾病或小分子任一入口反向查询证据网络，发现被单一数据库遗漏的连接。", C.blue],
      ["解释", "把表达、亚细胞定位、组装状态和结合位点放在同一身份框架中解释机制。", C.teal],
      ["优先级", "按证据强度、来源、定量值、化学可操作性与冲突状态筛选实验和 docking 候选。", C.coral],
      ["复现", "固定版本、哈希、QA 和审计队列，使论文、网站与后续 release 可以对账。", C.purple],
    ];
    sig.forEach((v, i) => {
      const y = 160 + i * 118;
      addText(s, `sig-num-${i}`, `0${i + 1}`, { left: 70, top: y, width: 72, height: 58 }, 42, v[2], true);
      addText(s, `sig-head-${i}`, v[0], { left: 166, top: y + 2, width: 180, height: 38 }, 30, C.navy, true);
      addText(s, `sig-body-${i}`, v[1], { left: 360, top: y, width: 820, height: 68 }, 23, C.ink, false, "left", "middle");
      if (i < 3) addRule(s, 166, y + 88, 1010);
    });
    addNotes(s, "科学价值不是单纯数量大，而是把多来源信息转成可比关系，并保留证据边界。对实验组，它是候选优先级工具；对计算组，它是可重放输入层。", [
      "Internal: D:\\finale\\00_交接说明\\01_项目与数据交接说明.md",
      "Li H et al. Nat Biotechnol 2024. https://doi.org/10.1038/s41587-023-01987-2",
    ]);
  }

  // 16. Limitations
  {
    const s = baseSlide(p, 16, "当前不足集中在功能、疾病、isoform与来源独立性", "14 · 局限性");
    const left = [
      ["功能分类仍不够成熟", "other_family_defined 2,526；unclassified 1,426。单一“家族”维度不足以表达功能。"],
      ["疾病本体尚未完全统一", "MONDO、OMIM、Orphanet、EFO 的粒度和同义关系仍需系统 crosswalk。"],
      ["isoform 层缺失", "关系层以 canonical UniProt accession 为主，尚未形成完整 isoform 级关系。"],
    ];
    const right = [
      ["来源独立性仍需追溯", "distinct source count 是数据库标签数；转载同一实验不能视作完全独立证据。"],
      ["结构与亚基状态依赖", "PDB 非对称单元不等于生理组装；构象、膜环境和激活状态可能改变寡聚体。"],
      ["授权限制影响再发布", "BRENDA、PDBbind 等来源在公开下载、API 或网站上线前必须逐项复核许可。"],
    ];
    [left, right].forEach((col, ci) => {
      col.forEach((v, i) => {
        const x = ci === 0 ? 68 : 664;
        const y = 160 + i * 158;
        addBox(s, `lim-accent-${ci}-${i}`, { left: x, top: y, width: 10, height: 110 }, ci === 0 ? C.coral : C.purple, "none", 4);
        addText(s, `lim-head-${ci}-${i}`, v[0], { left: x + 28, top: y, width: 500, height: 38 }, 27, C.navy, true);
        addText(s, `lim-body-${ci}-${i}`, v[1], { left: x + 28, top: y + 48, width: 500, height: 78 }, 20, C.ink);
      });
    });
    addNotes(s, "这些不足应在论文和网站中主动披露。特别是功能分类和来源独立性：当前数据完整不等于解释框架已经成熟。", [
      "Internal: D:\\finale\\00_交接说明\\04_已知限制与授权边界.md",
      "Internal protein master: D:\\finale\\01_正式数据_V6.2\\human_membrane_protein_master_v6_2.tsv",
    ]);
  }

  // 17. Roadmap
  {
    const s = baseSlide(p, 17, "下一阶段按成熟度增量补足，不再推倒重建", "15 · 正在补足");
    const phases = [
      { x: 60, w: 350, c: C.blue, h: "现在：注释深化", b: "• GO MF / BP\n• InterPro / Pfam\n• Reactome 通路\n• TCDB、GPCRdb、IUPHAR 正式层级" },
      { x: 465, w: 350, c: C.teal, h: "随后：身份与证据深化", b: "• 疾病本体 crosswalk\n• isoform 与复合物层\n• 原始实验谱系去重\n• 高风险冲突人工复核" },
      { x: 870, w: 350, c: C.coral, h: "发布：服务与应用", b: "• 授权矩阵与公开字段\n• 网站、下载与 API\n• 论文 Methods / Data Records\n• Docking 前处理与实验闭环" },
    ];
    // Direction arrows first.
    addBox(s, "road-arrow-1", { left: 410, top: 365, width: 55, height: 18 }, "#B8CDD5", "none", 8);
    addBox(s, "road-arrow-2", { left: 815, top: 365, width: 55, height: 18 }, "#B8CDD5", "none", 8);
    phases.forEach((v, i) => {
      addText(s, `road-step-${i}`, `0${i + 1}`, { left: v.x, top: 158, width: 80, height: 58 }, 46, v.c, true);
      addText(s, `road-head-${i}`, v.h, { left: v.x, top: 228, width: v.w, height: 58 }, 28, C.navy, true);
      addRule(s, v.x, 302, v.w, v.c, 5);
      addText(s, `road-body-${i}`, v.b, { left: v.x, top: 330, width: v.w, height: 220 }, 23, C.ink);
    });
    addBox(s, "road-bottom", { left: 150, top: 590, width: 980, height: 62 }, C.pale, "none", 14);
    addText(s, "road-bottom-text", "每次补足都生成新版本、差异报告、哈希与 QA；V6.2始终保持只读。", { left: 180, top: 607, width: 920, height: 30 }, 22, C.navy, true, "center");
    addNotes(s, "重点说明不是继续无限收集来源，而是把现有数据的功能本体、疾病 crosswalk、isoform、来源谱系与公开授权做深。", [
      "Internal: D:\\finale\\00_交接说明\\03_复现流程.md",
      "Internal: D:\\finale\\00_交接说明\\04_已知限制与授权边界.md",
    ]);
  }

  // 18. Docking and conclusion
  {
    const s = baseSlide(p, 18, "V6.2已经回答“有什么证据”，下一阶段回答“哪些值得验证”", "16 · 从数据库到验证");
    await addFigure(s, `${FIG}/M8_docking_prioritization.png`, "Structure-aware docking prioritization", { left: 48, top: 128, width: 820, height: 520 });
    addText(s, "dock-head", "三层工作量", { left: 918, top: 164, width: 260, height: 38 }, 30, C.navy, true);
    addMetric(s, 918, 220, 270, "6,019", "Pilot pairs", C.teal, "先验证环境、受体准备和 redocking 协议");
    addMetric(s, 918, 375, 270, "14,713", "Standard pairs", C.blue, "覆盖 1,773 个靶点的正式批次");
    addText(s, "dock-caution", "候选 PDB 只是入口\n仍需确定链、组装体、构象、微状态与 docking box", { left: 918, top: 540, width: 280, height: 90 }, 20, C.coral, true, "center");
    addBox(s, "closing-line", { left: 56, top: 655, width: 1168, height: 4 }, C.teal, "none", 2);
    addNotes(s, "结尾回到项目价值：数据库冻结不是终点，而是使 docking、实验验证、网站检索和后续版本更新拥有一致输入。建议先跑 6,019 对 pilot，再决定 standard/full。", [
      "Asset/internal data: D:\\finale\\02_展示图表_V6.2\\png\\M8_docking_prioritization.png",
      "Internal docking validation: D:\\finale\\03_Docking名单_V6.2\\qa\\DOCKING_HPC_SHORTLIST_V2_V62_VALIDATION.json",
      "Li H et al. Nat Biotechnol 2024. https://doi.org/10.1038/s41587-023-01987-2",
    ]);
  }

  // Export previews and deck.
  for (const [i, slide] of p.slides.items.entries()) {
    const png = await p.export({ slide, format: "png", scale: 1 });
    await fs.writeFile(path.join(TMP, "previews", `slide-${String(i + 1).padStart(2, "0")}.png`), new Uint8Array(await png.arrayBuffer()));
    const layout = await slide.export({ format: "layout" });
    await fs.writeFile(path.join(TMP, "previews", `slide-${String(i + 1).padStart(2, "0")}.layout.json`), await layout.text());
  }
  const montage = await p.export({ format: "webp", montage: true, scale: 1 });
  await fs.writeFile(path.join(TMP, "deck-montage.webp"), new Uint8Array(await montage.arrayBuffer()));
  const pptx = await PresentationFile.exportPptx(p);
  await pptx.save(OUT);
  console.log(JSON.stringify({ output: OUT, slides: p.slides.items.length }));
}

main().catch((err) => {
  console.error(err);
  process.exitCode = 1;
});
