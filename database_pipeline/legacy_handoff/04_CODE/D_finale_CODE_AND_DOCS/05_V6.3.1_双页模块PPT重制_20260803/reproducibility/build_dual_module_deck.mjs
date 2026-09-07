import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const INPUT = "starter_dual_module.pptx";
const OUTPUT = "MemPro_V6.3.1_数据库建立流程与科学意义_双页模块重制版.pptx";
const OUTDIR = "dual-module-final";
const FIGROOT = "D:/finale/05_V6.3.1_双页模块PPT重制_20260803/figures/png";

const COLORS = {
  navy: "#12324A",
  teal: "#0F8B8D",
  tealDark: "#0B6B6D",
  tealLight: "#DDF2F0",
  coral: "#E66B5B",
  coralLight: "#FBE8E4",
  blueLight: "#E8F1F7",
  amber: "#D99B2B",
  amberLight: "#FFF3D8",
  green: "#3B8C6E",
  greenLight: "#E5F3EC",
  gray: "#556575",
  pale: "#F8FBFA",
  white: "#FFFFFF",
  line: "#C8D7DC",
};

const deck = await PresentationFile.importPptx(await FileBlob.load(INPUT));
const inspected = await deck.inspect({
  kind: "slide,textbox,shape,image,notes",
  maxChars: 1000000,
});
const rows = inspected.ndjson.trim().split(/\n/).map((s) => JSON.parse(s));

function find(slide, name = null, kind = null) {
  const matches = rows.filter((r) => r.slide === slide && (!name || r.name === name) && (!kind || r.kind === kind));
  if (matches.length !== 1) throw new Error(`Target mismatch slide=${slide} name=${name} kind=${kind}: ${matches.length}`);
  return deck.resolve(matches[0].id);
}

function setText(slide, name, value, options = {}) {
  const shape = find(slide, name);
  shape.text = value;
  if (options.fontSize) shape.text.fontSize = options.fontSize;
  if (options.color) shape.text.color = options.color;
  if (options.bold !== undefined) shape.text.bold = options.bold;
  if (options.align) shape.text.alignment = options.align;
  if (options.valign) shape.text.verticalAlignment = options.valign;
  if (options.typeface) shape.text.typeface = options.typeface;
  return shape;
}

function setNotes(slide, text) {
  find(slide, null, "notes").setText(text);
}

function addShape(slide, { x, y, w, h, geometry = "rect", fill = COLORS.white, line = COLORS.line, lineWidth = 1, name }) {
  return slide.shapes.add({
    geometry,
    name,
    position: { left: x, top: y, width: w, height: h },
    fill: { type: "solid", color: fill },
    line: { style: "solid", fill: line, width: lineWidth },
  });
}

function addText(slide, text, { x, y, w, h, size = 20, color = COLORS.navy, bold = false, align = "left", valign = "top", fill = "#00000000", line = "#00000000", name }) {
  const box = addShape(slide, { x, y, w, h, geometry: "rect", fill, line, lineWidth: 0, name });
  box.text = text;
  box.text.fontSize = size;
  box.text.color = color;
  box.text.bold = bold;
  box.text.typeface = "Microsoft YaHei";
  box.text.alignment = align;
  box.text.verticalAlignment = valign;
  box.text.insets = { left: 6, right: 6, top: 4, bottom: 4 };
  return box;
}

function coverBody(slide, y = 128, h = 558) {
  return addShape(slide, { x: 24, y, w: 1232, h, geometry: "rect", fill: COLORS.pale, line: COLORS.pale, lineWidth: 0, name: `body-cover-${slide.id}` });
}

async function addFigure(slideNumber, path, alt) {
  const slide = deck.slides.items[slideNumber - 1];
  coverBody(slide, 126, 560);
  const image = slide.images.add({ blob: await fs.readFile(path), fit: "contain", alt, name: `main-figure-${slideNumber}` });
  image.position = { left: 36, top: 134, width: 1208, height: 536 };
}

function setupHeader(slide, kicker, title, page) {
  setText(slide, "section-kicker", kicker, { fontSize: 17, color: COLORS.tealDark, bold: true });
  setText(slide, "slide-title", title, { fontSize: 42, color: COLORS.navy, bold: true });
  setText(slide, "footer-left", "MemPro V6.3.1 · 人类膜蛋白—小分子互作数据库", { fontSize: 11, color: COLORS.gray });
  setText(slide, "footer-page", String(page).padStart(2, "0"), { fontSize: 11, color: COLORS.gray });
}

function fillSixBoxes(slideNumber, entries) {
  const slide = deck.slides.items[slideNumber - 1];
  const coords = [
    ["lim-head-0-0", "lim-body-0-0"],
    ["lim-head-0-1", "lim-body-0-1"],
    ["lim-head-0-2", "lim-body-0-2"],
    ["lim-head-1-0", "lim-body-1-0"],
    ["lim-head-1-1", "lim-body-1-1"],
    ["lim-head-1-2", "lim-body-1-2"],
  ];
  entries.forEach((entry, i) => {
    const [head, body] = coords[i];
    setText(slideNumber, head, entry[0], { fontSize: 21, color: COLORS.navy, bold: true });
    setText(slideNumber, body, entry[1], { fontSize: 15.5, color: COLORS.gray, bold: false });
  });
  return slide;
}

// 1 — title
setText(1, "deck-title", "MemPro V6.3.1", { fontSize: 52, color: COLORS.navy, bold: true });
setText(1, "deck-meta", "人类膜蛋白—小分子互作数据库 · 双页模块汇报版 · 2026-08-03", { fontSize: 20, color: COLORS.gray });
setText(1, "deck-tagline", "先讲清规则，再读一张大图：从来源、身份与证据到疾病、表达和 Docking", { fontSize: 24, color: COLORS.tealDark, bold: true });
setNotes(1, `本稿采用“双页模块制”：每个模块第一页解释术语、口径和读图方法，第二页只展示一张放大的主图。\n\n[Sources]\n- Internal: D:\\finale\\02_V6.3_candidate_20260803\\FREEZE_V6_3_CANDIDATE.json\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\05_qa\\V631_READINESS_VALIDATION.json`);

// 2 — background
setupHeader(2, "01 · 研究背景", "为什么需要一个可追溯的人类膜蛋白—小分子数据库", 2);
setText(2, "background-lead", "膜蛋白连接细胞内外，是药物靶点、疾病机制和结构生物学的交汇处；真正困难的不是“找到条目”，而是让身份、证据和版本彼此可核验。", { fontSize: 25, color: COLORS.navy, bold: true });
setText(2, "col-head-1", "生物学问题", { fontSize: 24, color: COLORS.tealDark, bold: true });
setText(2, "col-body-1", "• 哪些人类蛋白属于膜蛋白？\n• 哪些小分子与它们直接或功能性互作？\n• 互作发生在哪个结构、位点、组织和疾病背景？", { fontSize: 18, color: COLORS.gray });
setText(2, "col-head-2", "数据学难点", { fontSize: 24, color: COLORS.coral, bold: true });
setText(2, "col-body-2", "• 同一蛋白存在基因、canonical、isoform 与复合物层\n• 同一化合物存在母体、盐型、立体异构体和质子化形式\n• 数据库之间会重复转载、单位不一或结论相反", { fontSize: 18, color: COLORS.gray });
setText(2, "col-head-3", "MemPro 的回答", { fontSize: 24, color: COLORS.green, bold: true });
setText(2, "col-body-3", "• 冻结来源版本并保留原始 ID\n• 以结构身份与官方映射建立 canonical 实体\n• 对证据分层、保留冲突和负证据\n• 输出可查询、可复现、可进入 HPC 的发布层", { fontSize: 18, color: COLORS.gray });
setNotes(2, `这一页先建立叙事：膜蛋白数据库的价值不只在条目数量，而在跨来源身份统一、证据可追溯和计算可用性。\n\n[Sources]\n- Internal: D:\\finale\\02_V6.3_candidate_20260803\\README.md\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\README.md`);

// 3 — source prerequisites, rebuilt with a readable table
setupHeader(3, "02 · 数据来源与整合口径", "先明确：每类数据库提供什么，以及版本如何冻结", 3);
{
  const s = deck.slides.items[2];
  coverBody(s, 126, 560);
  const cols = [42, 205, 518, 850];
  const widths = [155, 305, 324, 386];
  ["模块", "冻结来源 / 版本", "主要信息贡献", "整合时保留的关键身份与溯源"].forEach((t, i) => {
    addShape(s, { x: cols[i], y: 145, w: widths[i], h: 43, fill: COLORS.navy, line: COLORS.navy, name: `source-header-${i}` });
    addText(s, t, { x: cols[i], y: 145, w: widths[i], h: 43, size: 17, color: COLORS.white, bold: true, align: i === 0 ? "center" : "left", valign: "middle", name: `source-header-text-${i}` });
  });
  const tableRows = [
    ["蛋白质", "UniProtKB 2026_02；HPA 25.1；HTP/UniTmp d.2.2；Membranome 快照", "序列、拓扑、膜定位、组织/细胞表达、家族与跨库链接", "UniProt accession、序列版本、基因 ID、ABC 与 E1–E3 证据"],
    ["结构", "PDBe 2026-07；OPM / PDBTM 2026-07；PDBbind 2020.R1（2024流程重处理）", "PDB 链、膜内朝向、生物学组装体、配体、接触残基和亲和力", "PDB、chain、assembly、ligand、residue 坐标与 SIFTS 映射"],
    ["互作", "ChEMBL 37；BindingDB 2026-07；PubChem 2026-07；BRENDA 2026.1", "直接结合、功能药理、酶—配体、定量值、阳性与阴性结果", "assay / AID / SID / CID / PubMed、原值、单位、关系符和条件"],
    ["疾病", "Open Targets 26.06；UniProtKB 2026_02；MONDO / DO / Uberon 冻结映射", "蛋白—疾病关系、遗传/临床/功能证据、疾病层级和解剖映射", "source disease ID、canonical MONDO、映射谓词、父子层级与来源"],
  ];
  tableRows.forEach((row, r) => {
    const y = 196 + r * 100;
    const bg = r % 2 === 0 ? COLORS.white : "#EFF6F5";
    row.forEach((t, c) => {
      addShape(s, { x: cols[c], y, w: widths[c], h: 92, fill: bg, line: COLORS.line, lineWidth: 1, name: `source-cell-${r}-${c}` });
      addText(s, t, { x: cols[c] + 4, y: y + 3, w: widths[c] - 8, h: 86, size: c === 0 ? 18 : 14.5, color: c === 0 ? COLORS.tealDark : COLORS.gray, bold: c === 0, align: c === 0 ? "center" : "left", valign: "middle", name: `source-cell-text-${r}-${c}` });
    });
  });
  addText(s, "统一原则：冻结版本/访问日期与许可；保留原始来源 ID；名称相同不能触发自动合并；受限来源的衍生数据单独经过再发布检查。", { x: 48, y: 604, w: 1184, h: 55, size: 17, color: COLORS.coral, bold: true, align: "center", valign: "middle", fill: COLORS.coralLight, line: COLORS.coral, name: "source-rule-banner" });
}
setNotes(3, `该页给出投稿时可直接引用的来源版本口径。版本、访问日期、许可和原始 ID 均进入 manifest；“整合”不代表抹掉来源差异。\n\n[Sources]\n- Internal: D:\\finale\\02_展示图表_V6.2\\data\\M0_source_versions_contributions.tsv\n- UniProtKB release 2026_02\n- Human Protein Atlas 25.1\n- ChEMBL 37\n- Open Targets Platform 26.06`);

// 4 — workflow figure, native shapes only
setupHeader(4, "03 · 数据库建立流程", "从原始来源到四张主表、复核层与可复现发布", 4);
{
  const s = deck.slides.items[3];
  coverBody(s, 126, 560);
  const arrowXs = [242, 482, 722, 962];
  arrowXs.forEach((x, i) => addShape(s, { x, y: 353, w: 44, h: 24, geometry: "rightArrow", fill: i === 3 ? COLORS.coral : COLORS.teal, line: i === 3 ? COLORS.coral : COLORS.teal, lineWidth: 0, name: `flow-arrow-${i}` }));
  const groups = [
    { x: 38, w: 204, title: "① 原始来源", fill: COLORS.blueLight },
    { x: 286, w: 196, title: "② 身份标准化", fill: COLORS.tealLight },
    { x: 526, w: 196, title: "③ Canonical 实体", fill: COLORS.greenLight },
    { x: 766, w: 196, title: "④ 证据与冲突", fill: COLORS.amberLight },
    { x: 1006, w: 236, title: "⑤ 发布与复现", fill: COLORS.coralLight },
  ];
  groups.forEach((g, i) => {
    addShape(s, { x: g.x, y: 150, w: g.w, h: 42, geometry: "roundRect", fill: i === 4 ? COLORS.coral : COLORS.navy, line: i === 4 ? COLORS.coral : COLORS.navy, name: `flow-head-${i}` });
    addText(s, g.title, { x: g.x, y: 150, w: g.w, h: 42, size: 18, color: COLORS.white, bold: true, align: "center", valign: "middle", name: `flow-head-text-${i}` });
    addShape(s, { x: g.x, y: 202, w: g.w, h: 405, geometry: "roundRect", fill: g.fill, line: COLORS.line, name: `flow-group-${i}` });
  });
  const sourceNodes = ["蛋白：UniProt / HPA / HTP", "结构：PDBe / OPM / PDBTM", "互作：ChEMBL / PubChem / BRENDA", "疾病：Open Targets / UniProtKB"];
  sourceNodes.forEach((t, i) => addText(s, t, { x: 50, y: 224 + i * 84, w: 180, h: 63, size: 15, color: COLORS.navy, bold: true, align: "center", valign: "middle", fill: COLORS.white, line: COLORS.line, name: `flow-source-${i}` }));
  const mapNodes = ["ID mapping\nUniProt / CID / InChIKey / MONDO", "Parent / form\n盐型 · 立体化学 · 电荷", "Cross-source dedup\n实验 / 结构 / 文献谱系键"];
  mapNodes.forEach((t, i) => addText(s, t, { x: 300, y: 232 + i * 112, w: 168, h: 82, size: 15, color: COLORS.navy, bold: true, align: "center", valign: "middle", fill: COLORS.white, line: COLORS.teal, name: `flow-map-${i}` }));
  const entityNodes = ["protein", "protein_isoform", "protein_complex", "compound parent / form", "canonical disease"];
  entityNodes.forEach((t, i) => addText(s, t, { x: 540, y: 220 + i * 71, w: 168, h: 52, size: 15, color: COLORS.navy, bold: true, align: "center", valign: "middle", fill: COLORS.white, line: COLORS.green, name: `flow-entity-${i}` }));
  const evidenceNodes = ["膜蛋白证据\nABC × E1/E2/E3", "互作证据\nBE1 / BE2 / BE3 / MECH", "阳性、阴性与冲突\n不互相覆盖", "主键 / 外键 / 来源对账"];
  evidenceNodes.forEach((t, i) => addText(s, t, { x: 780, y: 220 + i * 89, w: 168, h: 68, size: 15, color: COLORS.navy, bold: true, align: "center", valign: "middle", fill: COLORS.white, line: COLORS.amber, name: `flow-evidence-${i}` }));
  addText(s, "四张正式主表", { x: 1022, y: 226, w: 204, h: 50, size: 21, color: COLORS.coral, bold: true, align: "center", valign: "middle", fill: COLORS.white, line: COLORS.coral, name: "flow-release-main" });
  addText(s, "膜蛋白\n蛋白—疾病\n小分子—蛋白结合证据\n小分子及 parent/form", { x: 1022, y: 286, w: 204, h: 150, size: 17, color: COLORS.navy, bold: true, align: "center", valign: "middle", fill: COLORS.white, line: COLORS.line, name: "flow-release-tables" });
  addText(s, "复核队列 · 负证据表\nmanifest · 哈希 · QA\n下载 / 网站 / API / HPC", { x: 1022, y: 452, w: 204, h: 126, size: 16, color: COLORS.gray, bold: false, align: "center", valign: "middle", fill: COLORS.white, line: COLORS.coral, name: "flow-release-qa" });
  addText(s, "每一条发布记录都能回到：来源版本 → 原始记录 → 身份映射 → 证据等级 → 发布决策", { x: 130, y: 620, w: 1020, h: 42, size: 18, color: COLORS.tealDark, bold: true, align: "center", valign: "middle", name: "flow-bottom" });
}
setNotes(4, `流程图对应 integration 的完整含义：不仅把来源拼接起来，还包括身份映射、canonical parent/form、跨来源去重、证据分层、冲突保留、主外键质控和版本冻结。\n\n[Sources]\n- Internal: D:\\finale\\02_V6.3_candidate_20260803\\README.md\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\01_engineering\\requirements-scientific-lock.txt`);

// 5–6 protein
setupHeader(5, "04 · 蛋白质模块｜先讲规则", "ABC 是膜结合方式；E1–E3 是证据把握度：两条轴不能混为一谈", 5);
fillSixBoxes(5, [
  ["A · 严格整合膜蛋白", "具有跨膜区，或有结构/实验支持其穿越脂双层；用于核心统计与网站默认层。"],
  ["B · 直接嵌入膜", "没有经典跨膜螺旋，但存在单层嵌入、膜内区或脂质锚定；与 A 合成广义核心层。"],
  ["C · 外周膜相关", "通常不穿膜或嵌膜，但通过实验、复合物或稳定定位证据与膜持续关联。"],
  ["E1 · 直接强证据", "实验拓扑、脂锚、膜定位或膜结构直接支持。例如具有实验定位与明确跨膜区的受体。"],
  ["E2 · 可靠注释支持", "人工整理的拓扑/定位、可靠结构或多来源一致支持；不是“无实验就不可信”。"],
  ["E3 / E0 · 候选与排除", "E3 多为单一预测、by similarity 或弱定位；E0 为审计排除层。默认确认统计只用 E1+E2。"],
]);
setNotes(5, `ABC回答“蛋白怎样与膜相连”；E级回答“我们有多大把握确认这一点”。E2可以包含数据库人工注释和可靠结构支持，即使当前记录没有直接实验 ECO，也属于确认层。E3保留用于检索与后续复核，不进入默认核心统计。\n\n[Sources]\n- Internal: D:\\finale\\02_V6.3_candidate_20260803\\02_protein\\protein_master_v63_candidate.parquet\n- Internal: membrane evidence policy and V6.3 validation files`);
setupHeader(6, "05 · 蛋白质模块｜主图", "蛋白类别、证据等级与五轴交叉分类", 6);
await addFigure(6, `${FIGROOT}/M2_protein_classification_rebuilt.png`, "MemPro protein ABC and E1-E3 classification with five-axis functional annotation");
setNotes(6, `读图顺序：先看 ABC×E1/E2/E3/E0 的组成，再看膜角色，再看结构家族、分子功能、生物过程、膜角色和专科层级五轴覆盖。五轴交叉分类避免把家族误当成功能。\n\n[Sources]\n- Asset: D:\\finale\\05_V6.3.1_双页模块PPT重制_20260803\\figures\\png\\M2_protein_classification_rebuilt.png\n- Internal: frozen MemPro V6.3 protein master and annotation tables`);

// 7–8 expression
setupHeader(7, "06 · 表达与定位｜先讲规则", "先区分测量体系、零值与缺失，再解释热力图和表达广度", 7);
fillSixBoxes(7, [
  ["nTPM · 组织 RNA", "HPA normalized transcripts per million；主图以 nTPM≥1 作为“该组织检出”的描述阈值。"],
  ["nCPM · 单细胞 RNA", "normalized counts per million；用于细胞类型表达，不能直接与组织 nTPM 当作同一量纲比较。"],
  ["IHC · 蛋白染色", "Immunohistochemistry：抗体在组织切片中的蛋白染色强度和定位；是蛋白层证据，不等同于 RNA。"],
  ["0 ≠ 缺失", "mapped zero 表示已映射且在阈值下；unmapped/missing 表示没有可靠 HPA 映射或数据缺失，必须分开。"],
  ["log1p 转换", "log1p(x)=ln(1+x)：保留 0，同时压缩极高表达长尾；用于让分布与热图更可读。"],
  ["z-score 热力图", "z=(x−行均值)/行标准差。颜色表示同一蛋白类别内的相对高低，不是跨类别绝对表达量。"],
]);
setNotes(7, `表达页最重要的技术边界：RNA与IHC来自不同测量层；“未检出”不能自动解释成“不表达”；行z-score热力图只显示相对组织偏好。\n\n[Sources]\n- Human Protein Atlas 25.1 download documentation\n- Internal: V6.3 expression and localization matrices`);
setupHeader(8, "07 · 表达与定位｜主图", "HPA 映射状态、表达广度、组织偏好与亚细胞定位", 8);
await addFigure(8, `${FIGROOT}/M3_expression_localization_rebuilt.png`, "HPA mapping, tissue breadth, tissue heatmap, and subcellular localization");
setNotes(8, `面板A明确分开mapped detected、mapped zero和unmapped/missing；面板B比较各蛋白类别的组织检出广度；面板C为行z-score组织热图；面板D展示标准化亚细胞定位。\n\n[Sources]\n- Asset: D:\\finale\\05_V6.3.1_双页模块PPT重制_20260803\\figures\\png\\M3_expression_localization_rebuilt.png\n- Human Protein Atlas 25.1`);

// 9–10 disease
setupHeader(9, "08 · 疾病模块｜先讲规则", "疾病先统一身份，再建立层级、治疗领域与器官多标签", 9);
fillSixBoxes(9, [
  ["source disease", "Open Targets 或 UniProtKB 的原始疾病 ID、名称与记录永久保留，作为可追溯入口。"],
  ["canonical MONDO", "仅 MONDO 直接 ID 或官方 unique exact/equivalent 交叉引用可自动合并。"],
  ["不安全映射", "broad、narrow、related、同名无 xref、obsolete 与 1:N 映射只进关系/复核层，不强选 canonical。"],
  ["疾病层级", "疾病亚型不与大类合并身份；通过 parent–child 关系连接，统计可选择具体层级。"],
  ["治疗领域", "采用 Open Targets 26.06 的标准治疗领域，多标签；不再用自建关键词把疾病强塞入单一系统。"],
  ["解剖映射", "DO/MONDO → Uberon，区分 asserted 与 inferred；一个疾病可同时映射多个器官系统。"],
]);
setNotes(9, `该页解释为什么原先“其他/多系统”会堆积：关键词单标签规则不适合疾病。新版采用canonical MONDO、官方层级、Open Targets治疗领域和Uberon解剖多标签。\n\n[Sources]\n- MONDO ontology and official cross-reference mappings\n- Disease Ontology and Uberon\n- Open Targets Platform 26.06`);
setupHeader(10, "09 · 疾病模块｜主图", "疾病身份、本体层级、治疗领域与证据渠道", 10);
await addFigure(10, `${FIGROOT}/M4_disease_ontology_rebuilt.png`, "Disease canonicalization, therapeutic areas, anatomy mapping, and evidence channels");
setNotes(10, `面板A展示exact-only身份统一；B为Open Targets治疗领域多标签；C为DO/MONDO至Uberon的解剖多标签；D展示遗传、临床、体细胞、功能、动物、表达、文献和UniProt注释等证据渠道。\n\n[Sources]\n- Asset: D:\\finale\\05_V6.3.1_双页模块PPT重制_20260803\\figures\\png\\M4_disease_ontology_rebuilt.png\n- Open Targets 26.06; UniProtKB 2026_02; MONDO/DO/Uberon frozen mappings`);

// 11–12 compound
setupHeader(11, "10 · 小分子模块｜先讲规则", "化合物名称只用于检索；真正合并依赖结构身份与 parent/form 层级", 11);
fillSixBoxes(11, [
  ["SMILES", "可读的线性结构字符串；适合重建与计算，但同一结构可有多种写法，需 canonicalize。"],
  ["InChIKey", "由标准化结构生成的固定长度哈希；便于精确匹配和索引，但不能反向恢复完整结构。"],
  ["canonical parent", "跨来源统一的母体结构；用于蛋白—小分子主关系和跨库去重。"],
  ["exact form", "保留盐型、立体异构体、同位素、质子化/互变异构与电荷；有明确证据时关系可指向 form。"],
  ["UpSet 交集", "显示approved、clinical、probe、endogenous、natural等状态的交集；交集数字不是“获批药物总数”。"],
  ["Lipinski 与骨架", "四项规则只作成药性描述，不作为收录/对接硬门槛；Bemis–Murcko骨架衡量全量结构多样性。"],
]);
setNotes(11, `不建议把所有结构“只统一成SMILES”：SMILES是结构表达，InChIKey是身份索引，parent/form是数据库层级，三者承担不同职责。名称相同绝不自动合并。\n\n[Sources]\n- Internal: V6.3 compound parent/form master\n- RDKit 2026.03.5 definitions and calculations`);
setupHeader(12, "11 · 小分子模块｜主图", "生物状态交集、Lipinski描述、骨架覆盖与结构类别", 12);
await addFigure(12, `${FIGROOT}/M5_compound_identity_properties_rebuilt.png`, "Compound status intersections, Lipinski descriptors, scaffold diversity, and structural classes");
setNotes(12, `面板A是真正的状态UpSet；B把Lipinski四项分开统计，并明确仅为描述；C用Bemis–Murcko骨架秩和累计覆盖替代拥挤的Morgan降维图；D展示主要结构类别。\n\n[Sources]\n- Asset: D:\\finale\\05_V6.3.1_双页模块PPT重制_20260803\\figures\\png\\M5_compound_identity_properties_rebuilt.png\n- Internal: 899,591 core-QC compounds; RDKit 2026.03.5`);

// 13–14 evidence
setupHeader(13, "12 · 结合证据｜先讲规则", "证据等级、实验谱系、结构位点和正负冲突分别管理", 13);
fillSixBoxes(13, [
  ["BE1 · 结构/位点直证", "实验观察到的蛋白—小分子复合物，或具有可追溯的残基层结合位点。"],
  ["BE2 · 定量直接结合", "唯一映射到单一蛋白的 Kd、Ki 等直接结合测量；保留关系符、原值、单位与标准值。"],
  ["BE3 / MECH", "BE3 为靶点特异功能药理或酶—配体；MECH 为机制关系，但不冒充定量直接结合。"],
  ["数据库数 ≠ 实验数", "distinct_database_count 仅表示贡献数据库；同一论文/实验被转载不能算多个独立实验。"],
  ["谱系键", "结构键=PDB+chain+ligand；PubChem键=AID+SID+CID；实验/文献键结合 assay、target、compound 与 measurement。"],
  ["阴性与冲突", "inactive/negative 单独保存；与阳性同时出现时标记冲突并保留条件差异，绝不覆盖。"],
]);
setNotes(13, `来源谱系审核耗时的原因是数据库之间会转载同一篇论文、同一PDB或同一BioAssay。新版拆分数据库数、实验数、结构数和证据模态数，避免把“多个数据库”误称为“多个独立实验”。\n\n[Sources]\n- Internal: evidence_policy_v2.md\n- Internal: binding evidence lineage and conflict tables`);
setupHeader(14, "13 · 结合证据｜主图", "来源、BE1–BE3、谱系覆盖、结合位点与负证据", 14);
await addFigure(14, `${FIGROOT}/M6_binding_evidence_lineage_rebuilt.png`, "Binding evidence tiers, provenance lineage, binding-site availability, and negative evidence");
setNotes(14, `面板A按来源与BE等级展示；B报告结构、实验、PubChem与文献谱系覆盖；C仅统计具有PDB/链/配体/残基坐标溯源的位点；D展示负证据映射、复核与正负冲突。\n\n[Sources]\n- Asset: D:\\finale\\05_V6.3.1_双页模块PPT重制_20260803\\figures\\png\\M6_binding_evidence_lineage_rebuilt.png\n- Internal: V6.3 binding evidence, negative evidence, lineage and conflict tables`);

// 15–16 docking
setupHeader(15, "14 · Docking｜先讲规则", "Docking 名单按证据和准备价值分层，不用 Lipinski 作为单一淘汰器", 15);
fillSixBoxes(15, [
  ["小分子基础门槛", "有标准SMILES和InChIKey、核心结构、QC通过；MW 100–700、重原子6–60、|电荷|≤2、可旋转键≤20。"],
  ["蛋白基础门槛", "膜蛋白为ABC且E1/E2，进入网站默认层；必须能确定canonical蛋白或明确isoform/复合物。"],
  ["R0 · Redocking", "已有pair-specific实验PDB结合位点；先用于回对接验证方法，要求参考配体RMSD≤2 Å。"],
  ["D1 / D2 / D3", "D1高优先结构+强证据；D2扩展实验结构/定量证据；D3为已批准、临床、探针或内源性探索层。"],
  ["优先级分数", "基础分+数据库/BE/状态/效力/理化+谱系加分−冲突惩罚；它是准备优先级，不是预测亲和力。"],
  ["生产前四道门", "受体生物学组装体/链/状态；配体质子化/互变/立体；合理box；redocking验证通过。"],
]);
setNotes(15, `小分子筛选回到原有“广口径计算可用性”方案，不把Lipinski口服药规则当作硬门槛。R0先验证流程；D1–D3用于分配HPC资源。\n\n[Sources]\n- Internal: docking pilot eligibility and scoring policy\n- Internal: docking_pilot_manifest_v631.tsv`);
setupHeader(16, "15 · Docking｜主图", "6,019 个 pilot pair 的分层、信息完备度与谱系重排", 16);
await addFigure(16, `${FIGROOT}/M8_docking_readiness_rebuilt.png`, "Docking pilot tiers, readiness coverage, lineage score delta, and score distributions");
setNotes(16, `Pilot共6,019个蛋白—小分子pair、1,773个靶点；主图展示R0/D1/D2/D3分层、结构和位点准备信息、V6.3谱系重排前后变化及分数分布。名单完成不等于已经执行Docking。\n\n[Sources]\n- Asset: D:\\finale\\05_V6.3.1_双页模块PPT重制_20260803\\figures\\png\\M8_docking_readiness_rebuilt.png\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\04_docking_hpc\\docking_pilot_manifest_v631.tsv`);

// 17 — limitations
setupHeader(17, "16 · 当前边界", "哪些已经自动完成，哪些仍需人工、许可或真实计算环境", 17);
fillSixBoxes(17, [
  ["复合物直证原文", "候选已进入complex实体层；只有原文明确支持组分、计量和复合物特异结合时，才升级为正式直证。"],
  ["历史 isoform", "8个旧/删除accession可核查；原始来源未说明isoform的证据只能保留 gene_product_unspecified。"],
  ["受限来源许可", "许可未明确时可内部分析，但不能默认把受限原始/衍生记录公开再发布。"],
  ["第二台设备复现", "一键脚本、相对路径和环境锁已准备；仍需在另一台设备实际运行，验证无隐含本机依赖。"],
  ["受体/配体/box", "Docking名单不是可直接提交的坐标包；还需逐靶点准备生物组装体、微状态和搜索盒。"],
  ["真实 HPC Docking", "必须在目标集群完成调度、资源测试、redocking、失败重试、日志和结果聚合。"],
]);
setNotes(17, `这里区分“数据库自动流水线完成”和“科学/法律/计算执行门槛”。这些开放项不否定数据库主体完成，但不能在论文中写成已经完成。\n\n[Sources]\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\05_qa\\V631_READINESS_VALIDATION.json\n- Internal: review queues and HPC readiness documents`);

// 18 — conclusion
setupHeader(18, "17 · 结论与下一步", "MemPro 已形成可追溯的多实体数据库；下一阶段是人工收口、复现与 HPC 验证", 18);
setText(18, "road-head-0", "数据库主体", { fontSize: 24, color: COLORS.tealDark, bold: true });
setText(18, "road-body-0", "• protein / isoform / complex\n• compound parent / form\n• canonical disease / hierarchy\n• binding evidence / negative evidence\n• 表达、定位、结构与位点", { fontSize: 18, color: COLORS.gray });
setText(18, "road-head-1", "近期收口", { fontSize: 24, color: COLORS.amber, bold: true });
setText(18, "road-body-1", "• 复合物直证与8个旧accession核查\n• 高价值功能未决项人工复核\n• 疾病/isoform边界案例\n• 受限来源再发布许可\n• 第二台设备一键复现实测", { fontSize: 18, color: COLORS.gray });
setText(18, "road-head-2", "HPC 阶段", { fontSize: 24, color: COLORS.coral, bold: true });
setText(18, "road-body-2", "• 先用R0完成redocking验证\n• 分批准备受体、配体和box\n• 按D1→D2→D3提交数组任务\n• 记录失败原因与参数版本\n• 结果回写数据库并形成论文案例", { fontSize: 18, color: COLORS.gray });
setText(18, "road-bottom-text", "核心价值：数量可统计、身份可解释、证据可回溯、冲突不隐藏、流程可复现、计算可扩展。", { fontSize: 22, color: COLORS.navy, bold: true });
setNotes(18, `最终叙事：MemPro不是静态Excel集合，而是以实体身份、证据谱系和发布规则为核心的数据基础设施；Docking是下游验证模块，不应和数据库收集完成度混在一起。\n\n[Sources]\n- Internal: V6.3 candidate release and V6.3.1 readiness package\n- Internal: docking pilot manifest and reproducibility audit`);

// Ensure footer is current on all content slides.
for (let i = 2; i <= 18; i += 1) {
  try { setText(i, "footer-left", "MemPro V6.3.1 · 人类膜蛋白—小分子互作数据库", { fontSize: 11, color: COLORS.gray }); } catch {}
  try { setText(i, "footer-page", String(i).padStart(2, "0"), { fontSize: 11, color: COLORS.gray }); } catch {}
}

await fs.mkdir(`${OUTDIR}/slides`, { recursive: true });
await fs.mkdir(`${OUTDIR}/layouts`, { recursive: true });
for (let i = 0; i < deck.slides.items.length; i += 1) {
  const slide = deck.slides.items[i];
  const n = String(i + 1).padStart(2, "0");
  const png = await deck.export({ slide, format: "png", scale: 1.5 });
  await fs.writeFile(`${OUTDIR}/slides/slide-${n}.png`, new Uint8Array(await png.arrayBuffer()));
  const layout = await slide.export({ format: "layout" });
  await fs.writeFile(`${OUTDIR}/layouts/slide-${n}.layout.json`, await layout.text(), "utf8");
}
const montage = await deck.export({ format: "webp", montage: true, scale: 0.55 });
await fs.writeFile(`${OUTDIR}/montage.webp`, new Uint8Array(await montage.arrayBuffer()));
const finalInspect = await deck.inspect({ kind: "slide,textbox,shape,image,notes", maxChars: 1000000 });
await fs.writeFile(`${OUTDIR}/final-inspect.ndjson`, finalInspect.ndjson, "utf8");
const pptx = await PresentationFile.exportPptx(deck);
await pptx.save(OUTPUT);
console.log(JSON.stringify({ output: OUTPUT, slides: deck.slides.items.length, preview: `${OUTDIR}/montage.webp` }));
