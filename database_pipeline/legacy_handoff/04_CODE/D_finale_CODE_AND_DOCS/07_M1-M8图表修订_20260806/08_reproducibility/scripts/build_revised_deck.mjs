import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const sourcePptx = "C:/tmp/mempro_ppt_revision/source.pptx";
const revisionRoot = "D:/finale/07_M1-M8图表修订_20260806";
const outputPptx = `${revisionRoot}/07_ppt/MemPro_V6.3.1_M1-M8科学图表修订版_20260806.pptx`;
const qaDir = `${revisionRoot}/09_qa/ppt_render`;
const lecturePath = `${revisionRoot}/07_ppt/MemPro_V6.3.1_M1-M8修订版逐页讲稿.txt`;

const replacements = [
  ["im/7upwnm98", "M1_database_architecture_quality_revised.png", "M1 revised database architecture and quality-control figure"],
  ["im/8zutsju5", "M2_membrane_proteome_landscape_revised.png", "M2 revised membrane-proteome landscape with scoped Spearman analysis and AlphaFold coverage"],
  ["im/d4vqd4j2", "M3_expression_localization_atlas_revised.png", "M3 revised HPA expression and localization figure"],
  ["im/h47i1oj2", "M4_disease_association_landscape_revised.png", "M4 revised ontology-derived disease landscape"],
  ["im/dsv6903i", "M5_chemical_identity_space_revised.png", "M5 revised compound identity and PAGTN experiment figure"],
  ["im/0vud8v2d", "M6_binding_evidence_sites_revised.png", "M6 revised binding-evidence and membrane-side site figure"],
  ["im/6l0b25sv", "M7_negative_conflicts_quality_revised.png", "M7 revised negative-evidence and release-QA figure"],
  ["im/ul83e9gj", "M8_docking_prioritization_revised.png", "M8 revised docking prioritization figure"],
];

const commonSources = [
  "Local frozen data: D:/7.22/evidence_expansion_v2_working/releases/release_mempro_v6_2_20260730",
  "V6.3 candidate rules: D:/finale/02_V6.3_candidate_20260803",
  "Revision audit: D:/finale/07_M1-M8图表修订_20260806/01_source_audit/PANEL_SOURCE_FIELD_VERSION_UNIT.tsv",
];

const notes = [
`第1页｜标题
开场：本项目构建的是以人类膜蛋白为核心、连接小分子结合证据、疾病、表达定位、结构和Docking准备度的多实体数据库。今天重点不是逐个展示原始文件，而是说明数据如何被统一、哪些结果可以直接发布、哪些仍停留在复核层。

[Sources]\n- ${commonSources[0]}\n- ${commonSources[1]}`,

`第2页｜阅读方式
本报告按八个模块展开。每个模块先用一页定义统计对象、单位和限制，下一页只放一张放大的六面板科学图。读图时依次回答三件事：统计的实体是什么；证据来自哪里；图可以支持什么结论、不能支持什么结论。

[Sources]\n- ${commonSources[2]}`,

`第3页｜M1前置：数据库架构
四类来源分别提供蛋白身份与定位、三维结构、蛋白—小分子互作和蛋白—疾病关系。ID mapping不是按名称合并，而是把UniProt、CID、InChIKey、MONDO等标识指向明确实体。化合物同时保留canonical parent和exact form；蛋白区分canonical、isoform、complex和unresolved target。阳性、阴性、冲突和复核记录分层保存，只有通过主外键、数量、来源与哈希校验的记录进入冻结发布。

[Sources]\n- UniProtKB 2026_02; HPA 25.1; UniTmp/HTP d.2.2\n- PDBe, OPM, PDBTM, PDBbind 2020.R1 reprocessed\n- ChEMBL 37, BindingDB 2026-07, PubChem BioAssay, BRENDA 2026.1\n- Open Targets 26.06, MONDO, Disease Ontology, Uberon\n- ${commonSources[2]}`,

`第4页｜M1主图：从来源到冻结发布
A：四类数据源不是简单拼接，而是各自承担不同信息角色。B：整合顺序是原始记录、ID映射、parent/form、跨来源去重、证据分级、发布。C：阳性与阴性分表，冲突表示同一canonical pair在不同条件下同时存在阳性和阴性，并不自动判定哪一个错误。D：实体模型避免把同基因的isoform或多亚基复合物强塞成单蛋白。E：来源—模块矩阵显示每个数据库实际贡献的模块。F：冻结规模用于对账，不等于所有记录均为高置信默认层。

[Sources]\n- ${commonSources[0]}\n- ${commonSources[2]}\n- Database counts: 01_source_audit/DATABASE_EXACT_CONTRIBUTION_COUNTS.tsv`,

`第5页｜M2前置：膜蛋白分类与证据
A/B/C描述膜结合方式：A为跨膜或强结构支持的整合膜蛋白；B为单层嵌入、膜内区或脂质锚定；C为稳定外周膜相关。E1/E2/E3描述我们确认其膜蛋白身份的把握，不是功能等级。结构、人工拓扑和多模态一致证据优先；单一预测或by similarity通常更弱。功能采用结构家族、分子功能、生物过程、膜角色、专科层级五轴交叉分类。长度—跨膜次数相关性只在A类、E1/E2、canonical multipass、2–40个跨膜区、50–5000 aa范围内计算，并报告bootstrap区间，避免极端值主导。

[Sources]\n- UniProtKB 2026_02\n- HPA 25.1; UniTmp/HTP d.2.2; Membranome; OPM; PDBTM\n- 02_data/M2D_spearman_main_and_bootstrap.tsv`,

`第6页｜M2主图：膜蛋白景观
A：ABC与E1–E3是两个正交维度，柱高为唯一canonical蛋白数。B：五轴交叉分类减少“other”掩盖信息的问题。C：展示单跨膜、多跨膜和非跨膜膜结合架构。D：在预注册范围内，n=2,844，Spearman rho约0.445，bootstrap 95% CI约0.410–0.481；这是中等单调相关，不代表蛋白越长就必然有更多跨膜螺旋。E：实验结构覆盖与AlphaFoldDB标识并列；AlphaFold是预测覆盖，不能替代实验结构。F：组装体证据区分单体、二聚体等；状态依赖或研究条件差异应标冲突/多状态，而非强制唯一。

[Sources]\n- UniProtKB/HPA/UniTmp/Membranome/OPM/PDBTM\n- AlphaFold DB API, model version 6 metadata snapshot\n- 02_data/M2D_spearman_main_and_bootstrap.tsv\n- 05_alphafold/ALPHAFOLD_CONFIDENCE_QA.json`,

`第7页｜M3前置：HPA表达与定位
RNA与IHC不是同一测量层。组织RNA使用HPA consensus tissue RNA的nTPM；细胞类型RNA使用single-cell type RNA的nCPM；IHC是抗体免疫组织化学蛋白染色，反映组织切片中蛋白信号等级。log1p用于压缩长尾，row-wise z-score只表示同一蛋白在不同组织/细胞之间的相对高低。HPA mapped与unmapped/missing必须分开，表达广度为0不能直接解释为完全不表达。

[Sources]\n- Human Protein Atlas 25.1\n- V6.2 protein master expression/location fields\n- HPA mapping audit in M3 outputs`,

`第8页｜M3主图：表达与亚细胞定位
A：先报告HPA映射状态，决定后续0值能否解释。B：并列展示RNA与IHC字段的可用性，不把二者混成一个“表达量”。C：组织RNA热图是nTPM经log1p和蛋白内z-score后的相对模式。D：细胞类型热图同理，但单位为nCPM，来源是HPA单细胞类型层。E：亚细胞定位是多标签，一个蛋白可同时位于质膜、内质网或囊泡等。F：表达广度把mapped-zero、mapped-positive与missing/unmapped分开；这是覆盖率与广度图，不是因果或绝对丰度比较。

[Sources]\n- Human Protein Atlas 25.1\n- ${commonSources[0]}`,

`第9页｜M4前置：疾病身份与解剖映射
蛋白—疾病统计单位是唯一canonical protein—canonical disease pair。Open Targets和UniProtKB保留各自source disease ID；只有MONDO官方exact/equivalent映射才自动统一，broad/narrow/related只保留关系。治疗领域采用Open Targets多标签分类。器官与组织不再使用关键词硬分，而是通过MONDO/DO疾病本体和Uberon解剖实体建立多标签映射，因此同一疾病可属于多个系统。

[Sources]\n- Open Targets 26.06\n- UniProtKB 2026_02 disease annotations\n- MONDO, Disease Ontology, Uberon frozen mappings\n- 02_data/M4F_ontology_disease_anatomy_matrix.tsv`,

`第10页｜M4主图：疾病关联景观
A：Open Targets与UniProtKB贡献的是关系记录，来源数不等于独立实验数。B：按治疗领域比较very high/high/medium证据pair；等级来自上游证据渠道和整合评分，而非我们凭疾病名称判断。C：遗传、临床遗传、体细胞突变、功能实验、动物模型、表达和文献是不同渠道，不能互换。D：显示exact/equivalent身份统一及未安全映射的source-only实体。E：疾病系统为ontology-derived multi-label，不再把未命中关键词的疾病堆入“其他”。F：tissue location来自疾病—Uberon明确链接；空白表示没有安全映射，不等于疾病与该组织无关。

[Sources]\n- Open Targets 26.06; UniProtKB 2026_02\n- MONDO/DO/Uberon mappings\n- 02_data/M4F_ontology_disease_anatomy_matrix.tsv`,

`第11页｜M5前置：小分子身份与结构表示
canonical parent用于跨数据库归一，exact form保留盐型、立体化学、同位素、质子化和电荷。CID是PubChem记录ID，InChIKey是从标准化结构计算的固定长度结构键，两者功能不同；SMILES用于表达结构但不适合作为唯一跨库主键。M5C的类别由一套确定、互斥、可复算的描述符规则生成。PAGTN是路径增强图Transformer编码器，不是聚类算法；本实验用四项物化描述符监督训练编码器，再用HDBSCAN聚类，并与Morgan radius=2、1024-bit基线比较。

[Sources]\n- MemPro compound master v1.3\n- RDKit descriptors and Bemis–Murcko scaffolds\n- PAGTN paper: https://arxiv.org/abs/1905.12712\n- DGL-LifeSci PAGTN implementation`,

`第12页｜M5主图：化合物身份、结构类别与PAGTN实验
A：身份置信度说明canonical合并是否有唯一结构支持。B：Approved、Clinical、Endogenous、Natural、Probe可重叠，不能把单根柱误读为互斥总量。C：同一套规则把全部core QC化合物分入无碳、肽样、大环、多环、单环或无环等互斥结构类；这是结构描述，不是药效评级。D：10/50/90分位热图概括MW、XlogP、TPSA和可旋转键，用于Docking准备难度，不是数据库收录门槛。E：PAGTN二维图来自固定种子样本；颜色为HDBSCAN簇，灰色为噪声。二维坐标本身没有化学单位，聚类质量需结合与Morgan的silhouette、局部类别纯度和骨架纯度评价。F：Bemis–Murcko骨架频率直接回答骨架多样性，补充学习表示的可解释性。

[Sources]\n- ${commonSources[0]}/small_molecule_master_v1_3.tsv\n- 02_data/M5C_recomputed_structural_classes.tsv\n- 02_data/M5E_pagtn_cluster_summary.tsv\n- 04_pagtn/PAGTN_EXPERIMENT_QA.json`,

`第13页｜M6前置：结合证据与位点
BE等级描述证据类型，不是效力强弱。BE1为结构/位点直接证据；BE2为明确映射至具体人类蛋白的直接或竞争性定量结合；BE3为蛋白特异功能或药理结果。Kd是平衡解离常数，通常直接描述结合；Ki只有在竞争/结合模型明确时可作为结合参数；IC50和EC50受assay系统、底物、表达量与信号放大影响，只有实验上下文证明是binding assay时才可升为BE2。位点需要PDB、chain、ligand、接触残基和UniProt坐标映射。

[Sources]\n- ChEMBL 37; BindingDB 2026-07; PubChem BioAssay; BRENDA 2026.1\n- PDBe/SIFTS; PDBbind; BioLiP; UniProtKB sites\n- MemPro binding-evidence policy V6.2`,

`第14页｜M6主图：来源、证据等级和膜侧位点
A：不同来源贡献不同证据层，数据库行数不能等同于独立实验数。B：BE1/BE2/BE3按上下文分层，同一数值单位不会消除测量含义差异。C：位点可用性从残基接触、结构匹配、实验口袋、仅复合物PDB到UniProt人工位点逐级不同。D：9,036个可映射残基接触site按蛋白拓扑分为膜内、膜界面、胞质侧、非胞质侧、跨两侧或未解析；“非胞质侧”可能是细胞外或细胞器腔面。E：接触残基数是跨许多site的氨基酸接触累计，不是单个蛋白拥有几十万残基。F：能直接生成Docking box还需要链、配体、位点和结构上下文全部可用；AlphaFold局部pLDDT只表示模型局部置信度，不证明口袋真实。

[Sources]\n- PDBe/SIFTS + UniProt/UniTmp topology + OPM/PDBTM eligibility\n- 06_binding_site_side/binding_site_membrane_side_summary.tsv\n- 05_alphafold/alphafold_binding_site_local_plddt.tsv`,

`第15页｜M7前置：负证据、冲突与复核队列
负证据表示在特定assay条件、浓度和终点下未观察到活性，不是“永远不结合”。若同一canonical protein—compound pair同时有阳性和阴性，不互相覆盖，而是比较assay、浓度、构建体、终点、物种和文献。review queue保存身份或证据不足但尚不能安全删除的原始记录；它不是默认高置信数据库，也不是垃圾箱。冻结发布前必须通过主键、外键、来源数量、哈希、manifest和抽样验证。

[Sources]\n- PubChem BioAssay mapped and review layers\n- V6.2 validation report and release manifest\n- ${commonSources[0]}`,

`第16页｜M7主图：负证据、冲突和质量控制
A：负证据按已映射正式层、未映射复核层和精确重复分开。B：68,015个正负冲突pair只表示上下文不一致，不能自动判哪方错误。C：阳性身份、阴性身份、亚基状态和表达定位复核队列分别对应不同问题。D：缺失率描述字段可用性，不直接评价来源质量。E：QA gates中主键、外键、数量对账和hash/manifest均为通过/未通过门禁；blocking errors为0才允许冻结。F：默认高置信、阳性复核、映射负证据和阴性复核是并列发布层，用户应根据用途选择。

[Sources]\n- V6.2 frozen validation report\n- PubChem BioAssay negative-evidence reconciliation\n- ${commonSources[0]}`,

`第17页｜M8前置：Docking筛选与HPC准备
Docking优先级综合蛋白证据、互作证据、结构可用性、pair-specific位点、化合物身份、冲突与计算复杂度，不使用Lipinski作为数据库删除规则。候选应区分已有复合物可redocking、同靶点可cross-docking、预测结构探索和不可直接Docking。HPC输入需准备受体结构、质子化、缺失残基/金属/辅因子策略、配体3D与质子化/互变异构体、box坐标、参数、任务清单和可重启日志。

[Sources]\n- V6.3.1 docking readiness manifest\n- V6.3 candidate docking priority rules`,

`第18页｜M8主图：Docking优先级与工作量
A：按R0–D3或当前docking tier展示候选数量，层级越高表示越适合先做，并非生物学价值绝对更高。B：结构/身份准备度分别检查候选受体PDB、pair-specific site、SMILES和InChIKey。C：唯一蛋白、化合物和pair用于估算任务规模。D：site grade决定box能否自动生成；unresolved必须人工补位点或降低优先级。E：ranking score分位数用于批次切分，分数只在同一规则版本内可比。F：Pilot pairs、receptor jobs和conflict flags帮助设计HPC array；正式运行前仍需受体/配体前处理、参数验证和小规模redocking基准。

[Sources]\n- D:/finale/03_V6.3.1_readiness_20260803/04_docking_hpc/docking_pilot_manifest_v631.tsv\n- ${commonSources[1]}/07_docking/docking_priority_v6_3_candidate.tsv.gz`,
];

function setText(id, value) {
  const shape = presentation.resolve(id);
  if (!shape?.textFrame) throw new Error(`Text shape not found: ${id}`);
  shape.textFrame.setText(value);
}

await fs.mkdir(path.dirname(outputPptx), { recursive: true });
await fs.rm(qaDir, { recursive: true, force: true });
await fs.mkdir(path.join(qaDir, "slides"), { recursive: true });
await fs.mkdir(path.join(qaDir, "layouts"), { recursive: true });

const presentation = await PresentationFile.importPptx(await FileBlob.load(sourcePptx));

for (const [id, filename, alt] of replacements) {
  const image = presentation.resolve(id);
  if (!image?.replace) throw new Error(`Image anchor not found: ${id}`);
  const oldFrame = image.frame;
  const oldCrop = image.crop;
  const oldFit = image.fit;
  const oldGeometry = image.geometry;
  const oldBorderRadius = image.borderRadius;
  const oldRotation = image.rotation;
  const oldFlipHorizontal = image.flipHorizontal;
  const oldFlipVertical = image.flipVertical;
  const oldLockAspectRatio = image.lockAspectRatio;
  const bytes = await fs.readFile(`${revisionRoot}/03_figures/png/${filename}`);
  image.replace({ blob: bytes, contentType: "image/png", alt, fit: "contain" });
  image.frame = oldFrame;
  image.crop = oldCrop;
  image.fit = "contain";
  image.geometry = oldGeometry;
  image.borderRadius = oldBorderRadius;
  image.rotation = oldRotation;
  image.flipHorizontal = oldFlipHorizontal;
  image.flipVertical = oldFlipVertical;
  image.lockAspectRatio = oldLockAspectRatio;
}

// Only the M5 explanatory page had scientifically outdated wording after the approved revision.
setText("sh/xcryxg7y", "前置知识：parent/form、结构分类、描述符与PAGTN");
setText("sh/epobatgr", "描述符分位数");
setText("sh/zaxs3yhc", "MW、XlogP、TPSA和可旋转键用10/50/90分位概括，用于Docking准备分层，不作为数据库收录门槛。");
setText("sh/dofa18z6", "PAGTN编码与HDBSCAN聚类");
setText("sh/tsnip0ny", "PAGTN学习路径增强图表示，HDBSCAN负责聚类；固定种子抽样并与Morgan指纹比较，二维坐标没有直接化学单位。");

for (let i = 0; i < presentation.slides.items.length; i += 1) {
  const slide = presentation.slides.items[i];
  slide.speakerNotes.textFrame.setText(notes[i]);
  slide.speakerNotes.setVisible(true);
}

const pptx = await PresentationFile.exportPptx(presentation);
await fs.writeFile(outputPptx, new Uint8Array(await pptx.arrayBuffer()));
for (let i = 0; i < presentation.slides.items.length; i += 1) {
  const n = String(i + 1).padStart(2, "0");
  const slide = presentation.slides.items[i];
  const png = await presentation.export({ slide, format: "png", scale: 1 });
  const layout = await presentation.export({ slide, format: "layout" });
  await fs.writeFile(path.join(qaDir, "slides", `slide-${n}.png`), new Uint8Array(await png.arrayBuffer()));
  await fs.writeFile(path.join(qaDir, "layouts", `slide-${n}.layout.json`), new Uint8Array(await layout.arrayBuffer()));
}
const montage = await presentation.export({ format: "webp", montage: true, scale: 0.45 });
await fs.writeFile(path.join(qaDir, "montage.webp"), new Uint8Array(await montage.arrayBuffer()));
await fs.writeFile(lecturePath, notes.join("\n\n" + "=".repeat(88) + "\n\n"), "utf8");
await fs.writeFile(path.join(qaDir, "ppt_build_manifest.json"), JSON.stringify({ sourcePptx, outputPptx, slideCount: presentation.slides.items.length, replacements, generatedAt: new Date().toISOString() }, null, 2), "utf8");
console.log(JSON.stringify({ status: "PASS", outputPptx, lecturePath, slides: presentation.slides.items.length }));
