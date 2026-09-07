import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const INPUT = "template-starter.pptx";
const OUTPUT = "MemPro_V6.3.1_数据库建立流程与科学意义_M1-M8原图讲解终版.pptx";
const SCRIPT_OUTPUT = "MemPro_V6.3.1_M1-M8逐图讲稿.txt";

const deck = await PresentationFile.importPptx(await FileBlob.load(INPUT));
const inspected = await deck.inspect({
  kind: "slide,textbox,shape,image,notes",
  include: "id,slide,name,title,text,textPreview,bbox,bboxUnit,alt",
  maxChars: 1000000,
});
const rows = inspected.ndjson.trim().split(/\n/).filter(Boolean).map((line) => JSON.parse(line));

function hits(slide, name) {
  return rows.filter((row) => row.slide === slide && row.name === name);
}

function setText(slide, name, value, required = true) {
  const found = hits(slide, name);
  if (found.length !== 1) {
    if (required) throw new Error(`Expected one object named ${name} on slide ${slide}, found ${found.length}`);
    return false;
  }
  deck.resolve(found[0].id).text = value;
  return true;
}

function setNotes(slide, value) {
  const found = rows.filter((row) => row.slide === slide && row.kind === "notes");
  if (found.length !== 1) throw new Error(`Expected one notes object on slide ${slide}, found ${found.length}`);
  deck.resolve(found[0].id).setText(value);
}

function setBox(slide, col, row, head, body) {
  setText(slide, `lim-head-${col}-${row}`, head);
  setText(slide, `lim-body-${col}-${row}`, body);
}

// Slide 2: explain the repeated two-page teaching structure.
setText(2, "section-kicker", "01 · 汇报导读", false);
setText(2, "slide-title", "每个模块先讲清规则，再用一张原始主图集中展示结果", false);
setText(2, "background-lead", "前置页负责解释口径、指标和限制；后一页只放7月30日冻结的M1–M8原图，并放到最大。", false);
setText(2, "col-head-1", "先定义对象", false);
setText(2, "col-body-1", "• 蛋白、isoform、复合物如何区分\n• canonical化合物与exact form如何区分\n• disease source ID与canonical疾病如何区分", false);
setText(2, "col-head-2", "再定义指标", false);
setText(2, "col-body-2", "• ABC是膜结合方式，E1–E3是证据把握度\n• BE等级由实验上下文决定\n• log、z-score、ECDF和富集值只回答特定问题", false);
setText(2, "col-head-3", "最后读主图", false);
setText(2, "col-body-3", "• 按A–F面板依次讲\n• 先说图展示了什么，再说不能推出什么\n• 原图为V6.2冻结统计；V6.3.1规则与HPC清单另行说明", false);

// Slide 3: keep the detailed source/version table as the M1 prerequisite page.
setText(3, "section-kicker", "M1 · 数据库架构", false);
setText(3, "slide-title", "前置知识：数据来源、身份统一、跨库去重与冻结发布", false);
setText(3, "source-rule-banner", "核心规则：保留原始ID与版本；身份唯一才自动合并；名称相同不等于同一实体；冲突和复核记录不静默删除。", false);

// Slide 5: M2 prerequisites.
setText(5, "section-kicker", "M2 · 膜蛋白全景", false);
setText(5, "slide-title", "前置知识：ABC定义生物学范围，E1–E3定义证据把握度", false);
setBox(5, 0, 0, "A类 · 整合膜蛋白", "具有跨膜区，或结构/实验明确支持其穿越脂双层。类别说明膜结合方式，不代表证据等级。");
setBox(5, 0, 1, "B类 · 直接嵌入膜", "没有经典跨膜螺旋，但具有膜内区、单层嵌入或脂质锚定等直接膜嵌入方式。");
setBox(5, 0, 2, "C类 · 外周膜相关", "通常不穿膜；由实验、稳定复合物或可靠定位支持其与膜表面持续关联。");
setBox(5, 1, 0, "E1 / E2 / E3 / E0", "E1直接实验或结构；E2可靠人工注释/多来源；E3预测或较弱定位；E0未达发布标准。");
setBox(5, 1, 1, "结构等级", "S2为膜环境校正的实验PDB；S3为其他实验PDB；S4仅AlphaFold；S0尚无指定结构。");
setBox(5, 1, 2, "拓扑与组装体", "TM count是跨膜区数量；单体、二聚体等为生物学组装状态，unknown与state-dependent必须保留。");

// Slide 7: M3 prerequisites.
setText(7, "section-kicker", "M3 · 表达与定位", false);
setText(7, "slide-title", "前置知识：HPA映射状态、表达变换、IHC及相对热图", false);
setBox(7, 0, 0, "nTPM · 组织RNA", "HPA normalized transcripts per million；图中以nTPM≥1作为RNA检出阈值之一。");
setBox(7, 0, 1, "nCPM · 单细胞RNA", "normalized counts per million；用于比较细胞类型，但不能与组织nTPM直接互换。");
setBox(7, 0, 2, "IHC · 蛋白染色", "免疫组织化学显示组织切片中的蛋白染色；受抗体特异性、组织处理和判读影响。");
setBox(7, 1, 0, "0不等于缺失", "mapped zero表示有映射但未检出；unmapped/missing表示未映射或无完整数据，解释必须分开。");
setBox(7, 1, 1, "log1p与z-score", "log1p压缩长尾；行内z-score显示同一功能类的相对高低，不代表跨蛋白绝对表达量。");
setBox(7, 1, 2, "ECDF与RNA–IHC", "ECDF读累计比例；RNA与IHC一致性是检出层面对照，不代表两种实验完全等价。");

// Slide 9: M4 prerequisites.
setText(9, "section-kicker", "M4 · 疾病关联", false);
setText(9, "slide-title", "前置知识：疾病身份、证据渠道、网络筛选与器官一致性", false);
setBox(9, 0, 0, "唯一蛋白–疾病pair", "同一canonical UniProt与同一疾病只统计一次；来源记录和证据渠道仍分别保留。");
setBox(9, 0, 1, "器官系统分组", "原图采用疾病名称关键词分组，仅用于可视化；不是患病率，也不是最终本体分类。");
setBox(9, 0, 2, "疾病证据等级", "very high、high、medium综合遗传、临床、功能、模型、表达、文献及UniProt疾病注释。");
setBox(9, 1, 0, "证据渠道", "人类遗传、临床遗传、体细胞突变、功能实验、动物模型、表达、文献和UniProt注释。");
setBox(9, 1, 1, "网络与疾病负担", "网络只画高置信高连接节点；每蛋白疾病数反映关联广度，不等于疾病严重程度。");
setBox(9, 1, 2, "Jaccard一致性", "交集/并集衡量表达器官与疾病器官的重叠；是描述性关联，不能推出组织因果关系。");

// Slide 11: M5 prerequisites. The original panel B is not a true UpSet plot.
setText(11, "section-kicker", "M5 · 小分子身份与化学空间", false);
setText(11, "slide-title", "前置知识：parent/form层级、状态组合、ECDF、Morgan指纹与UMAP", false);
setBox(11, 0, 0, "Canonical parent / exact form", "parent用于跨库归一；exact form保留盐型、立体化学、同位素、质子化与电荷等具体形式。");
setBox(11, 0, 1, "状态组合柱状图", "原图B统计Approved、Clinical、Probe、Endogenous、Natural的组合；它不是矩阵式UpSet图。");
setBox(11, 0, 2, "计算结构类别", "按可计算结构特征分成有机小分子、脂质、糖类、核苷等；属于描述标签，不是药效等级。");
setBox(11, 1, 0, "物化ECDF", "每一点表示性质≤横轴数值的化合物比例；MW、XlogP、TPSA和可旋转键用于描述与Docking准备。");
setBox(11, 1, 1, "Morgan UMAP与靶点广度", "Morgan指纹编码局部子结构，UMAP做二维展示；距离仅描述样本结构邻近。靶点数反映promiscuity。");

// Slide 13: M6 prerequisites with assay-context correction.
setText(13, "section-kicker", "M6 · 结合证据与位点", false);
setText(13, "slide-title", "前置知识：BE等级看实验上下文，不按Kd、Ki、IC50名称机械分级", false);
setBox(13, 0, 0, "BE1 · 结构/位点直证", "实验复合物或可追溯残基接触直接显示蛋白与配体同处一个结合位点。");
setBox(13, 0, 1, "BE2 · 定量直接结合", "由直接或竞争结合assay支持。Kd通常属于此类；Ki或IC50只有在结合实验中才可进入BE2。");
setBox(13, 0, 2, "BE3 / MECH", "BE3为蛋白特异功能或药理结果，如功能IC50、EC50；MECH为机制关系，不能当作亲和力。");
setBox(13, 1, 0, "数据库数≠实验数", "同一论文、PDB或BioAssay被多个数据库转载，只能算多个贡献库，不能算多个独立实验。");
setBox(13, 1, 1, "定量值不可混同", "Kd、Ki、IC50、EC50即使都换算为nM，热力学与功能含义仍不同；曲线只作分层描述。");
setBox(13, 1, 2, "位点与富集", "位点需PDB/链/配体/残基溯源；富集热图为观察值相对期望值的Pearson残差，不是结合强度。");

// Slide 15: M7 prerequisites.
setText(15, "section-kicker", "M7 · 负证据、冲突与质量", false);
setText(15, "slide-title", "前置知识：阴性结果、身份解析、冲突标记与发布质量门", false);
setBox(15, 0, 0, "阳性与阴性分表", "active/结合证据与inactive/negative分别保存；阴性只代表特定条件下未检出。");
setBox(15, 0, 1, "正负冲突", "同一canonical pair同时出现阳性和阴性时不互相覆盖；需比较assay、浓度、构建体与终点。");
setBox(15, 0, 2, "冲突率与分层", "冲突率=冲突pair/阳性pair；按功能类和BE等级查看冲突集中在哪里。");
setBox(15, 1, 0, "负证据身份解析", "结构身份唯一且外键完整者进入正式负证据层；无法可靠解析者进入review queue。");
setBox(15, 1, 1, "缺失不等于错误", "missing、not applicable、unmapped和review含义不同；图中缺失率描述可用性，不直接评价来源质量。");
setBox(15, 1, 2, "发布QA", "主外键、来源数量、去重、manifest、文件哈希和抽样复核通过后，记录才进入冻结release。");

// Slide 17: M8 prerequisites.
setText(17, "section-kicker", "M8 · Docking优先级", false);
setText(17, "slide-title", "前置知识：候选漏斗、R0–D3、结构准备度和HPC工作量", false);
setBox(17, 0, 0, "Docking候选漏斗", "从全部阳性pair依次经过证据、冲突、结构和配体可用性过滤；筛选完成不等于已经运行Docking。");
setBox(17, 0, 1, "R0 · Redocking", "已有pair-specific实验复合物，用原配体回对接验证流程；通常要求参考构象RMSD≤2 Å。");
setBox(17, 0, 2, "D1 / D2 / D3", "D1优先级最高；D2为扩展结构/定量层；D3偏批准、临床、探针或内源性机制探索。");
setBox(17, 1, 0, "结构与配体门槛", "受体需确定链、构象、组装体和box；配体需标准结构、质子化/互变异构体和3D构象。");
setBox(17, 1, 1, "Pilot / Standard / Full", "三套规模用于逐级验证和控制算力。原图为V6.2冻结统计；正式HPC以V6.3.1 manifest为准。");

// Figure pages: add the untouched July 30 originals as the topmost, maximum-size content.
const FIGURES = {
  4: ["assets/M1_database_architecture_quality.png", "M1 database architecture, integration and quality control"],
  6: ["assets/M2_membrane_proteome_landscape.png", "M2 landscape of the human membrane proteome"],
  8: ["assets/M3_expression_localization_atlas.png", "M3 tissue, cell-type and subcellular localization atlas"],
  10: ["assets/M4_disease_association_landscape.png", "M4 disease-association landscape"],
  12: ["assets/M5_chemical_identity_space.png", "M5 chemical identity and small-molecule space"],
  14: ["assets/M6_binding_evidence_sites.png", "M6 binding evidence and binding-site landscape"],
  16: ["assets/M7_negative_conflicts_quality.png", "M7 negative evidence, conflicts and data reliability"],
  18: ["assets/M8_docking_prioritization.png", "M8 structure-aware docking prioritization"],
};
for (const [slideNumber, [path, alt]] of Object.entries(FIGURES)) {
  const slide = deck.slides.items[Number(slideNumber) - 1];
  const image = slide.images.add({
    blob: await fs.readFile(path),
    contentType: "image/png",
    alt,
    fit: "contain",
    name: `original-full-figure-${slideNumber}`,
  });
  image.position = { left: 35.9, top: 0, width: 1208.2, height: 720 };
}

const figureTalks = {
  4: `M1讲稿｜数据库架构、整合流程与质量控制

开场：这张图回答的不是“我们下载了多少文件”，而是不同来源怎样变成可追溯、可去重、可发布的统一数据库。

A 数据来源与输出模块：左侧按蛋白、结构、互作、表达和疾病等来源进入MemPro integration；右侧输出Protein、Disease、Binding & sites和Compound。这里强调每个来源贡献不同信息，不是简单把表拼在一起。

B 标准化流水线：Raw records依次经过ID mapping、Parent/form、跨来源去重、证据分级和Release。蛋白以UniProt accession为核心；化合物区分canonical parent和exact form；疾病保留source ID并进行安全canonical映射。

C 记录去向：增量阳性输入与负证据身份输入分别显示进入正式层、复核层、排除或重复的数量。Review不是删除，而是身份或证据尚不足以自动进入默认发布层。

D 关系模型：膜蛋白是主要锚点，并与基因–疾病、化合物parent/form、结合证据、结合位点及表达定位形成外键关系。这样可以从任一实体追溯到原始来源。

E 来源–模块覆盖热图：行是数据库，列是模块；颜色使用log10(count+1)压缩数量级，格内数字保留原始计数。深色表示贡献规模更大，不代表数据质量一定更高。

F V6.2冻结快照：用对数横轴并列展示膜蛋白、疾病关系、结合位点、阳性pair、化合物和负证据等规模。这里报告的是实体或关系数量，不能把不同类型数字直接相加。

结论：MemPro的核心贡献是身份统一、证据分层和可追溯发布，而不只是数据量大。

不能过度解释：这张图是2026-07-30的V6.2冻结统计；V6.3.1增加的规则和后续清单应以最新manifest为准。

[Sources]
- Internal figure: M1_database_architecture_quality.png, MemPro V6.2 frozen figure package, 2026-07-30
- Internal caption: CAPTIONS_M1_M4.md`,

  6: `M2讲稿｜人类膜蛋白库的组成、证据、拓扑、结构与组装体

开场：这张图以10,997个UniProt accession为统计单位，同时展示“属于哪类膜蛋白”和“我们有多大把握”两条不同轴。

A 膜类别×证据等级：A、B、C是膜结合方式；E1–E3是证据强度，E0为未达到正式证据层。堆叠条形图不能读成A一定比C可信，因为类别和证据是两个变量。

B 功能层级：展示数量最多的功能类别，包括酶、受体、转运体、离子通道及仍未充分分类的蛋白。other或unclassified说明功能本体仍可继续补充，不等于膜蛋白身份错误。

C 跨膜拓扑：按0、1、2–6、7、8–12和>12个跨膜区分组。0跨膜可以包括B类膜嵌入或C类外周膜相关蛋白。

D 序列长度与TM数：六边形密度图展示蛋白长度和跨膜区数量的关系，并报告Spearman相关。相关反映总体趋势，不能用于判定单个蛋白类别。

E 各功能类的结构覆盖：区分膜校正实验PDB、其他PDB、仅AlphaFold和无指定结构。颜色是每个功能类内部百分比，因此适合比较结构可用程度。

F 寡聚与组装状态：展示单体、二聚体、三聚体、四聚体、更大/未解析、状态依赖及未知。状态依赖可能源于配体、构象、膜环境或实验条件变化，不应强制压成一个数字。

结论：蛋白库不仅给出名单，还保留膜结合方式、证据、拓扑、功能、结构与组装状态。

[Sources]
- Internal figure: M2_membrane_proteome_landscape.png, MemPro V6.2 frozen figure package
- Internal caption: CAPTIONS_M1_M4.md`,

  8: `M3讲稿｜组织、细胞类型与亚细胞定位图谱

开场：这张图回答膜蛋白“在哪里表达”和“定位到哪个细胞区室”，同时把映射缺失和检测不到区分开。

A 解剖表达覆盖：人体轮廓是导航示意，旁边条形图才是定量结果。计数表示在相应器官系统中达到RNA检出阈值的唯一蛋白数。

B 组织表达热图：对组织RNA先做log1p，再在每个功能类内部计算z-score。红蓝表示该功能类内部相对高低，不是不同蛋白之间的绝对表达比较。

C 细胞类型热图：处理逻辑与B一致，但数据单位是单细胞表达。nCPM不能与组织nTPM直接比较绝对值。

D 表达广度ECDF：横轴是每个蛋白被检出的组织数，纵轴是累计比例。例如曲线上50%的位置可读中位表达广度。曲线越靠右，表示该类蛋白通常在更多组织检出。

E 亚细胞定位：条形图统计细胞质、核、质膜、囊泡等定位术语的唯一蛋白数；术语非互斥，一个蛋白可同时具有多个定位。

F RNA–IHC一致性：把蛋白–组织层面的RNA检出和免疫组化蛋白染色组成2×2矩阵。对角线是一致，非对角线是不一致；不一致可能来自抗体、转录后调控或检测灵敏度。

必须强调：mapped zero与unmapped/missing不是一回事；不能把所有0解释为完全不表达。

[Sources]
- Internal figure: M3_expression_localization_atlas.png, MemPro V6.2 frozen figure package
- HPA 25.1-derived internal expression/localization tables
- Internal caption: CAPTIONS_M1_M4.md`,

  10: `M4讲稿｜蛋白–疾病关联、证据渠道和器官系统图景

开场：这张图统计的是唯一蛋白–疾病pair，重点是关联范围和证据结构，不是疾病患病率。

A 按器官系统的疾病关联：人体图与条形图显示每个器官系统涉及多少唯一pair。原图使用透明的疾病名称关键词规则，因此“其他/多系统”可能同时包含真实跨系统疾病和未被关键词命中的疾病。

B 疾病系统×证据等级：堆叠条形图把very high、high和medium按器官系统展开。总条长代表pair数量，颜色构成代表证据等级分布。

C 证据架构热图：每一列是人类遗传、临床遗传、体细胞突变、功能实验、动物模型、表达、文献或UniProt疾病注释；格子表示该系统中具有相应证据的pair比例。

D 高置信网络：只选择高置信且连接度较高的蛋白和疾病以保证可读性。它用于举例说明枢纽关系，不是完整关系网；没画出的节点不代表不存在。

E 每个蛋白的疾病负担：互补累计曲线展示每个蛋白关联至少多少疾病。横纵轴均为对数，更容易看长尾。疾病数多表示关联广，不代表单个关系更可靠或疾病更严重。

F 表达–疾病器官一致性：用Jaccard指数比较表达器官蛋白集合和疾病器官蛋白集合。数值高表示集合重叠较多，但不能推出表达组织导致疾病。

方法限制：正式投稿应优先使用MONDO exact-only、Open Targets治疗领域和DO/Uberon多标签映射；原图关键词分组只适合作为历史探索图。

[Sources]
- Internal figure: M4_disease_association_landscape.png, MemPro V6.2 frozen figure package
- Open Targets 26.06 and UniProtKB disease-derived internal tables
- Internal caption: CAPTIONS_M1_M4.md`,

  12: `M5讲稿｜化合物身份、物化性质、结构空间和靶点广度

开场：这张图先解决“什么算同一个化合物”，再描述化合物具有什么结构和生物学状态。

A canonical parent与form层级：源记录先解析exact structure，再归一到canonical parent，同时保留盐型、立体异构体、同位素、质子化和电荷等form。Review表示结构无法安全自动合并。

B 生物状态组合：横条展示Approved、Clinical、Probe、Endogenous和Natural等标签组合的数量，并使用对数横轴。注意这是一张“组合柱状图”，不是标准UpSet矩阵；组合数字也不是获批药物总数。

C 计算结构类别：对QC通过的canonical化合物按结构规则分类。类别用于描述数据库组成，不等于成药性或药理作用分类。

D 物化性质ECDF：四个小图分别为MW、XlogP、TPSA和可旋转键。横轴给性质值，纵轴给≤该值的累计比例；它们影响配体准备、溶解性和构象搜索，但不能单独证明药效。

E Morgan指纹化学空间：以radius 2、1024-bit Morgan指纹编码局部子结构，再用UMAP投影到二维。相近点通常共享更多局部结构，但二维距离不是绝对相似度，颜色重叠也不代表类别错误。

F 靶点promiscuity：互补累计曲线显示一个化合物拥有至少多少支持靶点。长尾中的多靶点化合物可能是真实多药理，也可能受检测频率和非特异性影响。

结论：SMILES适合表达和计算，InChIKey适合检索与身份匹配；数据库必须同时保留parent和form层级。

[Sources]
- Internal figure: M5_chemical_identity_space.png, MemPro V6.2 frozen figure package
- Internal caption: CAPTIONS_M5_M8.md`,

  14: `M6讲稿｜结合证据、定量结果、位点可用性和类别富集

开场：这张图把“有关系”“有定量结合”“有结构位点”分开。最重要的规则是BE等级由实验上下文决定，而不是只看字段名。

A 来源×BE等级：每一条证据先取互斥的最高BE等级，再按数据库统计。BE1是结构/残基层直证；BE2是assay明确的定量直接或竞争结合；BE3是蛋白特异功能或药理证据。

需要纠正的口径：Kd通常直接描述平衡亲和力；Ki只有在明确竞争结合模型下才支持BE2；IC50若来自直接竞争结合assay也可以是BE2，而功能IC50和EC50通常属于BE3。因此不能写成Kd/Ki天然比IC50/EC50高一级。

B 来源支持数：统计每个唯一蛋白–canonical化合物pair来自多少数据库。这个字段只能叫distinct contributing database count；同一论文或实验被多库转载不能算多个独立实验。

C 定量值ECDF：把可标准化数值换算为nM并用对数横轴展示。不同BE层的曲线可以比较数据范围，但Kd、Ki、IC50和EC50不应被解释为热力学等价量。

D 结合位点准备度：区分PDB+残基接触、PDB结构匹配、实验口袋、仅复合物PDB和UniProt人工位点等。它统计位点实例，不覆盖所有阳性pair。

E 结合残基组成：从可解析位点描述中统计疏水、芳香、极性、正电、负电和特殊残基提及次数；它反映已解析位点的组成，不是全蛋白氨基酸频率。

F 蛋白类×化学类富集：颜色是Pearson残差，即观察计数相对独立模型期望值的偏离。正值表示组合出现多于期望，负值表示少于期望；它不是亲和力、效应量或因果关系。

[Sources]
- Internal figure: M6_binding_evidence_sites.png, MemPro V6.2 frozen figure package
- Internal caption: CAPTIONS_M5_M8.md`,

  16: `M7讲稿｜负证据、正负冲突、缺失与发布质量

开场：数据库不能只保留阳性结果。M7说明阴性结果如何映射、冲突如何保留，以及哪些记录仍需复核。

A 阳性与阴性pair覆盖：对主要功能类比较唯一阳性pair和负pair数量，横轴为对数。阴性量大常与大规模筛选有关，不代表该类蛋白更不容易结合。

B 正负冲突率：计算每个功能类中“同时具有阳性和至少一条阴性记录”的pair占阳性pair比例。冲突可能来自浓度、终点、构建体、细胞背景或检测灵敏度。

C 冲突×BE等级热图：显示冲突pair在功能类与BE1–BE3之间的分布；颜色用log10(count+1)，格内保留原始计数。高计数可能首先反映该类数据量大。

D 负证据身份解析：从10,246,948条源记录分到正式映射层、review queue和重复删除。未映射记录不进入默认负证据表，但不会被静默丢弃。

E 核心字段缺失率：对蛋白结构、化合物性质、pair定量值和位点字段报告missing/not available比例。高缺失说明该字段不适用于所有记录或来源覆盖不足，不一定是错误。

F 复核队列与验证门：并列显示阳性身份复核、负证据复核、正负冲突、亚基复核、表达定位复核和blocking QA error。Blocking error为0说明冻结门通过，不等于所有review项目都已人工解决。

结论：正、负、冲突和未解析是四种状态，应分别管理，不能用一个最终真假标签覆盖实验条件。

[Sources]
- Internal figure: M7_negative_conflicts_quality.png, MemPro V6.2 frozen figure package
- Internal caption: CAPTIONS_M5_M8.md`,

  18: `M8讲稿｜结构感知的Docking优先级与HPC工作量

开场：这张图展示怎样把150万级阳性pair缩减成可验证、可分批运行的Docking集合，而不是展示已经得到的Docking结果。

A 候选漏斗：从1,502,456个阳性pair到R0–D3 eligible、去冲突、Standard和6,019个Pilot。横轴为对数；每一级减少都来自明确质量门。

B Standard层级：R0用于方法验证；D1–D3是逐级降低优先级的前瞻或机制任务。柱高表示pair数量，不代表计算已完成。

C 靶点类别×Docking tier：热图显示受体、酶、转运体、离子通道等功能类在R0–D3中的pair数，颜色为log10(count+1)，格内为原始数量。

D 受体结构准备度：每个tier内部比较S1 pair-specific实验位点、S2膜环境实验结构和S3其他实验PDB的比例。结构等级高不代表蛋白准备已经自动完成。

E 配体可计算性：比较各tier的MW、XlogP、TPSA和可旋转键ECDF，用来估计质子化、构象采样和搜索空间难度；不是Lipinski硬性淘汰。

F HPC工作量：Pilot、Standard和Full分别统计pairs、targets和compounds，并使用对数纵轴。它帮助决定阵列任务数量、存储和预计核时。

执行前仍需完成：确定受体链和生物学组装体、处理缺失残基和辅因子、确定质子化状态、生成配体3D、定义box、用R0完成redocking，并在HPC上进行真实调度。

版本说明：原图来自V6.2冻结图包；正式执行应读取V6.3.1 docking manifest，原图用于解释筛选方法和统计结构。

[Sources]
- Internal figure: M8_docking_prioritization.png, MemPro V6.2 frozen figure package
- Internal caption: CAPTIONS_M5_M8.md
- Internal docking pilot/standard/full manifests`,
};

setNotes(1, "开场介绍项目名称和目标：建立可追溯的人类膜蛋白–小分子互作数据库。\n\n[Sources]\n- Internal: MemPro V6.3.1 project materials");
setNotes(2, "说明汇报采用固定双页结构：第一页建立术语和限制，第二页只讲一张原始主图。不要跳过前置页直接解释颜色和数值。\n\n[Sources]\n- Internal: MemPro presentation design and frozen figure captions");

const prereqNotes = {
  3: "M1前置页先介绍四类来源及统一规则。强调数据库整合不是名称拼接，必须保留来源ID、版本和映射路径。\n\n[Sources]\n- Internal: MemPro source registry and M1 caption",
  5: "M2前置页必须先把ABC和E等级分开。随后说明结构等级、跨膜次数和寡聚状态，避免把不同维度混成一个分类。\n\n[Sources]\n- Internal: MemPro protein master and M2 caption",
  7: "M3前置页先解释nTPM、nCPM和IHC，再解释log1p、z-score与ECDF。重点提醒0值和缺失状态必须区分。\n\n[Sources]\n- Internal: HPA-derived MemPro expression module and M3 caption",
  9: "M4前置页先声明原图器官分类是关键词可视化，而不是最终本体映射；正式分析保留MONDO/OMIM等原始ID。\n\n[Sources]\n- Internal: MemPro disease relation tables and M4 caption",
  11: "M5前置页明确原图B不是标准UpSet；Morgan UMAP是描述性投影，物化性质用于数据描述和Docking准备。\n\n[Sources]\n- Internal: MemPro compound parent/form tables and M5 caption",
  13: "M6前置页纠正证据分级口径：assay上下文优先于指标名称。Kd通常是直接结合；Ki和IC50必须查看实验类型；EC50通常为功能效力。\n\n[Sources]\n- Internal: MemPro binding evidence rules and M6 caption",
  15: "M7前置页区分阳性、阴性、冲突、未映射和缺失。Review queue既不是默认发布层，也不是删除区。\n\n[Sources]\n- Internal: MemPro negative evidence and QA tables",
  17: "M8前置页说明筛选清单不等于Docking结果；R0先验证协议，D1–D3再分批运行。\n\n[Sources]\n- Internal: MemPro docking manifests and M8 caption",
};
for (const [slide, note] of Object.entries(prereqNotes)) setNotes(Number(slide), note);
for (const [slide, note] of Object.entries(figureTalks)) setNotes(Number(slide), note);

const scriptText = [
  "MemPro V6.3.1｜M1–M8逐图讲稿",
  "说明：PPT中的主图均为2026-07-30冻结的原始M1–M8图；每幅图前一页用于解释术语、指标和限制。",
  "",
  ...Object.keys(figureTalks).sort((a, b) => Number(a) - Number(b)).map((key) => figureTalks[key]),
].join("\n\n============================================================\n\n");
await fs.writeFile(SCRIPT_OUTPUT, scriptText, "utf8");

await fs.mkdir("final/slides", { recursive: true });
await fs.mkdir("final/layouts", { recursive: true });
for (let index = 0; index < deck.slides.items.length; index += 1) {
  const slide = deck.slides.items[index];
  const number = String(index + 1).padStart(2, "0");
  const png = await deck.export({ slide, format: "png", scale: 1.5 });
  await fs.writeFile(`final/slides/slide-${number}.png`, new Uint8Array(await png.arrayBuffer()));
  const layout = await slide.export({ format: "layout" });
  await fs.writeFile(`final/layouts/slide-${number}.layout.json`, await layout.text(), "utf8");
}
const montage = await deck.export({ format: "webp", montage: true, scale: 0.55 });
await fs.writeFile("final/montage.webp", new Uint8Array(await montage.arrayBuffer()));
const finalInspect = await deck.inspect({ kind: "slide,textbox,shape,image,notes", maxChars: 1000000 });
await fs.writeFile("final/final-inspect.ndjson", finalInspect.ndjson, "utf8");

const pptx = await PresentationFile.exportPptx(deck);
await pptx.save(OUTPUT);
console.log(JSON.stringify({ output: OUTPUT, script: SCRIPT_OUTPUT, slides: deck.slides.items.length }));
