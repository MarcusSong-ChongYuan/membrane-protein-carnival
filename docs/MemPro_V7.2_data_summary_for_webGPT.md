# MemPro V7.2 数据统计总结与投稿作图审阅包

用途：将本文件直接上传给网页版 ChatGPT，请它审阅 MemPro V7.2 数据结构、统计边界和投稿级作图方案。本文档是统计摘要，不替代原始数据库；所有数字应以 V7.2 正式表和随附 Source Data 为准。

## 1. 请网页版 ChatGPT 完成的任务

请基于本文件和随后上传的 Source Data：

1. 审查数据库规模、分母、缺失值语义和统计口径是否一致；
2. 判断哪些分析可以作为 NAR Database Issue 主文结果，哪些应放补充材料；
3. 为每个模块提出非单一柱状图的统计图方案；
4. 检查是否存在过度解读、伪因果、分母错误、来源重复被误称为独立实验等问题；
5. 给出一套有叙事逻辑的 Graphical Abstract + 主图 + 补图建议，并明确每幅图的核心结论、输入字段、统计方法和潜在陷阱。

请不要把数据库行数直接解释为独立实验数，也不要把疾病—表达一致性解释为因果关系。

## 2. 数据库正式规模

| 模块 | 数量 | 统计单位 | 说明 |
|---|---:|---|---|
| Formal membrane proteins | 7,800 | canonical UniProt protein | V7.2 正式膜蛋白层，A/B/C 分类；不等于全部人类蛋白 |
| E1/E2 default proteins | 7,374 | protein | 默认确认层；E3 单独报告 |
| A/B/C membrane classes | A=5,680; B=613; C=1,507 | protein | A 严格整合膜蛋白；B 直接嵌入/脂质锚定；C 外周膜相关 |
| Compound registry | 646,670 | compound entity | 注册小分子主表，不代表都有膜蛋白互作 |
| Interaction-linked compounds | 240,055 | canonical compound | 至少进入正式蛋白—小分子互作关系 |
| Formal protein–compound pairs | 529,168 | unique protein–compound pair | 2,829 个膜蛋白有正式互作关系 |
| Public positive evidence | 942,455 | evidence record | 正证据记录；来源/实验谱系需分开解释 |
| Canonical diseases | 3,739 | canonical disease | MONDO/DO/来源实体标准化后的疾病实体 |
| Protein–disease relations | 6,753 | protein–disease relation | 涉及 2,347 个膜蛋白 |
| Expression measurements | 3,046,789 | expression record | 正常组织 RNA 等表达记录；不能直接称为蛋白表达 |
| Binding-site assertions | 55,610 | site assertion | 48,331 条膜侧为 unknown |
| Protein complexes | 4,584 | complex entity | 复合物上下文与 complex-specific binding 必须区分 |

### 重要分母边界

- 646,670 是小分子注册表分母；240,055 才是正式互作小分子分母。
- 529,168 是 pair 数；942,455 是证据记录数。一个 pair 可以有多个证据记录。
- 7,800 是正式膜蛋白层，不应写成“人类全部蛋白”。
- 3,739 是正式疾病实体；未映射解剖系统的疾病必须单独报告。

## 3. 蛋白模块

### 3.1 A/B/C 分类

| 类别 | 数量 | 含义 |
|---|---:|---|
| A | 5,680 | 具有跨膜区，或有较强结构/实验性膜蛋白证据 |
| B | 613 | 非典型跨膜，但直接嵌入膜、脂质锚定或膜内形式 |
| C | 1,507 | 外周膜相关、膜募集或膜稳定定位，但通常无跨膜段 |

### 3.2 E1/E2/E3 证据层

| 等级 | 含义 |
|---|---|
| E1 | 直接实验、结构或强膜拓扑证据 |
| E2 | 人工注释、拓扑/多来源支持，但实验强度相对较弱或不完整 |
| E3 | 单一预测、by-similarity、弱定位或外周膜关联；可保留，但不应和确认层混写 |

V7.2 的 A/B/C × E1/E2/E3 明细见：
`figures_nar_final/03_source_data/Figure2A_membrane_class_evidence.tsv`

### 3.3 五轴交叉分类

每个蛋白不再只有一个功能标签，而是并行保留：

- `structural_family`
- `molecular_function`
- `biological_process`
- `membrane_role`
- `specialist_classification`

这五个轴回答不同问题：同源结构、分子功能、生物过程、膜上的工作角色、专科数据库层级。应使用交叉热图、NPMI/log2 odds ratio、alluvial 或分类完整性矩阵，而不是用一个家族标签替代所有功能。

已知缺口：

- specialist classification 缺失/占位约 4,338 个蛋白；
- membrane role 未决约 1,722 个；
- molecular function 未分类约 1,260 个；
- biological process 未分类约 922 个。

这些是注释缺口，不等于膜身份错误。

## 4. 小分子与化学空间

- 646,670 个小分子注册实体；240,055 个进入正式互作 pair。
- 91,977 个 Bemis–Murcko scaffold group。
- scaffold Gini = 0.5505，提示化学空间具有明显长尾和骨架集中现象。
- compound status 是非互斥二元标签；一个化合物可以同时是 approved drug、clinical candidate、endogenous、natural product 或 chemical probe。
- 物化性质和 status 分组应使用 robust z-score 或 ridgeline/raincloud；不能把 Lipinski 规则当作数据库纳入标准。
- Morgan 分析：radius=2，2048 bits，top-20 邻居；98,603 个化合物、1,333,186 个邻居 pair；Spearman rho=0.5817，cluster-bootstrap 95% CI=[0.5566, 0.6066]，近似邻居抽样召回率=0.964。

建议主图优先：骨架 rank-abundance、状态×膜角色富集、靶点广度×角色熵、Tanimoto×target-set Jaccard。UMAP/PAGTN 只能作为探索性或补图，不能单凭二维距离合并结构。

## 5. 证据来源与证据谱系

主要来源字段包括：

- ChEMBL
- PubChem BioAssay
- BindingDB
- PDBe
- IUPHAR/BPS Guide to PHARMACOLOGY
- PDSP KiDatabase
- BRENDA、PDBbind、OPM/PDBTM 等结构或功能来源

必须严格区分：

- `distinct_database_count`：贡献数据库数；
- `putative_experiment_lineage_count`：按 assay、PubMed、PDB、AID/SID 等键推断的实验谱系数；
- 未完成原文级实验谱系审计前，不应称为“人工确认的独立实验数”。

建议图：UpSet、来源 Jaccard/overlap coefficient 三角热图、database count × lineage count hexbin、来源边际贡献 dot heatmap。

## 6. 结合位点与结构覆盖

总 binding-site assertions = 55,610：

| 膜侧状态 | 数量 | 比例 |
|---|---:|---:|
| unknown | 48,331 | 86.9% |
| extramembrane side unresolved | 4,390 | 7.9% |
| membrane interface or mixed | 1,328 | 2.4% |
| intramembrane | 773 | 1.4% |
| cytoplasmic side | 407 | 0.7% |
| non-cytoplasmic side | 378 | 0.7% |
| both sides or mixed | 3 | <0.1% |

注意：没有结合残基不代表没有活性证据；只有活性但无残基的记录可以进入互作表，并标记 `residue_not_reported`。来源残基与 PDB 坐标只能部分匹配时，应保留原始来源残基，并记录 `coordinate_mapping_fraction`，不强制删除关系。

## 7. 疾病、治疗领域与表达

疾病层：

- 3,739 个正式疾病；
- 2,926 个有解剖系统映射；813 个没有解剖映射；
- 3,080 个有治疗领域映射；659 个没有治疗领域映射；
- 6,753 个 protein–disease relations；证据等级：very high=1,411，high=3,196，medium=2,146。

疾病器官分类应使用 MONDO/DO/Uberon 的多标签映射，不用关键词强迫一个疾病只属于一个器官系统。

表达图必须：

- 只在 `mapped + measured` 中计算表达比例；
- missing、not measured、unmapped 单独显示；
- RNA 与 IHC 分开；
- RNA 信号不能直接写成蛋白表达。

### 解剖—表达一致性分析

当前 10,000 次置换结果：

- observed mean = 0.4866；
- null mean = 0.3783；
- z = 26.43；
- empirical P < 1e-4；
- 置换保持 protein degree、disease degree 和 anatomy-label count。

科学表述应为“anatomical-expression concordance”或“hypothesis-generating consistency”，不能写成疾病因果关系。

## 8. 已完成的统计分析

| 分析 | 结果 |
|---|---|
| protein role × molecular function | bias-corrected Cramer’s V = 0.6188；BH-FDR 处理单元格 |
| tissue Tau by membrane role | n=7,305；Kruskal–Wallis P=4.34e-92；epsilon-squared=0.0614；Dunn/BH post hoc |
| scaffold diversity | 240,055 compounds；91,977 groups；Gini=0.5505 |
| chemical similarity × target overlap | rho=0.5817；95% CI=[0.5566,0.6066]；ANN recall=0.964 |
| disease anatomy-expression concordance | z=26.43；10,000 degree-preserving permutations |
| disease chemical landscape | 55 个透明候选；921 个 disease-associated proteins 的 formal compound count 为 0 |

## 9. 推荐投稿叙事

### Graphical Abstract

多源蛋白、化合物、证据、结构、疾病/组织数据
→ canonical identity / deduplication / lineage / ontology
→ MemPro 多层知识结构
→ target preference、chemical preference、evidence-aware interpretation、disease-relevant underexplored targets。

### Figure 1：数据整合和可追溯发布

流程图、来源—证据模态—BE tier alluvial、模块覆盖、正式规模卡片。

### Figure 2：五轴膜蛋白组织结构

A/B/C×E1-E3 mosaic、role×function enrichment、四轴 alluvial、Tau 雨云图。

### Figure 3：小分子化学多样性与靶点偏好

骨架 rank-abundance、状态×膜角色富集、target breadth×role entropy、化学相似性×靶点集合 Jaccard。

### Figure 4：数据库覆盖不等于独立证据

来源 UpSet、Jaccard/overlap 热图、数据库数×谱系数、来源贡献、位点映射 ECDF、膜侧完整性。

### Figure 5：疾病与人体解剖背景

高清人体图、role×therapeutic area 置换热图、疾病证据 UpSet、疾病支持×化学探索度、解剖—表达置换检验。

## 10. 不应过度解读的内容

- E3 不是错误，但不应默认作为确认层统计；
- 多个数据库来源不等于多个独立实验；
- unknown binding-site membrane side 不等于没有结合位点；
- 没有 residue 不等于没有 protein–compound interaction；
- canonical 指向不等于 isoform 已精确确认；
- component-context-only 不等于 complex-specific binding；
- disease-expression concordance 不等于因果关系；
- 646,670 registry compounds 不等于 646,670 个正式膜蛋白配体。

## 11. 可上传的配套文件

建议同时上传：

1. 本摘要文件；
2. `D:\finale\figures_nar_final\07_legends_alttext\Figure1_SourceData.zip` 至 `Figure5_SourceData.zip`；
3. `D:\finale\figures_nar_final\08_QA\NAR_FIGURE_QA_FINAL.json`；
4. 如需看图，再上传 `D:\finale\figures_nar_final\04_main_figures\png` 中的 5 张主图。

当前正式图包：`D:\finale\MemPro_V7.2_NAR_figures_final_20260818.zip`。

## 12. 给网页版 ChatGPT 的审阅提示词

请把自己当作 NAR Database Issue 的数据库论文审稿人和统计图设计顾问。基于上述 MemPro V7.2 数据摘要及我上传的 Source Data：

1. 逐项检查分母、统计单位、缺失值和来源字段；
2. 指出任何数字冲突或可能误导读者的图；
3. 对 Figure 1–5 和 Graphical Abstract 给出保留/修改/移入补充材料的建议；
4. 每个建议必须写明：核心科学问题、推荐图型、输入字段、统计检验、效应量、缺失处理、可解释性风险；
5. 提出不依赖单纯柱状图、但能从现有数据得到稳健结论的新分析；
6. 不要提出需要重新收集大规模数据的分析，除非明确说明其必要性和替代方案。
