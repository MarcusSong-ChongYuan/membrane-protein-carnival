# M1–M8 图解与投稿式 Caption

## M1｜多来源整合的新增覆盖与冗余

Panel A 为蛋白—小分子 pair 的 UpSet 图，柱高表示指定来源组合中的唯一 pair 数，点阵表示来源组合。Panel B 以 Jaccard 系数衡量任意两来源的 pair 集合重叠。Panel C 按贪心顺序加入数据库，区分单个来源的新增 pair 与累计覆盖。Panel D 统计只由一个来源贡献的严格唯一 pair；ChEMBL、PubChem BioAssay、PDBe、BindingDB 和历史混合来源分别贡献 253,918、130,912、9,135、6,149 和 5,960 个严格唯一 pair。该图说明多来源整合既增加覆盖，也包含必须识别的转载或交叉收录。

## M2｜五轴交叉分类形成的膜蛋白注释空间

Panel A 将结构家族、分子功能、生物过程、具体膜角色和专科分类转换为多标签集合，以 Jaccard 距离进行 UMAP，点代表 7,800 个正式 A/B/C 膜蛋白。Panel B 展示最大的 HDBSCAN 社区及其主要注释；算法识别 23 个社区，并将 554 个证据弥散蛋白保留为 noise，避免强行分类。Panel C 比较大型社区由五个注释轴提供的相对支持。该图将原先单一 `other/unclassified` 标签改写为可交叉查询的生物学空间。

## M3｜正常组织表达与 RNA–IHC 一致性

Panel A 用人体器官图表示 HPA 中检测到 RNA 的膜蛋白数量。Panel B 对 7,496 个具组织 RNA 矩阵的膜蛋白进行表达 UMAP。Panel C 展示 150 个高变异膜蛋白在 20 个层次聚类组织中的标准化表达。Panel D 在大小写标准化后的同名组织中计算 `log1p(RNA)` 与 IHC 检出状态的 Spearman 相关；39 个组织满足至少 100 个匹配蛋白且 IHC 状态可变。RNA 与 IHC 是不同测量层，相关性不表示两者可互相替代。

## M4｜疾病器官负担与正常组织表达一致性

Panel A 将 MONDO/DO/Uberon 多标签疾病—解剖映射投射到人体器官。Panel B 比较每个器官中疾病相关膜蛋白的实际 RNA 检出比例与 500 次同规模随机抽样，显示多个器官存在高于随机期望的一致性；器官 z score 的中位数为 3.82，最高为 7.17。Panel C 展示不同疾病系统的膜角色构成。Panel D 展示 medium、high 和 very high 疾病证据在系统间的分布。此图的系统分类来自本体层级，不是关键词单选分类。

## M5｜化学空间与表示方法比较

Panel A 使用 ECFP4（radius=2，1,024 bits）对 15,000 个确定性抽样 canonical 小分子进行 UMAP，hexbin 表示局部密度，并叠加已批准、临床候选、内源性及天然产物状态。Panel B 展示 2,400 个化合物的 PAGTN pilot 表示，并以主要结构类别标色。Panel C 比较局部近邻保留率、聚类纯度和非噪声分配率。PAGTN pilot 由分子量、XlogP、TPSA 和可旋转键监督训练，因此仅表示物化性质学习的可行性，不是通用预训练模型。

## M6｜多靶点药理网络与膜侧化学环境

Panel A 将 220 个配体数最高的蛋白投影为共享配体网络；901 条边表示靶点间至少共享 3 个化合物，节点颜色为 Louvain 社区。Panel B/C 使用 7,279 条具有明确膜侧注释的位点—化合物记录，比较胞外、膜界面、膜内、胞质侧和非胞质侧的 XlogP 与 TPSA 分布。膜侧来自位点/结构注释，物化性质仅用于描述关联，不用于反向推断膜侧。

## M7｜证据组合、来源谱系与 BE 等级

Panel A 的 UpSet 图总结 942,455 条正式阳性证据的 7 类记录级组合，包括定量、功能、结构和位点信息。Panel B 展示不同数据库能够解析到数据库记录、文献实验、PubChem AID/SID/CID 或 PDB/链/配体结构键的数量。Panel C 展示 quantitative binding、bioassay、structure 和 functional pharmacology 等模态在 BE1/BE2/BE3 中的组成。数据库来源数量不等于独立实验数量；只有解析到相同原始实验键后才能识别转载。

## M8｜生物学价值、结构准备度与 Docking 缺口

Panel A 将疾病数量、正式配体数量和证据模态整合为透明的 biology score，将结构数量、位点数量和坐标可用性整合为 readiness score。Panel B/C 用分区和 Pareto 视角识别高生物学价值但结构准备不足的靶点。Panel D 将 529,168 个正式蛋白—小分子 pair 分为 G1 坐标就绪 14,503、G2 有位点但坐标不完整 1,269、G3 有结构无位点 3 和 G4 仅互作 513,393。该分层用于决定实验与计算资源，不改变数据库中的互作证据等级。

