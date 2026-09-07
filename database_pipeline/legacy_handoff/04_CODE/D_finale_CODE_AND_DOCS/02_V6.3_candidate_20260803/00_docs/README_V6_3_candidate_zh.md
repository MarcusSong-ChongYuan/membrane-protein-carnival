# MemPro V6.3 candidate

本目录是冻结的 V6.2 之上的增量候选版本。V6.2 文件没有被改写；V6.3 新增或更新的是实体层、疾病本体层、蛋白交叉分类、证据谱系、docking 排序和统计图。

## 本版新增

1. `gene -> canonical protein -> isoform -> complex` 四层靶点体系。相同基因名不触发蛋白或 isoform 合并；复合物不再强塞进单蛋白关系表。
2. 疾病采用 MONDO direct ID 或官方一对一 exact/equivalent 映射。broad、narrow、related、1:N、过时或非疾病实体均不自动合并。
3. 全部 10,997 个蛋白采用五轴交叉分类：结构家族、分子功能、生物过程、膜功能角色和专科分类。Reactome 作为独立通路层保留。
4. `independent_source_count` 的正式解释改为“贡献数据库数”。结构、PubChem assay、文献实验代理和证据模态分别建谱系键，互不替代。
5. docking 候选成员和 tier 沿用 V6.2；V6.3 只增加有上限的谱系加分及明确的正负证据冲突扣分。
6. M5E 改为全量 Bemis-Murcko 骨架多样性；899,591 个 core 且 QC 通过的 canonical 化合物全部处理，没有随机抽样或二维嵌入。

## 使用顺序

- 先读 `00_docs/MemPro_V6_3_candidate_review.xlsx` 和本文件。
- V6.2 完整冻结快照位于 `01_baseline_v6_2`。
- V6.3 的新增/替代表位于 `02_identity` 至 `07_docking`。
- 图片位于 `08_figures`，每张图同时提供 PNG、SVG 和 PDF，caption 为 Markdown。
- 所有质量控制文件位于 `09_qa`；总门控必须为 `PASS`。
- `10_scripts` 和 `11_source_snapshots` 用于复现本次增量构建。

## 重要边界

- 这是 `candidate_frozen`，不是已经对外发布的网站版本。
- 2,871,623 条记录构成机器筛选的谱系复核全集；没有声称这些记录已逐条人工审阅。固定分层抽样为 513 条。
- 复合物结合证据默认表目前为空，因为名称或来源标识尚不足以把候选升级为直接结合；候选没有删除。
- 任何对外发布前仍需完成各来源授权/再发布核查，尤其是受限来源的原始或衍生内容。
- V6.3 的 `family_defined_membrane_role_unresolved` 是诚实的未解状态；它不等同于数据缺失，也不应为图表美观而强行赋功能。

## 冻结规则

正式候选表、QA 报告、源码和来源快照均记录 SHA-256。任何修改都应产生新的候选版本和新 manifest，不允许原地覆盖。
