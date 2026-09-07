# MemPro V6.3 候选：基因—canonical 蛋白—isoform—复合物身份层

## 结论

本模块已完成并通过最终 QA。冻结的 V6.2 数据保持只读；新表以桥接方式叠加，不回写历史发布。

- 核心膜蛋白：10,997 条 canonical UniProt accession，数量不变。
- 基因实体：10,970 条。同一基因可连接多个 canonical 蛋白，但不会因此合并蛋白记录。
- 蛋白 isoform：20,211 条，均保存序列、长度、SHA-256 和 canonical 外键。
- 复合物：沿用 19,791 个审计实体；97,993 条组件断言已升级为 canonical/isoform 外键。
- 外部复合物蛋白：6,126 条，只用于解释复合物组成，不计入核心膜蛋白总数。
- 网站检索索引：支持 gene、canonical protein、protein isoform、protein complex 四种实体。

最终权威 QA 为 `qa/IDENTITY_LAYER_V0_3_VALIDATION.json`。早期 V0.1/V0.2 报告保留作开发审计，已被 V0.3 明确 supersede。

## 关键规则

1. canonical 蛋白主键为 `UNIPROT:<accession>`。
2. isoform 主键为 `UNIPROT_ISOFORM:<accession-n>`，并必须指向 canonical 蛋白。
3. 证据明确写 isoform accession 时才关联 isoform。
4. 只写 canonical accession 时只关联 canonical 蛋白，不推断 isoform。
5. Open Targets 和 HPA 的基因层记录标为 `gene_product_unspecified`；可保留 canonical 投影供兼容查询，但不能解释为特定 isoform 证据。
6. 相同基因名不会触发蛋白合并。
7. 复合物组件可以引用核心膜蛋白、核心 isoform、外部 canonical 蛋白或外部 isoform。
8. V6.2 的结合证据已经 canonical 化，历史 isoform 特异性无法从规范化表反推，因此统一标为 `not_recoverable_from_v62_normalized_record`。

## 主要文件

| 文件 | 用途 |
|---|---|
| `gene_entity_v0_1.tsv` | HGNC/Ensembl/NCBI Gene 统一后的基因实体 |
| `gene_identifier_v0_1.tsv` | 基因实体与各稳定 ID 的多值映射 |
| `canonical_protein_entity_v0_1.tsv.gz` | 10,997 条核心 canonical 膜蛋白 |
| `canonical_protein_gene_link_v0_1.tsv` | canonical 蛋白到基因的外键 |
| `protein_isoform_v0_2.tsv.gz` | isoform ID、序列、hash、canonical 外键和作用域 |
| `binding_evidence_target_resolution_v0_1.tsv.gz` | 3,003,306 条结合证据的目标层级桥接 |
| `binding_site_target_resolution_v0_1.tsv.gz` | 95,598 条位点证据的目标层级桥接 |
| `disease_relation_target_resolution_v0_1.tsv` | 蛋白级与 gene-product-unspecified 疾病关系区分 |
| `expression_target_resolution_v0_1.tsv` | HPA 基因层表达投影，禁止 isoform 推断 |
| `complex_target_components_v0_4.tsv.gz` | 带 canonical/isoform 外键的复合物组件表 |
| `external_complex_protein_entity_v0_2.tsv.gz` | 非核心但被复合物引用的外部蛋白实体 |
| `external_complex_protein_review_v0_2.tsv` | 8 个需审计的旧/删除 accession |
| `website_entity_search_index_v0_1.tsv.gz` | 四类实体的统一检索索引 |
| `target_resolution_policy_v0_1.tsv` | 机器可读的目标解析政策 |

## 数据来源与冻结信息

- UniProtKB REST API，Homo sapiens reference proteome `UP000005640`。
- UniProt release：`2026_02`，发布日期：`10-June-2026`。
- 批量 FASTA、响应头、显式 isoform 逐条请求清单和 SHA-256 均保存在 `raw/uniprot/`。
- 168 个复合物显式 isoform accession 已逐个核验，168/168 返回有效 FASTA。

## 复现

脚本均通过命令行参数接收输入路径，没有把 C 盘或 D 盘绝对路径写进生产逻辑。建议按以下顺序执行：

1. `fetch_explicit_isoforms.py`：核验复合物显式 isoform。
2. `build_identity_layer.py`：构建核心 gene/canonical/isoform 与证据桥接。
3. `extend_complex_identity_v2.py`（由 compatibility runner 调用）：补充外部复合物蛋白和统一检索索引。
4. `fetch_missing_external_canonical.py`：核验参考蛋白组之外的外部组件 accession。
5. `finalize_identity_layer_v3.py`：主外键、数量、序列和检索实体最终 QA。

运行环境只需 Python 3 标准库；不依赖 pandas。命令实例见 `REPRODUCE.md`。

## 发布边界

这是 V6.3 candidate 身份模块，不表示冻结完整 V6.3 数据库。建议在网站和论文中明确：

- 10,997 是核心膜蛋白统计口径。
- 外部复合物亚基属于上下文实体，不是新增膜蛋白。
- 结合证据中的 canonical 关联不能自动解释为“canonical isoform 特异性”。
- 8 个旧/删除外部 accession 保留复核，不静默删除，也不影响核心膜蛋白和正式结合证据统计。
