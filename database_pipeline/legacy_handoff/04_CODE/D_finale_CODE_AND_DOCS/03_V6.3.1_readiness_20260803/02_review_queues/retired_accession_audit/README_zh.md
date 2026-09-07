# 8个退役UniProt accession权威核查

核查日期：2026-08-03。V6.3冻结候选保持只读；本目录只提供下一增量版本的身份修正补丁。

## 判定规则

1. 从UniSave提取退役条目的最后版本、序列、蛋白名和基因信息。
2. 通过UniParc定位完全相同序列的现行人类UniProtKB交叉引用。
3. 只有唯一的active、reviewed人类条目才允许高置信替换。
4. 若完全相同序列对应多个基因，则保留`unresolved_gene_product_group`，不依据名称强选。
5. 来源文献与复合物名称冲突时，优先停止默认发布，而不是猜测正确文献。

## 逐条结论

|旧accession|结论|下一版本动作|
|---|---|---|
|D9YZV4|TPM1 isoform 1，与P09493-1序列完全一致|组件改为canonical P09493、isoform P09493-1|
|O60344|与现行reviewed P0DPD8完全一致|组件改为P0DPD8；原两条记录为预测复合物，不进入默认层|
|P01562|与现行reviewed IFNA1 P0DY56完全一致|组件改为P0DY56；作为非膜配体亚基保留|
|P04745|相同序列同时对应AMY1A、AMY1B、AMY1C|不得强选；保留未指定基因产物组，且相关记录均为预测层|
|P0CG05|与现行reviewed IGLC2 P0DOY2完全一致|九个免疫球蛋白复合物组件统一改为P0DOY2|
|P33765|与现行reviewed ADORA3 P0DMS8完全一致|组件改为P0DMS8；相关记录均为预测层|
|Q13748|相同序列同时对应TUBA3C、TUBA3D|不得强选；保留未指定基因产物组，且相关记录均为预测层|
|Q13791|历史40 aa ApoE3片段与P02649近似但并非完全一致；CORUM 1088所列PMID 16764853与该复合物无关|不自动映射；CORUM 1088退出默认层，等待CORUM更正或权威原文确认|

Q13791历史片段与P02649的最佳匹配位于前体序列135–174，存在1个残基差异。CORUM所列PMID 16764853对应EVI5与染色体乘客复合物研究；PMID 16953297的标题与PRNP–ApoE结合主题相符，但在CORUM正式更正前不能擅自替换文献ID。

## 影响

- 5条旧accession得到高置信身份解析。
- 2条变为明确的`unresolved_gene_product_group`，不再被描述为“尚未查到”。
- 1条变为`source_record_conflict`并建议从默认复合物层排除。
- 这8条的专项核查已经完成；更广泛的、原始论文未说明isoform的历史证据仍不能凭算法恢复。

## 主要文件

- `retired_accession_final_adjudication_v631.tsv`：首轮机器判定和最终结论。
- `complex_component_accession_patch_v631.tsv`：下一版本的增量修改说明。
- `current_uniprot_identity_reference.tsv`：现行UniProt条目名称和基因。
- `raw_uniprot/`：UniSave、UniParc和现行UniProtKB原始响应。
