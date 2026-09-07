# MemPro Docking重新评分（2026-08-11）

输入：docking_manifest.tsv（SHA256见SOURCE_MANIFEST.txt）。本次实际重评分421,926条记录；旧报告中的421,924不与本结果混用。

## 推荐决策
- A：P1，2,737 pairs，极严范围。
- B：P1+P2，10,440 pairs，高科学优先级。
- C：P1+P2，加P3且准备度为R1/R2，共13,606 pairs；推荐。
- D：C再加入批准药/临床候选/化学探针/内源配体且定量≤1 μM的P3，共14,404 pairs。

科学优先级P与Docking准备度R分开计算。评分没有使用Lipinski作为硬筛选。完整评分规则见Excel工作簿和rescore_summary.json。

## 文件
- all_pairs_rescored_421926.tsv：全部pair及新评分。
- A/B/C/D_*.tsv：四个可选范围。
- MemPro_Docking重新评分与范围审阅.xlsx：人类审阅工作簿。
- rescore_summary.json、scope_summary.json、scope_options.json：机器可读统计。
- SHA256SUMS.txt：文件完整性校验。

注意：选中pair不等于全部已经可运行。推荐C中10,735对仍缺pair特异结合残基，1,121对无PDB；下一阶段需构建结构—链—口袋分配。
