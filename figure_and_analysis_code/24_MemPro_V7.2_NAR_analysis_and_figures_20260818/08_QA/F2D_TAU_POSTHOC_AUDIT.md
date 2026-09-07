# Figure 2d Tau post-hoc audit

## 结论

Figure 2d 实际执行的不是 Dunn test，而是两两双侧 Mann–Whitney U 检验，并对全部 pairwise P 值进行 Benjamini–Hochberg FDR 校正。

因此不应在 caption、methods 或 manuscript 中写成：

`Dunn post-hoc test + BH-FDR`

应写成：

`pairwise two-sided Mann–Whitney U tests with Benjamini–Hochberg FDR correction`

## 代码位置

`01_scripts/01_compute_statistics.py`

函数：`tau_statistics(tau)`

实际代码逻辑：

1. 删除 `tau_consensus_tissue_rna` 或 `primary_membrane_role` 缺失记录；
2. 统计每个膜角色的蛋白数；
3. 仅保留 n ≥ 30 的角色进入总体 Kruskal–Wallis 检验和两两比较；
4. 每个角色对使用 `scipy.stats.mannwhitneyu(..., alternative="two-sided")`；
5. 将所有 pairwise P 值交给内部 `bh_fdr(pvals)` 函数；
6. 输出结果表。

## 结果文件

`03_source_data/F2D_tau_pairwise_statistics.tsv`

该表共 55 行，即 C(11, 2) = 55 个 pairwise membrane-role comparisons，字段为：

`role_a`, `role_b`, `n_a`, `n_b`, `median_a`, `median_b`, `mannwhitney_u`, `p`, `bh_fdr`

角色数为 11。`junction_or_adhesion` 只有 13 个蛋白，因此按 n ≥ 30 规则没有进入检验。

分组计数来自：

`03_source_data/F2D_tissue_specificity_tau.tsv`

该文件当前共有 7,496 条带 Tau 的蛋白记录。进入 n ≥ 30 组间检验的记录数为 7,483 条；13 条 junction_or_adhesion 记录被排除于正式组间检验之外。

## 总体检验

结果保存在：

`02_analysis_data/analysis_results.json`

当前值：

- Kruskal–Wallis H = 457.88552657421826
- P = 4.3424871562855384e-92
- groups_n_ge_30 = 11

Figure 2d 只在图中展示了 Kruskal–Wallis P，没有把 55 个 pairwise 结果直接画入主图；pairwise 结果在上述 TSV 中。

## 需要注意的版本对账问题

`06_captions_methods/FIGURE_PANEL_SOURCE_TRACEABILITY.tsv` 当前把 Figure 2d 分母写成 `7,318 proteins with Tau`，但实际 `F2D_tissue_specificity_tau.tsv` 有 7,496 条记录，绘图脚本也直接读取该文件。因此 7,318 是旧的或未同步的摘要数字，不能继续沿用；caption、source traceability 和正文应统一到实际分析文件的分母。

## 当前没有的内容

在整个 V7.2 投稿图包中没有找到：

- `posthoc_dunn` 或 `post_hoc_dunn` 的代码调用；
- `scikit_posthocs.posthoc_dunn` 调用；
- Dunn 专用结果 TSV/CSV；
- 任何以 Dunn 命名的 Figure 2d 结果表。

所以目前不能声称 Dunn test 已经完成。现有 pairwise 结果是 Mann–Whitney U + BH-FDR。
