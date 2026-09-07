# 01｜当前项目状态

## 数据库冻结结果

正式目录：01_V7.2_FROZEN_DATABASE  
原始位置：D:\finale\22_MemPro_V7.2_final_20260817

核心统计：
- canonical膜相关蛋白：7,800
- E1/E2确认核心：7,374
- E3候选：426
- 正式正证据：942,455
- protein—canonical compound pair：529,168
- 有正式互作的蛋白：2,829
- pair中的canonical小分子：240,055
- compound master总目录：646,700
- compound forms：649,792
- 位点：55,610
- canonical疾病：3,739
- 蛋白—疾病关系：6,753
- 复合物：4,584
- 正式复合物结合证据：24,954
- isoform：9,240
- processed form：8,454
- HPA表达测量：3,046,789
- 亚细胞定位：19,091
- v3.1_mixed_sources剩余：0
- 阻断性QA错误：0

验证报告SHA-256：
5DACD781AC4C8D3EF2A7A6D88270702725536CB6B0E1E871008444DDE04F7C8A

## 已完成

- 10,997个历史蛋白候选在线核实及确定性A/B/C重分。
- 最终7,800个A/B/C蛋白正式纳入；3,197个排除但保留审计。
- canonical蛋白、gene、isoform、processed form和state层。
- canonical compound与form层。
- 正向互作、pair、来源恢复、谱系、去重和外键。
- 位点S1/S2/S3。
- 疾病canonical、xref、hierarchy、治疗领域和解剖映射。
- HPA表达与定位语义。
- 统一复合物实体模块。
- V7.2 manifest、SHA-256、只读冻结。
- 网站数据层次与筛选规划。
- M1–M8既有图和PPT资料整理。

## 尚未完成

- M1–M8全部从V7.2重新计算并重新绘制的投稿终版。
- 网站前端、后端数据库导入和API实现。
- manuscript正式写作与投稿补充材料。
- 真实HPC正式Docking生产。
- 624条安全Docking任务的三随机种子方案尚未批准/运行。
- 13,709条被阻断任务仍需按grid、结构、受体、配体和关键原子原因分批救援。
- 人工抽样验证按用户要求暂时省略，因此论文不能宣称人工准确率≥95%。

## 版本关系

V6.x为早期扩展与标准化版本；V7系列完成公开层收口；V7.1.1是V7.2直接基线；V7.2解决E3默认展示、mixed source、位点分层、canonical fallback和统一复合物发布规则。旧版本只用于追溯，不应覆盖V7.2。
