# 06｜给新Codex账号的第一条提示词

请接手MemPro人类膜蛋白—小分子互作数据库项目。

先完整阅读：
1. 00_START_HERE/README_FIRST.md
2. 00_START_HERE/01_PROJECT_STATE.md
3. 00_START_HERE/02_END_TO_END_WORKFLOW.md
4. 00_START_HERE/03_DECISIONS_AND_RULES.md
5. 00_START_HERE/04_DOCKING_HANDOFF.md
6. 00_START_HERE/05_PATH_MAP.md
7. 02_WEB_PLANNING/MemPro_V7.2_WEB_DATABASE_PLAN_CN.md

然后核验：
- 01_V7.2_FROZEN_DATABASE/03_QA/V72_VALIDATION_REPORT.json
- 01_V7.2_FROZEN_DATABASE/02_metadata/V72_MANIFEST.tsv
- 08_MANIFESTS/HANDOFF_SHA256.tsv
- D盘目录链接是否可用
- ssh star是否可连接
- A100上是否有Vina进程

约束：
- MemPro V7.2是只读正式基线，不得直接修改。
- 新工作建立V7.2.1/V7.3 candidate。
- 不把E3算入默认确认核心。
- 不删除无residue互作。
- 不把数据库数量称作独立实验数。
- 不静默删除Docking口袋内未知/非标准关键原子。
- 不绕过失败的PRODUCTION_READY门控。

当前优先级：
1. 决定并实现624×3 seeds正式Docking策略，或重新设计可审计收敛门控。
2. 从V7.2重新计算并重绘投稿级M1–M8。
3. 设计网站MVP数据库导入、搜索、详情页和API。
4. 在固定V7.2统计上撰写manuscript。

任何数字必须从V7.2表重新计算，不得仅复述聊天记录。
