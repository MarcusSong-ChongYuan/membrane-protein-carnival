# MemPro 项目完整交接包｜请从这里开始

交接日期：2026-08-17  
当前正式数据库：MemPro V7.2  
数据库状态：FINAL_FROZEN / QA PASS  
当前Docking状态：正式生产尚未启动；安全可运行子集624条，最新收敛门控未通过。

## 阅读顺序

1. 00_START_HERE/01_PROJECT_STATE.md
2. 00_START_HERE/02_END_TO_END_WORKFLOW.md
3. 00_START_HERE/03_DECISIONS_AND_RULES.md
4. 00_START_HERE/04_DOCKING_HANDOFF.md
5. 00_START_HERE/05_PATH_MAP.md
6. 00_START_HERE/06_PROMPT_FOR_NEW_CODEX_AGENT.md
7. 02_WEB_PLANNING/MemPro_V7.2_WEB_DATABASE_PLAN_CN.md

## 目录说明

- 01_V7.2_FROZEN_DATABASE：完整V7.2正式冻结数据库，44个只读文件。
- 02_WEB_PLANNING：网站实体、层次、36表字段目录、完整筛选枚举和五轴标签。
- 03_FIGURES_AND_PPT：历次M1–M8、PPT、caption、QA和最终展示材料。
- 04_CODE：D:\7.22与D:\finale中的代码/Markdown镜像，当前构建、绘图与Docking修复代码。
- 05_DOCKING：A100当前管线、manifest、日志、benchmark结果和本地补丁。
- 06_QA_AND_RELEASE_HISTORY：旧交接说明、复现代码、质量控制、manifest和hash。
- 07_LOCAL_FULL_WORKSPACE_LINKS：指向D:\finale和D:\7.22的目录链接；只在当前电脑有效。
- 08_MANIFESTS：本交接包清单、哈希和外部路径目录。

## 最高优先级规则

- 不得直接修改01_V7.2_FROZEN_DATABASE。
- 新变化只能建立V7.2.1 candidate或V7.3 candidate。
- E1/E2是默认确认核心；E3保留可检索候选。
- 无结合残基不影响互作收录。
- ISOFORM_UNSPECIFIED保留canonical指向，不删除。
- 复合物统一发布，但component-context-only不得冒充complex-specific binding。
- distinct contributing databases不等于independent experiments。
- Docking不得静默删除口袋内UNK/UNX、金属、辅因子或非标准关键原子。
