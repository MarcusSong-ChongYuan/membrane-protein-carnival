# Codex Handoff — 2026-08-10

## 背景

分子对接 pipeline trial 阶段。7227 个 docking tasks（对应 7227 个 Vina configs），425 个受体 PDB 文件已转 PDBQT（UNK/UNX 已清理）。

当前阶段：计算 grid box→生成 Vina config→GPU docking 测试。

## 已完成的工作

### 1. V4.1 grid box 计算 (`04_compute_grid_boxes.py`)

V3.1 有 11 个代码级 bug（你上一轮 review 指出的），V4.1 已全部修复并在 A100 上重跑。

**修复清单：**
| # | Bug | 修复 |
|---|-----|------|
| 1 | SIFTS 未集成 | `load_sifts()` 加载 `sifts_mappings.tsv`，作为第 2 优先映射源 |
| 2 | 无 UniProt accession 验证 | 直接匹配成功时对比 manifest target_uniprot vs DBREF accessions |
| 3 | DBREF1/2 独立解析 | 按 (chain, pdb_begin, pdb_end) 配对：DBREF1 提供 accession，DBREF2 提供 db range |
| 4 | 非等长 segment 外推 | 仅 `pdb_len == db_len` 且无插入码的 segment 才映射 |
| 5 | 跨链 DBREF 借用 | 已删除。同源链匹配走独立 REVIEW 路径 |
| 6 | REVIEW 写 grid 坐标 | REVIEW/FAIL/NO_SITE_DATA → grid 字段全部留空 |
| 7 | mapped=requested 恒成立 | mapped=实际重映射的残基数，observed=PDB 中有坐标的残基数 |
| 8 | 未清旧 grid 字段 | 每个 task 开始时六字段清零 |
| 9 | config 生成无 confidence 检查 | `build_task_manifest.py` 按需检查 |
| 10 | match_residues 计数错误 | 只在 PDB 中找到匹配原子时才 mapped++ |
| 11 | 无 actual_chain 追踪 | 新增 `actual_chain_used` 列 |

**映射优先级：**
1. Direct match（残基编号即 PDB author numbering）
2. SIFTS 逐残基映射（sifts_mappings.tsv，如果存在）
3. DBREF segment-based（等长无插入码的 segment）
4. 同源链候选（REVIEW only，不写 grid）

### 2. V4.1 A100 运行结果（无 SIFTS）

```
Grid boxes with HIGH/MEDIUM: 4097/7227

QC breakdown:
  chain_mismatch                : 1803
  confidence_HIGH               : 3957
  confidence_MEDIUM             : 140
  confidence_REVIEW             : 425
  no_parseable                  : 187
  no_residue_observed           : 715
```

**vs V3.1 对比：**

| 指标 | V3.1 | V4.1 | 解读 |
|------|------|------|------|
| HIGH | 3617 | **3957** | +340 (DBREF1/2 修复 + 列位置修正) |
| MEDIUM | 2130 | **140** | -1990 (跨链/外推不再冒充 MEDIUM) |
| REVIEW | 394 | 425 | +31 |
| chain_mismatch | (混在其他) | **1803** | 新识别类别 |
| no_residue_obs | 899 | **715** | -184 |
| 有效 grid | 5747 | **4097** | 数据质量大幅提升 |

### 3. SIFTS 下载器 (`download_sifts_mappings.py` V2.0)

改用 **SIFTS XML 文件**（EBI FTP，非 PDBe REST API），逐残基 UniProt↔PDB 映射。支持断点续传、指数退避重试、原子写入、版本追踪。

**状态：在 Windows 上启动后中断。** 13871 个 PDB ID（看起来不只是 trial 的，是全量 pipeline 的），按 0.3s/个需要 ~70 分钟。

> ⚠ 需要决定：SIFTS 下载用 trial 的 425 个 PDB 还是全量 13871 个？

### 4. 其他已就绪

- `06_test_run.sh` — 5 任务测试脚本（pipefail、退出码检查、pose 数验证）
- `02_prepare_proteins.sh` — UNK/UNX 清理 + audit 集成
- receptors_pdbqt_clean/ — 425 个清理后的受体
- 7227 个 vina_configs/

---

## 当前问题（需要决策）

### A. 7R8D/7R8E — DBREF 力所不及

```
DOCK-00052: pdb=7r8d chain=B residues=TYR598;TYR604;VAL654
DBREF of 7r8d: PDB 5-300 → UniProt 20-315 (offset ~15)
Residues 598/604/654 >> PDB segment end (300)
```

DBREF 只覆盖有序坐标区段，ABCG1_HUMAN (666 aa) 的 C 端残基不在晶体结构里。**SIFTS 也解决不了——原子根本不存在于 PDB 中。** 这应该是正确行为：标记 FAIL，而不是给错误的 grid box。需要确认逻辑是否正确。

### B. chain_mismatch 1803 — 区分原因

已确认前几个：
- `7fdv ch=B` — PDB 只有 A/D
- `7r8d ch=D` — PDB 只有 A/B
- `3pyy ch=C` — PDB 只有 A/B

但尚未系统分析：
- 其中多少个 PDB 文件根本不存在？
- 多少个 manifest chain 写错（同 UniProt 的其他 chain 存在）→ 可通过同源链 fallback rescue？
- 多少个是因 PDB 是 NMR ensemble（chain 命名不同）？

### C. 尚未验证

1. **grid box 质量抽查** — 只随机看了 3 个 HIGH（2 direct + 1 dbref_segment），肉眼看起来合理，但没有系统验证（box center 是否真的在结合位点中心？margin 10Å 够不够？）
2. **6npe DBREF 映射是否正确** — DOCK-00118 用 dbref_segment 匹配了 7/7 residues，映射逻辑对但没逐残基核对
3. **Vina docking 实测** — `06_test_run.sh` 还沒跑。GPUs 已空闲（2×A100 80GB 全空）

### D. MEDIUM=140 怎么处理

V4.1 只有 140 个 MEDIUM（50-79% observed），比 V3.1 的 2130 更真实。它们：
- 部分残基在 PDB 中有坐标，部分没有（disordered loops）
- 用有坐标的残基算的 centroid 可能偏移
- 是否值得跑 docking？还是等 SIFTS 后再决定？

---

## 建议的下一步

1. **SIFTS 下载** — 只对 trial 的 ~400 PDB 下载（不是 13871），然后上传 A100 重跑 V4.1
2. **chain_mismatch 分类** — 跑一个脚本区分：PDB 缺失 / chain 错误 / 同源链存在
3. **grid box 抽查** — 从 HIGH/MEDIUM/REVIEW 各抽 5 个，PYMOL/ChimeraX 可视化验证
4. **先跑 5 个 HIGH docking 测试** — 即使 SIFTS 没到位，HIGH（100% observed）的 box 是可用的
5. **决定 MEDIUM+REVIEW 处理策略** — 哪些值得跑，哪些直接丢弃

---

## 文件清单

| 文件 | 位置 | 版本 |
|------|------|------|
| `04_compute_grid_boxes.py` | Windows `D:\...\docking_package\` + A100 `~/docking/docking_package_trial/` | V4.1 |
| `download_sifts_mappings.py` | Windows only | V2.0 |
| `docking_manifest_with_grid.tsv` | A100 `~/docking/docking_package_trial/` | V4.1 output |
| `grid_audit_report.tsv` | A100 `~/docking/docking_package_trial/` | V4.1 output |
| `grid_mapping_strategy.md` | Windows | Updated |
| 旧 grid | A100 `*_v3.1_backup.tsv` | V3.1 backup |
