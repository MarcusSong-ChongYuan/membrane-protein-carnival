# Grid Box 残基映射 — V4.1 实现

## 映射策略（按优先级，≥50% observed 即停止）

| 优先级 | 方法 | 条件 | 置信度 |
|--------|------|------|--------|
| 1 | 直接匹配 | 残基编号作为 PDB author numbering 直接查找 | ≥80% → HIGH, 50-79% → MEDIUM |
| 2 | SIFTS | sifts_mappings.tsv 逐残基 UniProt→PDB 映射 | ≥80% → HIGH, 50-79% → MEDIUM |
| 3 | DBREF segment | 等长、无插入码的 DBREF segment 线性映射 | ≥80% → HIGH, 50-79% → MEDIUM |
| 4 | 同源链候选 | 同一 UniProt 的其他 PDB chain 匹配 | 始终 REVIEW（不写 grid） |

## V4.1 修复的 11 个 bug（相对 V3.1）

| # | Bug | V3.1 行为 | V4.1 修复 |
|---|-----|----------|----------|
| 1 | SIFTS 未集成 | 代码中有 SIFTS 注释但未加载文件 | `load_sifts()` 解析 `sifts_mappings.tsv`，作为第 2 优先映射 |
| 2 | 无 accession 验证 | 直接匹配不检查 UniProt accession | 直接匹配成功时对比 manifest `target_uniprot` vs DBREF accessions |
| 3 | DBREF1/2 独立解析 | 各自当完整记录用，DBREF2 的 db_acc 读到空白 | DBREF1+DBREF2 按 (chain, pdb_begin, pdb_end) 配对；DBREF1 提供 accession，DBREF2 提供 db range |
| 4 | 非等长段外推 | 所有 DBREF segment 用比例映射/线性外推 | 仅等长（pdb_len == db_len）且无插入码的 segment 才映射 |
| 5 | 跨链 DBREF 借用 | 当前 chain 无 DBREF 时自动借其他 chain | 只查 manifest 指定的 chain；同源链匹配走尝试 #4（REVIEW only） |
| 6 | REVIEW 写 grid | 所有任务都算 box 并写 grid 坐标 | REVIEW/FAIL/NO_SITE_DATA → grid 字段留空 |
| 7 | mapped=requested | mapped 总是等于 requested（有数字就算 mapped） | mapped=实际通过映射方法重映射的残基数；observed=PDB 中有坐标的残基数 |
| 8 | 未清旧 grid | 依赖输出文件已有数据 | 每个 task 开始时清空所有 grid 字段 |
| 9 | config 生成无检查 | 不检查 confidence | 见 `build_task_manifest.py`（checks confidence 列） |
| 10 | try_match 计数 | matched 对找不到原子的残基也+1 | match_residues 只在 PDB 中找到匹配时才 mapped++ |
| 11 | 无实际 chain 追踪 | 不记录实际匹配的 chain | `actual_chain_used` 列记录（尤其同源链场景） |

## 硬规则（不可违反）

- **REVIEW / FAIL / NO_SITE_DATA → grid 字段为空**
- **跨链匹配永远 REVIEW，绝不自动写成 MEDIUM/HIGH**
- **DBREF 非等长 segment 绝不外推**
- **DBREF 跨链借用已删除**
- **mapped ≠ observed**（两个独立列）

## 输出列说明

| 列 | 含义 |
|----|------|
| `requested_residue_count` | manifest 中的残基总数 |
| `mapped_residue_count` | 通过当前映射方法成功重新编号的残基数 |
| `observed_residue_count` | 映射后在 PDB 中有原子坐标的残基数 |
| `grid_confidence` | HIGH / MEDIUM / REVIEW / FAIL / NO_SITE_DATA |
| `mapping_source` | direct / sifts / dbref_segment / homologous_chain / none |
| `actual_chain_used` | 实际用于匹配的 chain（同源链场景 ≠ manifest chain） |

## SIFTS 数据下载

`download_sifts_mappings.py` V2.0：
- 来源：EBI FTP SIFTS XML 文件（`{pdb_id}.xml.gz`）
- 逐残基 UniProt↔PDB 映射
- 断点续传（checkpoint 文件）
- 原子写入 + 版本追踪
- 3 次重试 + 指数退避
- 运行：Windows 有网环境 `python3 download_sifts_mappings.py`
- 输出：`sifts_mappings.tsv` → 上传 A100

## 运行顺序（A100）

```bash
# 1. 上传文件
#    - 04_compute_grid_boxes.py (V4.1)
#    - sifts_mappings.tsv (从 Windows 下载后上传)

# 2. 运行 grid 计算
python3 04_compute_grid_boxes.py
# 输出: docking_manifest_with_grid.tsv + grid_audit_report.tsv

# 3. 验证
# 检查 7R8D/7R8E/3PYY 之前失败的 task 是否被 rescue
grep -E '7R8D|7R8E|3PYY' grid_audit_report.tsv

# 4. 生成 configs（仅 HIGH + MEDIUM）
# 使用 build_task_manifest.py 或直接读取 docking_manifest_with_grid.tsv
```
