# Grid Box 残基映射问题 — 现状和策略

## 问题

Trial manifest 7227 个 task 中，881 个（12%）grid box 计算失败（`no_ca_match`），182 个（2.5%）无残基数据。原因：manifest 中部分残基用 UniProt 编号，PDB 文件用 author 编号，数字不一致。

例如：
```
DOCK-00052: pdb=7r8d chain=B  residues=TYR598;TYR604;VAL654 → 0/3 匹配
```
PDB 7r8d 的 DBREF 显示：PDB 5-300 → UniProt 20-315，偏移 ~15。但残基 598/604/654 远超此范围，可能是另一个 UniProt 区段。

## 约束条件

- A100 无外网，SIFTS 不可用
- Trial 已有 6164/7227（85%）直接匹配成功——说明大部分 manifest 残基编号本身就是 PDB numbering

## 当前策略（V3.1，待 A100 跑结果）

**三级 fallback，不破坏已有正确结果：**

1. **直接匹配**（优先）：用 manifest 残基编号直接查 PDB ATOM 行，按 chain→(resname, resnum, icode) 匹配
2. **DBREF 逐残基映射**（仅当直接匹配率 < 50%）：从 PDB 文件 `DBREF` 行提取逐残基 UniProt→PDB 映射，每个 DBREF segment 独立处理，不跨 segment 外推
3. **跨 chain 兜底**（仅当指定 chain 无匹配）：跨所有 chain 搜索，记录实际匹配的 chain

置信度分层：
- `HIGH`：直接匹配或 DBREF 映射，≥80% 残基有坐标
- `MEDIUM`：50-79% 覆盖，或跨 chain 匹配成功
- `REVIEW`：<50% 覆盖，或只有残基名匹配（无数值对应）
- `FAIL`：无 PDB/无残基/chain 不存在/0 残基有坐标

## 待确认

- DBREF 覆盖率：7865 个 chain 有 DBREF，6233/7227 task 可通过 DBREF 映射
- 但 DBREF 只覆盖 PDB 中有序坐标的区段，UniProt 编号 598 可能根本不在 PDB 结构里（是预测位点）
- 这种情况下标记为 REVIEW，等网络恢复后用 SIFTS 精确定位

## 其他已修复项

| 问题 | 状态 |
|------|------|
| parse_pdbqt.cpp:71 断言 | ✅ 已解决（删除 UNK/UNX 空白 atom type） |
| Task 数量对账 | ✅ 7227 tasks = 7227 configs |
| UNK/UNX 审计 | ✅ 记录坐标+距口袋距离，3 个原子在 11-14Å 需复核 |
| 测试脚本 pipefail | ✅ 检查退出码+输出文件+pose 数 |
| 废弃脚本 | ✅ 移入 deprecated/ |
| GPU | ❌ ESMFold 占满，等释放 |

## 开放问题

1. **SIFTS 映射**（长期）：有网络后做准确的 UniProt→PDB 逐残基映射
2. **mmCIF 保留**：PDB .cif 转 .pdb 时丢失 `label_seq_id`，应保留原始 mmCIF
3. **共晶配体 box**：优先使用原配体重原子中心，当前只用结合残基
4. **受体/配体质子化**：obabel 的 `-h` 不处理组氨酸互变异构、pH 依赖质子化
