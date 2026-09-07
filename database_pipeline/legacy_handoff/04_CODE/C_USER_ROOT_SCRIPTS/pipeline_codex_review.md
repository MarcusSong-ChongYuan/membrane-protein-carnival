# 150K Pair Docking Pipeline — Codex Review

## 背景
从 ~150K 实验结合数据对中筛选出分子对接任务。A100-80GB 服务器上运行 Vina-GPU。

## Pipeline 概览（6 步）

| 步骤 | 脚本 | 输入 | 输出 | 状态 |
|------|------|------|------|------|
| 1 | `01_download_pdbs.sh` | PDB ID 列表 | `pdb_structures/*.pdb` | trial 通过 |
| 2 | `02_prepare_proteins.sh` | PDB 文件 | `receptors_pdbqt_clean/*.pdbqt` | trial 通过 |
| 3 | `03_prepare_ligands.sh` | SMILES 列表 | `ligands_pdbqt/*.pdbqt` | trial 通过 |
| 4 | `04_compute_grid_boxes.py` | 残基注释 + PDB | `vina_configs/DOCK-*.conf` | trial 通过 |
| 5 | `05_run_docking.sh` | config + 受体 + 配体 | `results/*_out.pdbqt` | 解析已通，等 GPU |
| 6 | `06_test_run.sh` | - | 前5条测试 | 待跑 |

## 已知坑和修复（按时间顺序）

### 1. CRLF 行尾
- **症状**：所有 PDB 下载失败（URL 末尾带 `\r`）
- **修复**：`sed -i 's/\r$//'` 所有文本文件
- **涉及的脚本**：全部

### 2. conda activate 在非交互式脚本不生效
- **症状**：`CondaError: Run 'conda init' before 'conda activate'`
- **修复**：每个脚本前加 `source ~/miniconda3/etc/profile.d/conda.sh`
- **涉及的脚本**：01-06 全部

### 3. ADFR Suite 不可用 → 改用 obabel
- **症状**：`prepare_receptor` 命令不存在
- **修复**：用 `obabel "$pdb" -O "$out" -h` 替代
- **涉及的脚本**：02

### 4. obabel 给受体输出 torsion tree
- **症状**：受体 PDBQT 含 ROOT/ENDBRANCH/BRANCH/TORSDOF 标签
- **修复**：`grep -P '^(ATOM|HETATM|TER)' "$tmp" > "$out"`
- **涉及的脚本**：02

### 5. PDB .cif 格式（非废弃结构）
- **症状**：330 个 PDB 返回 404（仅提供 .cif 格式）
- **修复**：RCSB .pdb 失败后自动尝试 RCSB .cif → PDBe .pdb → PDBe .cif
- **涉及的脚本**：01

### 6. UNK/UNX 原子无 AutoDock 类型 → parse_pdbqt.cpp:71 断言崩溃
- **症状**：`assertion i-1 < str.size() failed`，`--randomize_only` 正常
- **根因**：obabel 对未知残基（UNK/UNX）无法分配 atom type，PDBQT 第 78-79 列为空白，Vina-GPU 解析器 `omit_whitespace(str, 78, 79)` 中 i 被递增越界
- **修复**：awk 删除 atom type 为空白的 ATOM/HETATM 行，输出到 `receptors_pdbqt_clean/`
- **涉及的脚本**：02（新增第二步清理）

### 7. Vina-GPU 必须在源码目录运行
- **症状**：`segfault`（找不到 OpenCL kernel 文件）
- **修复**：`cd ~/Vina-GPU && ./Vina-GPU --config "$conf"`，通过 symlink 引用数据目录
- **涉及的脚本**：05、06

### 8. Config 文件不兼容
- **症状**：Vina-GPU 不支持 `exhaustiveness`、`log =`
- **修复**：`fix_all_configs.py` 批量删除这些行，修正路径 `../receptors_pdbqt/` → `receptors_pdbqt/`
- **涉及的脚本**：04 生成 config，fix_all_configs.py 修正

### 9. UnicodeDecodeError in grid box 计算
- **症状**：CIF 转换的 PDB 含非 UTF-8 字节
- **修复**：`open(pdb_path, errors='ignore')`
- **涉及的脚本**：04

### 10. Vina-GPU 源码 bug
- **位置**：`parse_pdbqt.cpp:61`
- **错误**：`if(j < str.size()) j = str.size();` 应为 `if(j > str.size()) j = str.size();`
- **影响**：行长超 79 时错误扩大解析范围
- **状态**：未修复（不影响当前 79 字符行）

## 脚本文件位置
```
D:\finale\negative_upgrade\docking_package\
├── 01_download_pdbs.sh          # PDB + CIF 下载，三重 fallback
├── 02_prepare_proteins.sh       # obabel 转受体 + UNK/UNX 清理
├── 03_prepare_ligands.sh        # obabel SMILES → 3D PDBQT
├── 04_compute_grid_boxes.py     # 残基坐标 → grid box 参数
├── 05_run_docking.sh            # Vina-GPU 批量对接（resume-safe）
├── 06_test_run.sh               # 前 5 条快速测试
├── fix_all_configs.py           # 批量修正 421K config 文件
├── docking_manifest.tsv         # ~152K 任务主表
├── ligands.smi                  # 配体 SMILES
├── pdb_download_list.txt        # PDB ID 列表
├── vina_configs/                # 421523 个 DOCK-*.conf
├── README_DOCKING.md            # 原始说明
└── A100_docking_setup_guide.md  # 完整操作指南（同目录上级）
```

## Trial vs Full Package 差异

| 项目 | Trial (A100) | Full Package (本地) |
|------|-------------|-------------------|
| 任务数 | 7,227 | ~152,000 |
| 配体数 | 547 | ~10× |
| PDB 数 | 3,226 | ~10× |
| Config 数 | 7,227 | 421,523 |
| 脚本修复 | 全部同步 | 全部同步 |

Trial 是 full package 的子集，修复已全部双向同步。Trial 在 A100 上跑通了步骤 1-4，步骤 5 的解析断言已修复，只差 GPU 空闲即可完成验证。

## 当前阻塞
- GPU 被 ESMFold 进程占满（每卡 78GB/80GB），无法运行 docking
- 无网络，无法安装 CPU 版 Vina 替代验证
- 需要临时 kill 一个 ESMFold 进程或等待其自然结束
