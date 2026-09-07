# A100 Docking Setup Guide — 操作记录与坑

## 环境信息

| 项目 | 详情 |
|------|------|
| 服务器 | `csong@star` |
| GPU | 2× A100 80GB, CUDA 12.4, Driver 550.163.01 |
| 操作系统 | Ubuntu (Linux) |
| 用户权限 | 无 sudo |
| 工作目录 | `~/docking/docking_package_trial/` |
| 本地 Windows 目录 | `D:\finale\negative_upgrade\docking_package_trial\` |

## 环境变量（关键！）

每次 `conda activate docking` 会自动设置：

```bash
export LD_LIBRARY_PATH=/home/csong/miniconda3/envs/docking/lib:$LD_LIBRARY_PATH
```

已写入 `/home/csong/miniconda3/envs/docking/etc/conda/activate.d/env_vars.sh`

---

## 已完成：Vina-GPU 编译（最坑）

### Vina-GPU 源码位置
`~/Vina-GPU/`

### 编译成功的 Makefile 关键配置
```
BOOST_LIB_PATH=/home/csong/miniconda3/envs/docking
BOOST_INC_PATH=-I$(BOOST_LIB_PATH)/include
LIB1=-lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -lOpenCL
LIB_PATH=-L$(BOOST_LIB_PATH)/lib -L/usr/local/cuda/lib64
```

### 遇到的坑及解决方案

#### 坑1：conda 未安装
- 下载 Miniconda3：`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
- 安装（无需 sudo）：`bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3`

#### 坑2：conda ToS 未接受
- `conda tos accept --override-channels`

#### 坑3：adfr-suite 在 conda-forge 不可用
- 放弃 ADFR，改用 **meeko + obabel** 做蛋白/配体准备

#### 坑4：Boost 版本不兼容（最大坑）
- conda 只能装 Boost 1.85，Vina-GPU 源码需要 1.77 API
- 服务器网速极慢（~9KB/s），无法下载 Boost 1.77 源码
- **最终方案：用 Boost 1.85 + 手动打补丁**

#### 补丁清单：

**a) 头路径修复**：
```bash
sed -i "s|BOOST_INC_PATH=-I\$(BOOST_LIB_PATH) -I\$(BOOST_LIB_PATH)/boost|BOOST_INC_PATH=-I\$(BOOST_LIB_PATH)/include|" Makefile
```

**b) `boost/filesystem/convenience.hpp` 缺失修复**（Boost 1.85 已移除此文件）：
```bash
cat > /home/csong/miniconda3/envs/docking/include/boost/filesystem/convenience.hpp << 'EOF'
#ifndef BOOST_FILESYSTEM_CONVENIENCE_HPP
#define BOOST_FILESYSTEM_CONVENIENCE_HPP

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

namespace boost {
namespace filesystem {

inline std::string basename(const path& p) {
    return p.filename().string();
}

} // namespace filesystem
} // namespace boost

#endif
EOF
```

**c) `boost/progress.hpp` / `boost/timer.hpp` 被标为 `#error`**（Boost 1.85 deprecated → error）：
```bash
sed -i 's/-DOPENCL_3_0/-DBOOST_TIMER_ENABLE_DEPRECATED -DOPENCL_3_0/' Makefile
```

#### 最终编译命令
```bash
make clean 2>/dev/null; make
```

产物：`~/Vina-GPU/Vina-GPU`（922KB）

---

## 已完成：脚本修改

所有脚本在 `~/docking/docking_package_trial/` 下。

### 所有脚本 — `conda activate` 在 bash 脚本里不生效

脚本中用 `conda activate docking` 在非交互式 shell（`bash xxx.sh`）里会报错：
```
CondaError: Run 'conda init' before 'conda activate'
```

**修复**：在每个脚本里把 `conda activate docking` 前面加 source：
```bash
sed -i 's/conda activate docking/source ~\/miniconda3\/etc\/profile.d\/conda.sh\nconda activate docking/' 02_prepare_proteins.sh 03_prepare_ligands.sh 05_run_docking.sh
```

修改后脚本里变成：
```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate docking
```

### 02_prepare_proteins.sh — ADFR → obabel
原脚本用 `prepare_receptor`（ADFR，不可用），改为 OpenBabel：
```bash
obabel "$pdb_file" -O "$out" -h 2>/dev/null
```
`-h` = 加氢原子，输出 `receptors_pdbqt/*.pdbqt`

### 03_prepare_ligands.sh — 修复 obabel 语法
原脚本用的管道语法在 conda obabel 上不工作，改为：
```bash
obabel -:"$smiles" -opdbqt -O "$out" --gen3d 2>/dev/null
```
`-:"$smiles"` = 从 SMILES 字符串读取，`--gen3d` = 生成 3D 构象

### 05_run_docking.sh — 修正 Vina-GPU 路径和运行目录
原脚本用 `vina_gpu`（未安装的 PATH 命令），改为绝对路径。**必须从 `~/Vina-GPU/` 目录运行**（因为 OpenCL kernel 文件在 `./OpenCL/src/kernels/`），通过 symlink 引用数据目录。

必须在运行前创建 symlink：
```bash
cd ~/Vina-GPU
ln -sf ~/docking/docking_package_trial/receptors_pdbqt_clean receptors_pdbqt
ln -sf ~/docking/docking_package_trial/ligands_pdbqt ligands_pdbqt
ln -sf ~/docking/docking_package_trial/results results
```

**注意**：使用 `receptors_pdbqt_clean/` 而非 `receptors_pdbqt/`（UNK/UNX 已删除）。

### 01_download_pdbs.sh — CRLF 修复（坑）
Windows 传过来的文件有 CRLF（`\r\n`），导致 PDB ID 末尾带 `\r`，URL 无效，3226 个全失败。
```bash
sed -i 's/\r$//' pdb_download_list.txt 01_download_pdbs.sh docking_manifest.tsv ligands.smi *.sh *.py
```

---

## 已完成：下载结果

| 统计项 | 数值 |
|--------|------|
| 总 PDB ID 数 | 3,226 |
| 下载成功 | 2,896 |
| 下载失败 | 330 |

330 个 PDB ID 在当次下载中未成功获取（RCSB/PDBe 均返回 404），涉及 794 个 task。404 的原因可能是：URL 格式不兼容、结构当时尚未发布（如 10JT 于 2026 年 3 月才发布）、临时网络/API 错误、仅支持 `.cif` 格式而不支持旧版 `.pdb`、少部分才是真正撤销/替换的 ID。**不能统一标记为废弃结构，需重新审计。**

---

## 流水线 5 步（当前进度）

所有命令在 `~/docking/docking_package_trial/` 执行，确保 `conda activate docking` 生效。

| 步骤 | 脚本 | 作用 | 输入 | 输出 | 预计时间 | 状态 |
|------|------|------|------|------|----------|------|
| 1 | `bash 01_download_pdbs.sh` | 从 RCSB 下载蛋白结构 | PDB ID 列表 | `pdb_structures/*.pdb` | 30-60min | ✅ 完成（2896/3226） |
| 2 | `bash 02_prepare_proteins.sh` | PDB 加氢 → PDBQT | `pdb_structures/*.pdb` | `receptors_pdbqt/*.pdbqt` | ~30min | 📍 下一步 |
| 3 | `bash 03_prepare_ligands.sh` | SMILES → 3D PDBQT | `ligands.smi` | `ligands_pdbqt/*.pdbqt` | ~20min | ⏳ |
| 4 | `bash 04_compute_grids.sh` | 残基坐标 → Grid box | 残基列表 + PDB 文件 | 更新 `vina_configs/*.conf` | ~5min | ⏳ |
| 5 | `bash 05_run_docking.sh` | GPU 对接 | 蛋白+配体+config | `results/*_out.pdbqt` | 10-20h | ⏳ |

步骤 1 的后台运行方式：
```bash
nohup bash 01_download_pdbs.sh > logs/download.log 2>&1 &
```

步骤 5 同样用 nohup 跑：
```bash
nohup bash 05_run_docking.sh > logs/docking.log 2>&1 &
```

---

## 关键文件和目录结构

```
~/docking/docking_package_trial/
├── docking_manifest.tsv       # 7227 个 task 的主表（TSV，CRLF 已清理）
├── ligands.smi                # 547 个配体的 SMILES（CRLF 已清理）
├── pdb_download_list.txt      # 3226 个 PDB ID（CRLF 已清理）
├── vina_configs/              # 7227 个 .conf 文件
├── compute_grid_boxes.py      # 步骤4用：解析残基 → 计算包围盒 → 填入 config
├── pdb_structures/            # 步骤1产物：2896 个 .pdb 文件
├── receptors_pdbqt/           # 步骤2产物（原始，含 UNK/UNX）
├── receptors_pdbqt_clean/     # 步骤2产物（干净，UNK/UNX 已删除）← Vina-GPU 用这个
├── ligands_pdbqt/             # 步骤3产物
├── results/                   # 步骤5产物：对接结果
├── logs/                      # 日志
│   ├── receptor_cleanup.tsv   # UNK/UNX 删除记录
│   └── failed.txt             # 失败 task 列表
└── receptor_cleanup_logs/     # UNK/UNX 详细审计记录
```

## Vina Config 文件格式（示例）

**注意**：Vina-GPU 不支持 `exhaustiveness` 和 `log =` 参数，所有 config 已通过 `fix_all_configs.py` 移除。

**注意**：路径是相对于 `~/Vina-GPU/` 的（因为必须从该目录执行），需在 Vina-GPU 目录建立 symlink：
```bash
cd ~/Vina-GPU
ln -sf ~/docking/docking_package_trial/receptors_pdbqt_clean receptors_pdbqt
ln -sf ~/docking/docking_package_trial/ligands_pdbqt ligands_pdbqt
ln -sf ~/docking/docking_package_trial/results results
```

```
receptor = receptors_pdbqt/4wsq.pdbqt
ligand = ligands_pdbqt/HMPD-CMPD-0109985.pdbqt
out = results/DOCK-00001_out.pdbqt
center_x = 41.209
center_y = 33.564
center_z = 45.810
size_x = 36.3
size_y = 39.3
size_z = 45.6
num_modes = 9
```

---

## Conda 环境详情

```bash
conda activate docking
```

已安装的包：
- python=3.11
- openbabel (obabel 命令)
- curl (已自带)
- boost-cpp=1.85
- meeko (pip 安装)

编译依赖（编译 Vina-GPU 时需要，已在 conda 环境里）：
- boost 1.85（program_options, system, filesystem, thread）
- OpenCL（NVIDIA CUDA 自带）

---

## parse_pdbqt.cpp:71 断言崩溃（已解决）

### 症状
```
Reading input ... Vina-GPU: ./lib/parse_pdbqt.cpp:71: std::string omit_whitespace(const std::string&, sz, sz): Assertion `i-1 < str.size()' failed.
Aborted (core dumped)
```
`--randomize_only` 正常，正式 docking 必定崩溃。

### 根因
obabel 转换部分 PDB 结构为 PDBQT 时，遇到 `UNK`（未知残基）或 `UNX`（未知原子），无法分配 AutoDock 原子类型，导致 PDBQT 第 78-79 列（atom type 字段）为空白。

`parse_pdbqt_atom_string()` 调用 `omit_whitespace(str, 78, 79)` 提取原子类型时，因为 78-79 列全是空格，`i` 被递增到超出字符串长度，触发断言。

### 定位过程
1. 交叉测试：示例受体 + 我们的配体 ✅ / 我们的受体 + 示例配体 ❌ → 受体有问题
2. 逐行检查：`awk '/^(ATOM  |HETATM)/ && substr($0,78,2) ~ /^[[:space:]]*$/' 2gao.pdbqt` 发现 2 行 UNK/UNX 原子类型为空
3. 删除后立即通过解析

### 修复方案
**02_prepare_proteins.sh** 增加第二步清理：
```bash
# 删除 UNK/UNX 原子（atom type 为空白的 ATOM/HETATM 行）
awk '/^(ATOM  |HETATM)/ && substr($0,78,2) ~ /^[[:space:]]*$/ {next} {print}' "$in" > "$out"
```

清理输出到 `receptors_pdbqt_clean/`，原始文件保留在 `receptors_pdbqt/`。

### 注意事项
- UNK/UNX 可能是结晶添加剂（可删）、未识别小分子、修饰氨基酸、或结合口袋中的金属/辅因子
- 批量删除无法区分，对于结合位点附近的删除记录需要人工复核
- `logs/receptor_cleanup.tsv` 记录了每个受体被删除的原子数

### 源码层面的 bug
`parse_pdbqt.cpp:61` 存在逻辑错误：
```cpp
if(j < str.size()) j = str.size();  // 应为 j > str.size()
```
当前代码在行长超过 79 时会错误扩大解析范围。虽然行长为 79 时不会触发此 bug，但应当修复后重新编译。

---

## 待解决：330 个下载失败的 PDB → 需重新审计

### 现状
- 330 个 PDB ID 在 CRLF 修复后的第二轮下载中仍然失败
- 涉及 794 个 docking task（~11%）
- `logs/failed_ids.txt` 保存了完整列表

### 可能原因（非全部废弃）
1. **下载脚本限制**：只请求旧版 `.pdb` 格式，部分新结构仅提供 `.cif/mmCIF`
2. **发布时间差**：如 `10JT` 是 2026 年 3 月才发布的 KRAS G12C 结构，当时可能尚未可用；当前 RCSB 页面正常
3. **临时网络/API 抖动**：curl 的超时或 HTTP 状态码误判
4. **真正废弃/替换**：仅少数，需通过 RCSB Status API 逐个确认

### 不要做的
- ❌ 不要假定这些结构可被其他 PDB ID（如 4OBE、6Q0K）替换——即使蛋白相同，配体、链、口袋构象都可能不同，grid box 无法直接沿用
- ❌ 不要直接删除这 794 个 task

### 应该做的
- ✅ 用 `.cif` 格式重新尝试下载（`https://files.rcsb.org/download/${pdb}.cif` 或 `https://www.ebi.ac.uk/pdbe/entry-files/download/${pdb}.cif`）
- ✅ 对明确 `OBSOLETE`/`WITHDRAWN` 的条目，通过 RCSB REST API 查询官方 successor 映射
- ✅ 审计后更新 `01_download_pdbs.sh`，先 `http --head` 检查可用性，再区分 `404`（不存在）和 `5xx`（临时错误）

---

## 检查命令速查

```bash
# 激活环境
conda activate docking

# Vina-GPU 是否可用
~/Vina-GPU/Vina-GPU --help | head -5

# 单任务测试
~/Vina-GPU/Vina-GPU --config vina_configs/DOCK-00001.conf --gpu_id 0

# PDB 下载数量
ls pdb_structures/*.pdb | wc -l

# 蛋白受体数量
ls receptors_pdbqt/*.pdbqt | wc -l

# 配体数量
ls ligands_pdbqt/*.pdbqt | wc -l

# 对接完成数量
ls results/*_out.pdbqt | wc -l

# GPU 状态
nvidia-smi

# 后台进程
ps aux | grep bash
```
