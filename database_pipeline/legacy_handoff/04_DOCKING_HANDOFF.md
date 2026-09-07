# 04｜A100 Docking交接

## 连接

SSH别名：star  
主机：10.123.4.78  
用户：csong  
批次目录：
/home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2

## 当前真实状态

正式Docking没有开始。两次修复后benchmark均未通过pose稳定门槛。

14,333条原始任务经最新安全QC：
- PASS receptor/grid层：941
- oversized box：893
- grid QC blocked：12,413
- unsupported critical atom：86
- 加入clean receptor、ligand存在性和≥5 ligand atoms后正式可配置：624
- 阻断：13,709

已修复：
- Open Babel临时文件.pdbqt.tmp无法推断格式，未来使用显式-opdbqt。
- Boost动态库路径固定到docking conda环境。
- 5000/8000 benchmark输出不再互相覆盖。
- 每个任务使用确定性seed。
- 7W7T的BeF3位于口袋，被阻断，未删除。
- 4WNV转换超时部分文件已隔离。
- 缺失配体和缺失clean receptor不生成配置。
- 少于5原子的锌/离子配体不用于常规pose benchmark。

最新同seedbenchmark：
- 12/12有效
- median |score delta| = 0.2 kcal/mol
- median direct RMSD = 3.9374 Å
- stable fraction = 0.50
- 要求0.60
- status=FAIL
- 未创建PRODUCTION_READY
- 未启动正式624条

## 推荐下一步

推荐改为624 tasks × 3 deterministic seeds，thread=8000，共1,872次；保存每个seed结果，选择最佳分数pose，同时报告pose多模态。因为会将计算量扩大3倍，原线程最后要求用户确认后再启动。

禁止直接touch PRODUCTION_READY绕过失败门控。若改变门控标准，必须形成新validation脚本、版本和理由。

## 快速检查命令

ssh star
cd /home/csong/docking/handoff/MemPro_Docking_Dscope_14333_ProductionV3_20260813_upload2
nvidia-smi
cat metadata/thread_benchmark_validation_v4.txt
cat metadata/thread_benchmark_validation_v4.tsv
wc -l metadata/runnable_after_production_qc.tsv
head metadata/blocked_after_production_qc.tsv
ps -u csong -o pid,etime,%cpu,%mem,cmd

## 快照

05_DOCKING/MemPro_A100_state_20260817.tar.gz包含pipeline、metadata、logs、benchmark配置/日志/结果、manifest和当前修复脚本，不包含巨大的PDB/PDBQT实体文件。实体仍在A100原目录。
