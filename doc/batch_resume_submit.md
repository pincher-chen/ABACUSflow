# 批量续算与自动提交集成方案

## 概述

批量续算和自动提交系统现已完全集成，支持大规模续算作业的智能管理和自动提交。

## 🎯 完整工作流

```
1. 批量生成续算脚本
   ↓
2. 自动提交到Slurm
   ↓
3. 队列管理和资源控制
   ↓
4. 监控提交进度
```

## 📋 使用方法

### 方案1: 一键批量续算+自动提交（推荐）⭐

```bash
# 扫描所有需要续算的作业，生成脚本并自动提交
python abacus.py batch-resume batch_work/ --auto-submit
```

**特点**:
- ✅ 自动扫描所有失败/未完成的作业
- ✅ 生成续算脚本 (*_resume.sh)
- ✅ 自动提交到Slurm队列
- ✅ 使用队列管理系统
- ✅ 遵守资源限制（TOTAL_NODE, MAX_JOBS）
- ✅ 间隔提交，降低负载

### 方案2: 分步执行

#### 步骤1: 批量生成续算脚本

```bash
# 生成所有需要续算的脚本
python abacus.py batch-resume batch_work/
```

输出示例:
```
======================================================================
📦 Batch Resume - Found 100 jobs
======================================================================

📊 Summary:
  Total jobs:        100
  Need resume:       35 ✨
  Already complete:  60 ✅
  No stat.log:       5 ⚠️

🔄 Jobs needing resume:
    1. hmat_0
    2. hmat_12
    3. hmat_107
    ...

✨ Resume scripts generated:
   Success: 35
   Failed:  0
```

#### 步骤2: 自动提交续算脚本

```bash
# 使用自动提交系统提交所有续算脚本
python submit_jobs.py batch_work/ --resume
```

输出示例:
```
======================================================================
Slurm 自动提交器启动
======================================================================
工作目录: batch_work/
分区: deimos
节点限制: 10
任务限制: 1000
提交间隔: 0.5秒
模式: 续算脚本提交 (*_resume.sh)
======================================================================
扫描作业脚本: batch_work/
  包含续算脚本 (*_resume.sh)
找到 35 个脚本文件
队列中共有 35 个作业待提交
作业已提交: hmat_0 (ID: 12345)
作业已提交: hmat_12 (ID: 12346)
...
```

### 方案3: 手动提交（灵活控制）

#### 提交单个作业

```bash
cd batch_work/hmat_0
yhbatch Test_spin_resume.sh  # 或其他 *_resume.sh
```

#### 批量手动提交

```bash
cd batch_work/
for dir in */; do
    (cd "$dir" && [ -f *_resume.sh ] && yhbatch *_resume.sh)
done
```

## 🔧 高级选项

### 1. 查看哪些作业需要续算（不生成脚本）

```bash
python abacus.py batch-resume batch_work/ --dry-run
```

### 2. 从指定阶段开始续算

```bash
# 所有作业从 Relax 阶段开始续算
python abacus.py batch-resume batch_work/ --from-stage Relax --auto-submit
```

### 3. 自定义资源限制

编辑 `config/condor.ini`:

```ini
[ALLOW]
TOTAL_NODE = 20          # 增加节点限制
MAX_JOBS = 2000          # 增加任务数限制
INTERVAL_TIME = 1.0      # 调整提交间隔
CONFIG_RELOAD_INTERVAL = 600  # 配置热重载间隔
```

**注意**: 配置会自动重载，无需重启程序！

### 4. 不备份/不清理

```bash
# 不备份原始文件
python abacus.py batch-resume batch_work/ --no-backup --auto-submit

# 不清理失败的计算文件
python abacus.py batch-resume batch_work/ --no-clean --auto-submit
```

## 📊 实际案例

### 案例1: 100个作业的批量续算

```bash
# 1. 查看需要续算的作业数量
python abacus.py batch-resume batch_work/ --dry-run

# 输出:
# 📊 Summary:
#   Total jobs:        100
#   Need resume:       40 ✨
#   Already complete:  58 ✅
#   No stat.log:       2 ⚠️

# 2. 生成脚本并自动提交
python abacus.py batch-resume batch_work/ --auto-submit

# 3. 监控提交进度
tail -f logs/submit_*.log

# 4. 查看队列状态
squeue -u $USER
```

### 案例2: 部分作业从特定阶段续算

```bash
# 场景: Relax 阶段普遍失败，需要从 Relax 重新开始
python abacus.py batch-resume batch_work/ --from-stage Relax --auto-submit
```

### 案例3: 动态调整资源配额

```bash
# 1. 开始提交（初始配额：10节点）
python submit_jobs.py batch_work/ --resume

# 2. 运行中，配额增加到20节点
# 直接修改 config/condor.ini:
vim config/condor.ini
# TOTAL_NODE = 20

# 3. 系统会在10分钟内（CONFIG_RELOAD_INTERVAL）自动检测并应用新配置
# 无需停止或重启！
```

日志输出:
```
2025-12-26 10:10:00 - AutoSubmitter - INFO - 检测到配置文件更新，重新加载配置:
2025-12-26 10:10:00 - AutoSubmitter - INFO -   total_node: 10 -> 20
2025-12-26 10:10:00 - AutoSubmitter - INFO - ✅ 配置已更新并应用
2025-12-26 10:10:00 - AutoSubmitter - INFO -    当前节点限制: 20
```

## 🔍 状态监控

### 实时监控

```bash
# 终端1: 查看提交日志
tail -f logs/submit_*.log

# 终端2: 查看队列状态
watch -n 5 'squeue -u $USER'

# 终端3: 统计任务数
watch -n 10 'echo "Total: $(squeue -u $USER | tail -n +2 | wc -l)"; echo "Running: $(squeue -u $USER -t RUNNING | tail -n +2 | wc -l)"'
```

### 队列状态查询

```bash
# 查看所有作业
squeue -u $USER

# 查看运行中作业
squeue -u $USER -t RUNNING

# 统计节点使用
squeue -u $USER -h -o '%D' | awk '{sum+=$1} END {print "Using nodes:", sum}'
```

## ⚙️ 配置说明

### condor.ini 配置项

| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `TOTAL_NODE` | 总节点数限制 | 10-50 |
| `MAX_JOBS` | 最大任务数 | 500-5000 |
| `INTERVAL_TIME` | 提交间隔(秒) | 0.5-2.0 |
| `CHECK_TIME` | 队列检查间隔(秒) | 60-120 |
| `CONFIG_RELOAD_INTERVAL` | 配置重载间隔(秒) | 600 (10分钟) |

### 配置热重载

系统会定期检查配置文件是否有更新：

- 默认每10分钟检查一次
- 如果检测到变化，自动重新加载
- 无需停止正在运行的提交任务
- 适合根据集群负载动态调整配额

## 🎨 工作流对比

### 传统方式（手动）

```bash
# 1. 进入每个目录
cd batch_work/hmat_0
yhbatch hmat_0_resume.sh

cd ../hmat_12
yhbatch hmat_12_resume.sh

# ... 重复100次 😫
```

**缺点**:
- ❌ 手动操作，效率低
- ❌ 容易出错
- ❌ 无资源控制
- ❌ 可能压垮调度器

### 新方式（自动化）

```bash
# 一条命令完成所有操作
python abacus.py batch-resume batch_work/ --auto-submit
```

**优点**:
- ✅ 完全自动化
- ✅ 智能资源管理
- ✅ 间隔提交，保护调度器
- ✅ 详细日志记录
- ✅ 支持中断续传
- ✅ 配置热重载

## 🚨 常见问题

### Q1: 如何停止自动提交？

**A**: 按 `Ctrl+C` 安全停止，已提交的作业不受影响。

### Q2: 中断后如何继续？

**A**: 直接重新运行相同命令：

```bash
python submit_jobs.py batch_work/ --resume
```

系统会自动跳过已提交的作业。

### Q3: 如何调整提交速度？

**A**: 修改 `config/condor.ini`:

```ini
# 加快提交
INTERVAL_TIME = 0.1

# 放慢提交（更安全）
INTERVAL_TIME = 2.0
```

### Q4: 配额变化了怎么办？

**A**: 直接修改配置文件，系统会自动应用：

```bash
vim config/condor.ini
# 修改 TOTAL_NODE 和 MAX_JOBS
# 保存后，系统会在10分钟内自动检测并应用
```

### Q5: 如何只续算特定阶段？

**A**: 使用 `--from-stage` 选项：

```bash
python abacus.py batch-resume batch_work/ --from-stage Scf --auto-submit
```

## 📈 性能优化

### 大规模续算（1000+作业）

```bash
# 1. 增加资源限制
vim config/condor.ini
# TOTAL_NODE = 50
# MAX_JOBS = 5000

# 2. 调整提交间隔
# INTERVAL_TIME = 1.0
# CHECK_TIME = 120

# 3. 后台运行
nohup python abacus.py batch-resume batch_work/ --auto-submit > resume.out 2>&1 &

# 4. 监控进度
tail -f resume.out
tail -f logs/submit_*.log
```

## 🔗 相关文档

- [自动提交系统详细文档](auto_submit.md)
- [快速开始指南](../QUICK_START_SUBMIT.md)
- [设计方案总结](../DESIGN_SUMMARY.md)
- [续算工作流文档](restart_workflow.md)

## 📝 总结

批量续算与自动提交的集成提供了：

1. **完全自动化**: 一条命令完成所有操作
2. **智能管理**: 资源监控、队列管理、间隔提交
3. **灵活配置**: 支持热重载，动态调整
4. **可靠稳定**: 自动重试、详细日志、中断续传
5. **大规模支持**: 轻松管理1000+作业

通过这个集成方案，您可以：
- 🚀 大幅提升批量续算效率
- 🎯 精确控制资源使用
- 📊 实时监控提交进度
- 🔧 灵活调整配置参数
- ✅ 确保不超出集群限制

享受高效的批量续算体验！🎉


