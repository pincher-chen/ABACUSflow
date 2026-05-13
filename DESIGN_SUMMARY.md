# Slurm 自动提交系统 - 设计方案总结

## 概述

本方案实现了一个完整的大规模作业自动提交系统，支持10万+作业智能提交到Slurm调度系统。

## 🎯 设计目标

1. **大规模支持**: 支持10万+作业的提交管理
2. **资源控制**: 避免超出节点数和任务数限制
3. **负载友好**: 间隔提交，降低调度器负载
4. **智能管理**: 自动监控、等待、重试
5. **易于使用**: 简单命令，开箱即用

## 📐 系统架构

### 核心模块

```
┌──────────────────────────────────────────────────────────┐
│                    AutoSubmitter                         │
│                   (主控制器)                              │
└────────┬────────────────────────────┬────────────────────┘
         │                            │
    ┌────▼────────┐              ┌────▼─────────┐
    │  前端模块   │              │  后端模块     │
    └─────────────┘              └──────────────┘
         │                            │
    ┌────▼────────────────────────────▼──────────────┐
    │  submit_jobs.py (便捷入口)                      │
    │  python -m abacus.submit_manager (模块方式)    │
    └─────────────────────────────────────────────────┘
```

### 模块详解

#### 1. ConfigManager (配置管理器)

**功能**: 读取和管理配置参数

**配置项**:
- `PARTITION`: Slurm分区名称
- `TOTAL_NODE`: 总节点数限制 ⭐
- `MAX_JOBS`: 最大任务数限制 ⭐
- `INTERVAL_TIME`: 提交间隔时间(秒)
- `CHECK_TIME`: 队列检查时间(秒)

**代码位置**: `abacus/submit_manager.py` (L38-L82)

#### 2. SlurmMonitor (队列监控器)

**功能**: 
- 执行 `squeue -u $USER` 获取队列状态
- 统计总任务数、运行任务数、使用节点数
- 判断是否可以提交新作业

**核心逻辑**:
```python
# 检查节点限制
if used_nodes + job.nodes > TOTAL_NODE:
    return False, "节点资源不足"

# 检查任务数限制
if total_jobs >= MAX_JOBS:
    return False, "达到最大任务数限制"

return True, "资源充足"
```

**代码位置**: `abacus/submit_manager.py` (L85-L159)

#### 3. ScriptParser (脚本解析器)

**功能**: 从脚本中提取节点数和核心数信息

**支持格式**:
```bash
#SBATCH -N 2              # 方式1
#SBATCH --nodes=4         # 方式2
mpirun -np 16             # 方式3
```

**解析优先级**: SBATCH -N > SBATCH -n > mpirun -np > 默认值

**代码位置**: `abacus/submit_manager.py` (L162-L210)

#### 4. JobQueue (作业队列管理器)

**功能**: 
- 使用 `deque` 实现高效队列
- 管理待提交、已提交、失败作业
- 支持10万+作业

**内存占用**: 
- 单个JobInfo: ~200 bytes
- 10万作业: ~20 MB

**代码位置**: `abacus/submit_manager.py` (L213-L255)

#### 5. AutoSubmitter (自动提交引擎)

**功能**: 主控制流程

**核心流程**:

```
1. 扫描脚本
   ↓
2. 解析信息 → 生成JobInfo
   ↓
3. 加入队列
   ↓
4. 主循环:
   ├─ 获取下一个作业
   ├─ 检查资源 (SlurmMonitor)
   ├─ 资源充足? 
   │  ├─ 是: 提交作业 (sbatch)
   │  └─ 否: 放回队列，等待
   ├─ 间隔等待 (INTERVAL_TIME)
   └─ 继续循环
   ↓
5. 输出统计
```

**代码位置**: `abacus/submit_manager.py` (L258-L546)

## 🔧 关键技术实现

### 1. 资源限制机制

**节点数控制**:
```python
# 实时获取当前使用的节点数
used_nodes = sum(nodes for job in running_jobs)

# 预检查
if used_nodes + new_job.nodes <= TOTAL_NODE:
    submit_job(new_job)
else:
    # 资源不足，等待
    queue.push_back(new_job)
    sleep(CHECK_TIME)
```

**任务数控制**:
```python
# 获取队列中的总任务数
total_jobs = len(squeue_output)

if total_jobs < MAX_JOBS:
    submit_job(new_job)
else:
    # 达到上限，等待
    wait_and_retry()
```

### 2. 间隔提交机制

```python
for job in queue:
    submit_job(job)
    time.sleep(INTERVAL_TIME)  # 关键：避免短时间大量提交
```

**作用**:
- 降低调度器瞬时负载
- 避免触发系统保护机制
- 提高提交成功率

### 3. 自动重试机制

```python
if submit_failed:
    job.retry_count += 1
    if job.retry_count < max_retries:
        queue.add_job(job)  # 放回队列重试
    else:
        mark_failed(job)    # 最终失败
```

### 4. 脚本自动扫描

```python
# 递归查找所有.sh脚本
scripts = Path(work_dir).rglob("*.sh")

# 过滤resume脚本
scripts = [s for s in scripts if not s.name.endswith('_resume.sh')]

# 解析每个脚本
for script in scripts:
    nodes, cores = parse_script(script)
    jobs.append(JobInfo(...))
```

### 5. 作业提交命令

```bash
sbatch -N $node_num -p $partition --job-name $name $script_path
```

**参数说明**:
- `-N`: 节点数（从脚本解析或配置获取）
- `-p`: 分区（从condor.ini获取）
- `--job-name`: 作业名称（从目录名提取）

## 📊 性能优化

### 大规模场景优化

**10万作业场景**:

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| INTERVAL_TIME | 1.0 | 每秒1个，平衡速度和负载 |
| CHECK_TIME | 120 | 2分钟检查一次，减少开销 |
| TOTAL_NODE | 50-100 | 根据集群资源调整 |
| MAX_JOBS | 5000 | 避免队列过长 |

**预计提交时间**: 
- 理想速度: 1个/秒 = 10万个需要28小时
- 实际速度: 考虑等待时间，约30-40小时
- **支持中断续传**: 可以分多次提交

### 内存优化

- 使用 `deque` 而非 list，提高队列操作效率
- JobInfo 使用 `@dataclass`，减少内存占用
- 只在内存中保存必要信息，不加载完整脚本内容

### I/O 优化

- 批量扫描脚本，一次性加载
- 日志使用缓冲写入
- squeue 结果缓存，避免频繁调用

## 📝 配置说明

### config/condor.ini 完整配置

```ini
[ENV]
CONDA_PATH = /path/to/conda
CONDA_ENV = dftflow

[ABACUS]
ABACUS_DIR = /path/to/abacus/bin
ABACUS_EXE = abacus
PSEUDO_POTENTIAL_DIR = /path/to/pseudo
ORBITAL_DIR = /path/to/orbital

[STRU]
PATH = 
SUFFIX = *.poscar

[MODULE]
MODULES = module1 module2 module3

[ALLOW]
# Slurm 分区
PARTITION = deimos

# 默认节点数（脚本未指定时使用）
NODES = 1

# 每节点核心数
CORES_PER_NODE = 64

# 总节点数限制（重要！）
# 确保提交的所有作业使用的总节点数不超过此值
TOTAL_NODE = 10

# 最大任务数限制（重要！）
# 确保队列中的总任务数不超过此值
MAX_JOBS = 1000

# 提交间隔时间（秒）
# 每提交一个作业后等待的时间
# 建议值: 0.5-2.0，值越大越安全但提交越慢
INTERVAL_TIME = 0.5

# 队列检查时间（秒）
# 定期检查队列状态的时间间隔
# 建议值: 60-120，频繁检查会增加开销
CHECK_TIME = 80
```

### 参数调优建议

**快速提交（测试）**:
```ini
TOTAL_NODE = 5
MAX_JOBS = 100
INTERVAL_TIME = 0.1
CHECK_TIME = 30
```

**大规模生产**:
```ini
TOTAL_NODE = 50
MAX_JOBS = 5000
INTERVAL_TIME = 1.0
CHECK_TIME = 120
```

**资源受限**:
```ini
TOTAL_NODE = 2
MAX_JOBS = 10
INTERVAL_TIME = 2.0
CHECK_TIME = 30
```

## 🚀 使用流程

### 基本流程

```bash
# 1. 生成作业脚本
python abacus.py workflow InputPoscar/ work_cal7/

# 2. 检查配置
cat config/condor.ini

# 3. 自动提交
python submit_jobs.py work_cal7

# 4. 监控进度
tail -f logs/submit_*.log
squeue -u $USER
```

### 高级用法

```bash
# 使用模块方式
python -m abacus.submit_manager work_cal7

# 自定义配置
python -m abacus.submit_manager work_cal7 --config my_config.ini

# 调整重试次数
python -m abacus.submit_manager work_cal7 --max-retries 5

# 匹配特定脚本
python -m abacus.submit_manager work_cal7 --pattern "hmat_*.sh"
```

## 📂 文件清单

### 新增文件

| 文件 | 说明 | 行数 |
|------|------|------|
| `abacus/submit_manager.py` | 核心模块 | ~550 |
| `submit_jobs.py` | 便捷入口 | ~30 |
| `doc/auto_submit.md` | 详细文档 | ~800 |
| `QUICK_START_SUBMIT.md` | 快速指南 | ~300 |
| `test_submit.py` | 测试脚本 | ~150 |

### 修改文件

| 文件 | 修改内容 |
|------|----------|
| `config/condor.ini` | 添加 MAX_JOBS 配置项 |
| `README.md` | 添加自动提交系统说明 |

## 🔍 错误处理

### 常见错误及解决方案

| 错误 | 原因 | 解决方案 |
|------|------|----------|
| 节点资源不足 | 当前使用节点数达到上限 | 自动等待，无需处理 |
| 达到最大任务数 | 队列任务过多 | 自动等待，或增加 MAX_JOBS |
| 提交命令超时 | sbatch执行慢 | 自动重试，或检查网络 |
| 无法解析作业ID | sbatch输出格式异常 | 检查Slurm版本 |
| 脚本解析失败 | 脚本格式问题 | 手动添加SBATCH指令 |

### 重试机制

```python
# 默认重试3次
max_retries = 3

# 每次重试间隔
retry_interval = INTERVAL_TIME

# 最终失败后记录到日志
failed_jobs.append(job)
```

## 📈 监控和日志

### 日志系统

**位置**: `logs/submit_YYYYMMDD_HHMMSS.log`

**级别**:
- INFO: 重要操作
- DEBUG: 详细信息
- WARNING: 警告（重试）
- ERROR: 错误

**示例**:
```
2025-12-25 10:00:00 - AutoSubmitter - INFO - Slurm 自动提交器启动
2025-12-25 10:00:01 - AutoSubmitter - INFO - 找到 1000 个脚本文件
2025-12-25 10:00:02 - AutoSubmitter - INFO - 作业已提交: hmat_0 (ID: 12345)
2025-12-25 10:01:00 - AutoSubmitter - INFO - 队列状态 - 总任务:20, 运行中:15, 使用节点:8/10
```

### 实时监控

```bash
# 查看提交日志
tail -f logs/submit_*.log

# 查看队列状态
watch -n 5 'squeue -u $USER'

# 统计任务数
watch -n 10 'squeue -u $USER | tail -n +2 | wc -l'

# 统计节点使用
squeue -u $USER -h -o '%D' | awk '{sum+=$1} END {print sum}'
```

## 🎯 设计亮点

### 1. 智能资源管理
- ✅ 实时监控Slurm队列
- ✅ 预检查资源限制
- ✅ 自动等待资源释放

### 2. 高可扩展性
- ✅ 支持10万+作业
- ✅ 内存占用极小（20MB/10万作业）
- ✅ 使用高效数据结构（deque）

### 3. 负载友好
- ✅ 可配置提交间隔
- ✅ 避免瞬时大量提交
- ✅ 降低调度器压力

### 4. 健壮性
- ✅ 自动重试失败提交
- ✅ 详细日志记录
- ✅ 支持中断续传
- ✅ 优雅的错误处理

### 5. 易用性
- ✅ 简单命令行接口
- ✅ 合理的默认配置
- ✅ 详细的文档和示例
- ✅ 测试脚本验证

## 🔮 未来扩展

### 计划功能

1. **作业依赖管理**
   - 支持 `--dependency` 参数
   - 自动处理工作流依赖关系

2. **优先级控制**
   - 不同类型作业不同优先级
   - 支持插队和延迟提交

3. **Web 监控界面**
   - 实时可视化监控
   - 图形化配置管理

4. **邮件通知**
   - 作业完成通知
   - 错误告警

5. **断点续传优化**
   - 记录提交状态到数据库
   - 支持精确续传

### API 接口

```python
from abacus.submit_manager import AutoSubmitter, JobInfo

# 创建提交器
submitter = AutoSubmitter()

# 手动添加作业
job = JobInfo(
    script_path="test.sh",
    job_dir=".",
    job_name="my_job",
    nodes=2,
    cores=32,
    partition="deimos"
)
submitter.queue.add_job(job)

# 批量提交
submitter.run("work_dir", max_retries=3)
```

## 📖 参考文档

1. **详细设计文档**: [doc/auto_submit.md](doc/auto_submit.md)
2. **快速开始指南**: [QUICK_START_SUBMIT.md](QUICK_START_SUBMIT.md)
3. **主README**: [README.md](README.md)

## 🎉 总结

本设计方案实现了一个**完整、高效、易用**的大规模作业自动提交系统：

- ✅ **架构清晰**: 模块化设计，职责分明
- ✅ **功能完善**: 资源监控、智能提交、自动重试
- ✅ **性能优异**: 支持10万+作业，内存占用小
- ✅ **易于使用**: 简单命令，开箱即用
- ✅ **文档齐全**: 设计文档、使用指南、测试脚本
- ✅ **经过验证**: 所有模块通过测试

通过合理配置，可以高效管理大规模作业提交，显著提升研究效率！


