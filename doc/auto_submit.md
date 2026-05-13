# Slurm 自动提交系统设计方案

## 概述

本方案实现了大规模作业自动提交到Slurm系统的完整解决方案，支持10万+作业的智能提交管理。

## 核心特性

1. **智能资源监控** - 实时监控Slurm队列状态，避免超出节点和任务限制
2. **内存队列管理** - 高效的内存队列，支持海量作业管理
3. **间隔提交** - 可配置的提交间隔，降低调度系统负载
4. **自动重试** - 提交失败自动重试，提高成功率
5. **详细日志** - 完整的日志记录，便于追踪和调试

## 系统架构

```
┌─────────────────────────────────────────────────────────────┐
│                     AutoSubmitter                           │
│                   (自动提交主控制器)                          │
└─────────────────┬───────────────────────────────────────────┘
                  │
    ┌─────────────┼─────────────┬─────────────────┐
    │             │             │                 │
    ▼             ▼             ▼                 ▼
┌─────────┐  ┌─────────┐  ┌──────────┐    ┌──────────┐
│ Config  │  │  Slurm  │  │   Job    │    │  Script  │
│ Manager │  │ Monitor │  │  Queue   │    │  Parser  │
└─────────┘  └─────────┘  └──────────┘    └──────────┘
     │            │             │               │
     │            │             │               │
     ▼            ▼             ▼               ▼
  condor.ini   squeue      deque队列      解析.sh脚本
```

## 工作流程

```
1. 扫描脚本
   └─> 递归扫描work_dir下所有.sh脚本
       └─> 过滤_resume.sh脚本

2. 解析脚本信息
   └─> 提取节点数、核心数
       └─> 生成JobInfo对象

3. 加入队列
   └─> 添加到内存deque队列
       └─> 队列支持10万+作业

4. 监控资源
   └─> 执行squeue -u $USER
       └─> 统计当前使用的节点数和任务数

5. 检查限制
   └─> 检查是否超出TOTAL_NODE限制
       └─> 检查是否超出MAX_JOBS限制

6. 提交作业
   └─> 构建sbatch命令
       └─> sbatch -N $nodes -p $partition --job-name $name $script
           └─> 提取作业ID

7. 间隔等待
   └─> sleep(INTERVAL_TIME)
       └─> 继续下一个作业

8. 资源不足时
   └─> 作业放回队列
       └─> sleep(CHECK_TIME)
           └─> 重新检查资源
```

## 配置说明

### config/condor.ini 配置项

```ini
[ALLOW]
# Slurm分区名称
PARTITION = deimos

# 默认节点数（脚本中未指定时使用）
NODES = 1

# 每节点核心数
CORES_PER_NODE = 64

# 总节点数限制（最重要）
# 系统会确保提交的作业总节点数不超过此限制
TOTAL_NODE = 10

# 最大任务数限制
# 系统会确保队列中的总任务数不超过此限制
MAX_JOBS = 1000

# 提交间隔时间（秒）
# 每提交一个作业后等待的时间，避免短时间大量提交
INTERVAL_TIME = 0.5

# 队列检查时间（秒）
# 定期检查队列状态的时间间隔
CHECK_TIME = 80
```

### 配置项详解

| 参数 | 说明 | 推荐值 | 作用 |
|------|------|--------|------|
| `TOTAL_NODE` | 总节点数限制 | 10-100 | 防止占用过多集群资源 |
| `MAX_JOBS` | 最大任务数 | 1000-10000 | 防止队列过长影响调度 |
| `INTERVAL_TIME` | 提交间隔(秒) | 0.5-2.0 | 降低调度器负载 |
| `CHECK_TIME` | 检查间隔(秒) | 60-120 | 平衡响应速度和开销 |

## 使用方法

### 方法1: 使用便捷脚本

```bash
# 1. 生成作业脚本
python abacus.py workflow InputPoscar/ work_cal7/

# 2. 自动提交所有作业
python submit_jobs.py work_cal7
```

### 方法2: 使用模块方式

```bash
python -m abacus.submit_manager work_cal7
```

### 方法3: 在Python代码中使用

```python
from abacus.submit_manager import AutoSubmitter

# 创建提交器
submitter = AutoSubmitter(config_file="config/condor.ini")

# 运行提交流程
submitter.run("work_cal7", max_retries=3)
```

## 高级用法

### 1. 自定义配置文件

```bash
python submit_jobs.py work_cal7 --config my_config.ini
```

### 2. 调整重试次数

```bash
python -m abacus.submit_manager work_cal7 --max-retries 5
```

### 3. 自定义脚本匹配模式

```bash
python -m abacus.submit_manager work_cal7 --pattern "hmat_*.sh"
```

## 日志系统

### 日志位置

- **控制台输出**: 实时显示关键信息
- **文件日志**: `logs/submit_YYYYMMDD_HHMMSS.log`

### 日志级别

- **INFO**: 重要操作和状态信息
- **DEBUG**: 详细的调试信息
- **WARNING**: 警告信息（如重试）
- **ERROR**: 错误信息

### 日志示例

```
2025-12-25 10:00:00 - AutoSubmitter - INFO - ======================================================================
2025-12-25 10:00:00 - AutoSubmitter - INFO - Slurm 自动提交器启动
2025-12-25 10:00:00 - AutoSubmitter - INFO - ======================================================================
2025-12-25 10:00:00 - AutoSubmitter - INFO - 工作目录: work_cal7
2025-12-25 10:00:00 - AutoSubmitter - INFO - 分区: deimos
2025-12-25 10:00:00 - AutoSubmitter - INFO - 节点限制: 10
2025-12-25 10:00:00 - AutoSubmitter - INFO - 任务限制: 1000
2025-12-25 10:00:00 - AutoSubmitter - INFO - 提交间隔: 0.5秒
2025-12-25 10:00:00 - AutoSubmitter - INFO - 检查间隔: 80秒
2025-12-25 10:00:00 - AutoSubmitter - INFO - ======================================================================
2025-12-25 10:00:01 - AutoSubmitter - INFO - 扫描作业脚本: work_cal7
2025-12-25 10:00:01 - AutoSubmitter - INFO - 找到 100 个脚本文件
2025-12-25 10:00:01 - AutoSubmitter - INFO - 队列中共有 100 个作业待提交
2025-12-25 10:00:02 - AutoSubmitter - INFO - 作业已提交: hmat_0 (ID: 12345)
2025-12-25 10:00:03 - AutoSubmitter - INFO - 作业已提交: hmat_1 (ID: 12346)
...
2025-12-25 10:01:00 - AutoSubmitter - INFO - 队列状态 - 总任务:20, 运行中:15, 使用节点:15/10, 待提交:80
2025-12-25 10:01:00 - AutoSubmitter - INFO - 资源不足，等待中... (节点资源不足 (当前:15, 需要:1, 限制:10))
...
2025-12-25 10:05:00 - AutoSubmitter - INFO - ======================================================================
2025-12-25 10:05:00 - AutoSubmitter - INFO - 提交完成统计:
2025-12-25 10:05:00 - AutoSubmitter - INFO -   成功提交: 95
2025-12-25 10:05:00 - AutoSubmitter - INFO -   待提交: 0
2025-12-25 10:05:00 - AutoSubmitter - INFO -   失败: 5
2025-12-25 10:05:00 - AutoSubmitter - INFO - ======================================================================
```

## 资源限制机制

### 节点数限制

系统通过以下方式控制节点使用：

1. **实时监控**: 使用`squeue`命令获取当前使用的节点数
2. **预检查**: 提交前检查 `used_nodes + job_nodes <= TOTAL_NODE`
3. **等待机制**: 资源不足时，作业放回队列并等待

```python
# 伪代码
if used_nodes + job.nodes > TOTAL_NODE:
    # 资源不足，等待
    queue.push_back(job)
    sleep(CHECK_TIME)
else:
    # 可以提交
    submit_job(job)
```

### 任务数限制

```python
# 伪代码
if total_jobs >= MAX_JOBS:
    # 达到上限，等待
    queue.push_back(job)
    sleep(CHECK_TIME)
else:
    # 可以提交
    submit_job(job)
```

## 脚本解析

系统会自动从脚本中提取节点信息：

### 支持的格式

```bash
# 方式1: SBATCH指令
#SBATCH -N 2
#SBATCH --nodes=4

# 方式2: mpirun命令
mpirun -np 16 abacus

# 默认值
# 如果未找到，使用NODES=1, CORES=16
```

### 解析优先级

1. SBATCH -N / --nodes
2. SBATCH -n / --ntasks
3. mpirun -np
4. 默认值 (config/condor.ini)

## 错误处理

### 自动重试

提交失败会自动重试，最多重试3次（可配置）：

```python
# 重试逻辑
if submit_failed and job.retry_count < max_retries:
    queue.add_job(job)  # 放回队列重试
else:
    mark_as_failed(job)  # 标记为最终失败
```

### 常见错误

| 错误 | 原因 | 解决方案 |
|------|------|----------|
| `节点资源不足` | 当前使用节点数已达上限 | 等待作业完成后自动提交 |
| `达到最大任务数限制` | 队列中任务过多 | 等待任务执行完成 |
| `提交命令超时` | sbatch命令执行超时 | 检查网络和调度器状态 |
| `无法解析作业ID` | sbatch输出格式异常 | 检查Slurm版本兼容性 |

## 性能优化

### 大规模作业优化建议

对于10万+作业场景：

1. **调整提交间隔**
   ```ini
   INTERVAL_TIME = 1.0  # 增加到1秒，降低负载
   ```

2. **调整检查间隔**
   ```ini
   CHECK_TIME = 120  # 增加到2分钟，减少squeue调用
   ```

3. **分批提交**
   ```bash
   # 分多个目录
   python submit_jobs.py work_cal7/batch_01
   python submit_jobs.py work_cal7/batch_02
   ```

4. **增加节点限制**
   ```ini
   TOTAL_NODE = 50  # 如果有更多资源
   MAX_JOBS = 5000  # 增加任务数上限
   ```

### 内存使用

- 单个JobInfo对象: ~200 bytes
- 10万作业: ~20 MB
- 内存占用极小，可以轻松支持百万级作业

## 监控和管理

### 查看队列状态

```bash
# 查看所有作业
squeue -u $USER

# 查看运行中作业
squeue -u $USER -t RUNNING

# 查看等待中作业
squeue -u $USER -t PENDING

# 统计节点使用
squeue -u $USER -h -o '%D' | awk '{sum+=$1} END {print sum}'
```

### 取消作业

```bash
# 取消单个作业
scancel <job_id>

# 取消所有作业
scancel -u $USER

# 取消指定分区作业
scancel -u $USER -p deimos
```

### 暂停提交

按 `Ctrl+C` 可以安全停止提交器，已提交的作业不受影响。

## 实际案例

### 案例1: 10万个小作业

**场景**: 10万个单节点16核的结构优化计算

**配置**:
```ini
TOTAL_NODE = 50
MAX_JOBS = 5000
INTERVAL_TIME = 1.0
CHECK_TIME = 120
```

**预计时间**:
- 提交速度: 每秒1个
- 总提交时间: ~28小时（考虑等待时间）
- 可以中断后继续

### 案例2: 1000个大作业

**场景**: 1000个4节点64核的大规模计算

**配置**:
```ini
TOTAL_NODE = 100
MAX_JOBS = 200
INTERVAL_TIME = 0.5
CHECK_TIME = 60
```

**预计时间**:
- 提交速度: 每秒2个
- 总提交时间: ~10分钟

## 与现有系统集成

### 完整工作流

```bash
# 1. 准备结构文件
ls InputPoscar/*.poscar

# 2. 生成工作流脚本
python abacus.py workflow InputPoscar/ work_cal7/

# 3. 自动提交（新增）
python submit_jobs.py work_cal7

# 4. 监控作业状态
squeue -u $USER

# 5. 查看日志
tail -f logs/submit_*.log
```

## 故障排查

### 问题1: 提交速度过慢

**解决**:
- 减小 `INTERVAL_TIME`
- 增加 `TOTAL_NODE` 和 `MAX_JOBS`

### 问题2: 调度器报错

**解决**:
- 增加 `INTERVAL_TIME`
- 减小 `MAX_JOBS`
- 联系管理员调整用户限制

### 问题3: 脚本解析失败

**解决**:
- 检查脚本格式
- 手动在脚本中添加SBATCH指令：
  ```bash
  #SBATCH -N 2
  #SBATCH -p deimos
  ```

## 未来扩展

### 计划功能

1. **依赖关系管理**: 支持作业间依赖（--dependency）
2. **优先级控制**: 不同类型作业不同优先级
3. **Web界面**: 图形化监控和管理
4. **邮件通知**: 作业完成/失败通知
5. **断点续传**: 支持中断后从断点继续

### API接口

```python
from abacus.submit_manager import AutoSubmitter, JobInfo

# 创建提交器
submitter = AutoSubmitter()

# 手动添加作业
job = JobInfo(
    script_path="path/to/script.sh",
    job_dir="path/to/dir",
    job_name="my_job",
    nodes=2,
    cores=32,
    partition="deimos"
)
submitter.queue.add_job(job)

# 批量提交
submitter.run("work_dir")
```

## 总结

本自动提交系统具有以下优势：

✅ **可扩展性**: 支持10万+作业  
✅ **智能化**: 自动监控资源，智能调度  
✅ **可靠性**: 自动重试，详细日志  
✅ **易用性**: 简单命令，开箱即用  
✅ **可配置性**: 灵活配置，适应不同场景  
✅ **低侵入性**: 不修改现有脚本，无缝集成  

通过合理配置，可以高效管理大规模作业提交，显著提升工作效率。


