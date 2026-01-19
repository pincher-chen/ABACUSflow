# 自动提交系统 - 快速开始

## 一分钟上手

### 1. 生成作业脚本

```bash
python abacus.py workflow InputPoscar/ work_cal7/
```

### 2. 配置资源限制

编辑 `config/condor.ini`：

```ini
[ALLOW]
PARTITION = deimos        # 你的Slurm分区名
TOTAL_NODE = 10          # 最多使用多少节点
MAX_JOBS = 1000          # 队列中最多多少任务
INTERVAL_TIME = 0.5      # 每隔多少秒提交一个作业
CHECK_TIME = 80          # 每隔多少秒检查队列状态
```

### 3. 自动提交

```bash
python submit_jobs.py work_cal7
```

就这么简单！系统会自动：
- ✅ 扫描所有作业脚本
- ✅ 监控队列状态
- ✅ 智能提交作业
- ✅ 避免超出资源限制
- ✅ 失败自动重试

## 监控提交进度

### 实时查看日志

```bash
tail -f logs/submit_*.log
```

### 查看队列状态

```bash
squeue -u $USER
```

### 统计提交情况

日志会显示类似信息：

```
2025-12-25 10:00:00 - AutoSubmitter - INFO - ======================================================================
2025-12-25 10:00:00 - AutoSubmitter - INFO - Slurm 自动提交器启动
2025-12-25 10:00:00 - AutoSubmitter - INFO - 工作目录: work_cal7
2025-12-25 10:00:00 - AutoSubmitter - INFO - 队列中共有 1000 个作业待提交
2025-12-25 10:00:01 - AutoSubmitter - INFO - 作业已提交: hmat_0 (ID: 12345)
2025-12-25 10:00:02 - AutoSubmitter - INFO - 作业已提交: hmat_1 (ID: 12346)
...
2025-12-25 10:01:00 - AutoSubmitter - INFO - 队列状态 - 总任务:20, 运行中:15, 使用节点:8/10, 待提交:980
...
2025-12-25 10:05:00 - AutoSubmitter - INFO - 提交完成统计:
2025-12-25 10:05:00 - AutoSubmitter - INFO -   成功提交: 995
2025-12-25 10:05:00 - AutoSubmitter - INFO -   待提交: 0
2025-12-25 10:05:00 - AutoSubmitter - INFO -   失败: 5
```

## 常见场景

### 场景1: 小规模测试（10个作业）

```bash
# 1. 快速配置
cat > config/condor.ini << EOF
[ALLOW]
PARTITION = deimos
TOTAL_NODE = 5
MAX_JOBS = 100
INTERVAL_TIME = 0.5
CHECK_TIME = 60
EOF

# 2. 提交
python submit_jobs.py work_cal7
```

### 场景2: 大规模生产（10万个作业）

```bash
# 1. 优化配置
cat > config/condor.ini << EOF
[ALLOW]
PARTITION = deimos
TOTAL_NODE = 50           # 使用更多节点
MAX_JOBS = 5000           # 允许更多任务
INTERVAL_TIME = 1.0       # 增加间隔，降低负载
CHECK_TIME = 120          # 减少检查频率
EOF

# 2. 后台运行
nohup python submit_jobs.py work_cal7 > submit.out 2>&1 &

# 3. 监控
tail -f submit.out
tail -f logs/submit_*.log
```

### 场景3: 资源受限（节点数少）

```bash
# 配置小的资源限制
cat > config/condor.ini << EOF
[ALLOW]
PARTITION = deimos
TOTAL_NODE = 2            # 只使用2个节点
MAX_JOBS = 10             # 最多10个任务
INTERVAL_TIME = 2.0       # 慢慢提交
CHECK_TIME = 30           # 频繁检查
EOF

# 提交会自动等待资源释放
python submit_jobs.py work_cal7
```

## 停止提交

按 `Ctrl+C` 即可安全停止，已提交的作业不受影响。

## 中断后继续

直接重新运行即可：

```bash
python submit_jobs.py work_cal7
```

系统会：
1. 重新扫描所有脚本
2. 检查哪些已经在队列中
3. 只提交尚未提交的作业

## 取消所有作业

```bash
scancel -u $USER
```

## 常见问题

### Q1: 提交速度太慢？

**A**: 减小 `INTERVAL_TIME`，增加 `TOTAL_NODE` 和 `MAX_JOBS`

```ini
INTERVAL_TIME = 0.1      # 从0.5改为0.1
TOTAL_NODE = 50          # 从10改为50
```

### Q2: 调度器报错或卡死？

**A**: 增加 `INTERVAL_TIME`，减小 `MAX_JOBS`

```ini
INTERVAL_TIME = 2.0      # 从0.5改为2.0
MAX_JOBS = 500           # 从1000改为500
```

### Q3: 资源一直显示不足？

**A**: 检查配置是否过于保守

```bash
# 查看当前队列状态
squeue -u $USER

# 查看分区资源
sinfo -p deimos

# 可能需要增加TOTAL_NODE
```

### Q4: 如何只提交特定目录的作业？

```bash
# 只提交hmat_0到hmat_99
python submit_jobs.py work_cal7/hmat_[0-9]*

# 或使用通配符
python -m abacus.submit_manager work_cal7/batch_01
```

### Q5: 如何查看失败的作业？

查看日志文件末尾的失败作业列表：

```bash
tail -50 logs/submit_*.log
```

## 高级用法

### 自定义重试次数

```bash
python -m abacus.submit_manager work_cal7 --max-retries 5
```

### 使用自定义配置

```bash
python -m abacus.submit_manager work_cal7 --config my_config.ini
```

### 匹配特定脚本

```bash
python -m abacus.submit_manager work_cal7 --pattern "hmat_*.sh"
```

## 完整工作流示例

```bash
# 1. 准备结构文件
ls InputPoscar/*.poscar

# 2. 生成工作流脚本
python abacus.py workflow InputPoscar/ work_cal7/

# 3. 检查配置
cat config/condor.ini

# 4. 自动提交
python submit_jobs.py work_cal7

# 5. 监控进度
# 终端1: 查看提交日志
tail -f logs/submit_*.log

# 终端2: 查看队列
watch -n 5 'squeue -u $USER'

# 终端3: 统计运行状态
watch -n 10 'squeue -u $USER | tail -n +2 | wc -l'
```

## 更多信息

详细设计文档：[doc/auto_submit.md](doc/auto_submit.md)

系统架构、配置说明、故障排查等详细信息请参考完整文档。


