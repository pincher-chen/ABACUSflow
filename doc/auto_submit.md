# 自动提交系统

`submit_jobs.py` 实现了大批量作业自动提交到 Slurm 的功能，通过内存队列和实时资源监控
控制节点占用，避免超出集群限制。

## 使用方法

```bash
# 生成工作流脚本后，一键提交所有作业
python submit_jobs.py work_cal/

# 仅提交续算脚本（*_resume.sh）
python submit_jobs.py work_cal/ --resume

# 仅提交单阶段脚本（*_single.sh）
python submit_jobs.py work_cal/ --single
```

## 工作流程

```
扫描 work_dir 下所有 .sh 脚本
  → 解析每个脚本的节点数（#SBATCH -N 或默认值）
  → 加入内存队列
  → 循环：检查当前队列状态（squeue）
      如果 used_nodes + job.nodes <= TOTAL_NODE 且 total_jobs < MAX_JOBS：
          提交作业（sbatch），记录作业 ID
          等待 INTERVAL_TIME 秒
      否则：
          将作业放回队列，等待 CHECK_TIME 秒
  → 所有作业提交完成后退出
```

按 `Ctrl+C` 可以安全停止，已提交的作业不受影响。重新运行会跳过已提交的作业。

## 配置（config/condor.ini）

```ini
[ALLOW]
PARTITION      = deimos    # Slurm 分区
NODES          = 1         # 默认节点数（脚本未指定时使用）
CORES_PER_NODE = 64        # 每节点核心数
TOTAL_NODE     = 10        # 同时占用节点数上限（最重要）
MAX_JOBS       = 1000      # 队列中最大任务数
INTERVAL_TIME  = 0.5       # 每次提交后等待时间（秒）
CHECK_TIME     = 80        # 资源不足时重新检查的间隔（秒）
```

配置文件每隔 `CONFIG_RELOAD_INTERVAL`（默认 600 秒）自动重载，无需重启进程。

## 节点数解析优先级

1. `#SBATCH -N <n>` 或 `#SBATCH --nodes=<n>`
2. `mpirun -np <n>`（估算节点数）
3. `condor.ini` 中的 `NODES` 默认值

## 日志

提交日志写入 `logs/submit_YYYYMMDD_HHMMSS.log`，同时在控制台实时输出关键事件。

## 监控队列状态

```bash
# 查看当前作业
squeue -u $USER

# 统计占用节点数
squeue -u $USER -h -o '%D' | awk '{sum+=$1} END {print sum}'

# 实时刷新
watch -n 5 'squeue -u $USER'
```

## 常见问题

| 现象 | 原因 | 处理 |
|------|------|------|
| 长时间等待不提交 | used_nodes 超过 TOTAL_NODE | 等待运行中作业释放节点，或调高 TOTAL_NODE |
| 提交失败并重试 | sbatch 命令超时或报错 | 检查 Slurm 状态；会自动重试最多 3 次 |
| 脚本未被识别 | 后缀不是 `.sh` 或被过滤 | 检查 `--resume` / `--single` 选项是否匹配脚本名 |
