# 工作流续算

`resume` 命令用于从工作流中失败或中断的阶段继续计算，无需重跑整个流程。

## 基本用法

```bash
# 自动检测续算起点，生成续算脚本
python abacus.py resume work_cal/hmat_0/

# 预览续算计划，不生成脚本
python abacus.py resume work_cal/hmat_0/ --dry-run

# 强制从指定阶段开始
python abacus.py resume work_cal/hmat_0/ --from-stage Scf

# 不备份日志（覆盖原有日志）
python abacus.py resume work_cal/hmat_0/ --no-backup
```

生成的脚本为 `<job_name>_resume.sh`，用 `yhbatch` 提交即可。

## 阶段状态检测

`resume` 命令依次检查每个阶段目录，按以下优先级判断状态：

1. `error.txt` 存在 → failed
2. `converge.txt` 存在 → success
3. `ignore.txt` 存在 → ignored（未收敛但被允许继续）
4. `running.log` 存在 → incomplete
5. 其他 → not_started

`stat.log` 是工作流状态的权威记录，优先参考其内容。

## 续算策略

对失败阶段的处理分两种情况：

**清空重算**：该阶段只运行了一次就失败（目录中没有 `INPUT0`、`INPUT1` 等备份），
说明参数可能有问题，清空目录从头生成。建议先修改 template 或 `workflow.json` 再续算。

**继续迭代**：存在多次尝试的备份文件（`INPUT0`、`INPUT1`...）或已有部分优化结构
（`OUT.*/STRU_ION_D`），从上次结束处接着跑，利用已有结果。

## 日志备份

续算前自动备份 `stat.log` 和 `time.log`，命名规则：

```
stat.log  →  stat0.log（首次续算备份）
stat.log  →  stat1.log（二次续算备份）
...
```

续算脚本以追加方式写入日志，不覆盖历史记录。

## 文件命名规范

原始工作流每次尝试的备份：

```
INPUT → INPUT0, INPUT1, INPUT2, ...
```

续算模式的备份：

```
INPUT_R0_pre    首次续算开始前的状态
INPUT_R0_3      首次续算的第 4 次尝试（接续原始第 3 次之后）
INPUT_R1_pre    二次续算开始前的状态
INPUT_R1_5      二次续算的第 1 次尝试
```

## 典型场景

### 多次尝试未收敛

```bash
# 查看当前状态
python abacus.py resume work_cal/hmat_0/ --dry-run

# 修改参数（如放宽力收敛、调整 mixing）
vim config/template/Relax.yaml

# 生成续算脚本并提交
python abacus.py resume work_cal/hmat_0/
yhbatch work_cal/hmat_0/hmat_0_resume.sh
```

### 批量续算

```bash
# 扫描所有失败作业并生成续算脚本
python abacus.py batch-resume work_cal/

# 一键提交
python submit_jobs.py work_cal/ --resume
```

详见 [batch_resume_submit.md](batch_resume_submit.md)。

## 调试

```bash
# 查看生成的续算脚本内容
cat work_cal/hmat_0/hmat_0_resume.sh

# 手动测试单个阶段
cd work_cal/hmat_0/Relax/
python /path/to/abacus.py generate --work_dir . --stage Relax
yhrun -N 1 -n 16 -p deimos /path/to/abacus > test.log 2>&1
```

## 参数调整参考

```json
// config/workflow.json
{
  "Relax": {
    "kval": [0.04, 0.02, 0.015],
    "try_num": 3
  }
}
```

```yaml
# config/template/Relax.yaml
scf_nmax: 150
mixing_beta: 0.3
force_thr_ev: 0.02
```
