# 快速开始

## 基本流程

### 1. 生成作业脚本

```bash
# 多阶段工作流（Test_spin → Relax → Scf → Band/Dos）
python abacus.py workflow examples4stru/batch_poscar/ work_cal/

# 单阶段批量（以磁性共线计算为例）
python abacus.py single examples4stru/batch_mcif/ work_cal/ -t Scf-Spin2 -k 0.02 -spin 2
```

### 2. 配置资源限制

编辑 `config/condor.ini`：

```ini
[ALLOW]
PARTITION      = deimos   # Slurm 分区名
TOTAL_NODE     = 10       # 同时占用节点数上限
MAX_JOBS       = 1000     # 队列中最大任务数
INTERVAL_TIME  = 0.5      # 每次提交后等待时间（秒）
CHECK_TIME     = 80       # 资源不足时重新检查间隔（秒）
```

### 3. 自动提交

```bash
python submit_jobs.py work_cal/
```

系统自动扫描脚本、监控队列、等待资源释放、失败自动重试。按 `Ctrl+C` 安全停止，
已提交的作业不受影响。中断后直接重新运行，会跳过已提交的作业。

---

## 磁性体系

### 工作流模式：自动检测磁性

工作流包含 `Test_spin` 阶段，完成后自动调用 `abacus.py spin` 解析总磁矩并生成标记文件：

- `SPIN_ON`（`|mag| > 0.004`）：后续阶段使用 `nspin=2`
- `SPIN_OFF`（`|mag| <= 0.004`）：后续阶段使用 `nspin=1`

若阶段模板中显式写了 `nspin`，模板优先级更高。

查看磁性判断结果：

```bash
ls work_cal/hmat_0/SPIN_*
cat work_cal/hmat_0/SPIN_ON   # 文件内容为磁矩数值
```

手动强制开启/关闭自旋：

```bash
# 强制有自旋
rm -f work_cal/hmat_0/SPIN_OFF && echo 2.0 > work_cal/hmat_0/SPIN_ON

# 强制无自旋
rm -f work_cal/hmat_0/SPIN_ON && touch work_cal/hmat_0/SPIN_OFF
```

### 单阶段模式：直接指定自旋

手动调用 `generate` 或 `single` 时用 `-spin` 控制：

| 参数 | 含义 | STRU magmom 格式 |
|------|------|----------------|
| `-spin 1` | 非磁性（默认） | 不写 |
| `-spin 2` | 共线磁性 | `magmom <scalar>` |
| `-spin 4` | 非共线/SOC | `magmom <mx> <my> <mz>` |

```bash
# mcif 包含磁矩，共线计算
python abacus.py generate --work_dir work_cal/agm001014431 --stage Scf-Spin2 \
    --stru_file examples4stru/batch_mcif/agm001014431.mcif -spin 2 -k 0.02

# CIF 无磁矩，用 --guess-mag 给每个原子赋初始磁矩 2.0 μB
python abacus.py generate --work_dir work_cal/agm001003677 --stage Scf-Spin2 \
    --stru_file examples4stru/batch_cif/agm001003677.cif -spin 2 --guess-mag

# 批量单阶段
python abacus.py single examples4stru/batch_mcif/ work_cal/ -t Scf-Spin2 -k 0.02 -spin 2
```

mcif 中的磁矩字段（`_atom_site_moment_crystalaxis_x/y/z`）会被自动读取并经坐标变换写入
STRU。`spin=2` 时写标量，`spin=4` 时写三分量。详见 [doc/magnetic_calculation.md](doc/magnetic_calculation.md)。

---

## 监控与管理

```bash
# 查看提交日志
tail -f logs/submit_*.log

# 查看队列
squeue -u $USER

# 取消所有作业
scancel -u $USER
```

---

## 常见问题

**提交速度慢**：减小 `INTERVAL_TIME`，增大 `TOTAL_NODE` 和 `MAX_JOBS`。

**调度器报错**：增大 `INTERVAL_TIME`，减小 `MAX_JOBS`。

**资源一直不足**：`squeue -u $USER` 查看实际节点占用，对比 `TOTAL_NODE` 是否设置过小。

**磁性判断不生效**：检查 `work_cal/<job>/SPIN_*` 文件是否存在；若 `Test_spin/running*.log`
中没有 `total magnetism (Bohr mag/cell) = ...`，说明 Test_spin 未正常完成。

---

## 更多文档

- [磁性计算指南](doc/magnetic_calculation.md)
- [自动提交系统详细说明](doc/auto_submit.md)
- [批量续算](doc/batch_resume_submit.md)
- [续算工作流](doc/restart_workflow.md)
