# 自动提交系统 - 快速开始

## 一分钟上手

### 1. 生成作业脚本

```bash
python abacus.py workflow InputPoscar/ work_cal7/
```

**新增：磁性/自旋自动处理**

- 工作流默认包含 `Test_spin`（测试磁性）阶段：用于快速 SCF 测试并判断是否需要自旋极化。
- `Test_spin` 结束后，脚本会自动运行 `python abacus.py spin --work_dir .`，解析输出中的
  `total magnetism (Bohr mag/cell) = ...` 并在作业根目录（例如 `work_cal7/hmat_0/`）生成标记：
  - `SPIN_ON`：当 \(|mag| > 0.004\)（文件内容为最终磁矩数值）
  - `SPIN_OFF`：当 \(|mag| \le 0.004\)
- 后续阶段（`Coarse_relax/Relax/Scf/Band/Dos`）生成输入时会自动读取 `SPIN_ON/SPIN_OFF` 来设置 `nspin`
  （若对应阶段的 YAML 模板中显式写了 `nspin`，则以模板为准）。

**新增：CIF 自带磁矩（magmom）自动写入 STRU**

- 如果结构文件是 **CIF 格式**，且包含磁矩信息 loop（pymatgen/VASP 常见输出）：
  ```
  loop_
   _atom_site_moment_label
   _atom_site_moment_crystalaxis_x
   _atom_site_moment_crystalaxis_y
   _atom_site_moment_crystalaxis_z
    Fe0  0.00  0.00  2.50
    O1   0.00  0.00  0.05
    ...
  ```
- 在生成 `STRU` 时，**自动读取并写入**每个原子的初始磁矩到 `ATOMIC_POSITIONS` 块：
  - **`-spin 2`（共线磁性）**：写入 `magmom <mz>`（只用 z 分量）
  - **`-spin 4`（非共线磁性）**：写入 `magmom <mx> <my> <mz>`（完整矢量）
  - **`-spin 1`（无磁性）**：不写入 `magmom`
- 如果 CIF 没有磁矩信息，可使用 `--guess-mag` 猜测初始磁矩（默认 2.0 μB）
- 详见下文 **4.4 节**的手动使用示例

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

### 场景4: 磁性体系（查看/重跑/强制自旋开关）

#### 4.1 查看磁性判断结果

以 `hmat_0` 为例：

```bash
ls -l work_cal7/hmat_0/SPIN_*
cat work_cal7/hmat_0/SPIN_ON   # 若存在：表示需要自旋（文件内是磁矩数值）
```

#### 4.2 手动重跑磁性判断（需要时）

```bash
python abacus.py spin --work_dir work_cal7/hmat_0/Test_spin
```

**说明**：

- `abacus.py spin` 会从 `running*.log` 中解析 `total magnetism (Bohr mag/cell) = ...`
- 并在 `--work_dir` 的**父目录**生成 `SPIN_ON` / `SPIN_OFF`
  - `SPIN_ON`：当 \(|mag| > 0.004\)
  - `SPIN_OFF`：当 \(|mag| \le 0.004\)

#### 4.3 强制开启/关闭自旋（无需改代码）

- **方式A：改模板（推荐，最清晰）**
  - 在对应阶段模板中写死 `nspin`，例如编辑 `config/template/Scf.yaml`：
    - `nspin: 1`（强制无自旋）
    - `nspin: 2`（强制共线自旋）
    - `nspin: 4`（强制非共线；通常还需要 `noncolin: 1`，可参考 `doc/custom_parameters.md`）
- **方式B：直接写标记文件（快速但要自担风险）**
  - 强制无自旋：

```bash
rm -f work_cal7/hmat_0/SPIN_ON
touch work_cal7/hmat_0/SPIN_OFF
```

  - 强制有自旋：

```bash
rm -f work_cal7/hmat_0/SPIN_OFF
echo 2.0 > work_cal7/hmat_0/SPIN_ON
```

#### 4.4 手动生成输入时：指定自旋与初始磁矩（`-spin` / `--guess-mag`）

当你不走工作流脚本、而是手动调用 `abacus.py generate` 或 `abacus.py single` 时，可以用 `-spin` 显式控制 STRU 里是否写 `magmom ...`。

**参数说明：**
- **`-spin N`**（默认 1）：
  - `-spin 1`：不写 `magmom ...`（非磁性）
  - `-spin 2`：共线磁性，写 `magmom <value>`
  - `-spin 4`：非共线磁性，写 `magmom <mx> <my> <mz>`
- **`--guess-mag`**：仅当 CIF 无磁矩时生效，猜测初始磁矩为 2.0 μB

**快速示例：**

```bash
# CIF 有磁矩：自动读取并写入 STRU
python abacus.py generate --work_dir ./job --stage Scf \
  --stru_file agm001613165.cif -spin 2

# CIF 无磁矩：猜测初始磁矩
python abacus.py generate --work_dir ./job --stage Scf \
  --stru_file no_magmom.cif -spin 2 --guess-mag

# 单步批量计算模式
python abacus.py single cif_files/ work_out/ \
  -t Scf -k 0.02 -spin 2 --guess-mag
```

> 📖 **详细使用指南和更多示例**请参考下文 **场景5: 使用带磁矩的 CIF 文件**

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

### Q6: 为什么我已经有磁性，但后续阶段还是无自旋？

请检查以下几点：

1. `Test_spin` 是否确实跑完并生成 `SPIN_ON`（或至少没有生成 `SPIN_OFF`）：

```bash
ls -l work_cal7/<job>/SPIN_*
```

2. 阶段模板是否显式指定了 `nspin` 并覆盖了自动检测（模板优先级更高）。
3. 如果 `Test_spin` 输出里没有出现 `total magnetism (Bohr mag/cell) = ...`，`abacus.py spin` 会给出警告，
   此时请优先检查 `Test_spin/running*.log` 是否完整。

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

### （可选）手动生成输入时指定自旋

当你手动调用 `abacus.py generate`（不通过工作流脚本）时，可以用 `-spin` 参数控制磁性设置：

```bash
# 共线磁性（CIF 有磁矩会自动读取）
python abacus.py generate --work_dir <dir> --stage Scf -spin 2

# CIF 无磁矩时猜测初始磁矩
python abacus.py generate --work_dir <dir> --stage Scf -spin 2 --guess-mag
```

> 📖 完整说明和示例请参考 **场景4.4** 和 **场景5**

### 场景5: 使用带磁矩的 CIF 文件（Alexandria 数据库）

当你从 Alexandria 或其他数据库获得的 CIF 文件包含 DFT 计算得到的磁矩信息时：

#### 5.1 检查 CIF 是否包含磁矩

```bash
# 查看 CIF 文件中是否有磁矩信息
grep "_atom_site_moment" agm001613165.cif
```

如果有如下输出，说明包含磁矩：
```
loop_
 _atom_site_moment_label
 _atom_site_moment_crystalaxis_x
 _atom_site_moment_crystalaxis_y
 _atom_site_moment_crystalaxis_z
  Y0  0.00000000  0.00000000  -0.01700000
  Y1  0.00000000  0.00000000  -0.01700000
  C2  0.00000000  0.00000000  -0.00800000
  Fe3  0.00000000  0.00000000  0.16500000
  Cu4  0.00000000  0.00000000  0.00200000
```

#### 5.2 单个 CIF 文件测试

```bash
# 创建测试目录
mkdir -p test_magnetic
cd test_magnetic

# 共线磁性计算（自动读取 CIF 中的磁矩）
python ../abacus.py generate --work_dir . --stage Scf \
  --stru_file /path/to/agm001613165.cif -spin 2

# 检查生成的 STRU 文件
cat STRU | grep -A 2 "ATOMIC_POSITIONS"
# 应该看到每个原子行末尾有 magmom 值
```

生成的 STRU 格式示例：
```
ATOMIC_POSITIONS
Direct

Y
0.0
2
0.0 0.5 0.0 1 1 1 magmom -0.017
0.5 0.0 0.0 1 1 1 magmom -0.017
...
Fe
0.0
1
0.0 0.0 0.5 1 1 1 magmom 0.165
```

#### 5.3 批量处理带磁矩的 CIF 文件

```bash
# 将所有 CIF 文件放在一个目录
mkdir -p cif_magnetic/
cp /path/to/alexandria/convex_hull_res_mag/*/*.cif cif_magnetic/

# 使用 single 模式批量生成作业（共线磁性）
python abacus.py single cif_magnetic/ work_magnetic/ \
  -t Scf -k 0.02 -spin 2

# 自动提交
python submit_jobs.py work_magnetic/
```

#### 5.4 处理非共线磁性

如果 CIF 中的磁矩有明显的 x, y 分量（非共线）：

```bash
# 非共线磁性计算
python abacus.py generate --work_dir . --stage Scf \
  --stru_file noncollinear.cif -spin 4

# 生成的 STRU 会写入：magmom mx my mz
```

#### 5.5 CIF 无磁矩但需要猜测初始磁矩

```bash
# 对于过渡金属体系，CIF 没有磁矩信息，但想给初始磁矩
python abacus.py generate --work_dir . --stage Scf \
  --stru_file no_magmom.cif -spin 2 --guess-mag

# 这会给每个原子设置初始磁矩 2.0 μB
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


