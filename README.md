# ABACUS 高通量计算工作流管理系统

自动化 ABACUS 第一性原理计算的工作流管理系统，支持完整的多阶段计算流程、实时监控和资源统计。

## ✨ 主要特性

- 🔄 **自动化工作流**: Test_spin → Coarse_relax → Relax → Scf → Band → Dos
- 📊 **实时监控**: 支持 LCAO 和 PW 基组，实时显示离子步、电子步和收敛状态
- ⏱️ **资源统计**: 自动记录每个阶段的运行时间和核时消耗
- 🔗 **结构传递**: 自动使用上一阶段的优化结构
- 🛠️ **智能重试**: 自动调整参数并重试未收敛的计算
- 📝 **模块化设计**: 清晰的代码结构，易于维护和扩展

## 📂 文件结构

```
abacusflow/
├── abacus.py                # 主入口：CLI 命令集合
├── abacus/                  # 功能模块目录
│   ├── generator.py         # 输入文件生成核心逻辑
│   ├── stru_utils.py        # STRU 文件处理工具
│   ├── input_utils.py       # INPUT 文件处理工具
│   ├── monitor.py           # 实时监控脚本（LCAO + PW）
│   ├── kpoints_utils.py     # K点网格计算
│   ├── submit.py            # 工作流脚本生成器
│   ├── create_input.py      # ABACUS 输入对象
│   ├── potential.py         # 赝势管理
│   └── basis.py             # 轨道基组管理
├── config/                  # 配置文件
│   ├── condor.ini           # 集群配置
│   ├── workflow.json        # 工作流配置
│   └── template/            # 各阶段 YAML 模板
│       ├── Test_spin.yaml
│       ├── Coarse_relax.yaml
│       ├── Relax.yaml
│       ├── Scf.yaml
│       ├── Band.yaml
│       └── Dos.yaml
├── examples/InputPoscar/             # 输入结构文件（VASP POSCAR 格式）
└── examples/work_cal/                # 计算工作目录（自动生成）
```

## 🚀 快速开始

### 1. 环境要求

```bash
# 激活包含 ASE 的 conda 环境
conda activate dftflow

# 主要依赖
- Python 3.x
- ASE (Atomic Simulation Environment)
- Click (CLI 工具)
- NumPy
```

### 2. 生成工作流

```bash
# 方式1：使用新格式（推荐）
python abacus.py workflow InputPoscar/ work_cal/

# 方式2：使用旧格式（向后兼容）
python abacus.py InputPoscar/ work_cal/ --workflow
```

生成的脚本包含：
- 完整的 6 阶段工作流
- 自动结构传递
- 实时监控
- 时间和资源统计
- 智能错误处理和重试

### 3. 提交作业

```bash
# 进入工作目录
cd work_cal/structure_name/

# 提交作业
yhbatch structure_name.sh

# 或直接运行（测试）
bash structure_name.sh
```

### 4. 查看进度

作业运行时会自动显示：
```
======================================================================
[10:30:15] 🔄 ION STEP 1
======================================================================
  ELEC  1  E=  -1289.089017 eV  ΔDens=2.207e-01  mag=0.648104
  ELEC  2  E=  -1288.848889 eV  ΔE=    0.240128 eV  ΔDens=1.157e-01
  ELEC  5  E=  -1288.203222 eV  ΔE=   -0.000417 eV  ΔDens=3.607e-03
  ELEC 10  E=  -1288.203658 eV  ΔE=    0.000004 eV  ΔDens=4.191e-05  🔸
  ✅ SCF converged!
  🔧 Max Force: 0.045678 eV/Å

[INFO] Test_spin completed in 450s (0.1250h)
[INFO] Test_spin consumed 2.00 Core-hours (1nodes * 16cores * 0.1250h)
```

### 5. 查看统计

作业完成后：
```bash
# 查看详细时间统计
cat work_cal/structure_name/time.log

# 查看状态日志
cat work_cal/structure_name/stat.log

# 查看总结
python abacus.py summary --root work_cal/structure_name/
```

## 🔧 CLI 命令详解

### workflow - 生成工作流

```bash
python abacus.py workflow <stru_path> <work_dir>

# 示例
python abacus.py workflow InputPoscar/ work_cal/
```

### generate - 生成输入文件

```bash
python abacus.py generate --work_dir <dir> --stage <stage>

# 示例
python abacus.py generate --work_dir . --stage Scf
```

### errors - 检查错误

```bash
python abacus.py errors --work_dir <dir>

# 检测的错误类型：
# - ABACUS_NOTICE_ERROR
# - INPUT_ERROR
# - OPERATION_FAILED
# - SCF_NOT_CONVERGED
# - RELAX_NOT_CONVERGED
# - SEGMENTATION_FAULT
# - OUT_OF_MEMORY
```

### converge - 检查收敛

```bash
python abacus.py converge --work_dir <dir>

# 正确识别以下状态：
# - SCF converged
# - Relaxation converged
# - Relaxation converged but SCF unconverged (标记为不可靠)
```

### spin - 判断磁性

```bash
python abacus.py spin --work_dir <dir>

# 自动判断：
# - |mag| > 0.004 → 创建 SPIN_ON 文件
# - |mag| ≤ 0.004 → 创建 SPIN_OFF 文件
```

### update - 更新参数

```bash
python abacus.py update --work_dir <dir> --try_num <n> --stage <stage>

# 自动调整：
# - K点密度（根据 workflow.json 配置）
# - mixing_beta（SCF 未收敛时）
# - scf_nmax（SCF 未收敛时）
# - force_thr_ev（Relax 未收敛时）
# - relax_nmax（Relax 未收敛时）
```

### summary - 结果汇总

```bash
python abacus.py summary --root <dir>

# 显示所有阶段的状态和统计信息
```

## 📊 工作流阶段详解

### 结构传递链

```
STRU.vasp (原始结构)
  ↓
Test_spin (测试磁性)
  ↓
Coarse_relax (粗优化) → OUT.Coarse_relax/STRU_ION_D
  ↓
Relax (精细优化，使用 Coarse_relax 优化结构) → OUT.Relax/STRU_ION_D
  ↓
Scf (自洽计算，使用 Relax 优化结构) → Scf/STRU
  ↓
Band/Dos (能带/态密度，使用 Scf 结构 + 电荷密度)
```

### 各阶段配置

所有阶段参数在 `config/workflow.json` 和 `config/template/*.yaml` 中配置：

1. **Test_spin**: 快速测试，判断是否需要自旋极化
2. **Coarse_relax**: 粗网格优化，快速获得近似结构
3. **Relax**: 精细优化，获得最终结构
4. **Scf**: 自洽场计算，生成电荷密度
5. **Band**: 能带结构（需要 Scf 的电荷密度）
6. **Dos**: 态密度（需要 Scf 的电荷密度）

## 📈 时间和资源统计

### 实时统计

每个阶段完成后自动显示：
```
[INFO] Relax started at 2025-12-24 18:30:15
[INFO] Relax completed in 3600s (1.0000h)
[INFO] Relax consumed 16.00 Core-hours (1nodes * 16cores * 1.0000h)
```

### 最终汇总

所有阶段完成后显示总结表格：
```
======================================================================
                    WORKFLOW SUMMARY
======================================================================

[INFO] Individual Stage Statistics:
----------------------------------------------------------------------
Stage           Duration(s)  Nodes    Cores    Core-hours      Status         
----------------------------------------------------------------------
Test_spin       450s (0.13h) 1        16       2.00            success        
Coarse_relax    2850s (0.79h)1        16       12.67           success        
Relax           3600s (1.00h)1        16       16.00           success        
Scf             1800s (0.50h)1        16       8.00            success        
Band            900s (0.25h) 1        16       4.00            success        
Dos             600s (0.17h) 1        16       2.67            success        
----------------------------------------------------------------------

[INFO] Total Workflow Statistics:
  Total Duration    : 10200s (2.8333h)
  Total Core-hours  : 45.34
  Workflow Start    : 2025-12-24 18:30:15
  Workflow End      : 2025-12-24 21:20:15

======================================================================
```

### time.log 格式

```
# Stage statistics
# Format: Stage | Duration(s) | Nodes | Cores | Core-hours | Status
Test_spin|450|1|16|2.00|success
Coarse_relax|2850|1|16|12.67|success
...
```

## 🎯 实时监控说明

### 监控信息

监控脚本自动显示：
- **ION STEP**: 离子步编号
- **ELEC**: 电子步编号
- **E**: E_KohnSham 能量（eV）
- **ΔE**: 能量变化（eV）
- **ΔDens**: 密度误差（收敛判据）
- **mag**: 磁矩（仅当 |mag| > 1e-4 时显示）
- **收敛标记**: ✅ converged、🔸 approaching

### 收敛判据

- `ΔDens < scf_thr (1e-7)` → ✅ converged
- `ΔDens < 1e-4` → 🔸 approaching
- `ΔDens > 1e-4` → 继续迭代

### 格式支持

- **LCAO**: 局域轨道基组（详细信息在 `OUT.*/running_*.log`）
- **PW**: 平面波基组（详细信息在标准输出）

## 🔄 智能重试机制

### 多次尝试

每个阶段支持多次尝试（在 `workflow.json` 中配置 `try_num`）：
```json
"Relax": {
    "try_num": 2,
    "kval": [0.04, 0.02],
    ...
}
```

### 自动调整策略

1. **K点密度**: 使用 `kval` 列表中的下一个值
2. **SCF 参数**: 
   - `mixing_beta` × 0.7
   - `scf_nmax` + 50
3. **Relax 参数**:
   - `force_thr_ev` × 1.5
   - `relax_nmax` + 20
4. **结构**: 使用上一次迭代的 `STRU_ION_D`

### 错误处理

- **error.txt 存在**: 立即退出（真正的错误）
- **ignore.txt 存在**: 允许未收敛继续（需要在 workflow.json 中设置）
- **converge.txt 存在**: 成功，进入下一阶段

## ⚙️ 配置文件说明

### condor.ini

```ini
[ABACUS]
ABACUS_BIN = /path/to/abacus
PSEUDO_POTENTIAL_DIR = /path/to/pseudopotentials
ORBITAL_DIR = /path/to/orbitals

[PARAMETERS]
KSPACING = 0.13

[ALLOW]
PARTITION = deimos
```

### workflow.json

```json
{
  "Test_spin": {
    "template": "Test_spin",
    "node": 1,
    "core": 16,
    "try_num": 2,
    "kval": [0.15],
    "ktype": "Gamma"
  },
  ...
}
```

### template/*.yaml

```yaml
calculation: "relax"
relax_nmax: 200
force_thr_ev: 0.01
stress_thr: 0.5
relax_method: "cg"
out_stru: 1
...
```

## 📚 最近更新

### v2.0 (2025-12-24)

- ✅ **代码模块化**: abacus.py 从 1063行 精简到 522行（51% 精简）
- ✅ **结构自动传递**: 实现完整的结构优化链
- ✅ **实时监控**: 支持 LCAO 和 PW 格式的实时显示
- ✅ **时间统计**: 自动记录每阶段用时和核时消耗
- ✅ **优化提示**: 更详细的信息输出
- ✅ **Bug 修复**: 修复多个路径和逻辑错误

### v1.x

- ✅ 修复收敛检查逻辑
- ✅ 清理冗余代码（~700MB）
- ✅ 优化代码结构

## 🐛 故障排查

### 常见问题

1. **ImportError: No module named 'ase'**
   ```bash
   conda activate dftflow
   ```

2. **找不到 abacus 可执行文件**
   - 检查 `config/condor.ini` 中的 `ABACUS_BIN` 路径

3. **实时监控不工作**
   ```bash
   # 检查 stdbuf 是否可用
   command -v stdbuf
   ```

4. **stat.log 路径错误**
   - 已在最新版本中修复

### 调试模式

```bash
# 查看详细输出
bash -x structure_name.sh

# 单独测试某个阶段
cd work_cal/structure_name/Scf
python /path/to/abacus.py generate --work_dir . --stage Scf
```


