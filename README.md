# abacusflow

ABACUS 第一性原理计算的工作流自动化工具，支持多阶段计算、批量提交和磁性材料处理。

## 环境要求

```bash
conda activate dftflow  # 包含 Python 3.x, ASE, Click, NumPy, pymatgen
```

## 文件结构

```
abacusflow/
├── abacus.py              # CLI 主入口
├── submit_jobs.py         # 批量提交脚本
├── abacus/
│   ├── generator.py       # 输入文件生成
│   ├── stru_utils.py      # STRU 文件处理（含磁矩读写）
│   ├── kpoints_utils.py   # K 点网格计算
│   ├── input_utils.py     # INPUT 文件处理
│   ├── submit.py          # 提交脚本生成
│   ├── submit_manager.py  # 自动提交管理
│   ├── potential.py       # 赝势管理
│   └── basis.py           # 轨道基组管理
├── config/
│   ├── condor.ini         # 集群和路径配置
│   ├── workflow.json      # 工作流阶段配置
│   └── template/          # 各阶段 YAML 参数模板
├── examples4stru/         # 输入结构示例
│   ├── batch_poscar/      # VASP POSCAR 格式
│   ├── batch_cif/         # CIF 格式
│   └── batch_mcif/        # 磁结构 magCIF 格式
└── doc/                   # 详细文档
```

## 快速开始

### 工作流模式（多阶段）

```bash
# 生成完整工作流脚本（Test_spin → Relax → Scf → Band/Dos）
python abacus.py workflow examples4stru/batch_poscar/ work_cal/

# 进入某材料目录提交
cd work_cal/<structure>/
yhbatch <structure>.sh
```

### 单阶段模式（批量生成+提交）

```bash
# 为某目录下所有结构文件生成指定阶段的输入文件和提交脚本
python abacus.py single examples4stru/batch_mcif/ work_cal/ -t Scf-Spin2 -k 0.02 -spin 2

# 批量提交
python submit_jobs.py work_cal/ --single
```

### 生成单个材料的输入文件

```bash
python abacus.py generate \
    --work_dir work_cal/agm001014431 \
    --stage Scf-Spin2 \
    --stru_file examples4stru/batch_mcif/agm001014431.mcif \
    --spin 2 \
    --kval 0.02
```

## 磁性计算

支持三种自旋模式：

| 参数 | nspin | 适用场景 | STRU magmom 格式 |
|------|-------|---------|----------------|
| `-spin 1` | 1 | 非磁性 | 无 |
| `-spin 2` | 2 | 共线磁性（FM/AFM） | `magmom <scalar>` |
| `-spin 4` | 4 | 非共线/SOC | `magmom <mx> <my> <mz>` |

输入结构文件应为 `.mcif` 格式（包含 `_atom_site_moment_crystalaxis_x/y/z` 字段）。
详见 [doc/magnetic_calculation.md](doc/magnetic_calculation.md)。

## 可用计算模板

| 模板名 | 用途 |
|--------|------|
| `Test_spin` | 快速测试磁性 |
| `Coarse_relax` / `Relax` / `Fine_relax` | 结构优化 |
| `Scf` | 非磁性自洽 |
| `Scf-Spin2` | 共线磁性自洽 (nspin=2) |
| `Scf-Soc` | 非共线/SOC 自洽 (nspin=4) |
| `Band` / `Dos` | 能带/态密度 |

## 主要 CLI 命令

```bash
# 生成工作流
python abacus.py workflow <stru_path> <work_dir>

# 生成输入文件
python abacus.py generate --work_dir <dir> --stage <stage> [--spin N] [--kval 0.02]

# 批量单阶段
python abacus.py single <input_dir> <output_dir> -t <template> -k <kval> [-spin N]

# 批量提交
python submit_jobs.py <work_dir> [--single]

# 查看错误
python abacus.py errors --work_dir <dir>

# 检查收敛
python abacus.py converge --work_dir <dir>

# 判断磁性（写 SPIN_ON / SPIN_OFF 文件）
python abacus.py spin --work_dir <dir>
```

## 配置

`config/condor.ini` 配置集群参数：

```ini
[ABACUS]
ABACUS_BIN     = /path/to/abacus
PSEUDO_DIR     = /path/to/pseudopotentials
ORBITAL_DIR    = /path/to/orbitals

[ALLOW]
PARTITION      = your_partition
TOTAL_NODE     = 10
MAX_JOBS       = 1000
```

`config/workflow.json` 配置每个阶段的节点数、K 点密度、重试次数等。
`config/template/*.yaml` 配置各阶段的 ABACUS INPUT 参数。

## 文档

- [磁性计算指南](doc/magnetic_calculation.md)
- [自动提交系统](doc/auto_submit.md)
- [自定义参数](doc/custom_parameters.md)
- [批量提交](doc/batch_resume_submit.md)
