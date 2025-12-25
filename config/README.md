# ABACUSflow 配置文件说明

## 配置文件位置

配置文件位于 `config/condor.ini`，使用 INI 格式。

## 配置节说明

### [ABACUS]
ABACUS 相关配置：
- `ABACUS_DIR`: ABACUS 可执行文件所在目录
- `ABACUS_EXE`: ABACUS 可执行文件名（默认：abacus）
- `PSEUDO_POTENTIAL_DIR`: 赝势文件目录（绝对路径）
- `ORBITAL_DIR`: 轨道基组文件目录（绝对路径）

### [STRU]
结构文件相关配置：
- `PATH`: 结构文件搜索路径（可选）
- `SUFFIX`: 结构文件后缀（默认：*.vasp）

### [MODULE]
模块加载配置：
- `MODULES`: 需要加载的模块列表，用空格分隔
  - 例如：`intel/oneapi2023.2_noimpi mpi/mpich/4.1.2-icc-oneapi2023.2-ch4 elpa/2024.05.001-icc-oneapi2023.2-mpich-4.1.2-ch4`

### [ALLOW]
作业提交相关配置：
- `PARTITION`: 计算分区名称（如：deimos, mars）
- `NODES`: 每个作业使用的节点数（默认：1）
- `CORES_PER_NODE`: 每个节点的核心数（默认：64）
- `TOTAL_NODE`: 总节点数限制（用于作业调度）
- `INTERVAL_TIME`: 作业提交间隔时间（秒）
- `CHECK_TIME`: 作业检查间隔时间（秒）

### [PARAMETERS]
计算参数配置：
- `KSPACING`: K点间距（默认：0.13）
- `ECUTWFC`: 平面波截断能（默认：60，单位：Ry）
- `NSPIN`: 自旋设置（1=无自旋，2=自旋极化，默认：2）
- `CALCULATION`: 计算类型（scf, relax, nscf 等，默认：scf）

## 使用示例

### 修改配置文件

编辑 `config/condor.ini`：

```ini
[ABACUS]
ABACUS_DIR = /path/to/abacus/bin
PSEUDO_POTENTIAL_DIR = /path/to/pseudopotentials
ORBITAL_DIR = /path/to/orbitals

[ALLOW]
PARTITION = deimos
NODES = 1
CORES_PER_NODE = 64
```

### 在代码中使用配置

```python
from config import get

# 读取配置
abacus_dir = get('ABACUS', 'ABACUS_DIR')
partition = get('ALLOW', 'PARTITION', 'deimos')  # 带默认值
```

## 注意事项

1. 路径使用绝对路径，避免相对路径带来的问题
2. 修改配置后，无需重启程序，下次运行时会自动读取新配置
3. 如果配置文件不存在或某个配置项缺失，会使用代码中的默认值


