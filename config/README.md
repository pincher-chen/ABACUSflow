# ABACUSflow 配置文件说明

## 配置文件位置

配置文件使用 INI 格式。查找顺序如下（优先级从高到低）：

1. 环境变量 `ABACUSFLOW_CONFIG`（指向一个 ini 文件）
2. 环境变量 `ABACUSFLOW_CONFIG_DIR`（指向配置目录；既可以是 `.../config`，
   也可以是项目根目录 `.../proj`，程序会再找其下的 `config/condor.ini`）
3. 从当前工作目录开始，向上逐级查找的 `condor.ini` 或 `config/condor.ini`
   （因此可以把配置放在作业目录的父目录，便于一个项目共用一份配置）
4. 环境变量 `ABACUSFLOW_HOME` 下的 `config/condor.ini`
5. 程序自带的 `config/condor.ini`（位于安装目录内）

**无需修改安装目录里的文件**即可适配新集群：在作业目录放一份 `condor.ini`，
或直接用环境变量覆盖单项。

### 环境变量覆盖

任意配置项都可以用环境变量临时覆盖，格式为 `ABACUSFLOW_<大写键名>`：

```bash
export ABACUSFLOW_PSEUDO_POTENTIAL_DIR=/data/pp/sg15oncv
export ABACUSFLOW_ORBITAL_DIR=/data/orb/SG15std
export ABACUSFLOW_ABACUS_BIN=/opt/abacus/bin/abacus
export ABACUSFLOW_PARTITION=normal
export ABACUSFLOW_MODULES="gcc/12.3 openmpi/4.1"
export ABACUSFLOW_SCHEDULER=slurm          # 或 yh
export ABACUSFLOW_RUN_CMD="mpirun -np {cores} {cmd}"
```

环境变量优先级最高。必填项（`PSEUDO_POTENTIAL_DIR`、`ORBITAL_DIR`）缺失时
会明确报错并提示需要配置哪一项，不会静默使用其他集群的路径。

### [SCHEDULER] 调度器（新增）

- `TYPE`: `slurm`（默认）或 `yh`（天河系）
- `SUBMIT_CMD`: 提交命令，留空按 TYPE 推断（`sbatch` / `yhbatch`）
- `RUN_CMD`: 运行 MPI 的命令模板，支持 `{nodes}` `{cores}` `{partition}` `{cmd}`
  - 留空按 TYPE 推断（`srun` / `yhrun`）
  - 可替换为 `mpirun -np {cores} {cmd}` 或 Apptainer 形式：
    `apptainer exec --cleanenv --bind $PWD:/data --pwd /data $SIF {cmd}`
- `QUEUE_CMD`: 队列查询命令，默认 `squeue`

### [ENV] 环境初始化（新增）

- `PYTHON`: 生成脚本中调用的 Python，留空则用 `python`
- `SETUP`: 自定义初始化片段，原样插入生成脚本（`\n` 表示换行）
- `CONDA_PATH` / `CONDA_ENV`: 留空表示不使用 conda
- 生成的作业脚本**只包含你实际配置过的内容**，不再写入任何集群专有的
  module 初始化命令或 conda 环境名。

## 配置节说明

### [ABACUS]
ABACUS 相关配置：
- `ABACUS_BIN`: ABACUS 可执行文件路径（推荐；容器内可直接写 `abacus`）
- `ABACUS_DIR`: ABACUS 可执行文件所在目录（与 `ABACUS_EXE` 配合）
- `ABACUS_EXE`: ABACUS 可执行文件名（默认：abacus）
- `PSEUDO_POTENTIAL_DIR`: 赝势文件目录（绝对路径）
- `ORBITAL_DIR`: 轨道基组文件目录（绝对路径）

### [STRU]
结构文件相关配置：
- `PATH`: 结构文件搜索路径（可选）
- `SUFFIX`: 结构文件匹配模式，空格分隔多个（默认：*.vasp *.poscar *.cif *.mcif）

### [MODULE]
模块加载配置：
- `MODULES`: 需要加载的模块列表，用空格分隔
  - 例如：`intel/oneapi2023.2_noimpi mpi/mpich/4.1.2-icc-oneapi2023.2-ch4 elpa/2024.05.001-icc-oneapi2023.2-mpich-4.1.2-ch4`

### [ALLOW]
作业提交相关配置：
- `PARTITION`: 计算分区名称（如：normal, debug；留空则不传 -p）
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

[SCHEDULER]
TYPE = slurm

[ALLOW]
PARTITION = normal
NODES = 1
CORES_PER_NODE = 64
```

### 在代码中使用配置

```python
from config import get

# 读取配置
abacus_dir = get('ABACUS', 'ABACUS_DIR')
partition = get('ALLOW', 'PARTITION')  # 未配置时返回 None
```

## 注意事项

1. 路径使用绝对路径，避免相对路径带来的问题
2. 修改配置后，无需重启程序，下次运行时会自动读取新配置
3. 如果配置文件不存在或某个配置项缺失，会使用代码中的默认值


