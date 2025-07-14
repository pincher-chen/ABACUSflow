# ABACUSflow

ABACUSflow是一个专为[ABACUS](https://abacus.ustc.edu.cn/)第一性原理计算软件设计的自动化工作流工具包。它提供了高通量计算、结果后处理、数据格式转换等功能，支持与ASE（Atomic Simulation Environment）的无缝集成。

## ✨ 主要特性

- **🚀 高通量计算**: 批量处理结构文件，自动生成ABACUS计算输入
- **🔄 ASE集成**: 完整的ABACUS与ASE接口，支持结构优化和属性计算
- **📊 后处理分析**: 能带结构、态密度、磁性等电子结构属性分析
- **🔧 格式转换**: 支持VASP、CIF、POSCAR等多种格式转换
- **📈 数据可视化**: 内置能带图、DOS图等可视化功能
- **🧹 清理工具**: 自动清理失败计算和临时文件

## 📁 项目结构

```
ABACUSflow/
├── abacus/                    # 核心ABACUS接口模块
│   ├── abacus_out.py         # ABACUS计算器类
│   ├── basis_set.py          # 基组管理
│   └── example/              # 使用示例
├── abacus_helper/            # 辅助工具和后处理
│   ├── core/                 # 核心处理模块
│   │   ├── bandstructure.py  # 能带结构分析
│   │   ├── dos.py           # 态密度分析
│   │   ├── structure.py     # 结构处理
│   │   └── optimize.py      # 优化分析
│   ├── utils/               # 实用工具
│   └── main.py              # 主程序入口
├── postprocess/             # 后处理分析工具
├── InputPoscar/            # 示例输入结构
├── OrbSG15std/             # SG15轨道文件库
├── PotSG15/                # SG15赝势文件库
├── calc.py                 # 批量计算脚本
├── abacusHT.py            # 高通量计算脚本
└── clean_null.py          # 清理工具
```

## 🛠️ 安装要求

### 依赖软件

- Python 3.7+
- [ABACUS](https://abacus.ustc.edu.cn/) 软件包
- [ASE](https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment)

### Python依赖包

```bash
pip install -r abacus_helper/requirements.txt
```

主要依赖：

- `numpy`, `scipy` - 数值计算
- `matplotlib`, `plotly` - 数据可视化
- `pymatgen` - 材料科学工具
- `pandas` - 数据处理
- `ase` - 原子模拟环境

## 🚀 快速开始

### 1. 高通量批量计算

```python
# 使用abacusHT.py进行批量计算
python abacusHT.py /path/to/structures /path/to/calculations
```

### 2. ASE接口使用示例

```python
from ase.io import read
from abacus.abacus_out import Abacus

# 读取结构文件
atoms = read('structure.cif')

# 设置ABACUS计算器
calc = Abacus(
    label='calculation',
    pseudo_dir='./PotSG15/',
    basis_dir='./OrbSG15std/',
    ecutwfc=60,
    calculation='scf',
    nspin=2
)

atoms.set_calculator(calc)
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
```

### 3. 后处理分析

```python
from abacus_helper.main import *

# 获取能带结构数据
band_data = get_bandstructure(
    band_dir='./Band/',
    stru_filename='structure', 
    scf_log_filepath='./SCF/OUT.ABACUS/running_scf.log'
)

# 获取态密度数据
dos_data = get_density_of_states('./DOS/')

# 获取优化后的结构
final_structure = get_optimized_stru_filepath('./RELAX/')
```

## 📋 支持的计算类型

### 基础计算

- **SCF**: 自洽场计算
- **NSCF**: 非自洽场计算（能带）
- **RELAX**: 离子位置优化
- **CELL_RELAX**: 晶格参数优化

### 电子结构属性

- **能带结构**: 沿高对称k点路径的能带计算
- **态密度**: 总态密度和投影态密度
- **磁性**: 磁矩和自旋极化计算

### 数据格式支持

- **输入格式**: VASP POSCAR, CIF, XYZ
- **输出格式**: 能带数据、DOS数据、结构优化轨迹

## 🔧 命令行工具

### abacus_helper 主程序

```bash
# 获取结构属性数据
python main.py -t sp -s ./SCF/structure_file

# 获取能带数据  
python main.py -t band -d ./Band/ -l ./SCF/OUT.ABACUS/running_scf.log -s structure

# 获取态密度数据
python main.py -t dos -d ./DOS/

# 获取磁性数据
python main.py -t mag -d ./SCF/

# 生成k点路径
python main.py -t kpath -s structure_file -n 20
```

### 批量计算脚本

```bash
# 高通量计算（带参数）
python abacusHT.py /path/to/structures /path/to/calculations

# 固定路径批量计算
python calc.py
```

### 清理工具

```bash
# 清理失败的计算
python clean_null.py failed_list.txt
```

## 📊 数据处理功能

### 能带结构分析

- 自动识别能隙
- 价带顶和导带底定位
- 金属性判断
- 能带可视化

### 态密度分析

- 总态密度计算
- 原子和轨道投影态密度
- 自旋极化DOS
- DOS积分和可视化

### 结构优化分析

- 优化轨迹跟踪
- 能量收敛分析
- 最终结构提取

## 🎯 应用场景

- **材料筛选**: 高通量计算大量候选材料
- **电子结构研究**: 详细的能带和DOS分析
- **结构优化**: 自动化的几何优化流程
- **数据挖掘**: 批量提取和整理计算结果
- **工作流集成**: 与其他计算工具的接口

## 📚 示例和教程

项目包含多个使用示例：

- `abacus/example/property/ase_abacus.py` - ASE接口使用示例
- `abacus/example/neb/` - 过渡态计算示例
- `abacus_helper/test/` - 后处理分析示例
- `postprocess/test/` - 完整的计算流程示例

## 🤝 贡献

欢迎提交Issues和Pull Requests来改进项目！

## 📄 许可证

本项目采用开源许可证，具体请查看LICENSE文件。

## 🔗 相关链接

- [ABACUS官网](https://abacus.ustc.edu.cn/)
- [ASE文档](https://wiki.fysik.dtu.dk/ase/)
- [Pymatgen文档](https://pymatgen.org/)

---

如有问题或建议，请通过GitHub Issues联系我们。
