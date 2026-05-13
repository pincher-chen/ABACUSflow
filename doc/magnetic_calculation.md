# 磁性计算使用指南

本文说明如何用 abacusflow 为磁性材料准备输入文件，重点覆盖 mcif 的生成方式
（以 Alexandria 数据库为例）以及 `spin=2`（共线）和 `spin=4`（非共线/SOC）两种
计算模式下 STRU 的生成规则。

---

## 1. 磁矩信息的来源

VASP `ISPIN=2` 计算输出的是**每原子一个标量磁矩**，例如 Alexandria 数据库中的
JSON 字段：

```json
{
  "label": "Cr",
  "properties": {
    "magmom": -0.907,
    "charge": 10.863
  }
}
```

这个标量 `magmom` 是 VASP 内部量化轴（z 方向）上的投影，包含符号。
非共线（SOC，`LSORBIT=T`）计算才会输出三分量向量磁矩。

---

## 2. 生成 mcif 文件

### 2.1 为什么用 mcif

CIF 文件不包含磁矩字段。磁结构标准格式 **magCIF（.mcif）** 使用
`_atom_site_moment_crystalaxis_x/y/z` 字段存储每原子的磁矩，
pymatgen 的 `CifWriter(write_magmoms=True)` 可直接写出。

### 2.2 从 Alexandria 数据库批量生成

Alexandria 提供 `json.bz2` 格式的材料数据库。使用随附脚本：

```bash
# 处理单个 bz2 文件，筛选磁性材料，按类别输出 mcif
python analyze_convex_hull_magnetic_with_mag.py \
    --input alexandria_000.json.bz2 \
    --outdir ./convex_hull_mag_cifs \
    --strict-threshold 0.02 \
    --loose-threshold 0.1
```

输出目录结构：

```
convex_hull_mag_cifs/
├── ferromagnetic/       # |Σm| > strict_threshold
├── ferrimagnetic/       # 部分原子异号
├── antiferromagnetic/   # Σm ≈ 0，但各原子 |m| > threshold
└── weak_magnetic/       # 小磁矩
```

每个子目录下为 `<mat_id>.mcif` 文件。

### 2.3 mcif 内容示例

```
# generated using pymatgen
data_K2MnCoS4
_cell_length_a   13.13922997
...
loop_
 _atom_site_moment_label
 _atom_site_moment_crystalaxis_x
 _atom_site_moment_crystalaxis_y
 _atom_site_moment_crystalaxis_z
  Mn4   1.29418119   1.15440695  -0.69015242
  Co6  -1.02154035  -0.91121189   0.54476031
  S8   -0.05809436  -0.05182005   0.03098018
  ...

# === abacusflow: original scalar magmoms (uB, VASP ISPIN=2) ===
# scalar_magmom Mn4  2.250000
# scalar_magmom Co6 -1.776000
# scalar_magmom S8  -0.101000
```

> **注意**：mcif 的 crystalaxis 三分量是**晶轴系**投影，不是 Cartesian 坐标。
> 文件末尾的 `# scalar_magmom` 注释由 `write_mcif_with_magmom` 追加，
> 保留了原始 VASP 标量（含符号），供 `spin=2` 下游优先使用。

### 2.4 手动查看某材料原始磁矩

```bash
python read_and_search_bz2.py \
    --in convex_hull.json.bz2 \
    --search agm019989792
```

---

## 3. 生成 STRU 文件

### 3.1 `spin=2`（共线磁性，nspin=2）

模板 `Scf-Spin2.yaml` 设置 `nspin: 2`，STRU 中每原子写一行：

```
magmom <signed_value>
```

**磁矩值的取法**（优先级由高到低）：

1. 若 mcif 末尾有 `# scalar_magmom` 注释（abacusflow 生成的 mcif），
   直接使用该标量，**完整保留原始 VASP 符号和量值**。
2. 若无注释，从 crystalaxis 三分量反变换到 Cartesian，
   取 `sign(m_z) × |m|` 作为标量（符号依赖晶格取向，可能与 VASP 不一致）。

示例 STRU 片段（`agm019989792`，K₂MnCoS₄）：

```
Mn
0.0000000000
2
0.2500000000  0.2500000000  0.2500000000  1 1 1  magmom  2.250000
0.0000000000  0.0000000000  0.0000000000  1 1 1  magmom  2.250000
Co
0.0000000000
2
0.7500000000  0.7500000000  0.7500000000  1 1 1  magmom -1.776000
0.5000000000  0.5000000000  0.5000000000  1 1 1  magmom -1.776000
S
0.0000000000
8
0.0955030000  0.7973170000  0.7117300000  1 1 1  magmom -0.101000
...
```

#### 生成命令

```bash
python abacus.py generate \
    --work_dir <output_dir> \
    --stage Scf-Spin2 \
    --stru_file <mat_id>.mcif \
    --spin 2 \
    --kval 0.02
```

#### 批量生成（`single` 模式）

```bash
python abacus.py single \
    <mcif_input_dir> <output_root_dir> \
    -t Scf-Spin2 -k 0.02 -spin 2
```

### 3.2 `spin=4`（非共线磁性/SOC，nspin=4）

模板 `Scf-Soc.yaml` 设置 `nspin: 4`、`noncolin: 1`、`lspinorb: 1`。
STRU 中每原子写三分量 Cartesian 磁矩：

```
magmom <mx> <my> <mz>
```

**磁矩值的取法**：
从 mcif 的 crystalaxis 三分量经矩阵反变换得到 Cartesian 向量，
直接写入三分量。非共线计算要求 mcif 来源本身是非共线（SOC）数据，
从 VASP ISPIN=2 标量换算的结果不能正确表示非共线磁矩方向。

#### 生成命令

```bash
python abacus.py generate \
    --work_dir <output_dir> \
    --stage Scf-Soc \
    --stru_file <mat_id>.mcif \
    --spin 4 \
    --kval 0.02
```

### 3.3 `spin=1`（无磁性，nspin=1）

不写 `magmom`，STRU 原子行无磁矩字段。

---

## 4. 坐标变换说明

pymatgen `CifWriter` 的写入约定：

```
m_crys = m_cart @ inv(unit_M)
unit_M = M / |row_i(M)|        # M 的每行为晶格向量 a, b, c（单位归一化）
```

逆变换（crystalaxis → Cartesian）：

```python
m_cart = m_crys @ unit_M
```

代码实现见 `abacus/stru_utils.py: _crystal_axis_to_cartesian_magmom`。

**注意**：对于从 VASP ISPIN=2 标量生成的 mcif，
pymatgen 先把标量 `m` 当作 `(0, 0, m)` 在原始 Cartesian 系下转换，
但读取端（ASE）使用 CIF 标准化晶格（a 沿 x，b 在 xy 平面），
两个 Cartesian 系差一个旋转，导致逆变换后的 `m_z` 符号可能翻转。
这就是为什么必须保留 `# scalar_magmom` 注释——在 `spin=2` 模式下
优先从注释读取原始值，避免坐标系旋转导致的符号错误。

---

## 5. 批量提交示例

以 Alexandria `convex_hull_mag_cifs/` 下四类磁性材料的 spin=2 计算为例：

```bash
#!/bin/bash
set -euo pipefail

ABACUSFLOW=/XYFS01/nscc-gz_pinchen_1/Storage01/sf_box/abacusflow
INPUT_ROOT=/path/to/convex_hull_mag_cifs
OUT_ROOT=/path/to/convex_hull_mag_spin2

CLASSES=(antiferromagnetic ferrimagnetic ferromagnetic weak_magnetic)

# Stage 1: 生成输入文件和提交脚本
for cls in "${CLASSES[@]}"; do
    python "${ABACUSFLOW}/abacus.py" single \
        "${INPUT_ROOT}/${cls}" "${OUT_ROOT}/${cls}" \
        -t Scf-Spin2 -k 0.02 -spin 2
done

# Stage 2: 提交到队列
for cls in "${CLASSES[@]}"; do
    python "${ABACUSFLOW}/submit_jobs.py" "${OUT_ROOT}/${cls}" --single
done
```

---

## 6. 常见问题

**Q: STRU 里的 `magmom` 值符号跟 VASP 结果相反怎么办？**

原因是 mcif 是旧版本生成的（不含 `# scalar_magmom` 注释）。
重新用新版 `analyze_convex_hull_magnetic_with_mag.py` 生成 mcif 即可。
旧 mcif 的问题只影响 `spin=2` 的初始磁矩符号；
对于能量、|m|、磁性类型分类，结果不受影响（能量在时间反演下简并）。

**Q: spin=2 的 `magmom` 初始值和 VASP 最终值一致就够了吗？**

是的。`spin=2` 的 `magmom` 仅用于 SCF 迭代的初始猜测，
ABACUS 会自洽到正确的基态；初始值的绝对量值和符号 pattern 对即可。

**Q: 什么情况下需要 spin=4？**

当材料涉及强自旋轨道耦合效应（拓扑绝缘体、重元素体系等）或
原始数据已含非共线磁结构（如 magndata 数据库）时需要 `spin=4`。
Alexandria 数据来自 VASP ISPIN=2，本身是共线计算，通常用 `spin=2`。
