# 自定义参数使用说明

## 📋 概述

现在 `config/template/*.yaml` 文件支持添加**任意 ABACUS 参数**，不再受代码预定义列表限制。

## ✨ 功能说明

### 修改前
- ❌ 只能使用代码中预定义的参数
- ❌ 添加新参数会抛出 `TypeError: Parameter not defined` 异常
- ❌ 需要修改 `create_input.py` 的参数列表才能添加新参数

### 修改后
- ✅ 支持所有 ABACUS 官方参数
- ✅ 支持自定义参数（用于测试新功能）
- ✅ 自动识别并添加到 INPUT 文件
- ✅ 添加时会显示提示信息：`[INFO] Custom parameter added: param_name = value`

## 📝 使用方法

### 1. 在 YAML 模板中添加参数

编辑 `config/template/Scf.yaml`（或其他阶段的 yaml）：

```yaml
# 标准参数
calculation: scf
basis_type: lcao
ecutwfc: 100
scf_thr: 1.0e-7
scf_nmax: 100

# ========================================
# 新增自定义参数（ABACUS 新版本功能）
# ========================================

# DeePKS 相关
deepks_out_labels: 1
deepks_scf: 1
deepks_bandgap: 1

# 输出控制
out_proj_band: 1        # 输出投影能带
out_wfc_lcao: 1         # 输出 LCAO 波函数
out_bandgap: 1          # 输出带隙信息

# 对角化参数
pw_diag_thr: 1.0e-2     # PW 对角化阈值
diago_david_ndim: 4     # Davidson 维度

# 电荷外推
chg_extrap: atomic      # 电荷外推方法

# DFT+U 参数
dft_plus_u: 1
orbital_corr: -1 3 -1   # 对第2种元素的 d 轨道应用 U
hubbard_u: 0 5.0 0      # U 值（eV）

# 任何您需要的 ABACUS 参数...
```

### 2. 参数会自动生成到 INPUT 文件

运行工作流后，自定义参数会自动写入 INPUT 文件：

```
INPUT_PARAMETERS
# Created by Atomic Simulation Enviroment
calculation         scf
basis_type          lcao
ecutwfc             100
scf_thr             1.0e-7

# 自定义参数会出现在文件末尾
deepks_out_labels   1
out_proj_band       1
pw_diag_thr         0.01
chg_extrap          atomic
```

## 🎯 常用自定义参数示例

### DeePKS 机器学习
```yaml
deepks_out_labels: 1
deepks_scf: 1
deepks_model: model.ptg
```

### 输出控制
```yaml
out_proj_band: 1       # 投影能带
out_bandgap: 1         # 带隙信息
out_wfc_lcao: 1        # LCAO 波函数
out_mat_hs: 1          # 哈密顿和重叠矩阵
out_mat_r: 1           # 位置算符矩阵
printe: 100            # 每 100 步打印能量
```

### DFT+U 计算
```yaml
dft_plus_u: 1
orbital_corr: -1 3 -1  # 对第2种元素的d轨道
hubbard_u: 0 5.0 0     # U值（eV）
```

### 杂化泛函
```yaml
exx_hybrid_alpha: 0.25     # PBE0: 25% Fock交换
exx_separate_loop: 0       # 是否分离 EXX 循环
```

### 范德华修正
```yaml
vdw_method: d3_bj         # DFT-D3(BJ) 修正
vdw_s6: 1.0               # s6 参数
vdw_s8: 0.722             # s8 参数
```

### 磁性计算
```yaml
nspin: 4                  # 非共线磁性
noncolin: 1               # 非共线计算
lspinorb: 1               # 自旋轨道耦合
starting_magnetization: 0.5 1.0  # 初始磁矩
```

## 🔍 验证参数是否生效

### 方法 1：查看生成脚本的输出
```bash
python abacus.py generate --work_dir Test_spin --stage Test_spin
```

如果看到：
```
[INFO] Custom parameter added: deepks_out_labels = 1
```
说明自定义参数已被识别。

### 方法 2：检查生成的 INPUT 文件
```bash
cat Test_spin/INPUT | tail -20
```

确认自定义参数出现在文件末尾。

## ⚠️ 注意事项

1. **参数名称必须正确**
   - 参数名要与 ABACUS 官方文档一致
   - 拼写错误会导致 ABACUS 忽略该参数或报错

2. **参数值格式要正确**
   - 数值类型：`ecutwfc: 100`
   - 字符串：`calculation: scf`
   - 布尔值：`symmetry: 1` 或 `symmetry: 0`
   - 列表：`orbital_corr: -1 3 -1`（用空格分隔）

3. **兼容性**
   - 确保您的 ABACUS 版本支持该参数
   - 建议查阅 ABACUS 官方文档确认参数可用性

4. **参数优先级**
   - 代码中硬编码的参数会覆盖 YAML 中的设置
   - 例如：`generator.py` 中对 `nspin` 的处理

## 📚 参考资源

- **ABACUS 官方文档**: http://abacus.ustc.edu.cn/
- **INPUT 参数列表**: 查看 ABACUS 源码 `docs/input-main.md`
- **本项目预定义参数**: `abacus/create_input.py` 第 37-281 行

## 🔧 代码修改说明

修改了 `abacus/create_input.py` 的 3 处：

1. **添加 `custom_params` 字典**（第 316 行）
   ```python
   self.custom_params = {}  # 用于存储任意自定义参数
   ```

2. **修改 `set()` 方法**（第 401-404 行）
   ```python
   else:
       # 允许任意自定义参数
       self.custom_params[key] = kwargs[key]
       print(f"[INFO] Custom parameter added: {key} = {kwargs[key]}")
   ```

3. **修改 `write_input_input()` 方法**（第 496-499 行）
   ```python
   # 写入自定义参数
   for key, val in self.custom_params.items():
       if val is not None:
           params = str(key) + ' ' * (20 - len(key)) + str(val)
           input_file.write('%s\n' % params)
   ```

## ✅ 测试结果

已通过完整测试：
- ✅ 预定义参数正常工作
- ✅ 自定义参数成功添加
- ✅ 从 YAML 加载自定义参数
- ✅ INPUT 文件正确生成

---

**更新日期**: 2025-12-30  
**修改人**: AI Assistant  
**版本**: v1.0


