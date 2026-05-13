# NBANDS 计算指南

## 🐛 问题描述

在 Band/Dos 计算中，出现 `NLOCAL < NBANDS` 错误，导致计算失败。

### 错误示例
```
ABACUS_NOTICE_ERROR
NLOCAL < NBANDS
```

**原因**：在 LCAO 基组下，设置的 `nbands` 超过了基函数总数 `NLOCAL`。

---

## 📚 理论基础

### 1. LCAO vs PW 的关键区别

| 基组类型 | 基函数数量 | NBANDS 上限 | 典型值 |
|---------|-----------|------------|--------|
| **PW** (平面波) | 成千上万 | 几乎无限制 | 可设置 100+ |
| **LCAO** (局域轨道) | 由 .orb 文件定义 | **NLOCAL** | 通常 20-100 |

**关键约束**：
```
NBANDS ≤ NLOCAL (LCAO 的硬性物理上限)
```

### 2. nspin=2 的正确理解

**常见误区**：
> "有自旋时，每个能带容纳1个电子，需要更多能带" ❌

**正确理解**：
- `nspin=2` 时，ABACUS 的 `nbands` 是指**每个自旋通道**（Spin Up 和 Spin Down）各自拥有的能带数
- **不需要**因为开启自旋就把 nbands 翻倍！

---

## ✅ 正确的计算方法

### 计算公式

```python
# 1. 计算占据能带数 (Occupied bands)
if nspin == 1:
    occ_bands = total_electrons / 2  # 每个能带容纳2个电子
else:  # nspin == 2
    occ_bands = total_electrons      # 每个能带容纳1个电子

# 2. 增加冗余量 (Buffer) 用于观察空带
buffer = 20  # 对于 LCAO，+10 到 +20 足够

nbands = int(occ_bands + buffer)

# 3. LCAO 的硬性约束
nbands = min(nbands, NLOCAL)
```

### 计算示例

#### 示例 1: Ir1C1 体系

**体系信息**：
- 元素：Ir (17 价电子) + C (4 价电子)
- 总电子数：21
- NLOCAL：40（基函数总数）
- nspin：2（有自旋）

**计算过程**：
```
1. occ_bands = 21 (nspin=2, 每带1电子)
2. nbands = 21 + 20 = 41
3. 检查约束：41 > 40 (NLOCAL)
4. 调整：nbands = 40 ✅
```

#### 示例 2: 无自旋体系

**体系信息**：
- 总电子数：20
- NLOCAL：50
- nspin：1（无自旋）

**计算过程**：
```
1. occ_bands = 20 / 2 = 10 (nspin=1, 每带2电子)
2. nbands = 10 + 20 = 30
3. 检查约束：30 < 50 (NLOCAL) ✅
4. 最终：nbands = 30
```

---

## 🔧 代码修复说明

### 修复前的错误逻辑

```python
if default_nspin == 1:
    nbands = int(total_ne / 2 * 2.0 + 20)  # 正确
else:
    nbands = int(max(total_ne * 2.0, total_ne + 40))  # ❌ 错误！
    # 示例：21 * 2 = 42 或 21 + 40 = 61 → 取 max = 61
    # 61 > NLOCAL=40 → 崩溃！
```

**问题**：
1. 对 nspin=2 的误解，使用了 `total_ne * 2.0`
2. 冗余量过大（+40）
3. 没有检查 NLOCAL 约束

### 修复后的正确逻辑

```python
# 1. 计算占据能带数
if default_nspin == 1:
    occ_bands = total_ne / 2  # 每带2电子
else:
    occ_bands = total_ne      # 每带1电子 ✅

# 2. 增加合理的冗余量
buffer = 20  # ✅ 不过大
nbands = int(occ_bands + buffer)

# 3. 确保 nbands 至少是最小值
nbands = max(nbands, 20)

# 4. LCAO 硬性约束
if nlocal is not None:
    nbands = min(nbands, nlocal)  # ✅ 关键修复！
```

---

## 🎯 使用建议

### 1. 让 ABACUS 自动计算（推荐）

对于 `Test_spin`, `Coarse_relax`, `Relax`, `Scf` 阶段：
- **不设置** nbands
- 让 ABACUS 自动计算，最安全

### 2. 手动设置（Band/Dos）

如果自动计算仍有问题，可以在 yaml 中手动设置：

```yaml
# config/template/Band.yaml
calculation: nscf
init_chg: 1
nbands: 40  # 手动设置，确保 <= NLOCAL
```

### 3. 查看 NLOCAL 值

```bash
# 从 Scf 的 running.log 中查看
grep "NLOCAL" batch_work2/hmat_103/Scf/OUT.*/running*.log
# 输出：NLOCAL = 40
```

---

## ⚠️ 常见错误排查

### 错误 1: `NLOCAL < NBANDS`

**原因**：nbands 设置过大

**解决方法**：
1. 检查 NLOCAL 值：`grep NLOCAL Scf/OUT.*/running*.log`
2. 手动设置 `nbands = NLOCAL` 或更小

### 错误 2: Band 能带图不完整

**原因**：nbands 设置过小，看不到高能级空态

**解决方法**：
1. 在 NLOCAL 约束下，适当增加 nbands
2. 如果需要更多能带，考虑使用更大的基组（如 DZP → TZP）

### 错误 3: 无法读取 NLOCAL

**现象**：
```
[WARNING] Could not read NLOCAL from Scf/OUT.*/running*.log
```

**原因**：Scf 阶段未完成或日志文件不存在

**解决方法**：
1. 确保 Scf 阶段已成功完成
2. 手动在 yaml 中设置保守的 nbands 值（如 30）

---

## 📖 参考资源

1. **ABACUS 官方文档**：http://abacus.ustc.edu.cn/
2. **INPUT 参数说明**：查看 ABACUS 源码 `docs/input-main.md`
3. **LCAO 基组说明**：`docs/basis-set.md`

---

## 📝 修改历史

| 日期 | 版本 | 修改内容 |
|------|------|---------|
| 2025-12-30 | v2.0 | 根据社区反馈重写 nbands 计算逻辑 |
| 2025-12-30 | v1.0 | 初次修复，添加 NLOCAL 检查 |

---

**贡献者**: AI Assistant  
**审核**: 用户专业建议  
**状态**: ✅ 已测试通过


