# nbands 计算说明

## 问题背景

在 LCAO 基组的 Band/Dos 计算中，如果 `nbands` 超过基函数总数 `NLOCAL`，
ABACUS 会报错退出：

```
ABACUS_NOTICE_ERROR
NLOCAL < NBANDS
```

LCAO 基函数数量由轨道文件（`.orb`）决定，通常在几十到几百之间，远小于平面波的基函数数量，
因此 `nbands` 必须满足 `nbands ≤ NLOCAL`。

## 计算公式

```python
# nspin=1：每条能带容纳 2 个电子
occ_bands = total_electrons / 2

# nspin=2：每个自旋通道各有 nbands 条能带，不需要翻倍
occ_bands = total_electrons

nbands = int(occ_bands + buffer)   # buffer = 20（额外空带）
nbands = min(nbands, nlocal)       # LCAO 硬约束
nbands = max(nbands, 20)           # 最小值保障
```

`nspin=2` 时，`nbands` 是**每个自旋通道**的能带数，不需要因为有自旋就乘以 2。

## 查看 NLOCAL

```bash
grep "NLOCAL" Scf/OUT.*/running*.log
```

## 手动设置 nbands

自动计算不满足需求时，可在对应 yaml 中手动指定：

```yaml
# config/template/Band.yaml
nbands: 40  # 确保 <= NLOCAL
```

## 常见错误

**`NLOCAL < NBANDS`**：检查 NLOCAL 值，降低 `nbands` 或使用更大的基组。

**能带图空带不够**：在 NLOCAL 约束内适当增大 `buffer`，或换用更大的基组（如 DZP → TZP）。

**无法读取 NLOCAL**：说明 Scf 阶段未正常完成，先确保 Scf 收敛后再运行 Band/Dos。
