# 自定义参数

`config/template/*.yaml` 支持写入任意 ABACUS 参数，不受代码内预定义列表限制。
自定义参数会被识别并写入 INPUT 文件末尾，并在生成时打印提示：

```
[INFO] Custom parameter added: deepks_out_labels = 1
```

## 使用方法

直接在对应阶段的 yaml 文件中添加参数即可：

```yaml
# config/template/Scf.yaml
calculation: scf
basis_type: lcao
ecutwfc: 100

# 以下为自定义参数示例
deepks_out_labels: 1
out_proj_band: 1
vdw_method: d3_bj
```

## 常用示例

```yaml
# DeePKS
deepks_out_labels: 1
deepks_scf: 1
deepks_model: model.ptg

# 输出控制
out_proj_band: 1       # 投影能带
out_bandgap: 1         # 带隙信息
out_wfc_lcao: 1        # LCAO 波函数
out_mat_hs: 1          # 哈密顿和重叠矩阵
out_mat_r: 1           # 位置算符矩阵

# DFT+U
dft_plus_u: 1
orbital_corr: -1 3 -1  # 对第 2 种元素的 d 轨道
hubbard_u: 0 5.0 0

# 杂化泛函
exx_hybrid_alpha: 0.25

# 范德华修正
vdw_method: d3_bj

# 磁性（非共线/SOC）
nspin: 4
noncolin: 1
lspinorb: 1
```

## 注意事项

- 参数名称须与 ABACUS 官方文档一致，拼写错误会导致 ABACUS 忽略或报错
- 代码中硬编码的参数（如 `nspin`）优先级高于 yaml 中的设置，
  建议通过 CLI 参数（`--spin`）控制，而不是在 yaml 中修改 `nspin`
- 确认目标 ABACUS 版本支持该参数
