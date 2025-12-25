#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K点计算工具模块
统一的 K点网格生成逻辑，供 generate 和 update 命令共享使用
"""

import numpy as np
from math import pi, ceil


def calculate_kpoints_from_cell(cell, kspacing, ktype='Gamma'):
    """
    从晶胞矩阵计算 K点网格
    
    Parameters:
    -----------
    cell : array-like or ase.cell.Cell
        晶胞矩阵，3x3 数组，每行是一个晶格矢量（单位：Å）
    kspacing : float
        倒空间 K点间距参数（Å⁻¹），与 VASP 的 KSPACING 一致
        越小越密集，典型值：0.02-0.10
    ktype : str
        K点类型，'Gamma' 或 'Line'
    
    Returns:
    --------
    kpoints : list of int
        [k1, k2, k3] K点网格
    """
    if ktype == 'Line':
        # 能带计算，kspacing 表示点数
        npoints = int(kspacing) if kspacing >= 1 else 20
        return [npoints]
    
    # Gamma 网格：计算倒易晶格
    # ASE Cell.reciprocal() 返回的是物理学约定（包含 2π）
    if hasattr(cell, 'reciprocal'):
        # ASE Cell 对象
        reciprocal_cell = cell.reciprocal()
    else:
        # 普通数组，手动计算倒易晶格（物理学约定，包含 2π）
        # b_i = 2π × (a_j × a_k) / V
        cell = np.array(cell)
        volume = np.abs(np.dot(cell[0], np.cross(cell[1], cell[2])))
        reciprocal_cell = np.zeros((3, 3))
        reciprocal_cell[0] = 2 * pi * np.cross(cell[1], cell[2]) / volume
        reciprocal_cell[1] = 2 * pi * np.cross(cell[2], cell[0]) / volume
        reciprocal_cell[2] = 2 * pi * np.cross(cell[0], cell[1]) / volume
    
    # 计算 K点数（与 VASP 相同的公式）
    # k_i = floor(|b_i| / (2π × kspacing))
    # 其中 |b_i| 是倒易晶格矢量的模（Å⁻¹）
    Kpoints = []
    for i in range(3):
        b_norm = np.linalg.norm(reciprocal_cell[i])
        if b_norm > 0:
            # 使用 floor 和 round，与 VASP 一致
            k = int(np.floor(np.round(b_norm / (2 * pi * kspacing))))
            k = max(1, k)  # 至少为 1
        else:
            k = 1
        Kpoints.append(k)
    
    # 检测真空层/层状结构
    cell_lengths = [np.linalg.norm(cell[i]) for i in range(3)]
    
    for i in range(3):
        # 条件：
        # 1. 晶格常数 > 15 Å
        # 2. K点数 > 50
        # 3. 该方向远大于其他方向的平均值（> 5倍）
        other_lengths = [cell_lengths[j] for j in range(3) if j != i]
        avg_other = sum(other_lengths) / len(other_lengths) if other_lengths else 1.0
        
        if cell_lengths[i] > 15.0 and Kpoints[i] > 50 and cell_lengths[i] > 5 * avg_other:
            print(f"[WARNING] Direction {i}: lattice = {cell_lengths[i]:.2f} Å, K-points = {Kpoints[i]}")
            print(f"[WARNING] Detected vacuum/layered direction, setting K-point {i} = 1")
            Kpoints[i] = 1
    
    return Kpoints


def calculate_kpoints_from_stru(stru_file, kspacing, ktype='Gamma'):
    """
    从 STRU 文件计算 K点网格
    
    Parameters:
    -----------
    stru_file : str
        STRU 文件路径
    kspacing : float
        倒空间 K点间距参数（Å⁻¹），与 VASP 的 KSPACING 一致
    ktype : str
        K点类型
    
    Returns:
    --------
    kpoints : list of int
        [k1, k2, k3] K点网格
    """
    # 读取 STRU 文件
    with open(stru_file, 'r') as f:
        lines = f.readlines()
    
    # 解析晶格信息
    lattice_constant = 1.0
    cell_vectors = []
    
    for i, line in enumerate(lines):
        if 'LATTICE_CONSTANT' in line:
            if i + 1 < len(lines):
                try:
                    lattice_constant = float(lines[i + 1].strip())
                except:
                    pass
        elif 'LATTICE_VECTORS' in line:
            for j in range(1, 4):
                if i + j < len(lines):
                    try:
                        vector = [float(x) for x in lines[i + j].split()[:3]]
                        cell_vectors.append(vector)
                    except:
                        pass
            break
    
    if len(cell_vectors) < 3:
        return [1, 1, 1]
    
    # 构建晶胞矩阵（单位：Å）
    cell = np.array(cell_vectors) * lattice_constant
    
    return calculate_kpoints_from_cell(cell, kspacing, ktype)


def update_kpt_file(kpt_file, kpoints):
    """
    更新 KPT 文件
    
    Parameters:
    -----------
    kpt_file : str
        KPT 文件路径
    kpoints : list of int
        新的 K点网格 [k1, k2, k3]
    """
    with open(kpt_file, 'r') as f:
        lines = f.readlines()
    
    # 修改第4行（K点数）
    if len(lines) > 3:
        old_kpts = lines[3].split()
        if len(old_kpts) >= 6:
            # 保持后3个参数不变（shift）
            lines[3] = f"{kpoints[0]}  {kpoints[1]}  {kpoints[2]}  {old_kpts[3]}  {old_kpts[4]}  {old_kpts[5]}\n"
        else:
            lines[3] = f"{kpoints[0]}  {kpoints[1]}  {kpoints[2]}  0  0  0\n"
    
    with open(kpt_file, 'w') as f:
        f.writelines(lines)
    
    print(f"[INFO] K-points updated: {kpoints[0]} {kpoints[1]} {kpoints[2]}")


if __name__ == '__main__':
    # 测试
    import sys
    if len(sys.argv) > 2:
        stru_file = sys.argv[1]
        kval = float(sys.argv[2])
        kpoints = calculate_kpoints_from_stru(stru_file, kval)
        print(f"K-points for kval={kval}: {kpoints}")

