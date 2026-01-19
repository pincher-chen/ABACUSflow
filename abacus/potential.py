# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 20:02:39 2018

@author: Shen Zhen-Xiong

Updated: 2026-01 - Added dynamic pseudopotential selection for SG15 ONCV
"""

import re
from pathlib import Path
from typing import Optional, Tuple, List

# ============================================================================
# 动态赝势选择函数（用于 sg15oncv 等新格式目录）
# ============================================================================

# 版本号解析正则表达式
_VER_RE = re.compile(r"-(\d+)\.(\d+)\.upf$", re.IGNORECASE)


def parse_version(fname: str) -> Tuple[int, int]:
    """
    解析赝势文件名中的版本号
    
    Args:
        fname: 文件名，如 "Ge_ONCV_PBE-1.2.upf"
    
    Returns:
        (major, minor) 版本元组，如 (1, 2)
        解析失败返回 (-1, -1)
    """
    m = _VER_RE.search(fname)
    if not m:
        return (-1, -1)
    return (int(m.group(1)), int(m.group(2)))


def pick_upf(pseudo_dir: str, element: str, need_fr: bool = False) -> str:
    """
    为指定元素选择最佳赝势文件
    
    规则：
    - SOC 计算 (need_fr=True): 必须选 FR 版本，选最高版本号
    - 非 SOC 计算 (need_fr=False): 选非 FR 版本，选最高版本号
    
    Args:
        pseudo_dir: 赝势目录路径
        element: 元素符号（如 "Ge"）
        need_fr: 是否需要 FR（全相对论）赝势用于 SOC 计算
    
    Returns:
        选中的赝势文件名（不含路径）
    
    Raises:
        FileNotFoundError: 赝势目录不存在
        RuntimeError: 找不到合适的赝势
    """
    d = Path(pseudo_dir)
    if not d.is_dir():
        raise FileNotFoundError(f"pseudo_dir not found: {pseudo_dir}")

    # 收集候选文件
    fr_candidates: List[str] = []
    nonfr_candidates: List[str] = []
    
    for p in d.iterdir():
        if not p.is_file():
            continue
        name = p.name
        
        # 严格匹配元素前缀: {Element}_ONCV_PBE
        if not name.startswith(f"{element}_ONCV_PBE"):
            continue
        
        if name.lower().endswith(".upf"):
            if "_FR-" in name:
                fr_candidates.append(name)
            else:
                nonfr_candidates.append(name)

    # 按版本号降序排序（最高版本在前）
    fr_candidates.sort(key=parse_version, reverse=True)
    nonfr_candidates.sort(key=parse_version, reverse=True)

    if need_fr:
        # SOC 计算：必须使用 FR 赝势
        if fr_candidates:
            return fr_candidates[0]
        raise RuntimeError(
            f"[SOC ERROR] No FR pseudopotential for element '{element}' in {pseudo_dir}\n"
            f"SOC calculation requires FR (fully relativistic) pseudopotentials.\n"
            f"Available non-FR files: {nonfr_candidates[:3] if nonfr_candidates else 'None'}"
        )
    else:
        # 非 SOC 计算：优先使用非 FR 赝势
        if nonfr_candidates:
            return nonfr_candidates[0]
        # 回退：允许使用 FR 但警告
        if fr_candidates:
            print(f"[WARNING] No non-FR pseudopotential for '{element}', using FR: {fr_candidates[0]}")
            return fr_candidates[0]
        raise RuntimeError(
            f"No pseudopotential for element '{element}' in {pseudo_dir}"
        )


def input_need_fr(input_text: str) -> bool:
    """
    从 INPUT 文件内容判断是否需要 FR 赝势（SOC 计算）
    
    Args:
        input_text: INPUT 文件的文本内容
    
    Returns:
        True 如果 lspinorb=1（需要 SOC），否则 False
    """
    m = re.search(r"^\s*lspinorb\s+(\d+)\s*$", input_text, flags=re.IGNORECASE | re.MULTILINE)
    if not m:
        return False
    return int(m.group(1)) == 1


def get_dynamic_potential(pseudo_dir: str, atoms_list: List[str], need_fr: bool = False) -> List[str]:
    """
    为元素列表动态选择赝势文件
    
    Args:
        pseudo_dir: 赝势目录路径
        atoms_list: 元素符号列表，如 ['Ge', 'In', 'Te']
        need_fr: 是否需要 FR 赝势
    
    Returns:
        赝势文件路径列表（完整路径）
    """
    potentials = []
    for elem in atoms_list:
        upf_name = pick_upf(pseudo_dir, elem, need_fr)
        potentials.append(str(Path(pseudo_dir) / upf_name))
    return potentials


def is_dynamic_pseudo_dir(pseudo_dir: str) -> bool:
    """
    判断赝势目录是否需要使用动态选择逻辑
    
    检测目录中是否有 {Elem}_ONCV_PBE-*.upf 格式的文件
    
    Args:
        pseudo_dir: 赝势目录路径
    
    Returns:
        True 如果目录包含版本化的 ONCV 赝势文件
    """
    d = Path(pseudo_dir)
    if not d.is_dir():
        return False
    
    # 检查是否有 _ONCV_PBE- 格式的文件
    for p in d.iterdir():
        if p.is_file() and "_ONCV_PBE" in p.name and p.name.endswith(".upf"):
            return True
    return False


# ============================================================================
# 静态赝势字典（向后兼容）
# ============================================================================

PotentialDict = {'PotLDA': {'H': 'H.LDA.UPF',
                            'Li': 'Li.LDA.UPF',
                            'Be': 'Be.LDA.UPF',
                            'B': 'B.LDA.UPF',
                            'C': 'C.pz-vbc.UPF',
                            'N': 'N.LDA.UPF',
                            'O': 'O.LDA.UPF',
                            'F': 'F.LDA.UPF',
                            'Na': 'Na.LDA.UPF',
                            'Mg': 'Mg.LDA.UPF',
                            'Al': 'Al.pz-vbc.UPF',
                            'Si': 'Si.LDA.UPF',
                            'P': 'P.LDA.UPF',
                            'S': 'S.LDA.UPF',
                            'Cl': 'Cl.pz-bhs.UPF',
                            'K': 'K.LDA.UPF',
                            'Ca': 'Ca.pz-n-vbc.UPF',
                            'Ti': 'Ti.LDA.UPF',
                            'Mn': 'Mn.LDA.UPF',
                            'Fe': 'Fe.LDA.UPF',
                            'Cu': 'Cu.LDA.UPF',
                            'Zn': 'Zn.LDA.UPF',
                            'Ga': 'Ga.LDA.UPF',
                            'Ge': 'Ge.pz-bhs.UPF',
                            'As': 'As.pz-bhs.UPF',
                            'Se': 'Se.pz-bhs.UPF',
                            'Br': 'Br.LDA.UPF',
                            'In': 'In.pz-bhs.UPF',
                            'Sn': 'Sn.pz-bhs.UPF',
                            'I': 'I.LDA.UPF'},
                 'PotPBE': {'H': 'H.PBE.UPF',
                            'Li': 'Li.PBE.UPF',
                            'Be': 'Be.PBE.UPF',
                            'B': 'B.PBE.UPF',
                            'C': 'C.PBE.UPF',
                            'N': 'N.PBE.UPF',
                            'O': 'O.PBE.UPF',
                            'F': 'F.PBE.UPF',
                            'Na': 'Na.PBE.UPF',
                            'Mg': 'Mg.PBE.UPF',
                            'Al': 'Al.PBE.UPF',
                            'Si': 'Si.PBE.UPF',
                            'P': 'P.PBE.UPF',
                            'S': 'S.PBE.UPF',
                            'Cl': 'Cl.PBE.UPF',
                            'K': 'K.PBE.UPF',
                            'Ti': 'Ti.PBE.UPF',
                            'Mn': 'Mn.PBE.UPF',
                            'Fe': 'Fe.PBE.UPF',
                            'Cu': 'Cu.PBE.UPF',
                            'Zn': 'Zn.pbe.UPF',
                            'Ga': 'Ga.PBE.UPF',
                            'Br': 'Br.PBE.UPF',
                            'I': 'I.PBE.UPF'},
                 'PotSG15': {'Ag': 'Ag_ONCV_PBE-1.0.upf',
                                'Al': 'Al_ONCV_PBE-1.0.upf',
                                'Ar': 'Ar_ONCV_PBE-1.0.upf',
                                'As': 'As_ONCV_PBE-1.0.upf',
                                'Au': 'Au_ONCV_PBE-1.0.upf',
                                'B': 'B_ONCV_PBE-1.0.upf',
                                'Ba': 'Ba_ONCV_PBE-1.0.upf',
                                'Be': 'Be_ONCV_PBE-1.0.upf',
                                'Bi': 'Bi_ONCV_PBE-1.0.upf',
                                'Br': 'Br_ONCV_PBE-1.0.upf',
                                'C': 'C_ONCV_PBE-1.0.upf',
                                'Ca': 'Ca_ONCV_PBE-1.0.upf',
                                'Cd': 'Cd_ONCV_PBE-1.0.upf',
                                'Cl': 'Cl_ONCV_PBE-1.0.upf',
                                'Co': 'Co_ONCV_PBE-1.0.upf',
                                'Cr': 'Cr_ONCV_PBE-1.0.upf',
                                'Cs': 'Cs_ONCV_PBE-1.0.upf',
                                'Cu': 'Cu_ONCV_PBE-1.0.upf',
                                'F': 'F_ONCV_PBE-1.0.upf',
                                'Fe': 'Fe_ONCV_PBE-1.0.upf',
                                'Ga': 'Ga_ONCV_PBE-1.0.upf',
                                'Ge': 'Ge_ONCV_PBE-1.0.upf',
                                'H': 'H_ONCV_PBE-1.0.upf',
                                'He': 'He_ONCV_PBE-1.0.upf',
                                'Hf': 'Hf_ONCV_PBE-1.0.upf',
                                'Hg': 'Hg_ONCV_PBE-1.0.upf',
                                'I': 'I_ONCV_PBE-1.0.upf',
                                'In': 'In_ONCV_PBE-1.0.upf',
                                'Ir': 'Ir_ONCV_PBE-1.0.upf',
                                'K': 'K_ONCV_PBE-1.0.upf',
                                'Kr': 'Kr_ONCV_PBE-1.0.upf',
                                'La': 'La_ONCV_PBE-1.0.upf',
                                'Li': 'Li_ONCV_PBE-1.0.upf',
                                'Mg': 'Mg_ONCV_PBE-1.0.upf',
                                'Mn': 'Mn_ONCV_PBE-1.0.upf',
                                'Mo': 'Mo_ONCV_PBE-1.0.upf',
                                'N': 'N_ONCV_PBE-1.0.upf',
                                'Na': 'Na_ONCV_PBE-1.0.upf',
                                'Nb': 'Nb_ONCV_PBE-1.0.upf',
                                'Ne': 'Ne_ONCV_PBE-1.0.upf',
                                'Ni': 'Ni_ONCV_PBE-1.0.upf',
                                'O': 'O_ONCV_PBE-1.0.upf',
                                'Os': 'Os_ONCV_PBE-1.0.upf',
                                'P': 'P_ONCV_PBE-1.0.upf',
                                'Pb': 'Pb_ONCV_PBE-1.0.upf',
                                'Pd': 'Pd_ONCV_PBE-1.0.upf',
                                'Pt': 'Pt_ONCV_PBE-1.0.upf',
                                'Rb': 'Rb_ONCV_PBE-1.0.upf',
                                'Re': 'Re_ONCV_PBE-1.0.upf',
                                'Rh': 'Rh_ONCV_PBE-1.0.upf',
                                'Ru': 'Ru_ONCV_PBE-1.0.upf',
                                'S': 'S_ONCV_PBE-1.0.upf',
                                'Sb': 'Sb_ONCV_PBE-1.0.upf',
                                'Sc': 'Sc_ONCV_PBE-1.0.upf',
                                'Se': 'Se_ONCV_PBE-1.0.upf',
                                'Si': 'Si_ONCV_PBE-1.0.upf',
                                'Sn': 'Sn_ONCV_PBE-1.0.upf',
                                'Sr': 'Sr_ONCV_PBE-1.0.upf',
                                'Ta': 'Ta_ONCV_PBE-1.0.upf',
                                'Tc': 'Tc_ONCV_PBE-1.0.upf',
                                'Te': 'Te_ONCV_PBE-1.0.upf',
                                'Ti': 'Ti_ONCV_PBE-1.0.upf',
                                'Tl': 'Tl_ONCV_PBE-1.0.upf',
                                'V': 'V_ONCV_PBE-1.0.upf',
                                'W': 'W_ONCV_PBE-1.0.upf',
                                'Xe': 'Xe_ONCV_PBE-1.0.upf',
                                'Y': 'Y_ONCV_PBE-1.0.upf',
                                'Zn': 'Zn_ONCV_PBE-1.0.upf',
                                'Zr': 'Zr_ONCV_PBE-1.0.upf'}}
