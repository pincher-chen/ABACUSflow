#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STRU 文件生成工具模块
包含 ABACUS 结构文件的读取、生成和处理功能
"""

from os.path import join, basename, exists
import shutil
try:
    from ase import Atoms
except ImportError:
    Atoms = None

try:
    from abacus.potential import PotentialDict, is_dynamic_pseudo_dir, pick_upf
    from abacus.basis import BasisDict
except ImportError:
    PotentialDict = {}
    BasisDict = {}
    is_dynamic_pseudo_dir = lambda x: False
    pick_upf = None

# 常见元素的价电子数字典
# 注意：这里的价电子数是指 SG15 赝势中包含的电子数
# 对于某些元素（如 Mg, Ca, Sr, Ba），赝势包含了更多的电子层
VALENCE_ELECTRONS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8,
    'Na': 9, 'Mg': 10, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8,  # Na/Mg 修正为包含内层
    'K': 9, 'Ca': 10, 'Sc': 11, 'Ti': 12, 'V': 13, 'Cr': 14, 'Mn': 15, 'Fe': 16, 'Co': 17, 'Ni': 18, 'Cu': 19, 'Zn': 20,  # K/Ca 修正
    'Ga': 13, 'Ge': 14, 'As': 15, 'Se': 16, 'Br': 17, 'Kr': 18,
    'Rb': 9, 'Sr': 10, 'Y': 11, 'Zr': 12, 'Nb': 13, 'Mo': 14, 'Tc': 15, 'Ru': 16, 'Rh': 17, 'Pd': 18, 'Ag': 19, 'Cd': 20,  # Rb/Sr 修正
    'In': 13, 'Sn': 14, 'Sb': 15, 'Te': 16, 'I': 17, 'Xe': 18,
    'Cs': 9, 'Ba': 10, 'La': 11, 'Hf': 12, 'Ta': 13, 'W': 14, 'Re': 15, 'Os': 16, 'Ir': 17, 'Pt': 18, 'Au': 19, 'Hg': 20,  # Cs/Ba 修正
    'Tl': 13, 'Pb': 14, 'Bi': 15, 'Po': 16, 'At': 17, 'Rn': 18
}

def potential_list():
    """返回可用的赝势列表"""
    return list(PotentialDict.keys())

def basis_list():
    """返回可用的轨道基组列表"""
    return list(BasisDict.keys())

def get_supported_elements(potential_name='PotSG15', basis_name='SG15std'):
    """
    获取指定赝势和基组支持的元素列表
    
    Args:
        potential_name: 赝势名称（默认 'PotSG15'）
        basis_name: 基组名称（默认 'SG15std'）
    
    Returns:
        set: 两者都支持的元素集合
    """
    supported_pot = set()
    supported_basis = set()
    
    if potential_name in PotentialDict:
        supported_pot = set(PotentialDict[potential_name].keys())
    
    if basis_name in BasisDict:
        supported_basis = set(BasisDict[basis_name].keys())
    
    # 返回两者的交集（LCAO 需要同时有赝势和基组）
    if supported_pot and supported_basis:
        return supported_pot & supported_basis
    elif supported_pot:
        # PW 只需要赝势
        return supported_pot
    else:
        return set()

def check_elements_supported(atoms, potential_name='PotSG15', basis_name='SG15std'):
    """
    检查结构中的所有元素是否被支持
    
    Args:
        atoms: ASE Atoms 对象
        potential_name: 赝势名称
        basis_name: 基组名称
    
    Returns:
        (bool, list): (是否全部支持, 不支持的元素列表)
    """
    if atoms is None:
        return False, []
    
    # 获取结构中的所有元素
    elements = set(atoms.get_chemical_symbols())
    
    # 获取支持的元素
    supported = get_supported_elements(potential_name, basis_name)
    
    # 找出不支持的元素
    unsupported = elements - supported
    
    return (len(unsupported) == 0, sorted(list(unsupported)))

def judge_exist_stru(stru=None):
    """判断结构对象是否存在"""
    return stru is not None

def get_total_valence_electrons(atoms):
    """计算体系总价电子数"""
    total_ne = 0
    for symbol in atoms.get_chemical_symbols():
        ne = VALENCE_ELECTRONS.get(symbol)
        if ne is None:
            print(f"[WARNING] Valence electrons for {symbol} not found, using default 8")
            ne = 8
        total_ne += ne
    return total_ne

def read_ase_stru(stru=None, coordinates_type="Cartesian"):
    """读取ASE结构并提取原子信息"""
    if not judge_exist_stru(stru):
        return None
    
    atoms_list = []
    atoms_position = []
    atoms_masses = []
    atoms_magnetism = []
    atoms_all = stru.get_chemical_symbols()
    
    # 按元素排序
    for atoms_all_name in atoms_all:
        if atoms_all_name not in atoms_list:
            atoms_list.append(atoms_all_name)
    
    for atoms_list_name in atoms_list:
        atoms_position.append([])
        atoms_masses.append([])
        atoms_magnetism.append([])
    
    # 获取位置、质量、磁矩
    for i in range(len(atoms_list)):
        for j in range(len(atoms_all)):
            if atoms_all[j] == atoms_list[i]:
                if coordinates_type == 'Cartesian':
                    atoms_position[i].append(list(stru.get_positions()[j]))
                elif coordinates_type == 'Direct':
                    atoms_position[i].append(list(stru.get_scaled_positions()[j]))
                else:
                    raise ValueError("coordinates_type must be 'Cartesian' or 'Direct'")
                atoms_masses[i] = stru.get_masses()[j]
                atoms_magnetism[i] = stru.get_initial_magnetic_moments()[j]
    
    return atoms_list, atoms_masses, atoms_position, atoms_magnetism

def set_potential(atoms_list=None, pseudo_dir="./", potential_name=None, need_fr=False):
    """
    设置赝势文件
    
    Args:
        atoms_list: 元素符号列表
        pseudo_dir: 赝势目录路径
        potential_name: 赝势名称（字典键名）或文件列表
        need_fr: 是否需要 FR（全相对论）赝势用于 SOC 计算
    
    Returns:
        赝势文件路径列表
    """
    if atoms_list is None:
        print("[ERROR] Please set atoms_list")
        return None
    
    potential = []
    
    # 检查是否是动态赝势目录（包含版本化的 ONCV 赝势）
    if is_dynamic_pseudo_dir(pseudo_dir) and pick_upf is not None:
        # 使用动态选择逻辑
        print(f"[INFO] Using dynamic pseudopotential selection from {pseudo_dir}")
        print(f"[INFO] SOC mode (FR required): {need_fr}")
        
        for elem in atoms_list:
            try:
                upf_name = pick_upf(pseudo_dir, elem, need_fr)
                potential.append(join(pseudo_dir, upf_name))
                print(f"[INFO] {elem}: {upf_name}")
            except RuntimeError as e:
                print(f"[ERROR] {e}")
                raise
        
        return potential
    
    # 原有的静态字典逻辑（向后兼容）
    if potential_name is None:
        for atoms_list_name in atoms_list:
            potential.append(join(pseudo_dir, PotentialDict['PotLDA'][atoms_list_name]))
    elif type(potential_name) == str:
        if potential_name in potential_list():
            for atoms_list_name in atoms_list:
                potential.append(join(pseudo_dir, PotentialDict[potential_name][atoms_list_name]))
        else:
            raise ValueError(f"potential_name '{potential_name}' not found")
    elif type(potential_name) == list:
        ele_name = {}
        for i in potential_name:
            with open(join(pseudo_dir, i), 'r') as f:
                lines = f.readlines()
            for line in lines:
                line = line.replace('=', ' = ').replace('"', ' " ')
                data = line.split()
                if len(data) >= 4 and data[0] == 'element':
                    ele_name[data[3]] = i
                    break
                elif len(data) == 2 and data[1] == 'Element':
                    ele_name[data[0]] = i
                    break
        for atoms_list_name in atoms_list:
            potential.append(join(pseudo_dir, ele_name[atoms_list_name]))
    else:
        raise ValueError("Invalid potential_name type")
    
    return potential

def set_basis(atoms_list=None, basis_dir="./", basis_name=None):
    """设置轨道文件"""
    if atoms_list is None:
        print("[ERROR] Please set atoms_list")
        return None
    
    basis = []
    if basis_name is None:
        for atoms_list_name in atoms_list:
            basis.append(join(basis_dir, BasisDict['LDAmin'][atoms_list_name]))
    elif type(basis_name) == str:
        if basis_name in basis_list():
            for atoms_list_name in atoms_list:
                basis.append(join(basis_dir, BasisDict[basis_name][atoms_list_name]))
        else:
            raise ValueError(f"basis_name '{basis_name}' not found")
    elif type(basis_name) == list:
        ele_name = {}
        for i in basis_name:
            with open(join(basis_dir, i), 'r') as f:
                lines = f.readlines()
            for line in lines:
                data = line.split()
                if len(data) >= 2 and data[0] == 'Element':
                    ele_name[data[1]] = i
                    break
        for atoms_list_name in atoms_list:
            basis.append(join(basis_dir, ele_name[atoms_list_name]))
    else:
        raise ValueError("Invalid basis_name type")
    
    return basis

def write_input_stru(stru=None, pseudo_dir='./', potential_name=None,
                     basis_dir='./', basis_name=None, fix=1,
                     filename='STRU', directory='./', coordinates_type='Cartesian',
                     spin=1, copy_files=True, need_fr=False, **kwargs):
    """
    生成 STRU 文件
    
    Args:
        stru: ASE Atoms 对象
        pseudo_dir: 赝势目录
        potential_name: 赝势名称
        basis_dir: 轨道基组目录
        basis_name: 轨道基组名称
        fix: 固定原子标记
        filename: 输出文件名
        directory: 输出目录
        coordinates_type: 坐标类型 ('Cartesian' 或 'Direct')
        spin: 自旋设置
        copy_files: 是否复制赝势/轨道文件
        need_fr: 是否需要 FR 赝势（SOC 计算）
    """
    if not judge_exist_stru(stru):
        return "No input structure!"
    
    atoms_list, atoms_masses, atoms_position, atoms_magnetism = read_ase_stru(stru, coordinates_type)
    
    potential = set_potential(atoms_list, pseudo_dir, potential_name, need_fr=need_fr)
    basis = set_basis(atoms_list, basis_dir, basis_name)
    
    # 如果 spin=2，设置初始磁矩
    if spin == 2:
        for i in range(len(atoms_list)):
            atoms_magnetism[i] = 1.0
    
    # 写入 STRU 文件
    with open(join(directory, filename), 'w') as f:
        f.write('ATOMIC_SPECIES\n')
        for i in range(len(atoms_list)):
            if copy_files:
                if not exists(join(directory, basename(potential[i]))):
                    shutil.copyfile(potential[i], join(directory, basename(potential[i])))
                pot_file = basename(potential[i])
            else:
                pot_file = basename(potential[i])
            
            temp1 = ' ' * (4 - len(atoms_list[i]))
            temp2 = ' ' * (14 - len(str(atoms_masses[i])))
            atomic_species = atoms_list[i] + temp1 + str(atoms_masses[i]) + temp2 + pot_file
            f.write(atomic_species + '\n')
        
        f.write('\nNUMERICAL_ORBITAL\n')
        for i in range(len(atoms_list)):
            if copy_files:
                if not exists(join(directory, basename(basis[i]))):
                    shutil.copyfile(basis[i], join(directory, basename(basis[i])))
                basis_file = basename(basis[i])
            else:
                basis_file = basename(basis[i])
            f.write(basis_file + '\n')
        
        f.write('\nLATTICE_CONSTANT\n1.0\n\nLATTICE_VECTORS\n')
        for i in range(3):
            for j in range(3):
                f.write(f"{stru.get_cell()[i][j]:0<12f}   ")
            f.write('\n')
        
        f.write(f'\nATOMIC_POSITIONS\n{coordinates_type}\n\n')
        for i in range(len(atoms_list)):
            f.write(f"{atoms_list[i]}\n{atoms_magnetism[i]:0<12f}\n{len(atoms_position[i])}\n")
            for j in range(len(atoms_position[i])):
                f.write(f"{atoms_position[i][j][0]:0<12f}   ")
                f.write(f"{atoms_position[i][j][1]:0<12f}   ")
                f.write(f"{atoms_position[i][j][2]:0<12f}   ")
                f.write((str(fix) + '   ') * 3 + '\n')
    
    return {'pseudo_dir': pseudo_dir, 'basis_dir': basis_dir,
            'potential_name': potential, 'basis_name': basis}


