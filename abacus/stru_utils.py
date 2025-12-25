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
    from abacus.potential import PotentialDict
    from abacus.basis import BasisDict
except ImportError:
    PotentialDict = {}
    BasisDict = {}

# 常见元素的价电子数字典
VALENCE_ELECTRONS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8,
    'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8,
    'K': 1, 'Ca': 2, 'Sc': 11, 'Ti': 12, 'V': 13, 'Cr': 14, 'Mn': 15, 'Fe': 16, 'Co': 17, 'Ni': 18, 'Cu': 19, 'Zn': 20,
    'Ga': 13, 'Ge': 14, 'As': 15, 'Se': 16, 'Br': 17, 'Kr': 18,
    'Rb': 1, 'Sr': 2, 'Y': 11, 'Zr': 12, 'Nb': 13, 'Mo': 14, 'Tc': 15, 'Ru': 16, 'Rh': 17, 'Pd': 18, 'Ag': 19, 'Cd': 20,
    'In': 13, 'Sn': 14, 'Sb': 15, 'Te': 16, 'I': 17, 'Xe': 18,
    'Cs': 1, 'Ba': 2, 'La': 11, 'Hf': 12, 'Ta': 13, 'W': 14, 'Re': 15, 'Os': 16, 'Ir': 17, 'Pt': 18, 'Au': 19, 'Hg': 20,
    'Tl': 13, 'Pb': 14, 'Bi': 15, 'Po': 16, 'At': 17, 'Rn': 18
}

def potential_list():
    """返回可用的赝势列表"""
    return list(PotentialDict.keys())

def basis_list():
    """返回可用的轨道基组列表"""
    return list(BasisDict.keys())

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

def set_potential(atoms_list=None, pseudo_dir="./", potential_name=None):
    """设置赝势文件"""
    if atoms_list is None:
        print("[ERROR] Please set atoms_list")
        return None
    
    potential = []
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
                     spin=1, copy_files=True, **kwargs):
    """生成 STRU 文件"""
    if not judge_exist_stru(stru):
        return "No input structure!"
    
    atoms_list, atoms_masses, atoms_position, atoms_magnetism = read_ase_stru(stru, coordinates_type)
    
    potential = set_potential(atoms_list, pseudo_dir, potential_name)
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

