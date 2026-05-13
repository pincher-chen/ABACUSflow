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
    import numpy as np
except ImportError:
    np = None

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

def _crystal_axis_to_cartesian_magmom(m_crys, lattice_matrix):
    """
    把 magCIF 中 ``_atom_site_moment_crystalaxis_x/y/z`` 三分量反变换回笛卡尔坐标系。

    pymatgen 的写入约定（见 ``Magmom.get_moment_relative_to_crystal_axes``）::

        unit_m = M / |row_i(M)|                # 每行除以该行的模长
        m_crys = m_cart @ inv(unit_m)

    因此逆变换为::

        m_cart = m_crys @ unit_m

    其中 ``M`` 的行向量为晶格向量 a, b, c（与 ASE/pymatgen 一致）。

    Args:
        m_crys: 长度 3 的可迭代对象，单位 μB
        lattice_matrix: 3x3 数组，行向量为 a, b, c。单位无所谓（norm 会被约掉）。

    Returns:
        长度 3 的 numpy 数组 (mx, my, mz)，单位 μB
    """
    if np is None:
        raise RuntimeError("numpy is required for crystalaxis→Cartesian transform")
    M = np.asarray(lattice_matrix, dtype=float)
    norms = np.linalg.norm(M, axis=1)
    if not np.all(norms > 0):
        raise ValueError(f"Degenerate lattice rows: {norms}")
    unit_m = M / norms[:, None]
    m = np.asarray(m_crys, dtype=float).reshape(3)
    return m @ unit_m


def read_magmoms_from_cif(cif_file, lattice=None):
    """
    从 CIF/MCIF 文件读取磁矩信息。

    Args:
        cif_file: CIF/MCIF 文件路径
        lattice: 可选；3x3 数组，行向量为 a, b, c（如 ``ase.Atoms.get_cell()`` 的返回）。
            **强烈建议传入**，否则会直接返回 magCIF 中 ``_atom_site_moment_crystalaxis_*``
            字面值，对非正交晶格会导致坐标系错配。
            传入后：内部把每个原子的磁矩从晶轴系反变换到笛卡尔系，
            再按 Cartesian 大小判断共线性。

    Returns:
        list: 磁矩列表
            - 共线（Cartesian 空间 |mx|,|my| 均小于阈值）: ``[mz1, mz2, ...]``（标量列表）
            - 非共线: ``[[mx1,my1,mz1], [mx2,my2,mz2], ...]``
            - 文件中无磁矩 loop: ``None``

        当 ``lattice is None`` 时返回值的语义是晶轴系（crystalaxis）三分量；
        当 ``lattice`` 传入时返回值的语义是笛卡尔系。
    """
    def _cif_float(s):
        """CIF 数值字段可能带不确定度括号，如 '3.96(5)'，需要剥掉再转 float。"""
        s = s.strip()
        if '(' in s:
            s = s.split('(', 1)[0]
        return float(s)

    def _parse_loop(lines, start_idx):
        """
        解析 CIF 的 loop_:
        - 返回 (headers, rows, next_idx)
        """
        headers = []
        rows = []
        i = start_idx
        # headers
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue
            if s.startswith('_'):
                headers.append(s.split()[0])
                i += 1
                continue
            break
        # data rows
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith('#') or s.startswith('loop_') or s.startswith('_'):
                break
            rows.append(s.split())
            i += 1
        return headers, rows, i

    try:
        with open(cif_file, 'r') as f:
            lines = f.read().splitlines()

        # 优先解析 abacusflow 自定义注释 `# scalar_magmom <label> <value>`
        # 这一段由上游 write_mcif_with_magmom 写入，保留 VASP ISPIN=2 标量符号。
        # 在 pymatgen magCIF round-trip 中，crystalaxis 三分量恢复出的 Cartesian
        # m_z 与原始量化轴方向不一定一致 (受晶格 Cartesian 标准化旋转影响)，
        # 因此 sign(m_z) 不能可靠还原原标量符号；直接读 scalar 是最稳的做法。
        scalar_map = {}
        for ln in lines:
            s = ln.strip()
            if s.startswith('#') and 'scalar_magmom' in s:
                parts = s.lstrip('#').split()
                if len(parts) >= 3 and parts[0] == 'scalar_magmom':
                    try:
                        scalar_map[parts[1]] = float(parts[2])
                    except ValueError:
                        continue

        # 先解析 atom_site loop（为了得到 atom label 的顺序）
        atom_labels = []
        i = 0
        while i < len(lines):
            s = lines[i].strip()
            if s == 'loop_':
                headers, rows, next_i = _parse_loop(lines, i + 1)
                # 兼容 CIF1（_atom_site_label）和 CIF2/MCIF（_atom_site.label）
                _atom_label_key = next(
                    (h for h in headers if h.replace('.', '_') == '_atom_site_label'), None
                )
                if _atom_label_key:
                    label_col = headers.index(_atom_label_key)
                    for r in rows:
                        if label_col < len(r):
                            atom_labels.append(r[label_col])
                i = next_i
            else:
                i += 1

        # 再解析 moment loop：label -> (mx, my, mz)
        moment_map = {}
        i = 0
        while i < len(lines):
            s = lines[i].strip()
            if s == 'loop_':
                headers, rows, next_i = _parse_loop(lines, i + 1)
                # 兼容 CIF1（下划线）和 CIF2/MCIF（点号）两种命名风格
                _mh = {h.replace('.', '_') for h in headers}
                _hmap = {h.replace('.', '_'): h for h in headers}
                if ('_atom_site_moment_label' in _mh and
                    '_atom_site_moment_crystalaxis_x' in _mh and
                    '_atom_site_moment_crystalaxis_y' in _mh and
                    '_atom_site_moment_crystalaxis_z' in _mh):
                    c_label = headers.index(_hmap['_atom_site_moment_label'])
                    c_x = headers.index(_hmap['_atom_site_moment_crystalaxis_x'])
                    c_y = headers.index(_hmap['_atom_site_moment_crystalaxis_y'])
                    c_z = headers.index(_hmap['_atom_site_moment_crystalaxis_z'])
                    for r in rows:
                        try:
                            lab = r[c_label]
                            mx = _cif_float(r[c_x])
                            my = _cif_float(r[c_y])
                            mz = _cif_float(r[c_z])
                            moment_map[lab] = (mx, my, mz)
                        except Exception as _e:
                            print(f"[WARNING] Failed to parse magnetic moment row {r}: {_e}")
                            continue
                i = next_i
            else:
                i += 1

        if not moment_map:
            return None

        # 若上游写了 scalar_magmom 提示且 label 能对得上，直接返回标量列表，
        # 完全保留原始 VASP ISPIN=2 的带符号磁矩，跳过 crystalaxis 反变换。
        # 这是 spin=2 的最准确路径；对 spin=4 该数据无方向信息，仍会回退到 crystalaxis。
        if scalar_map and atom_labels and all(lab in scalar_map for lab in atom_labels):
            scalar_list = [float(scalar_map[lab]) for lab in atom_labels]
            print(f"[INFO] Using scalar_magmom hints from CIF "
                  f"(N={len(scalar_list)}, preserves original ISPIN=2 sign)")
            return scalar_list

        # 将 moment_map 按 atom_labels 顺序展开；如果拿不到 atom_labels，则按 moment_map 的出现顺序退化
        magmoms_vec = []
        if atom_labels:
            for lab in atom_labels:
                if lab in moment_map:
                    magmoms_vec.append(list(moment_map[lab]))
                else:
                    magmoms_vec.append([0.0, 0.0, 0.0])
        else:
            magmoms_vec = [list(v) for v in moment_map.values()]

        # 若提供了 lattice，则把 crystalaxis 三分量反变换到 Cartesian 后再判共线。
        # 这是必要的：对于非正交晶格（三斜/单斜等），一个纯笛卡尔 z 方向磁矩
        # 在晶轴系下三个分量都非零，旧逻辑会把"共线沿 z"误判为"非共线"。
        if lattice is not None:
            try:
                magmoms_vec = [
                    list(_crystal_axis_to_cartesian_magmom(v, lattice))
                    for v in magmoms_vec
                ]
            except Exception as _e:
                print(f"[WARNING] crystalaxis→Cartesian transform failed, "
                      f"will fall back to literal crystalaxis values: {_e}")

        # 判断是否共线：所有原子的 |mx|, |my| < 阈值时视为沿 z 轴的共线磁矩。
        # 1e-4 μB 对应数值噪声/写入精度，足够区分真实的非共线分量。
        thresh = 1e-4
        is_collinear = all(abs(m[0]) < thresh and abs(m[1]) < thresh for m in magmoms_vec)
        if is_collinear:
            return [m[2] for m in magmoms_vec]
        return magmoms_vec

    except Exception as e:
        print(f"[WARNING] Failed to read magnetic moments from CIF: {e}")
        return None


def read_ase_stru(stru=None, coordinates_type="Cartesian", cif_magmoms=None):
    """
    读取ASE结构并提取原子信息
    
    Args:
        stru: ASE Atoms 对象
        coordinates_type: 坐标类型 ('Cartesian' 或 'Direct')
        cif_magmoms: 从 CIF 读取的磁矩列表（如果有）
    
    Returns:
        atoms_list, atoms_masses, atoms_position, atoms_magnetism, atoms_indices
    """
    if not judge_exist_stru(stru):
        return None
    
    atoms_list = []
    atoms_position = []
    atoms_masses = []
    atoms_magnetism = []
    atoms_indices = []
    atoms_all = stru.get_chemical_symbols()
    
    # 按元素排序
    for atoms_all_name in atoms_all:
        if atoms_all_name not in atoms_list:
            atoms_list.append(atoms_all_name)
    
    for atoms_list_name in atoms_list:
        atoms_position.append([])
        atoms_masses.append([])
        atoms_magnetism.append([])
        atoms_indices.append([])
    
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
                atoms_indices[i].append(j)
                
                # 如果提供了 CIF 磁矩，使用它
                if cif_magmoms is not None and j < len(cif_magmoms):
                    atoms_magnetism[i] = cif_magmoms[j]
                else:
                    atoms_magnetism[i] = stru.get_initial_magnetic_moments()[j]
    
    return atoms_list, atoms_masses, atoms_position, atoms_magnetism, atoms_indices

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
                     spin=1, copy_files=True, need_fr=False, cif_file=None, guess_mag=False, **kwargs):
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
        spin: 自旋设置 (1=无磁性, 2=共线磁性, 4=非共线磁性)
        copy_files: 是否复制赝势/轨道文件
        need_fr: 是否需要 FR 赝势（SOC 计算）
        cif_file: CIF 文件路径（用于读取磁矩信息）
        guess_mag: 是否猜测初始磁矩（当CIF中没有磁矩时）
    """
    if not judge_exist_stru(stru):
        return "No input structure!"
    
    # 尝试从 CIF 读取磁矩
    # 关键：把 ASE 晶格（行向量 a, b, c, 单位 Å）传给 read_magmoms_from_cif，
    # 以便把 magCIF 中 _atom_site_moment_crystalaxis_* 反变换到 Cartesian。
    # 这样下面写入 STRU 时 magmom 总是 Cartesian 系下的 μB，与 ABACUS 约定一致。
    cif_magmoms = None
    if cif_file and exists(cif_file):
        lat_for_cif = None
        try:
            if np is not None:
                lat_for_cif = np.array(stru.get_cell(), dtype=float)
        except Exception as _e:
            print(f"[WARNING] failed to obtain lattice from ASE Atoms: {_e}")
        cif_magmoms = read_magmoms_from_cif(cif_file, lattice=lat_for_cif)
        if cif_magmoms:
            print(f"[INFO] Read {len(cif_magmoms)} magnetic moments from CIF file (Cartesian, μB)")
            # 判断是共线还是非共线（基于 Cartesian 分量）
            if isinstance(cif_magmoms[0], list):
                # 计算最大面内分量，给用户一个量级感
                max_xy = 0.0
                for m in cif_magmoms:
                    max_xy = max(max_xy, abs(m[0]), abs(m[1]))
                print(f"[INFO] Non-collinear magnetism detected from CIF "
                      f"(max in-plane |m| = {max_xy:.4f} μB)")
                if spin == 2:
                    # spin=2 时写入 sign(m_z)*|m|（带符号模长），
                    # 因此 |m| 与原始数据严格守恒；只是方向被强制沿 z。
                    print(f"[INFO] spin=2 requested: writing sign(m_z)*|m| per atom "
                          f"(|m| preserved, direction projected onto z). "
                          f"Use spin=4 if mx/my orientation is needed.")
            else:
                print(f"[INFO] Collinear magnetism detected from CIF (along z)")
    
    atoms_list, atoms_masses, atoms_position, atoms_magnetism, atoms_indices = read_ase_stru(
        stru, coordinates_type, cif_magmoms
    )
    
    potential = set_potential(atoms_list, pseudo_dir, potential_name, need_fr=need_fr)
    basis = set_basis(atoms_list, basis_dir, basis_name)
    
    # 处理磁矩设置
    # 优先级: CIF 磁矩 > guess_mag > 默认值
    if cif_magmoms is None:
        # CIF 中没有磁矩信息
        if spin >= 2 and guess_mag:
            # 使用猜测的磁矩
            print(f"[INFO] Using guess_mag: setting initial magnetic moment to 2.0")
            for i in range(len(atoms_list)):
                if spin == 2:
                    atoms_magnetism[i] = 2.0  # 共线磁性
                elif spin == 4:
                    atoms_magnetism[i] = [0.0, 0.0, 2.0]  # 非共线磁性，沿z轴
        elif spin >= 2:
            # 不使用 guess_mag，但 spin > 1，给一个默认小磁矩
            for i in range(len(atoms_list)):
                if spin == 2:
                    atoms_magnetism[i] = 0.0
                elif spin == 4:
                    atoms_magnetism[i] = [0.0, 0.0, 0.0]
    
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
        
        # LATTICE_CONSTANT: 转换因子从 Ångström 到 Bohr
        # ASE 的 get_cell() 返回 Ångström，ABACUS 期望 Bohr
        # 1 Ångström = 1.889726 Bohr
        f.write('\nLATTICE_CONSTANT\n1.889726\n\nLATTICE_VECTORS\n')
        for i in range(3):
            for j in range(3):
                f.write(f"{stru.get_cell()[i][j]:0<12f}   ")
            f.write('\n')
        
        f.write(f'\nATOMIC_POSITIONS\n{coordinates_type}\n\n')
        for i in range(len(atoms_list)):
            if spin in (2, 4):
                species_mag = 0.0
            else:
                species_mag = atoms_magnetism[i] if isinstance(atoms_magnetism[i], (int, float)) else 0.0
            f.write(f"{atoms_list[i]}\n{species_mag:0<12f}\n{len(atoms_position[i])}\n")
            for j in range(len(atoms_position[i])):
                f.write(f"{atoms_position[i][j][0]:0<12f}   ")
                f.write(f"{atoms_position[i][j][1]:0<12f}   ")
                f.write(f"{atoms_position[i][j][2]:0<12f}   ")
                f.write((str(fix) + '   ') * 3)
                
                # 添加磁矩信息（如果 spin >= 2）
                original_atom_index = atoms_indices[i][j] if atoms_indices and atoms_indices[i] else None

                # 注意: cif_magmoms 元素已是 Cartesian 系 (μB)
                #   - 共线 (沿 z): 元素是标量 float
                #   - 非共线:       元素是长度 3 的 list [mx, my, mz]
                # spin 由 INPUT 中的 nspin 决定: 1 不写 magmom; 2 写一个标量;
                # 4 写三个分量 (mx my mz) - 这与 ABACUS STRU 官方约定一致。
                if spin == 1:
                    # 非磁性: 不写 magmom，留给 INPUT 默认行为
                    pass
                elif spin == 2 and cif_magmoms is not None and original_atom_index is not None:
                    if original_atom_index < len(cif_magmoms):
                        v = cif_magmoms[original_atom_index]
                        if isinstance(v, list):
                            # spin=2 共线: 取带符号模长 sign(m_z) * |m|。
                            # 这样磁矩大小 |m| 与原始 ISPIN=2 标量严格一致，
                            # 不受 pymatgen (c 沿 z) vs ASE (a 沿 x) 笛卡尔标准化
                            # 差异影响——后者会让 m_z 投影偏小（如 |m|=0.907 但
                            # ASE 系下 m_z=-0.876）。
                            mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5
                            mag_val = mag if v[2] >= 0.0 else -mag
                        else:
                            mag_val = float(v)
                        f.write(f' magmom {mag_val:.6f}')
                elif spin == 4 and cif_magmoms is not None and original_atom_index is not None:
                    if original_atom_index < len(cif_magmoms):
                        v = cif_magmoms[original_atom_index]
                        if isinstance(v, list):
                            mx, my, mz = float(v[0]), float(v[1]), float(v[2])
                        else:
                            mx, my, mz = 0.0, 0.0, float(v)
                        f.write(f' magmom {mx:.6f} {my:.6f} {mz:.6f}')
                elif spin == 2 and guess_mag:
                    f.write(f' magmom 2.0')
                elif spin == 4 and guess_mag:
                    f.write(f' magmom 0.0 0.0 2.0')
                
                f.write('\n')
    
    return {'pseudo_dir': pseudo_dir, 'basis_dir': basis_dir,
            'potential_name': potential, 'basis_name': basis}


