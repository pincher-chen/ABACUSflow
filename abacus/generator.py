#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
输入文件生成器模块
包含 INPUT、KPT、STRU 文件生成的核心逻辑
"""

import os
import re
from pathlib import Path

try:
    from ase.io import read
except ImportError:
    read = None

# 导入本地模块
from abacus.stru_utils import (
    write_input_stru, 
    get_total_valence_electrons,
    check_elements_supported,
    get_supported_elements
)


def get_nlocal_from_scf(work_dir):
    """
    从 Scf 阶段的 running.log 中读取 NLOCAL 值
    
    Args:
        work_dir: 当前工作目录 (Path 对象)
    
    Returns:
        NLOCAL 值 (int)，如果找不到返回 None
    """
    scf_dir = work_dir.parent / 'Scf'
    if not scf_dir.exists():
        return None
    
    # 查找 running.log
    log_files = list(scf_dir.glob('OUT.*/running*.log'))
    if not log_files:
        log_files = list(scf_dir.glob('running*.log'))
    
    if not log_files:
        return None
    
    # 读取最新的 log 文件
    log_file = max(log_files, key=lambda p: p.stat().st_mtime)
    
    try:
        content = log_file.read_text()
        # 搜索 NLOCAL = xx 的行
        match = re.search(r'NLOCAL\s*=\s*(\d+)', content)
        if match:
            return int(match.group(1))
    except Exception:
        pass
    
    return None

def generate_input_files(work_dir, stage, stru_file, config_get, incar_template, workflow, click_echo=print):
    """
    生成 INPUT、KPT、STRU 文件的核心逻辑
    
    参数:
        work_dir: 工作目录（Path 对象）
        stage: 计算阶段名称
        stru_file: 结构文件路径        config_get: 配置读取函数
        incar_template: INPUT 参数模板字典
        workflow: 工作流配置字典
        click_echo: 输出函数（默认 print，可传入 click.echo）
    """
    # 检查 ASE 是否可用
    if read is None:
        raise RuntimeError("ASE not available. Please run in conda environment: conda activate dftflow")
    
    # 延迟导入 AbacusInput
    try:
        from abacus.create_input import AbacusInput
    except ImportError:
        raise RuntimeError("AbacusInput not available. Please check abacus/create_input.py")
    
    # 读取配置
    POTPATH = config_get('ABACUS', 'PSEUDO_POTENTIAL_DIR', 
                  '/XYFS01/nscc-gz_pinchen_1/sf_box/abacus-develop-LTSv3.10.0/PotSG15')
    ORBPATH = config_get('ABACUS', 'ORBITAL_DIR',
                  '/XYFS01/nscc-gz_pinchen_1/sf_box/abacus-develop-LTSv3.10.0/OrbSG15std')
    KSPACING = float(config_get('PARAMETERS', 'KSPACING', '0.13'))
    
    # 切换到工作目录
    os.chdir(work_dir)
    
    # 确定使用哪个结构文件
    # Relax 使用 Coarse_relax 的优化结构（STRU_ION_D）
    # Scf 使用 Relax 的优化结构（STRU_ION_D）
    # Band/Dos 使用 Scf 的 STRU
    use_optimized_stru = False
    optimized_stru_path = None
    
    if stage.lower() == 'relax':
        # 查找 Coarse_relax 的优化结构
        coarse_out_dir = work_dir.parent / 'Coarse_relax' / 'OUT.Coarse_relax'
        stru_ion_d = coarse_out_dir / 'STRU_ION_D'
        if stru_ion_d.exists():
            use_optimized_stru = True
            optimized_stru_path = stru_ion_d
            click_echo(f"[INFO] Will use Coarse_relax optimized structure: {stru_ion_d.name}")
        else:
            click_echo(f"[WARNING] Coarse_relax STRU_ION_D not found, using original: {stru_file}")
    
    elif stage.lower() == 'scf':
        # 查找 Relax 的优化结构
        relax_out_dir = work_dir.parent / 'Relax' / 'OUT.Relax'
        stru_ion_d = relax_out_dir / 'STRU_ION_D'
        if stru_ion_d.exists():
            use_optimized_stru = True
            optimized_stru_path = stru_ion_d
            click_echo(f"[INFO] Will use Relax optimized structure: {relax_out_dir.name}/STRU_ION_D")
        else:
            click_echo(f"[WARNING] Relax STRU_ION_D not found, using original: {stru_file}")
    
    elif stage.lower() in ['band', 'dos']:
        # 查找 Scf 的 STRU 文件
        scf_stru = work_dir.parent / 'Scf' / 'STRU'
        if scf_stru.exists():
            use_optimized_stru = True
            optimized_stru_path = scf_stru
            click_echo(f"[INFO] Will use Scf STRU file: {scf_stru}")
        else:
            click_echo(f"[WARNING] Scf STRU not found, using original: {stru_file}")
    
    # 读取结构文件（用于获取元素类型等信息）
    # 根据文件扩展名选择正确的格式
    stru_file_lower = str(stru_file).lower()
    if stru_file_lower.endswith('.cif'):
        stru = read(stru_file, format='cif')
    elif stru_file_lower.endswith('.xyz'):
        stru = read(stru_file, format='xyz')
    elif stru_file_lower.endswith(('.vasp', '.poscar', '.contcar')):
        stru = read(stru_file, format='vasp')
    else:
        # 尝试自动检测格式
        stru = read(stru_file)
    ntype = len(set(stru.get_chemical_symbols()))
    
    # ========================================
    # 关键检查：验证所有元素是否被支持
    # ========================================
    potential_name = 'PotSG15'  # 默认使用 SG15 赝势
    basis_name = 'SG15std'      # 默认使用 SG15 标准基组
    
    is_supported, unsupported_elements = check_elements_supported(
        stru, potential_name, basis_name
    )
    
    if not is_supported:
        supported_elements = get_supported_elements(potential_name, basis_name)
        
        click_echo("")
        click_echo("=" * 70)
        click_echo("❌ ERROR: Unsupported Elements Detected")
        click_echo("=" * 70)
        click_echo(f"Structure contains unsupported elements: {', '.join(unsupported_elements)}")
        click_echo("")
        click_echo(f"Currently using:")
        click_echo(f"  - Pseudopotential: {potential_name} (69 elements)")
        click_echo(f"  - Orbital basis:   {basis_name} (68 elements)")
        click_echo("")
        click_echo(f"Supported elements:")
        click_echo(f"  {', '.join(sorted(list(supported_elements)))}")
        click_echo("")
        click_echo("=" * 70)
        click_echo("🔧 Solutions:")
        click_echo("=" * 70)
        click_echo("1. Remove unsupported elements from your structure")
        click_echo("2. Replace with similar supported elements")
        click_echo("3. Use a different pseudopotential library (if available)")
        click_echo("")
        click_echo("📚 References:")
        click_echo("  - SG15 pseudopotentials: http://www.quantum-simulation.org/potentials/sg15_oncv/")
        click_echo("  - ABACUS documentation: http://abacus.ustc.edu.cn/")
        click_echo("=" * 70)
        
        raise ValueError(
            f"Unsupported elements: {', '.join(unsupported_elements)}. "
            f"Only {len(supported_elements)} elements are supported by {potential_name}."
        )
    
    click_echo(f"[INFO] All elements are supported: {', '.join(sorted(set(stru.get_chemical_symbols())))}")
    
    # 处理需要电荷密度的阶段（Band/Dos）
    if stage.lower() in ['band', 'dos']:
        # 查找 Scf 的电荷密度文件
        scf_out_dir = work_dir.parent / 'Scf' / 'OUT.Scf'
        if scf_out_dir.exists():
            # 创建 OUT.{stage} 目录
            out_dir = work_dir / f'OUT.{stage}'
            out_dir.mkdir(exist_ok=True)
            
            # 查找并链接电荷密度文件
            chg_files = list(scf_out_dir.glob('SPIN*_CHG.cube')) + list(scf_out_dir.glob('SPIN*_CHG'))
            
            if chg_files:
                for chg_file in chg_files:
                    link_name = out_dir / chg_file.name
                    if link_name.exists() or link_name.is_symlink():
                        link_name.unlink()
                    
                    # 创建相对路径的符号链接
                    rel_target = os.path.relpath(chg_file, out_dir)
                    os.symlink(rel_target, link_name)
                    click_echo(f"[INFO] Linked charge density: {chg_file.name}")
            else:
                click_echo(f"[WARNING] No charge density files found in {scf_out_dir}")
        else:
            click_echo(f"[WARNING] Scf output directory not found: {scf_out_dir}")
    
    # 生成 INPUT 文件
    input_obj = AbacusInput()
    
    # 获取模板参数
    stage_params = incar_template.get(stage, {})
    
    # 基础设置
    input_obj.set(ntype=ntype,
                  pseudo_dir=POTPATH,
                  orbital_dir=ORBPATH,
                  suffix=stage,
                  gamma_only=0)
    
    # 检查自旋设置
    # 优先使用模板中的 nspin 值
    template_nspin = stage_params.get('nspin', None)
    
    # 如果模板中没有设置，则根据 SPIN_ON/SPIN_OFF 文件决定
    default_nspin = 2
    if stage.lower() != 'test_spin':
        spin_on_file = work_dir.parent / 'SPIN_ON'
        spin_off_file = work_dir.parent / 'SPIN_OFF'
        
        if spin_off_file.exists():
            default_nspin = 1
        elif spin_on_file.exists():
            default_nspin = 2
    
    # 如果模板中有 nspin 设置，使用模板的值（用户明确指定优先）
    if template_nspin is not None:
        default_nspin = template_nspin
        click_echo(f"[INFO] Using nspin={template_nspin} from template")
    
    # 计算 nbands
    # 注意：对于 Test_spin, Coarse_relax, Relax, Scf 不设置 nbands
    # 让 ABACUS 使用默认值更保险，只对 Band/Dos 手动设置
    if stage.lower() in ['band', 'dos']:
        # Band/Dos 需要更多能带以获得完整的能带结构/态密度
        total_ne = get_total_valence_electrons(stru)
        
        # 尝试从 Scf 阶段读取 NLOCAL（基函数总数）
        nlocal = get_nlocal_from_scf(work_dir)
        
        # 1. 计算占据能带数 (Occupied bands)
        # 重要说明：nspin=2 时，NBANDS 是每个自旋通道的能带数，不需要翻倍！
        # - nspin=1: 每个能带容纳 2 个电子 → occ_bands = total_ne / 2
        # - nspin=2: 每个能带容纳 1 个电子 → occ_bands = total_ne
        if default_nspin == 1:
            occ_bands = total_ne / 2
        else:
            occ_bands = total_ne
        
        # 2. 增加冗余量 (Buffer) 用于观察空带
        # 对于 LCAO 体系，冗余量不宜过大；+10 到 +20 足够
        buffer = 20
        nbands = int(occ_bands + buffer)
        
        # 3. 确保 nbands 至少是某个最小值（对于很小的体系）
        nbands = max(nbands, 20)
        
        # 4. 【关键】LCAO 的硬性约束：nbands 必须 <= NLOCAL
        # NLOCAL 是基函数总数，这是 LCAO 的物理上限
        if nlocal is not None:
            if nbands > nlocal:
                click_echo(f"[WARNING] Calculated nbands={nbands} exceeds NLOCAL={nlocal} (LCAO limit)")
                click_echo(f"[INFO] Adjusting nbands to {nlocal} (maximum for LCAO basis)")
                nbands = nlocal
            else:
                click_echo(f"[INFO] nbands={nbands} is valid (NLOCAL={nlocal}, occupied={int(occ_bands)})")
        else:
            click_echo(f"[WARNING] Could not read NLOCAL from Scf/OUT.*/running*.log")
            click_echo(f"[WARNING] Using nbands={nbands} (may fail if > NLOCAL for LCAO basis)")
            click_echo(f"[WARNING] If you see 'NLOCAL < NBANDS' error, manually set nbands in yaml")
        
        input_obj.set(nbands=nbands)
        click_echo(f"[INFO] Set nbands={nbands} for {stage} (electrons={total_ne}, nspin={default_nspin})")
    else:
        # Test_spin, Coarse_relax, Relax, Scf 不设置 nbands
        # 让 ABACUS 自动计算，更安全
        click_echo(f"[INFO] Using ABACUS default nbands for {stage}")
    
    # 应用模板参数
    if stage_params:
        input_obj.set(**stage_params)
    else:
        # 默认参数
        input_obj.set(calculation='scf',
                      nspin=default_nspin,
                      ecutwfc=100,
                      scf_thr=1e-7,
                      scf_nmax=100,
                      ks_solver='dav',
                      basis_type='lcao',
                      smearing_method='gaussian',
                      smearing_sigma=0.02,
                      mixing_type='broyden',
                      mixing_beta=0.4,
                      symmetry=1)
    
    # 对于 Band/Dos，设置 init_chg=1（从文件读取电荷密度）
    if stage.lower() in ['band', 'dos']:
        input_obj.set(init_chg=1)
        click_echo(f"[INFO] Set init_chg=1 for {stage} calculation")
    
    # 设置 nspin（只有在模板中没有明确指定时才使用 default_nspin）
    # 如果模板中已经有 nspin，上面 input_obj.set(**stage_params) 已经设置过了
    # 这里只处理模板中没有 nspin 的情况
    if template_nspin is None and stage.lower() != 'test_spin':
        input_obj.set(nspin=default_nspin)
        click_echo(f"[INFO] Using nspin={default_nspin} (auto-detected from SPIN_ON/SPIN_OFF or default)")
    
    # Band 计算设置 symmetry=0
    stage_workflow_config = workflow.get(stage, {})
    ktype_check = stage_workflow_config.get('ktype', 'Gamma')
    if ktype_check in ['Line', 'L'] or stage.lower() == 'band':
        input_obj.set(symmetry=0)
    
    input_obj.write_input_input()
    
    # 生成 KPT 文件
    if 'kspacing' not in stage_params:
        kspacing_list = stage_workflow_config.get('kval', None)
        ktype = stage_workflow_config.get('ktype', 'Gamma')
        
        if kspacing_list and len(kspacing_list) > 0:
            kspacing = kspacing_list[0]
        else:
            kspacing = KSPACING
        
        if ktype == 'Line' or ktype == 'L':
            npoints = int(kspacing) if kspacing >= 1 else 20
            kpt_content = f"""K_POINTS
6
Line
0.5 0.0 0.5 {npoints}
0.0 0.0 0.0 {npoints}
0.5 0.5 0.5 {npoints}
0.5 0.25 0.75 {npoints}
0.375 0.375 0.75 {npoints}
0.0 0.0 0.0 1
"""
            with open('KLINES', 'w') as f:
                f.write(kpt_content)
        else:
            # 导入 K点计算模块
            from abacus.kpoints_utils import calculate_kpoints_from_cell
            Kpoints = calculate_kpoints_from_cell(stru.cell, kspacing, ktype='Gamma')
            Kpoints += [0, 0, 0]
            input_obj.set(knumber=0, kmode='Gamma', kpts=Kpoints)
            input_obj.write_input_kpt()
    
    # 生成 STRU 文件
    # 对于 Relax/Band/Dos，使用优化后的结构
    if use_optimized_stru and optimized_stru_path:
        import shutil
        shutil.copy2(optimized_stru_path, work_dir / 'STRU')
        click_echo(f"[INFO] Copied optimized STRU from previous stage")
    else:
        # 检测是否需要 FR 赝势（SOC 计算）
        need_fr = False
        lspinorb_value = stage_params.get('lspinorb', 0)
        if isinstance(lspinorb_value, str):
            lspinorb_value = int(lspinorb_value)
        if lspinorb_value == 1:
            need_fr = True
            click_echo(f"[INFO] SOC calculation detected (lspinorb=1), will use FR pseudopotentials")
        
        # 其他阶段正常生成 STRU
        # 注意：potential_name='PotSG15' 用于静态字典回退
        # 如果是动态目录，会自动使用 pick_upf 选择最佳版本
        write_input_stru(stru=stru,
                        pseudo_dir=POTPATH,
                        basis_dir=ORBPATH,
                        potential_name='PotSG15',
                        basis_name='SG15std',
                        coordinates_type='Direct',
                        spin=default_nspin,
                        filename='STRU',
                        copy_files=False,
                        need_fr=need_fr)

