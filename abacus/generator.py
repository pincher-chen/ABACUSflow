#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
输入文件生成器模块
包含 INPUT、KPT、STRU 文件生成的核心逻辑
"""

import os
from pathlib import Path

try:
    from ase.io import read
except ImportError:
    read = None

# 导入本地模块
from abacus.stru_utils import write_input_stru, get_total_valence_electrons

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
    stru = read(stru_file, format='vasp')
    ntype = len(set(stru.get_chemical_symbols()))
    
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
    default_nspin = 2
    if stage.lower() != 'test_spin':
        spin_on_file = work_dir.parent / 'SPIN_ON'
        spin_off_file = work_dir.parent / 'SPIN_OFF'
        
        if spin_off_file.exists():
            default_nspin = 1
        elif spin_on_file.exists():
            default_nspin = 2
    
    # 计算 nbands
    total_ne = get_total_valence_electrons(stru)
    if default_nspin == 1:
        nbands = int(max(total_ne / 2 * 1.2, total_ne / 2 + 4))
    else:
        nbands = int(max(total_ne * 1.2, total_ne + 4))
    input_obj.set(nbands=nbands)
    
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
    
    # 覆盖 nspin
    if stage.lower() != 'test_spin':
        input_obj.set(nspin=default_nspin)
    
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
        # 其他阶段正常生成 STRU
        write_input_stru(stru=stru,
                        pseudo_dir=POTPATH,
                        basis_dir=ORBPATH,
                        potential_name='PotSG15',
                        basis_name='SG15std',
                        coordinates_type='Direct',
                        spin=default_nspin,
                        filename='STRU',
                        copy_files=False)

