#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ABACUS Workflow Manager and CLI Tool
主入口文件 - 只包含 CLI 命令定义，具体功能从模块导入

主要命令：
1. workflow: 生成完整计算工作流
2. single: 单步计算模式（不使用workflow，只执行一个step）
3. generate: 生成输入文件  
4. errors: 检查计算错误
5. converge: 检查收敛状态
6. spin: 判断是否需要自旋
7. update: 更新输入参数
8. summary: 汇总计算结果

Usage:
    python abacus.py workflow InputPoscar/ work_cal/
    python abacus.py single batch_poscar/ batch_work_single/ -t Scf -k 0.02
    python abacus.py converge --work_dir .
    python abacus.py errors --work_dir .
"""

import os
import sys
import click
import re
import shutil
import json
from pathlib import Path

# 导入功能模块
from abacus.stru_utils import get_total_valence_electrons
from abacus.input_utils import parse_input_file, write_input_file
from abacus.generator import generate_input_files

# 导入配置
try:
    from config import CONDOR, get, INCAR_TEMPLATE, WORKFLOW
except ImportError:
    import configparser
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config', 'condor.ini')
    CONDOR = configparser.ConfigParser()
    if os.path.exists(config_path):
        CONDOR.read(config_path)
    def get(section, key, fallback=None):
        try:
            return CONDOR.get(section, key, fallback=fallback)
        except (configparser.NoSectionError, configparser.NoOptionError):
            return fallback
    INCAR_TEMPLATE = {}
    WORKFLOW = {}
    # 加载 workflow.json
    try:
        workflow_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config', 'workflow.json')
        if os.path.exists(workflow_path):
            with open(workflow_path, 'r') as f:
                WORKFLOW = json.load(f)
    except:
        pass

# ============================================================================
# 辅助函数
# ============================================================================

def get_atom_count_from_structure(file_path):
    """
    从结构文件中获取原子数（支持多种格式：VASP, CIF, XYZ等）
    
    Args:
        file_path: 结构文件路径
    
    Returns:
        int: 原子总数
    """
    file_path = str(file_path)
    file_lower = file_path.lower()
    
    try:
        from ase.io import read
        
        # 根据文件扩展名选择格式
        if file_lower.endswith('.mcif'):
            # mcif 是磁性 CIF 格式；优先用 mcif，不支持则退回 cif
            try:
                atoms = read(file_path, format='mcif')
            except Exception:
                atoms = read(file_path, format='cif')
        elif file_lower.endswith('.cif'):
            atoms = read(file_path, format='cif')
        elif file_lower.endswith('.xyz'):
            atoms = read(file_path, format='xyz')
        elif file_lower.endswith(('.vasp', '.poscar', '.contcar')):
            atoms = read(file_path, format='vasp')
        else:
            # 尝试自动检测格式
            atoms = read(file_path)
        
        return len(atoms)
    except Exception as e:
        print(f"[WARNING] Failed to read with ASE: {e}")
        
        # 手动解析POSCAR格式（作为后备）
        if not file_lower.endswith('.cif'):
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                
                # 第6行或第7行是原子数（取决于是否有元素符号行）
                for i in range(5, 8):
                    if i < len(lines):
                        try:
                            counts = [int(x) for x in lines[i].split()]
                            if counts:
                                return sum(counts)
                        except ValueError:
                            continue
            except:
                pass
        
        print(f"[ERROR] Cannot determine atom count for {file_path}")
        return 0


# 保持向后兼容
def get_atom_count_from_vasp(file_path):
    """向后兼容的别名"""
    return get_atom_count_from_structure(file_path)


def auto_select_resources(atom_count):
    """
    根据原子数自动选择计算资源，进程数不超过 condor.ini 中的 CORES_PER_NODE。
    
    规则：
    - 1~3个原子: 16核, 1节点
    - 4~8个原子: 32核, 1节点
    - 9~20个原子: 64核, 1节点
    - 21~60个原子: 64核, 1节点  (上限由 CORES_PER_NODE 决定)
    - 超过60个原子: 128核, 2节点
    
    Args:
        atom_count: 原子数
    
    Returns:
        tuple: (node_num, proc_num)
    """
    cores_per_node = int(get('ALLOW', 'CORES_PER_NODE', fallback=64))

    if atom_count <= 3:
        nodes, procs = 1, 16
    elif atom_count <= 8:
        nodes, procs = 1, 32
    elif atom_count <= 20:
        nodes, procs = 1, 64
    elif atom_count <= 60:
        nodes, procs = 1, 88
    else:
        nodes, procs = 2, 128

    # 单节点时，进程数不得超过该节点的核心数上限
    if nodes == 1:
        procs = min(procs, cores_per_node)

    return nodes, procs


# ============================================================================
# CLI 主程序
# ============================================================================

@click.group()
def cli():
    """ABACUS Workflow Manager and CLI Tool"""
    pass

# ============================================================================
# workflow 命令
# ============================================================================

@cli.command()
@click.argument('stru_path', type=click.Path(exists=True))
@click.argument('work_dir', type=click.Path())
def workflow(stru_path, work_dir):
    """生成完整ABACUS计算工作流"""
    click.echo("[INFO] Running in WORKFLOW mode")
    
    try:
        from abacus.submit import AbacusFlowManager
    except ImportError:
        click.echo("[ERROR] Cannot import submit module")
        sys.exit(1)
    
    manager = AbacusFlowManager(work_dir)
    scripts = manager.generate_scripts(stru_path=stru_path)
    
    if scripts:
        click.echo(f"[INFO] Generated workflow scripts: {len(scripts)} files")
        click.echo(f"[INFO] Example: yhbatch {scripts[0]}")
    else:
        click.echo("[ERROR] No scripts generated")

# ============================================================================
# generate 命令
# ============================================================================

@cli.command()
@click.option('--work_dir', default='.', help='Working directory')
@click.option('--stage', default='Scf', help='Calculation stage')
@click.option('--stru_file', help='Structure file path')
@click.option('--spin', '-spin', type=int, default=1, show_default=True,
              help='Spin setting (1=non-magnetic, 2=collinear, 4=non-collinear)')
@click.option('--guess-mag', is_flag=True, help='Guess initial magnetic moment (2.0) if not in CIF')
@click.option('--kval', '-k', type=float, default=None,
              help='K 点密度参数 (Å⁻¹)；不指定则沿用 workflow.json 中 stage.kval[0]')
def generate(work_dir, stage, stru_file, spin, guess_mag, kval):
    """生成 ABACUS 输入文件（INPUT, KPT, STRU）"""
    work_dir = Path(work_dir).absolute()
    
    # 查找结构文件
    if not stru_file:
        for candidate in ['STRU.vasp', '../STRU.vasp', 'POSCAR', 'CONTCAR',
                          'structure.mcif', 'structure.cif']:
            if (work_dir / candidate).exists():
                stru_file = str(work_dir / candidate)
                break
    
    if not stru_file or not Path(stru_file).exists():
        click.echo("[ERROR] No structure file found", err=True)
        sys.exit(1)
    
    click.echo(f"[INFO] Generating input files for stage: {stage}")
    click.echo(f"[INFO] Structure file: {stru_file}")
    click.echo(f"[INFO] Work directory: {work_dir}")
    
    if spin is not None:
        click.echo(f"[INFO] Spin setting: {spin}")
    if guess_mag:
        click.echo(f"[INFO] Guess magnetic moment: enabled")
    if kval is not None:
        click.echo(f"[INFO] K-spacing override: {kval} Å⁻¹")
    
    try:
        generate_input_files(work_dir, stage, stru_file, get, INCAR_TEMPLATE, WORKFLOW, click.echo, 
                           spin=spin, guess_mag=guess_mag, kspacing=kval)
        click.echo("[INFO] Input files generated successfully")
    except Exception as e:
        click.echo(f"[ERROR] Failed to generate: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)

# ============================================================================
# errors 命令
# ============================================================================

@cli.command()
@click.option('--work_dir', default='.', help='Working directory')
def errors(work_dir):
    """检查 ABACUS 计算错误"""
    work_dir = Path(work_dir).absolute()
    
    # 查找日志
    log_files = []
    if (work_dir / 'running.log').exists():
        log_files.append(work_dir / 'running.log')
    
    if not log_files:
        for out_dir in work_dir.glob('OUT.*'):
            log_files.extend(out_dir.glob('running*.log'))
            log_files.extend(out_dir.glob('warning.log'))
    
    if not log_files:
        click.echo("[WARNING] No log file found")
        return
    
    log_file = max(log_files, key=lambda p: p.stat().st_mtime)
    click.echo(f"[INFO] Checking errors in: {log_file}")
    
    content = log_file.read_text()
    errors_found = []
    error_details = []
    
    # 检查各类错误
    if 'NOTICE' in content and '!!!!!!!!!!!!' in content:
        notice_pattern = r'!+\s*NOTICE\s*!+\s*(.*?)\s*!+\s*NOTICE\s*!+'
        for match in re.findall(notice_pattern, content, re.DOTALL):
            for line in match.split('\n'):
                line = line.strip()
                if (line and 'CHECK IN FILE' not in line and
                    '!!!!!!!!!!!!' not in line and
                    not line.startswith('---') and
                    'TIME STATISTICS' not in line):
                    error_details.append(line)
        if error_details:
            errors_found.append('ABACUS_NOTICE_ERROR')
    
    if 'Something wrong' in content:
        if 'ABACUS_NOTICE_ERROR' not in errors_found:
            errors_found.append('INPUT_ERROR')
    
    # 检查严重的操作失败（排除警告信息）
    # 注意：不要把警告信息当作错误
    critical_failure_patterns = [
        'Failed to allocate',  # 内存分配失败
        'Failed to read',      # 文件读取失败
        'Failed to write',     # 文件写入失败
        'Failed to open',      # 文件打开失败
        'Cannot allocate',     # 无法分配内存
        'Cannot open file',    # 无法打开文件
        'Cannot find file',    # 找不到文件
    ]
    
    # 只检查严重错误，忽略普通警告
    if any(p in content for p in critical_failure_patterns):
        if 'ABACUS_NOTICE_ERROR' not in errors_found:
            errors_found.append('OPERATION_FAILED')
    
    # 检查 ERROR 关键词（但要更谨慎）
    # 排除一些不是真正错误的情况
    if 'ERROR' in content or 'Error' in content:
        # 排除常见的非错误情况
        if not any(exclude in content for exclude in [
            'error in',  # 可能只是变量名
            'ERROR_THRESHOLD',  # 可能是参数名
        ]):
            if not errors_found:
                errors_found.append('GENERIC_ERROR')
    
    if 'convergence has NOT been achieved' in content:
        errors_found.append('SCF_NOT_CONVERGED')
    
    if 'relax is not converged' in content:
        errors_found.append('RELAX_NOT_CONVERGED')
    
    if 'segmentation fault' in content.lower():
        errors_found.append('SEGMENTATION_FAULT')
    
    if 'out of memory' in content.lower() or 'cannot allocate' in content.lower():
        errors_found.append('OUT_OF_MEMORY')
    
    if errors_found:
        click.echo(f"[ERROR] Found errors: {', '.join(errors_found)}")
        if error_details:
            click.echo("[ERROR] Error details:")
            for detail in error_details[:5]:
                click.echo(f"  - {detail}")
        
        error_file = work_dir / 'error.txt'
        error_content = '\n'.join(errors_found)
        if error_details:
            error_content += '\n\nDetails:\n' + '\n'.join(error_details)
        error_file.write_text(error_content)
    else:
        click.echo("[INFO] No obvious errors detected")

# ============================================================================
# converge 命令
# ============================================================================

@cli.command()
@click.option('--work_dir', default='.', help='Working directory')
def converge(work_dir):
    """检查 ABACUS 计算是否收敛"""
    work_dir = Path(work_dir).absolute()
    
    # 查找日志
    log_files = []
    for out_dir in work_dir.glob('OUT.*'):
        log_files.extend(out_dir.glob('running*.log'))
    
    if not log_files:
        log_files = list(work_dir.glob('running*.log'))
    
    if not log_files:
        click.echo("[WARNING] No log file found")
        return
    
    log_file = max(log_files, key=lambda p: p.stat().st_mtime)
    click.echo(f"[INFO] Checking convergence in: {log_file}")
    
    content = log_file.read_text()
    is_finished = False
    is_converged = False
    
    if 'Total  Time' in content or 'Total time' in content:
        is_finished = True
        click.echo("[INFO] Calculation finished")
    
    is_relax_calc = ('Relax' in str(log_file) or 
                     'relax' in content.lower()[:5000] or
                     'Ion relaxation' in content or
                     'Relaxation is not converged' in content)
    
    scf_converged = 'charge density convergence is achieved' in content
    if scf_converged:
        click.echo("[INFO] SCF (electronic) converged")
    
    if is_relax_calc:
        # 重要：需要检查最后的收敛状态，而不是整个文件中是否出现过
        # 因为 Relax 计算会多次输出 "Relaxation is not converged yet"
        # 只有最后一次的状态才是最终结果
        
        # 查找所有 Relaxation 相关的行
        relax_lines = [line for line in content.split('\n') if 'Relaxation is' in line]
        
        if relax_lines:
            # 取最后一行判断最终状态
            last_relax_status = relax_lines[-1]
            
            if 'Relaxation is converged, but the SCF is unconverged' in last_relax_status:
                click.echo("[WARNING] Relax (ionic) converged, but SCF unconverged - results unreliable")
                is_converged = False
            elif 'Relaxation is not converged yet' in last_relax_status:
                click.echo("[WARNING] Relax (ionic) NOT converged")
                is_converged = False
            elif 'Relaxation is converged!' in last_relax_status:
                click.echo("[INFO] Relax (ionic) converged")
                is_converged = True
            else:
                click.echo("[WARNING] Relax convergence status unclear")
                is_converged = False
        else:
            # 兼容旧版本 ABACUS 的输出格式
            if ('Ion relaxation is converged' in content or 
                'RELAX CONVERGED' in content):
                click.echo("[INFO] Relax (ionic) converged")
                is_converged = True
            else:
                click.echo("[WARNING] Relax convergence status unclear")
                is_converged = False
    else:
        if scf_converged:
            click.echo("[INFO] SCF calculation converged")
            is_converged = True
    
    if is_finished and is_converged:
        (work_dir / 'converge.txt').touch()
        click.echo("[INFO] Created converge.txt")
    elif is_finished:
        click.echo("[WARNING] Calculation finished but did not converge")
    else:
        click.echo("[WARNING] Calculation not finished")

# ============================================================================
# spin 命令
# ============================================================================

@cli.command()
@click.option('--work_dir', default='.', help='Working directory')
def spin(work_dir):
    """判断是否需要自旋"""
    work_dir = Path(work_dir).absolute()
    
    # 查找日志
    log_files = []
    for out_dir in work_dir.glob('OUT.*'):
        log_files.extend(out_dir.glob('running*.log'))
    
    if not log_files:
        log_files = list(work_dir.glob('running*.log'))
    
    if not log_files:
        click.echo("[WARNING] No log file found")
        return
    
    log_file = max(log_files, key=lambda p: p.stat().st_mtime)
    click.echo(f"[INFO] Checking magnetic moment in: {log_file}")
    
    content = log_file.read_text()
    
    mag_pattern = r'total magnetism \(Bohr mag/cell\)\s*=\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)'
    mag_matches = re.findall(mag_pattern, content)
    
    if mag_matches:
        final_mag = float(mag_matches[-1])
        click.echo(f"[INFO] Final magnetic moment: {final_mag}")
        
        if abs(final_mag) > 0.004:
            spin_on_file = work_dir.parent / 'SPIN_ON'
            spin_on_file.write_text(str(final_mag))
            
            spin_off_file = work_dir.parent / 'SPIN_OFF'
            if spin_off_file.exists():
                spin_off_file.unlink()
            
            click.echo(f"[INFO] Spin needed: |mag| = {abs(final_mag)} > 0.004")
            click.echo(f"[INFO] Created {spin_on_file}")
        else:
            spin_off_file = work_dir.parent / 'SPIN_OFF'
            spin_off_file.touch()
            
            spin_on_file = work_dir.parent / 'SPIN_ON'
            if spin_on_file.exists():
                spin_on_file.unlink()
            
            click.echo(f"[INFO] Spin not needed: |mag| = {abs(final_mag)} <= 0.004")
            click.echo(f"[INFO] Created {spin_off_file}")
    else:
        click.echo("[WARNING] Could not extract magnetic moment")

# ============================================================================
# update 命令
# ============================================================================

@cli.command()
@click.option('--work_dir', default='.', help='Working directory')
@click.option('--try_num', type=int, help='Current retry number')
@click.option('--stage', help='Current calculation stage')
def update(work_dir, try_num, stage):
    """更新输入文件（用于重试）"""
    work_dir = Path(work_dir).absolute()
    
    click.echo(f"[INFO] Updating input files in: {work_dir}")
    click.echo(f"[INFO] Try number: {try_num}, Stage: {stage}")
    
    # 备份文件
    if try_num is not None:
        for fname in ['INPUT', 'KPT', 'STRU']:
            f = work_dir / fname
            if f.exists():
                shutil.copy2(f, work_dir / f'{fname}{try_num}')
                click.echo(f"[INFO] Backed up {fname} -> {fname}{try_num}")
        
        # Relax: 更新 STRU（使用上一次迭代的优化结构）
        if stage and ('relax' in stage.lower()):
            for out_dir in work_dir.glob('OUT.*'):
                stru_ion_d = out_dir / 'STRU_ION_D'
                if stru_ion_d.exists():
                    shutil.copy2(stru_ion_d, work_dir / 'STRU')
                    click.echo(f"[INFO] Updated STRU from last iteration's optimized geometry:")
                    click.echo(f"       Copied {out_dir.name}/STRU_ION_D -> STRU")
                    break
    
    # 读取 workflow.json
    if stage:
        workflow_file = Path(__file__).parent / 'config' / 'workflow.json'
        if workflow_file.exists():
            workflow = json.loads(workflow_file.read_text())
            stage_config = workflow.get(stage, {})
        else:
            stage_config = {}
    else:
        stage_config = {}
    
    # 更新 KPT
    if try_num is not None and 'kval' in stage_config:
        kpt_file = work_dir / 'KPT'
        stru_file = work_dir / 'STRU'
        
        if kpt_file.exists() and stru_file.exists():
            kspacing_list = stage_config['kval']
            ktype = stage_config.get('ktype', 'Gamma')
            next_try_num = try_num + 1
            
            if next_try_num < len(kspacing_list):
                new_kspacing = kspacing_list[next_try_num]
                click.echo(f"[INFO] Updating K-points with kspacing={new_kspacing} Å⁻¹")
                
                try:
                    from abacus.kpoints_utils import calculate_kpoints_from_stru, update_kpt_file
                    kpoints = calculate_kpoints_from_stru(str(stru_file), new_kspacing, ktype)
                    update_kpt_file(str(kpt_file), kpoints)
                except Exception as e:
                    click.echo(f"[WARNING] Failed to update K-points: {e}")
    
    # 调整 INPUT 参数
    input_file = work_dir / 'INPUT'
    error_file = work_dir / 'error.txt'
    converge_file = work_dir / 'converge.txt'
    
    if not input_file.exists():
        click.echo("[WARNING] INPUT file not found")
        return
    
    input_content = input_file.read_text()
    input_params = parse_input_file(input_content)
    updated = False
    
    is_converged = converge_file.exists()
    
    error_types = []
    if error_file.exists():
        error_content = error_file.read_text()
        error_types = error_content.split('\n')[0].split(', ') if error_content else []
        click.echo(f"[INFO] Detected errors: {error_types}")
    
    if not is_converged:
        # SCF 未收敛
        if 'SCF_NOT_CONVERGED' in error_types:
            click.echo("[INFO] Adjusting for SCF convergence")
            if 'mixing_beta' in input_params:
                old_val = float(input_params['mixing_beta'])
                new_val = max(0.1, old_val * 0.7)
                input_params['mixing_beta'] = f"{new_val:.2f}"
                click.echo(f"  - mixing_beta: {old_val} -> {new_val}")
                updated = True
            
            if 'scf_nmax' in input_params:
                old_val = int(input_params['scf_nmax'])
                new_val = min(300, old_val + 50)
                input_params['scf_nmax'] = str(new_val)
                click.echo(f"  - scf_nmax: {old_val} -> {new_val}")
                updated = True
        
        # Relax 未收敛
        if 'RELAX_NOT_CONVERGED' in error_types:
            click.echo("[INFO] Adjusting for Relax convergence")
            if 'force_thr_ev' in input_params:
                old_val = float(input_params['force_thr_ev'])
                new_val = min(0.1, old_val * 1.5)
                input_params['force_thr_ev'] = f"{new_val:.3f}"
                click.echo(f"  - force_thr_ev: {old_val} -> {new_val}")
                updated = True
            
            if 'relax_nmax' in input_params:
                old_val = int(input_params['relax_nmax'])
                new_val = old_val + 20
                input_params['relax_nmax'] = str(new_val)
                click.echo(f"  - relax_nmax: {old_val} -> {new_val}")
                updated = True
    
    if updated:
        new_content = write_input_file(input_params, input_content)
        input_file.write_text(new_content)
        click.echo("[INFO] INPUT file updated")
    else:
        click.echo("[INFO] No parameter updates needed")

# ============================================================================
# summary 命令
# ============================================================================

@cli.command()
@click.option('--root', required=True, help='Root directory')
def summary(root):
    """汇总计算结果"""
    root = Path(root).absolute()
    
    click.echo(f"[INFO] Summary for: {root}")
    
    stat_log = root / 'stat.log'
    if stat_log.exists():
        click.echo("\n=== Calculation Status ===")
        click.echo(stat_log.read_text())
    else:
        click.echo("[WARNING] No stat.log found")
    
    stages = ['Test_spin', 'Coarse_relax', 'Relax', 'Scf', 'Band', 'Dos']
    click.echo("\n=== Stage Status ===")
    for stage in stages:
        stage_dir = root / stage
        if stage_dir.exists():
            if (stage_dir / 'converge.txt').exists():
                click.echo(f"{stage}: ✓ Converged")
            else:
                click.echo(f"{stage}: ✗ Not converged")
        else:
            click.echo(f"{stage}: - Not run")

# ============================================================================
# single 命令 - 单步计算（不使用workflow）
# ============================================================================

@cli.command()
@click.argument('stru_path', type=click.Path(exists=True))
@click.argument('work_dir', type=click.Path())
@click.option('--template', '-t', required=True, help='使用的模板名称（如 Scf, Relax, Band 等）')
@click.option('--kval', '-k', type=float, default=0.02, help='K点密度参数（Å⁻¹，默认0.02）')
@click.option('--node-num', '-n', type=int, default=None, help='节点数（不指定则自动根据原子数设置）')
@click.option('--proc-num', '-p', type=int, default=None, help='进程数（不指定则自动根据原子数设置）')
@click.option('--auto-resource/--no-auto-resource', default=True, help='是否自动根据原子数设置资源（默认开启）')
@click.option('--omp-threads', type=int, default=1, help='OMP线程数（默认1）')
@click.option('--spin', '-spin', type=int, default=1, show_default=True,
              help='Spin setting (1=non-magnetic, 2=collinear, 4=non-collinear)')
@click.option('--guess-mag', is_flag=True, help='Guess initial magnetic moment (2.0) if not in CIF')
@click.option('--dry-run', is_flag=True, help='只生成脚本，不提交')
@click.option('--filter', '-f', 'name_filter', default=None, multiple=True,
              help='只处理文件名包含指定字符串的结构（可多次使用，取并集）。支持 glob 通配符，如 "Fe*"。')
def single(stru_path, work_dir, template, kval, node_num, proc_num, auto_resource, omp_threads, spin, guess_mag, dry_run, name_filter):
    """
    单步计算模式 - 不使用workflow，只执行一个step
    
    不需要retry等复杂逻辑，只进行一次计算。
    
    示例:
      # 使用Scf模板，kval=0.02，自动设置资源
      python abacus.py single batch_poscar/ batch_work_single/ -t Scf -k 0.02
      
      # 手动指定节点和进程数
      python abacus.py single batch_poscar/ batch_work_single/ -t Relax -k 0.04 -n 1 -p 64
      
      # 不自动设置资源，必须手动指定
      python abacus.py single batch_poscar/ work_out/ -t Scf --no-auto-resource -n 2 -p 128
      
      # 只处理特定结构（名称包含 FeI2 的文件）
      python abacus.py single cif/ work_out/ -t Scf-Soc-PW -f "FeI2"
      
      # 多个 filter 取并集
      python abacus.py single cif/ work_out/ -t Scf-Soc-PW -f "FeI2" -f "MnO"
      
      # 支持 glob 通配符
      python abacus.py single cif/ work_out/ -t Scf-Soc-PW -f "3.15.*" -f "0.1.*"
    """
    from datetime import datetime
    
    click.echo("=" * 70)
    click.echo("🚀 SINGLE STEP CALCULATION MODE")
    click.echo("=" * 70)
    click.echo(f"[INFO] Template: {template}")
    click.echo(f"[INFO] K-spacing: {kval} Å⁻¹")
    if spin is not None:
        click.echo(f"[INFO] Spin setting: {spin}")
    if guess_mag:
        click.echo(f"[INFO] Guess magnetic moment: enabled")
    
    # 检查模板是否存在
    if template not in INCAR_TEMPLATE:
        available_templates = list(INCAR_TEMPLATE.keys())
        click.echo(f"[ERROR] Template '{template}' not found!", err=True)
        click.echo(f"[INFO] Available templates: {', '.join(available_templates)}", err=True)
        sys.exit(1)
    
    # 准备结构文件列表
    stru_path = Path(stru_path)
    if stru_path.is_file():
        structure_files = [stru_path]
    else:
        # 目录模式，查找所有结构文件
        suffixes = get('STRU', 'SUFFIX', '*.vasp *.poscar').split()
        structure_files = []
        for suffix in suffixes:
            structure_files.extend(stru_path.glob(suffix))
        # 去重并排序
        structure_files = sorted(set(structure_files))
    
    if not structure_files:
        click.echo(f"[ERROR] No structure files found in {stru_path}", err=True)
        sys.exit(1)

    # --filter 过滤：名称匹配任意一个 pattern 即保留（glob 通配符或子字符串）
    if name_filter:
        import fnmatch
        total_before = len(structure_files)
        def _matches(f):
            stem = f.stem
            for pat in name_filter:
                # 先尝试 glob 匹配，再尝试子字符串包含
                if fnmatch.fnmatch(stem, pat) or pat in stem:
                    return True
            return False
        structure_files = [f for f in structure_files if _matches(f)]
        click.echo(f"[INFO] Filter '{', '.join(name_filter)}': {total_before} → {len(structure_files)} structure(s)")
        if not structure_files:
            click.echo(f"[ERROR] No structures matched the filter.", err=True)
            sys.exit(1)

    click.echo(f"[INFO] Found {len(structure_files)} structure file(s)")
    
    # 创建工作目录
    work_dir = Path(work_dir).absolute()
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # 获取配置
    conda_path = get('ENV', 'CONDA_PATH', '/XYFS01/nscc-gz_pinchen_1/sf_install/miniconda3')
    conda_env = get('ENV', 'CONDA_ENV', 'dftflow')
    abacus_dir = get('ABACUS', 'ABACUS_DIR', '/XYFS01/nscc-gz_pinchen_1/sf_box/abacus-develop-LTSv3.10.0/bin')
    abacus_exe = get('ABACUS', 'ABACUS_EXE', 'abacus')
    abacus_bin = f"{abacus_dir}/{abacus_exe}"
    modules_str = get('MODULE', 'MODULES', 'intel/oneapi2023.2_noimpi mpi/mpich/4.1.2-icc-oneapi2023.2-ch4')
    partition = get('ALLOW', 'PARTITION', 'deimos')
    
    abacus_py = Path(__file__).absolute()
    
    generated_scripts = []
    skipped_files = []

    # 为每个结构文件生成脚本
    for stru_file in structure_files:
        job_name = stru_file.stem
        job_dir = work_dir / job_name
        job_dir.mkdir(parents=True, exist_ok=True)
        
        # 复制结构文件（保持原始扩展名或使用通用名称）
        stru_ext = stru_file.suffix.lower()
        if stru_ext == '.cif':
            stru_dest = job_dir / 'structure.cif'
        elif stru_ext in ['.vasp', '.poscar', '.contcar']:
            stru_dest = job_dir / 'STRU.vasp'
        else:
            stru_dest = job_dir / f'structure{stru_ext}'
        
        if not stru_dest.exists():
            shutil.copy2(stru_file, stru_dest)
        
        # 获取原子数（支持多种格式）
        atom_count = get_atom_count_from_structure(str(stru_file))
        click.echo(f"[INFO] {job_name}: {atom_count} atoms")

        # 元素支持预检查：不支持的结构直接跳过，避免提交必然失败的作业
        try:
            from ase.io import read as _ase_read
            _stru_ext = str(stru_file).lower()
            if _stru_ext.endswith('.mcif') or _stru_ext.endswith('.cif'):
                _atoms_tmp = _ase_read(str(stru_file), format='cif')
            else:
                _atoms_tmp = _ase_read(str(stru_file))
            from abacus.stru_utils import check_elements_supported
            _ok, _bad = check_elements_supported(_atoms_tmp, 'PotSG15', 'SG15std')
            if not _ok:
                click.echo(f"[SKIP] {job_name}: unsupported elements {_bad}, skipping.", err=True)
                skipped_files.append(f"{job_name} (unsupported elements: {_bad})")
                continue
        except Exception as _e:
            click.echo(f"[WARNING] {job_name}: element pre-check failed ({_e}), proceeding anyway.")

        # 确定计算资源
        if auto_resource and (node_num is None or proc_num is None):
            auto_nodes, auto_procs = auto_select_resources(atom_count)
            final_nodes = node_num if node_num is not None else auto_nodes
            final_procs = proc_num if proc_num is not None else auto_procs
            click.echo(f"[INFO] {job_name}: Auto-selected {final_nodes} node(s), {final_procs} processes")
        else:
            if node_num is None or proc_num is None:
                click.echo(f"[ERROR] --no-auto-resource requires both --node-num and --proc-num", err=True)
                sys.exit(1)
            final_nodes = node_num
            final_procs = proc_num
        
        # 生成脚本
        script_path = job_dir / f"{job_name}_single.sh"
        script_content = _generate_single_script(
            job_name=job_name,
            job_dir=job_dir,
            template=template,
            kval=kval,
            node_num=final_nodes,
            proc_num=final_procs,
            omp_threads=omp_threads,
            abacus_py=abacus_py,
            abacus_bin=abacus_bin,
            conda_path=conda_path,
            conda_env=conda_env,
            modules_str=modules_str,
            partition=partition,
            stru_file=stru_dest,  # 传入实际的结构文件路径
            spin=spin,
            guess_mag=guess_mag
        )
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        
        generated_scripts.append(str(script_path))
        click.echo(f"[INFO] Generated: {script_path}")
    
    click.echo("=" * 70)
    click.echo(f"✅ Generated {len(generated_scripts)} scripts")
    if skipped_files:
        click.echo(f"⚠️  Skipped {len(skipped_files)} structure(s):")
        for s in skipped_files:
            click.echo(f"   - {s}")
    click.echo("=" * 70)
    
    if generated_scripts:
        click.echo("\n📝 Next steps:")
        click.echo(f"   cd {work_dir}")
        click.echo(f"\n   # 提交单个作业:")
        example_script = Path(generated_scripts[0]).name
        example_dir = Path(generated_scripts[0]).parent.name
        click.echo(f"   cd {example_dir} && yhbatch {example_script}")
        click.echo(f"\n   # 批量提交所有作业:")
        click.echo(f"   for dir in */; do")
        click.echo(f"       (cd \"$dir\" && [ -f *_single.sh ] && yhbatch *_single.sh)")
        click.echo(f"   done")


def _generate_single_script(job_name, job_dir, template, kval, node_num, proc_num,
                             omp_threads, abacus_py, abacus_bin, conda_path, 
                             conda_env, modules_str, partition, stru_file, spin=None, guess_mag=False):
    """
    生成单步计算的脚本内容
    """
    from datetime import datetime
    
    # stru_file 已经从外部传入
    
    # 构建 generate 命令的参数（路径用双引号包裹，防止特殊字符如括号引发 bash 语法错误）
    generate_args = f'--work_dir "." --stage {template} --stru_file "{stru_file}" --kval {kval}'
    if spin is not None:
        generate_args += f" --spin {spin}"
    if guess_mag:
        generate_args += " --guess-mag"
    
    script_lines = [
        "#!/bin/bash",
        f"# Single-step calculation script for {job_name}",
        f"# Template: {template}",
        f"# Generated at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"# K-spacing: {kval} Å⁻¹",
        f"# Resources: {node_num} node(s), {proc_num} processes, OMP={omp_threads}",
        "",
        "# ===== 环境初始化 =====",
        f"export OMP_NUM_THREADS={omp_threads}",
        "",
        f"source {conda_path}/etc/profile.d/conda.sh",
        f"conda activate {conda_env}",
        "",
        "# 加载必要的 module",
        f"# MODULES: {modules_str}",
        "export MODULESHOME=/usr/share/modules",
        "export MODULEPATH=/APP/u22/x86/modulepath/Compilers:/APP/u22/x86/modulepath/application",
        "export MODULES_CMD=/usr/lib/x86_64-linux-gnu/modulecmd.tcl",
        "ml() { module ml \"$@\"; }",
        "module() { _module_raw \"$@\" 2>&1; }",
        "_module_raw() { eval `/usr/bin/tclsh8.6 /usr/lib/x86_64-linux-gnu/modulecmd.tcl bash \"$@\"`; }",
        "",
        f"for mod in {modules_str}; do",
        "  module load $mod",
        "done",
        "",
        f'cd "{job_dir}" || exit 1',
        "",
        "echo '[...]SINGLE STEP CALCULATION START!'",
        f"echo '[...]Template: {template}'",
        f"echo '[...]K-spacing: {kval}'",
        'echo "[...]Working directory: $(pwd)"',
        "",
        "# 记录开始时间",
        "START_TIME=$(date +%s)",
        'echo "[INFO] Started at $(date \'+%Y-%m-%d %H:%M:%S\')"',
        "",
        "# ===== 生成输入文件 =====",
        f"echo '[...]Generating INPUT files with template: {template}'",
        f"python {abacus_py} generate {generate_args}",
        "GENERATE_EXIT=$?",
        "if [ $GENERATE_EXIT -ne 0 ]; then",
        "    echo '[ERROR] generate failed (exit code '$GENERATE_EXIT'), aborting calculation.'",
        "    echo 'generate_failed' > status.txt",
        "    exit $GENERATE_EXIT",
        "fi",
        "",
        "# KPT 已经由 generate 命令使用 --kval 直接生成（公式已修复，与 dftflow bit-exact）",
        "",
        "# ===== 运行ABACUS计算 =====",
        f"echo '[...]Running ABACUS on {node_num} node(s), {proc_num} processes'",
        f"yhrun -N {node_num} -n {proc_num} -p {partition} {abacus_bin} > running.log 2>&1",
        "RUN_EXIT_CODE=$?",
        "",
        "# ===== 检查计算结果 =====",
        "echo '[...]Checking calculation result...'",
        f'python "{abacus_py}" errors --work_dir "."',
        f'python "{abacus_py}" converge --work_dir "."',
        "",
        "# 记录结束时间",
        "END_TIME=$(date +%s)",
        "DURATION=$((END_TIME - START_TIME))",
        'HOURS=$(awk "BEGIN {printf \\"%.2f\\", $DURATION/3600}")',
        f"CORE_HOURS=$(awk \"BEGIN {{printf \\\"%.2f\\\", {node_num}*{proc_num}*$HOURS}}\")",
        "",
        'echo ""',
        "echo '======================================================================'",
        "echo '                    CALCULATION SUMMARY'",
        "echo '======================================================================'",
        f'echo "[INFO] Job Name      : {job_name}"',
        f'echo "[INFO] Template      : {template}"',
        f'echo "[INFO] K-spacing     : {kval} Å⁻¹"',
        f'echo "[INFO] Resources     : {node_num} node(s), {proc_num} processes"',
        'echo "[INFO] Duration      : ${DURATION}s (${HOURS}h)"',
        'echo "[INFO] Core-hours    : ${CORE_HOURS}"',
        "",
        "if [ -f 'converge.txt' ]; then",
        "    echo '[✓] Calculation CONVERGED successfully!'",
        "    echo 'success' > status.txt",
        "elif [ -f 'error.txt' ]; then",
        "    echo '[✗] Calculation FAILED with errors'",
        "    cat error.txt",
        "    echo 'failed' > status.txt",
        "else",
        "    echo '[?] Calculation status unclear'",
        "    echo 'unknown' > status.txt",
        "fi",
        "",
        "echo '======================================================================'",
        'echo "[...]SINGLE STEP CALCULATION FINISHED!"',
        "",
    ]
    
    return '\n'.join(script_lines)


# ============================================================================
# resume 命令
# ============================================================================

@cli.command()
@click.argument('work_dir', type=click.Path(exists=True))
@click.option('--from-stage', help='强制从指定阶段开始（忽略自动检测）')
@click.option('--dry-run', is_flag=True, help='仅显示续算计划，不执行')
@click.option('--no-backup', is_flag=True, help='不备份 stat.log 和 time.log')
@click.option('--no-clean', is_flag=True, help='不清理 error.txt')
def resume(work_dir, from_stage, dry_run, no_backup, no_clean):
    """
    智能续算工作流
    
    自动判断失败阶段的处理策略：
    - 首次尝试就失败：清空重算
    - 多次尝试未收敛：继续迭代
    
    示例:
      python abacus.py resume work_cal/hmat_0/
      python abacus.py resume work_cal/hmat_0/ --from-stage Relax
      python abacus.py resume work_cal/hmat_0/ --dry-run
    """
    from abacus.resume_utils import (
        detect_stage_status,
        parse_stat_log,
        determine_resume_point,
        format_status_icon,
        detect_resume_number,
        detect_previous_try_count,
        should_clean_and_restart
    )
    
    work_dir = Path(work_dir).absolute()
    
    # 检查是否是有效的工作目录
    if not (work_dir / 'STRU.vasp').exists():
        click.echo("[ERROR] Invalid work directory (no STRU.vasp found)", err=True)
        sys.exit(1)
    
    # 加载工作流配置
    workflow_stages = list(WORKFLOW.keys())
    
    # 解析当前状态
    click.echo("=" * 70)
    click.echo("📊 Current Workflow Status")
    click.echo("=" * 70)
    
    stat_status = parse_stat_log(work_dir)
    stage_details = []
    
    for stage in workflow_stages:
        stage_dir = work_dir / stage
        file_status = detect_stage_status(stage_dir, stage)
        log_status = stat_status.get(stage, 'not_started')
        
        status_icon = format_status_icon(file_status)
        
        click.echo(f"  {status_icon} {stage:15s} - File: {file_status:12s} | Log: {log_status}")
        stage_details.append((stage, file_status, log_status))
    
    click.echo("=" * 70)
    
    # 确定续算起点
    if from_stage:
        if from_stage not in workflow_stages:
            click.echo(f"[ERROR] Invalid stage: {from_stage}", err=True)
            sys.exit(1)
        resume_stage = from_stage
        reason = f"User specified: --from-stage {from_stage}"
    else:
        resume_index, reason = determine_resume_point(work_dir, workflow_stages)
        if resume_index is None:
            click.echo(f"\n[INFO] {reason}")
            click.echo("[INFO] No need to resume, workflow completed.")
            sys.exit(0)
        resume_stage = workflow_stages[resume_index]
    
    # 显示续算计划
    resume_index = workflow_stages.index(resume_stage)
    remaining_stages = workflow_stages[resume_index:]
    
    click.echo(f"\n🎯 Resume Plan: {reason}")
    click.echo(f"   Starting from: {resume_stage}")
    click.echo(f"   Stages to run: {', '.join(remaining_stages)}")
    
    # 分析续算策略
    click.echo(f"\n📋 Resume Strategy Analysis:")
    for stage in remaining_stages:
        stage_dir = work_dir / stage
        if stage_dir.exists():
            should_clean, clean_reason = should_clean_and_restart(stage_dir, stage)
            resume_num = detect_resume_number(stage_dir)
            prev_tries = detect_previous_try_count(stage_dir)
            
            click.echo(f"\n  📁 {stage}:")
            click.echo(f"     Resume attempt: #{resume_num}")
            click.echo(f"     Previous tries: {prev_tries}")
            click.echo(f"     Action: {'🧹 Clean & Restart' if should_clean else '🔄 Continue'}")
            click.echo(f"     Reason: {clean_reason}")
    
    # 检查是否需要备份和清理
    if not no_backup:
        stat_log = work_dir / 'stat.log'
        time_log = work_dir / 'time.log'
        if stat_log.exists() or time_log.exists():
            click.echo(f"\n📦 Backup Plan:")
            if stat_log.exists():
                click.echo(f"   stat.log will be backed up")
            if time_log.exists():
                click.echo(f"   time.log will be backed up")
    
    if not no_clean:
        error_file = work_dir / resume_stage / 'error.txt'
        if error_file.exists():
            click.echo(f"\n🧹 Cleanup Plan:")
            click.echo(f"   {resume_stage}/error.txt will be removed")
            click.echo(f"   '{resume_stage} failed' will be removed from stat.log")
    
    if dry_run:
        click.echo("\n[DRY-RUN] No actions taken.")
        sys.exit(0)
    
    # 生成续算脚本
    click.echo("\n🔧 Generating resume script...")
    
    try:
        from abacus.submit import AbacusFlowManager
        manager = AbacusFlowManager(work_dir.parent)
        
        script_path = manager.generate_resume_script(
            work_dir=work_dir,
            start_stage=resume_stage,
            clean_errors=(not no_clean),
            no_backup=no_backup
        )
        
        click.echo(f"\n✅ Resume script generated successfully!")
        click.echo(f"   Script: {script_path}")
        click.echo(f"\n📝 Next steps:")
        click.echo(f"   cd {work_dir}")
        click.echo(f"   yhbatch {Path(script_path).name}")
        
    except Exception as e:
        click.echo(f"[ERROR] Failed to generate resume script: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)

# ============================================================================
# batch-resume 命令：批量续算
# ============================================================================

@cli.command('batch-analysis')
@click.argument('batch_dir', type=click.Path(exists=True))
@click.option('--show-all', is_flag=True, help='显示所有作业（包括成功的）')
@click.option('--output', '-o', type=click.Path(), help='输出到文件')
def batch_analysis(batch_dir, show_all, output):
    """
    批量分析作业状态和错误信息
    
    显示每个失败作业的详细错误信息，包括失败阶段和具体错误内容。
    
    示例:
      # 分析 batch_work3 目录
      python abacus.py batch-analysis batch_work3/
      
      # 显示所有作业（包括成功的）
      python abacus.py batch-analysis batch_work3/ --show-all
      
      # 输出到文件
      python abacus.py batch-analysis batch_work3/ -o analysis.txt
    """
    from abacus.resume_utils import (
        parse_stat_log,
        get_running_jobs_from_slurm
    )
    
    batch_dir = Path(batch_dir).absolute()
    
    # 扫描所有子目录
    job_dirs = []
    for item in batch_dir.iterdir():
        if item.is_dir() and (item / 'STRU.vasp').exists():
            job_dirs.append(item)
    
    if not job_dirs:
        click.echo(f"[ERROR] No valid job directories found in {batch_dir}", err=True)
        sys.exit(1)
    
    click.echo("=" * 70)
    click.echo(f"📊 Batch Analysis - {batch_dir.name}")
    click.echo("=" * 70)
    
    # 获取正在运行的作业
    click.echo("\n🔍 Checking Slurm queue...")
    running_jobs = get_running_jobs_from_slurm(batch_dir)
    
    # 统计信息
    needs_resume = []
    already_done = []
    no_stat_log = []
    currently_running = []
    
    # 错误信息字典
    error_details = {}
    
    # 分析每个作业
    for job_dir in sorted(job_dirs):
        job_name = job_dir.name
        
        # 检查是否正在运行
        if job_name in running_jobs:
            currently_running.append(job_name)
            continue
        
        # 检查 stat.log
        stat_log = job_dir / 'stat.log'
        if not stat_log.exists():
            no_stat_log.append(job_name)
            continue
        
        # 解析状态
        stat_status = parse_stat_log(job_dir)
        stat_content = stat_log.read_text()
        workflow_completed = 'Workflow completed' in stat_content
        
        # 判断是否需要续算
        has_failed = any(status in ['failed', 'running'] for status in stat_status.values())
        all_stages_ok = all(status in ['success', 'ignored'] for status in stat_status.values() if status != 'not_started')
        
        if workflow_completed or (all_stages_ok and not has_failed and stat_status):
            already_done.append(job_name)
        else:
            needs_resume.append(job_name)
            
            # 分析错误信息
            error_info = analyze_job_errors(job_dir, stat_status)
            if error_info:
                error_details[job_name] = error_info
    
    # 显示统计
    click.echo(f"\n📊 Summary:")
    click.echo(f"  Total jobs:        {len(job_dirs)}")
    click.echo(f"  Need resume:       {len(needs_resume)} ✨")
    click.echo(f"  Already complete:  {len(already_done)} ✅")
    click.echo(f"  Currently running: {len(currently_running)} 🏃")
    click.echo(f"  No stat.log:       {len(no_stat_log)} ⚠️")
    
    if currently_running:
        click.echo(f"\n🏃 Currently running:")
        for i, job_name in enumerate(currently_running[:10], 1):
            click.echo(f"  {i:3d}. {job_name}")
        if len(currently_running) > 10:
            click.echo(f"  ... and {len(currently_running) - 10} more")
    
    if needs_resume:
        click.echo(f"\n🔄 Jobs needing resume:")
        for i, job_name in enumerate(needs_resume[:20], 1):
            click.echo(f"  {i:3d}. {job_name}")
        if len(needs_resume) > 20:
            click.echo(f"  ... and {len(needs_resume) - 20} more")
    
    # 显示详细错误信息
    if error_details:
        click.echo("\n" + "=" * 70)
        click.echo("❌ Error Details")
        click.echo("=" * 70)
        
        for job_name in sorted(error_details.keys()):
            errors = error_details[job_name]
            click.echo(f"\n📁 Job: {job_name}")
            
            for error in errors:
                stage = error['stage']
                error_type = error.get('error_type', 'UNKNOWN')
                details = error.get('details', [])
                
                click.echo(f"  ├─ Failed Stage: {stage}")
                click.echo(f"  ├─ Error Type:   {error_type}")
                
                if details:
                    click.echo(f"  └─ Details:")
                    for detail in details[:3]:  # 只显示前3行
                        click.echo(f"       {detail}")
                    if len(details) > 3:
                        click.echo(f"       ... ({len(details) - 3} more lines)")
                else:
                    click.echo(f"  └─ No detailed error message found")
    
    if show_all and already_done:
        click.echo("\n" + "=" * 70)
        click.echo("✅ Completed Jobs")
        click.echo("=" * 70)
        for i, job_name in enumerate(already_done[:20], 1):
            click.echo(f"  {i:3d}. {job_name}")
        if len(already_done) > 20:
            click.echo(f"  ... and {len(already_done) - 20} more")
    
    # 输出到文件
    if output:
        with open(output, 'w') as f:
            f.write(f"Batch Analysis Report - {batch_dir.name}\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Total jobs: {len(job_dirs)}\n")
            f.write(f"Need resume: {len(needs_resume)}\n")
            f.write(f"Already complete: {len(already_done)}\n")
            f.write(f"Currently running: {len(currently_running)}\n\n")
            
            if error_details:
                f.write("Error Details:\n")
                f.write("=" * 70 + "\n\n")
                for job_name, errors in sorted(error_details.items()):
                    f.write(f"Job: {job_name}\n")
                    for error in errors:
                        f.write(f"  Stage: {error['stage']}\n")
                        f.write(f"  Error: {error.get('error_type', 'UNKNOWN')}\n")
                        if error.get('details'):
                            f.write(f"  Details: {', '.join(error['details'][:3])}\n")
                        f.write("\n")
        
        click.echo(f"\n✅ Analysis report saved to: {output}")

def analyze_job_errors(job_dir, stat_status):
    """
    分析作业的错误信息
    
    Args:
        job_dir: 作业目录 (Path)
        stat_status: stat.log 解析结果
    
    Returns:
        list: 错误信息列表
    """
    errors = []
    
    # 工作流阶段
    stages = ['Test_spin', 'Coarse_relax', 'Relax', 'Scf', 'Band', 'Dos']
    
    for stage in stages:
        stage_dir = job_dir / stage
        
        # 检查是否失败
        if stat_status.get(stage) == 'failed' or (stage_dir.exists() and (stage_dir / 'error.txt').exists()):
            error_file = stage_dir / 'error.txt'
            
            error_info = {
                'stage': stage,
                'error_type': 'UNKNOWN',
                'details': []
            }
            
            if error_file.exists():
                try:
                    content = error_file.read_text().strip()
                    lines = content.split('\n')
                    
                    # 第一行通常是错误类型
                    if lines:
                        error_info['error_type'] = lines[0]
                    
                    # 其他行是详细信息
                    if len(lines) > 1:
                        # 跳过空行和分隔符
                        details = [line.strip() for line in lines[1:] 
                                  if line.strip() and not line.strip().startswith('Details:')]
                        error_info['details'] = details
                
                except Exception as e:
                    error_info['details'] = [f"Failed to read error.txt: {e}"]
            
            errors.append(error_info)
    
    return errors

@cli.command('batch-resume')
@click.argument('batch_dir', type=click.Path(exists=True))
@click.option('--from-stage', default=None, help='指定从哪个阶段开始续算')
@click.option('--dry-run', is_flag=True, help='只显示需要续算的作业，不生成脚本')
@click.option('--no-backup', is_flag=True, help='不备份原始文件')
@click.option('--no-clean', is_flag=True, help='不清理失败的计算文件')
@click.option('--auto-submit', is_flag=True, help='生成脚本后自动提交到Slurm队列')
@click.option('--parallel', '-j', default=1, type=int, help='并行处理的作业数（默认1）')
def batch_resume(batch_dir, from_stage, dry_run, no_backup, no_clean, auto_submit, parallel):
    """
    批量续算多个作业
    
    自动扫描指定目录下的所有作业目录，对每个作业执行续算操作。
    
    示例:
      # 生成续算脚本
      python abacus.py batch-resume batch_work/
      
      # 生成脚本并自动提交（推荐）⭐
      python abacus.py batch-resume batch_work/ --auto-submit
      
      # 只查看需要续算的作业，不生成脚本
      python abacus.py batch-resume batch_work/ --dry-run
      
      # 从指定阶段开始续算
      python abacus.py batch-resume batch_work/ --from-stage Relax --auto-submit
    """
    from abacus.resume_utils import (
        detect_stage_status,
        parse_stat_log,
        determine_resume_point,
        get_running_jobs_from_slurm
    )
    
    batch_dir = Path(batch_dir).absolute()
    
    # 扫描所有子目录（包含STRU.vasp的目录）
    job_dirs = []
    for item in batch_dir.iterdir():
        if item.is_dir() and (item / 'STRU.vasp').exists():
            job_dirs.append(item)
    
    if not job_dirs:
        click.echo(f"[ERROR] No valid job directories found in {batch_dir}", err=True)
        click.echo("[INFO] A valid job directory should contain STRU.vasp", err=True)
        sys.exit(1)
    
    click.echo("=" * 70)
    click.echo(f"📦 Batch Resume - Found {len(job_dirs)} jobs")
    click.echo("=" * 70)
    
    # 获取正在运行的作业列表
    click.echo("\n🔍 Checking Slurm queue for running jobs...")
    running_jobs = get_running_jobs_from_slurm(batch_dir)
    if running_jobs:
        click.echo(f"   Found {len(running_jobs)} running jobs in queue")
    else:
        click.echo(f"   No running jobs detected")
    
    # 统计信息
    needs_resume = []
    already_done = []
    no_stat_log = []
    currently_running = []
    
    # 分析每个作业的状态
    for job_dir in sorted(job_dirs):
        job_name = job_dir.name
        
        # 检查是否正在运行
        if job_name in running_jobs:
            currently_running.append(job_name)
            continue
        
        # 检查是否有stat.log
        stat_log = job_dir / 'stat.log'
        if not stat_log.exists():
            no_stat_log.append(job_name)
            continue
        
        # 解析状态
        stat_status = parse_stat_log(job_dir)
        
        # 检查是否有 "Workflow completed" 标记
        stat_log = job_dir / 'stat.log'
        workflow_completed = False
        if stat_log.exists():
            stat_content = stat_log.read_text()
            workflow_completed = 'Workflow completed' in stat_content
        
        # 判断是否需要续算
        has_failed = any(status in ['failed', 'running'] for status in stat_status.values())
        # 修复 bug: 状态应该是 'success' 或 'ignored'，而不是 'done'
        all_stages_ok = all(status in ['success', 'ignored'] for status in stat_status.values() if status != 'not_started')
        
        # 如果有 "Workflow completed" 标记或所有阶段都成功，则认为已完成
        if workflow_completed or (all_stages_ok and not has_failed and stat_status):
            already_done.append(job_name)
        else:
            needs_resume.append(job_name)
    
    # 显示统计
    click.echo(f"\n📊 Summary:")
    click.echo(f"  Total jobs:        {len(job_dirs)}")
    click.echo(f"  Need resume:       {len(needs_resume)} ✨")
    click.echo(f"  Already complete:  {len(already_done)} ✅")
    click.echo(f"  Currently running: {len(currently_running)} 🏃")
    click.echo(f"  No stat.log:       {len(no_stat_log)} ⚠️")
    
    if currently_running:
        click.echo(f"\n🏃 Jobs currently running (will not be resumed):")
        for i, job_name in enumerate(currently_running[:20], 1):
            click.echo(f"  {i:3d}. {job_name}")
        if len(currently_running) > 20:
            click.echo(f"  ... and {len(currently_running) - 20} more")
    
    if needs_resume:
        click.echo(f"\n🔄 Jobs needing resume:")
        for i, job_name in enumerate(needs_resume[:20], 1):  # 只显示前20个
            click.echo(f"  {i:3d}. {job_name}")
        if len(needs_resume) > 20:
            click.echo(f"  ... and {len(needs_resume) - 20} more")
    
    if no_stat_log:
        click.echo(f"\n⚠️  Jobs without stat.log (may not have started):")
        for job_name in no_stat_log[:10]:
            click.echo(f"  - {job_name}")
        if len(no_stat_log) > 10:
            click.echo(f"  ... and {len(no_stat_log) - 10} more")
    
    if dry_run:
        click.echo("\n[INFO] Dry-run mode: No scripts will be generated")
        return
    
    if not needs_resume:
        click.echo("\n✅ All jobs are complete! Nothing to resume.")
        return
    
    # 生成续算脚本
    click.echo(f"\n🚀 Generating resume scripts...")
    
    success_count = 0
    failed_count = 0
    
    from abacus.submit import AbacusFlowManager
    from config import CONDOR, WORKFLOW
    
    for job_name in needs_resume:
        job_dir = batch_dir / job_name
        
        try:
            # 使用resume逻辑生成续算脚本
            workflow_stages = list(WORKFLOW.keys())
            resume_stage_idx, reason = determine_resume_point(
                job_dir,
                workflow_stages
            )
            
            # 检查是否需要续算
            if resume_stage_idx is None:
                # 所有阶段都完成，不需要续算
                click.echo(f"  ⏭️  {job_name}: {reason}")
                continue
            
            # 将索引转换为阶段名称
            resume_stage = workflow_stages[resume_stage_idx]
            
            # 如果用户指定了起始阶段，使用用户指定的
            if from_stage:
                if from_stage in WORKFLOW:
                    resume_stage = from_stage
                    reason = f"用户指定从 {from_stage} 开始"
                else:
                    click.echo(f"  ❌ {job_name}: 无效的阶段 {from_stage}")
                    failed_count += 1
                    continue
            
            # 生成续算脚本
            manager = AbacusFlowManager(work_dir=job_dir.parent)
            
            script_path = manager.generate_resume_script(
                work_dir=job_dir,
                start_stage=resume_stage,
                clean_errors=not no_clean,
                no_backup=no_backup
            )
            
            if script_path and Path(script_path).exists():
                click.echo(f"  ✅ {job_name}: {Path(script_path).name}")
                success_count += 1
            else:
                click.echo(f"  ⚠️  {job_name}: 脚本生成失败")
                failed_count += 1
            
        except Exception as e:
            click.echo(f"  ❌ {job_name}: 生成失败 - {str(e)[:50]}")
            failed_count += 1
            import traceback
            if '--debug' in sys.argv:
                traceback.print_exc()
    
    # 最终统计
    click.echo("\n" + "=" * 70)
    click.echo(f"✨ Resume scripts generated:")
    click.echo(f"   Success: {success_count}")
    click.echo(f"   Failed:  {failed_count}")
    click.echo("=" * 70)
    
    if success_count > 0:
        if auto_submit:
            # 自动提交模式
            click.echo(f"\n🚀 Auto-submit mode: Submitting resume scripts...")
            
            try:
                from abacus.submit_manager import AutoSubmitter
                
                submitter = AutoSubmitter(config_file="config/condor.ini")
                submitter.run(str(batch_dir), max_retries=3, include_resume=True)
                
            except KeyboardInterrupt:
                click.echo("\n\n用户中断，正在停止...")
            except Exception as e:
                click.echo(f"\n[ERROR] 自动提交失败: {e}", err=True)
                click.echo("[INFO] 您可以手动提交：python submit_jobs.py batch_work/ --resume")
        else:
            # 手动提交说明
            click.echo(f"\n📝 Next steps:")
            click.echo(f"   cd {batch_dir}")
            click.echo(f"\n   # 方式1: 使用自动提交系统（推荐）⭐")
            click.echo(f"   python submit_jobs.py {Path(batch_dir).name}/ --resume")
            click.echo(f"\n   # 方式2: 手动提交单个作业")
            if needs_resume:
                example_job = needs_resume[0]
                click.echo(f"   cd {example_job} && yhbatch *_resume.sh")
            click.echo(f"\n   # 方式3: 批量手动提交")
            click.echo(f"   for dir in */; do")
            click.echo(f"       (cd \"$dir\" && [ -f *_resume.sh ] && yhbatch *_resume.sh)")
            click.echo(f"   done")

# ============================================================================
# 主入口（支持旧格式）
# ============================================================================

if __name__ == '__main__':
    # 兼容旧格式: python abacus.py InputPoscar/ work_cal/ --workflow
    args = sys.argv[1:]
    
    if len(args) >= 2 and not args[0].startswith('-'):
        known_commands = ['workflow', 'generate', 'errors', 'converge', 'spin', 'update', 'summary', 'resume', 'batch-resume', 'batch-analysis', 'single']
        if args[0] not in known_commands:
            if '--workflow' in args:
                # 转换为新格式
                sys.argv = [sys.argv[0], 'workflow', args[0], args[1]]
                click.echo("[INFO] Using legacy format (auto-converted)")
            else:
                click.echo("[ERROR] Command format changed!")
                click.echo("[INFO] Old: python abacus.py InputPoscar/ work_cal/ --workflow")
                click.echo("[INFO] New: python abacus.py workflow InputPoscar/ work_cal/")
                sys.exit(1)
    
    cli()
