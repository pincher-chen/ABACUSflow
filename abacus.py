#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ABACUS Workflow Manager and CLI Tool
主入口文件 - 只包含 CLI 命令定义，具体功能从模块导入

主要命令：
1. workflow: 生成完整计算工作流
2. generate: 生成输入文件  
3. errors: 检查计算错误
4. converge: 检查收敛状态
5. spin: 判断是否需要自旋
6. update: 更新输入参数
7. summary: 汇总计算结果

Usage:
    python abacus.py workflow InputPoscar/ work_cal/
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
def generate(work_dir, stage, stru_file):
    """生成 ABACUS 输入文件（INPUT, KPT, STRU）"""
    work_dir = Path(work_dir).absolute()
    
    # 查找结构文件
    if not stru_file:
        for candidate in ['STRU.vasp', '../STRU.vasp', 'POSCAR', 'CONTCAR']:
            if (work_dir / candidate).exists():
                stru_file = str(work_dir / candidate)
                break
    
    if not stru_file or not Path(stru_file).exists():
        click.echo("[ERROR] No structure file found", err=True)
        sys.exit(1)
    
    click.echo(f"[INFO] Generating input files for stage: {stage}")
    click.echo(f"[INFO] Structure file: {stru_file}")
    click.echo(f"[INFO] Work directory: {work_dir}")
    
    try:
        generate_input_files(work_dir, stage, stru_file, get, INCAR_TEMPLATE, WORKFLOW, click.echo)
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
    
    failure_patterns = ['Failed to', 'failed to', 'Cannot', 'cannot']
    if any(p in content for p in failure_patterns):
        if 'ABACUS_NOTICE_ERROR' not in errors_found:
            errors_found.append('OPERATION_FAILED')
    
    if 'ERROR' in content or 'Error' in content:
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
        if 'Relaxation is converged, but the SCF is unconverged' in content:
            click.echo("[WARNING] Relax (ionic) converged, but SCF unconverged - results unreliable")
            is_converged = False
        elif 'Relaxation is not converged yet' in content:
            click.echo("[WARNING] Relax (ionic) NOT converged")
            is_converged = False
        elif ('Relaxation is converged!' in content or 
              'Ion relaxation is converged' in content or 
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
# 主入口（支持旧格式）
# ============================================================================

if __name__ == '__main__':
    # 兼容旧格式: python abacus.py InputPoscar/ work_cal/ --workflow
    args = sys.argv[1:]
    
    if len(args) >= 2 and not args[0].startswith('-'):
        known_commands = ['workflow', 'generate', 'errors', 'converge', 'spin', 'update', 'summary']
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
