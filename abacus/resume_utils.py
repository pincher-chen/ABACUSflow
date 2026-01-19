#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
续算功能辅助函数模块

提供工作流续算所需的各种检测和处理功能：
- 状态检测
- stat.log 解析
- 续算起点确定
- 日志备份
- 错误清理
"""

import re
import shutil
import subprocess
from pathlib import Path


def detect_stage_status(stage_dir, stage_name):
    """
    检测阶段状态
    
    Args:
        stage_dir: 阶段目录路径 (Path对象)
        stage_name: 阶段名称
    
    Returns:
        状态字符串: 'success', 'ignored', 'failed', 'incomplete', 'not_started'
    """
    if not stage_dir.exists():
        return 'not_started'
    
    # 优先级1: error.txt 存在 → 失败
    if (stage_dir / 'error.txt').exists():
        return 'failed'
    
    # 优先级2: converge.txt 存在 → 成功
    if (stage_dir / 'converge.txt').exists():
        return 'success'
    
    # 优先级3: ignore.txt 存在但无 converge.txt → 被忽略
    if (stage_dir / 'ignore.txt').exists():
        return 'ignored'
    
    # 优先级4: 有运行痕迹但无结论 → 未完成
    if (stage_dir / 'running.log').exists():
        return 'incomplete'
    
    return 'not_started'


def parse_stat_log(work_dir):
    """
    解析 stat.log，返回阶段状态映射
    
    Args:
        work_dir: 工作目录路径 (Path对象)
    
    Returns:
        字典: {stage_name: status}
        例如: {'Test_spin': 'success', 'Relax': 'failed'}
    """
    stat_log = work_dir / 'stat.log'
    if not stat_log.exists():
        return {}
    
    status_map = {}
    try:
        for line in stat_log.read_text().strip().split('\n'):
            if not line or line.startswith('#'):
                continue
            
            # 解析格式: "Test_spin success" 或 "Coarse_relax success (ignored)"
            parts = line.split()
            if len(parts) >= 2:
                stage = parts[0]
                if 'success' in line:
                    status_map[stage] = 'ignored' if 'ignored' in line else 'success'
                elif 'failed' in line:
                    status_map[stage] = 'failed'
    except Exception as e:
        print(f"[WARNING] Failed to parse stat.log: {e}")
    
    return status_map


def determine_resume_point(work_dir, workflow_stages):
    """
    确定续算起点
    
    Args:
        work_dir: 工作目录路径 (Path对象)
        workflow_stages: 工作流阶段列表
    
    Returns:
        (resume_stage_index, reason) 元组
        如果所有阶段都完成，返回 (None, reason)
    """
    stat_status = parse_stat_log(work_dir)
    
    for i, stage in enumerate(workflow_stages):
        stage_dir = work_dir / stage
        stage_status = detect_stage_status(stage_dir, stage)
        stat_log_status = stat_status.get(stage, 'not_started')
        
        # 情况1: stat.log 记录失败 → 从这里开始
        if stat_log_status == 'failed':
            return i, f"{stage} failed in previous run"
        
        # 情况2: 阶段未开始 → 从这里开始
        if stage_status == 'not_started':
            return i, f"{stage} not started yet"
        
        # 情况3: 阶段未完成 → 从这里开始
        if stage_status == 'incomplete':
            return i, f"{stage} incomplete"
        
        # 情况4: 失败但 stat.log 未记录 → 从这里开始
        if stage_status == 'failed':
            return i, f"{stage} has error.txt"
        
        # 情况5: 成功或被忽略 → 继续检查下一个
        if stage_status in ['success', 'ignored']:
            continue
    
    # 所有阶段都完成了
    return None, "All stages completed"


def backup_logs(work_dir):
    """
    备份 stat.log 和 time.log
    自动递增编号: stat0.log, stat1.log, ...
    
    Args:
        work_dir: 工作目录路径 (Path对象)
    
    Returns:
        (stat_backup, time_backup) 元组，如果文件不存在则为 None
    """
    backups = []
    
    for log_name in ['stat.log', 'time.log']:
        log_file = work_dir / log_name
        
        if not log_file.exists():
            backups.append(None)
            continue
        
        # 查找可用的备份编号
        backup_num = 0
        base_name = log_name.replace('.log', '')
        
        while True:
            backup_name = work_dir / f'{base_name}{backup_num}.log'
            if not backup_name.exists():
                break
            backup_num += 1
        
        # 执行备份
        try:
            shutil.copy2(log_file, backup_name)
            print(f"[BACKUP] {log_name} -> {backup_name.name}")
            backups.append(backup_name.name)
        except Exception as e:
            print(f"[WARNING] Failed to backup {log_name}: {e}")
            backups.append(None)
    
    return tuple(backups)


def clean_stage_errors(work_dir, start_stage):
    """
    清理续算起点阶段的 error.txt
    
    Args:
        work_dir: 工作目录路径 (Path对象)
        start_stage: 续算起点阶段名称
    """
    stage_dir = work_dir / start_stage
    error_file = stage_dir / 'error.txt'
    
    if error_file.exists():
        try:
            error_file.unlink()
            print(f"[CLEAN] Removed {start_stage}/error.txt")
        except Exception as e:
            print(f"[WARNING] Failed to remove {start_stage}/error.txt: {e}")
        
        # 同时清理 stat.log 中该阶段的失败记录
        stat_log = work_dir / 'stat.log'
        if stat_log.exists():
            try:
                lines = stat_log.read_text().strip().split('\n')
                # 移除该阶段的 failed 记录
                new_lines = [line for line in lines 
                           if not (line.startswith(f'{start_stage} ') and 'failed' in line)]
                stat_log.write_text('\n'.join(new_lines) + '\n')
                print(f"[CLEAN] Removed '{start_stage} failed' from stat.log")
            except Exception as e:
                print(f"[WARNING] Failed to clean stat.log: {e}")


def detect_resume_number(stage_dir):
    """
    检测当前是第几次续算
    
    Args:
        stage_dir: 阶段目录路径 (Path对象)
    
    Returns:
        resume_num: 0表示首次续算，1表示第2次续算，以此类推
    """
    if not stage_dir.exists():
        return 0
    
    # 查找现有的 *_R* 文件
    max_resume = -1
    
    for f in stage_dir.glob('*_R*_*'):
        if f.is_file():
            # 文件名格式: INPUT_R0_0, KPT_R1_2 等
            match = re.match(r'.*_R(\d+)_', f.stem)
            if match:
                resume_num = int(match.group(1))
                max_resume = max(max_resume, resume_num)
    
    return max_resume + 1


def detect_previous_try_count(stage_dir):
    """
    检测之前失败时尝试了多少次
    
    Args:
        stage_dir: 阶段目录路径 (Path对象)
    
    Returns:
        try_count: 0表示首次尝试就失败，1表示尝试了2次，以此类推
    """
    if not stage_dir.exists():
        return 0
    
    # 查找备份文件: INPUT0, INPUT1, INPUT2...
    count = 0
    while (stage_dir / f'INPUT{count}').exists():
        count += 1
    
    return count


def should_clean_and_restart(stage_dir, stage_name):
    """
    判断是否应该清空重算
    
    Args:
        stage_dir: 阶段目录路径 (Path对象)
        stage_name: 阶段名称
    
    Returns:
        (should_clean, reason) 元组
    """
    if not stage_dir.exists():
        return False, "Directory does not exist"
    
    try_count = detect_previous_try_count(stage_dir)
    
    # 情况A：首次尝试就失败 → 清空重算
    if try_count == 0:
        return True, "First attempt failed, likely input parameter error"
    
    # 情况B：多次尝试未收敛 → 检查是否有可用的中间结果
    # 对于 Relax 类型的计算，检查是否有 STRU_ION_D
    if 'relax' in stage_name.lower():
        out_dirs = list(stage_dir.glob('OUT.*'))
        has_stru_ion = False
        for out_dir in out_dirs:
            if (out_dir / 'STRU_ION_D').exists():
                has_stru_ion = True
                break
        
        if not has_stru_ion:
            return True, f"Multiple attempts ({try_count}) but no STRU_ION_D found"
    
    return False, f"Continue from attempt {try_count}"


def verify_stage_completion(stage_dir, stage_name):
    """
    验证阶段是否真正完成（更严格的检查）
    
    Args:
        stage_dir: 阶段目录路径 (Path对象)
        stage_name: 阶段名称
    
    Returns:
        (is_complete, checks_dict) 元组
    """
    checks = {
        "converge_marker": (stage_dir / "converge.txt").exists(),
        "no_error": not (stage_dir / "error.txt").exists(),
        "has_output_dir": len(list(stage_dir.glob("OUT.*"))) > 0,
        "has_log": len(list(stage_dir.glob("OUT.*/running*.log"))) > 0 or 
                   (stage_dir / "running.log").exists()
    }
    
    # 阶段特定检查
    if "relax" in stage_name.lower():
        checks["has_stru_ion"] = any((d / "STRU_ION_D").exists() 
                                     for d in stage_dir.glob("OUT.*"))
    elif stage_name == "Scf":
        checks["has_charge"] = any((d / "SPIN1_CHG.cube").exists() 
                                   for d in stage_dir.glob("OUT.*"))
    
    is_complete = all(checks.values())
    return is_complete, checks


def format_status_icon(status):
    """
    根据状态返回合适的图标
    
    Args:
        status: 状态字符串
    
    Returns:
        图标字符串
    """
    icons = {
        'success': '✅',
        'ignored': '⚠️ ',
        'failed': '❌',
        'incomplete': '🔄',
        'not_started': '⭕'
    }
    return icons.get(status, '❓')


def get_running_jobs_from_slurm(batch_dir):
    """
    从 Slurm 队列中获取正在运行的作业列表
    
    Args:
        batch_dir: 批处理目录路径 (Path对象)
    
    Returns:
        set: 正在运行的作业名称集合（从 WorkDir 中提取）
    """
    running_jobs = set()
    
    try:
        # 尝试使用 squeue 命令（Slurm 标准命令）
        result = subprocess.run(
            ['squeue', '-u', subprocess.getoutput('whoami'), '-h', '-o', '%i %T %Z'],
            capture_output=True,
            text=True,
            timeout=5
        )
        
        if result.returncode == 0:
            # 解析输出: JOBID STATE WORK_DIR
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 3:
                    state = parts[1]
                    work_dir = parts[2]
                    
                    # 只关注 RUNNING, PENDING, CONFIGURING 状态的作业
                    if state in ['RUNNING', 'PENDING', 'CONFIGURING', 'R', 'PD', 'CF']:
                        # 从 WorkDir 中提取作业名称
                        work_path = Path(work_dir)
                        if batch_dir.resolve() in work_path.resolve().parents:
                            job_name = work_path.name
                            running_jobs.add(job_name)
        else:
            # squeue 失败，尝试 yhq 命令
            result = subprocess.run(
                ['yhq'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                # 解析 yhq 输出
                for line in result.stdout.strip().split('\n'):
                    if not line or 'JOBID' in line:
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 3:
                        job_name = parts[2]  # NAME 列
                        state = parts[4]     # ST 列
                        
                        # R = RUNNING, PD = PENDING
                        if state in ['R', 'PD', 'CF']:
                            running_jobs.add(job_name)
    
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
        # 如果命令失败，返回空集合（不影响主流程）
        pass
    
    return running_jobs


