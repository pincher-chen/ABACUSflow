#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import shutil
from datetime import datetime
from pathlib import Path
from config import WORKFLOW, CONDOR, get
from abacus.resume_utils import (
    backup_logs, 
    clean_stage_errors,
    detect_resume_number,
    detect_previous_try_count,
    should_clean_and_restart
)

class AbacusFlowManager:
    def __init__(self, work_dir):
        self.work_dir = Path(work_dir).absolute()
        self.flow_stages = list(WORKFLOW.keys())
        
    def generate_scripts(self, stru_path=None, suffix='*.vasp'):
        """生成提交脚本"""
        # 读取结构文件列表
        if stru_path:
            search_path = Path(stru_path)
        else:
            # 回退到配置文件或默认值
            conf_path = get('STRU', 'PATH', '')
            if conf_path:
                search_path = Path(conf_path)
            else:
                search_path = Path('InputPoscar')
                if not search_path.exists():
                    search_path = Path('.')
        
        # 确定后缀
        stru_suffix = get('STRU', 'SUFFIX', suffix)
                
        # 处理单个文件还是目录
        if search_path.is_file():
            structure_files = [search_path]
        else:
            structure_files = list(search_path.glob(stru_suffix))
        
        if not structure_files:
            print(f"No structure files found in {search_path} with suffix {stru_suffix}")
            return
            
        print(f"Found {len(structure_files)} structure files.")
        
        # 为每个结构文件生成提交脚本
        generated_scripts = []
        for stru_file in structure_files:
            script = self._create_job_script(stru_file)
            generated_scripts.append(script)
        
        return generated_scripts
            
    def _create_job_script(self, stru_file):
        """为单个结构创建完整的流程控制脚本"""
        job_name = stru_file.stem
        
        # 使用绝对路径，但保留符号链接（不解析到物理路径）
        job_dir = (self.work_dir / job_name).absolute()
        
        # 确保计算目录存在
        try:
            # 先确保父目录存在
            if not job_dir.parent.exists():
                job_dir.parent.mkdir(parents=True, exist_ok=True)
            
            # 再创建作业目录
            if not job_dir.exists():
                job_dir.mkdir(parents=False, exist_ok=True)
                
        except Exception as e:
            print(f"[ERROR] 创建目录失败 {job_dir}: {e}")
            print(f"[DEBUG] work_dir={self.work_dir}, exists={self.work_dir.exists()}")
            print(f"[DEBUG] job_name={job_name}")
            print(f"[DEBUG] job_dir.parent={job_dir.parent}, exists={job_dir.parent.exists()}")
            raise
        
        # 再次验证目录确实已创建
        if not job_dir.exists():
            error_msg = f"目录创建后仍不存在: {job_dir}"
            print(f"[ERROR] {error_msg}")
            # 尝试列出父目录内容
            if job_dir.parent.exists():
                print(f"[DEBUG] 父目录内容: {list(job_dir.parent.iterdir())[:10]}")
            raise RuntimeError(error_msg)
        
        # 复制结构文件
        stru_file_dest = job_dir / 'STRU.vasp'
        if not stru_file_dest.exists():
            try:
                shutil.copy2(stru_file, stru_file_dest)
            except Exception as e:
                print(f"[ERROR] 复制结构文件失败 {stru_file} -> {stru_file_dest}: {e}")
                print(f"[DEBUG] 源文件存在: {stru_file.exists()}")
                print(f"[DEBUG] 目标目录存在: {job_dir.exists()}")
                print(f"[DEBUG] 目标目录内容: {list(job_dir.iterdir()) if job_dir.exists() else 'N/A'}")
                raise
        
        # 脚本生成在计算目录内：work_cal/job_name/job_name.sh
        script_name = job_dir / f"{job_name}.sh"
        
        # 生成完整工作流脚本（使用统一的脚本生成方法）
        script_content = self._create_workflow_script(
            job_dir=job_dir,
            job_name=job_name,
            stages=self.flow_stages,  # 所有阶段
            is_resume=False
        )
        
        # 写入脚本文件
        with open(script_name, 'w') as f:
            f.write(script_content)
        
        print(f"Generated script: {script_name}")
        os.chmod(script_name, 0o755)
        return str(script_name)
    
    def generate_resume_script(self, work_dir, start_stage, clean_errors=True, no_backup=False):
        """
        生成续算脚本
        
        Args:
            work_dir: 工作目录 (如 work_cal/hmat_0/)
            start_stage: 从哪个阶段开始
            clean_errors: 是否清理起点阶段的 error.txt
            no_backup: 是否跳过日志备份
        
        Returns:
            生成的脚本路径
        """
        work_dir = Path(work_dir).absolute()
        
        # 1. 备份日志
        backup_info = (None, None)
        if not no_backup:
            backup_info = backup_logs(work_dir)
        
        # 2. 清理错误标记
        if clean_errors:
            clean_stage_errors(work_dir, start_stage)
        
        # 3. 确定续算阶段列表
        if start_stage not in self.flow_stages:
            raise ValueError(f"Invalid stage: {start_stage}")
        
        start_index = self.flow_stages.index(start_stage)
        resume_stages = self.flow_stages[start_index:]
        
        # 4. 生成脚本
        job_name = work_dir.name
        script_name = work_dir / f"{job_name}_resume.sh"
        
        script_content = self._create_workflow_script(
            job_dir=work_dir,
            job_name=job_name,
            stages=resume_stages,
            is_resume=True,
            backup_info=backup_info
        )
        
        with open(script_name, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_name, 0o755)
        print(f"Generated resume script: {script_name}")
        
        return str(script_name)
    
    def _create_workflow_script(self, job_dir, job_name, stages, is_resume=False, backup_info=None):
        """
        生成工作流脚本的核心逻辑（完整或部分）
        
        Args:
            job_dir: 作业目录
            job_name: 作业名称
            stages: 要包含的阶段列表
            is_resume: 是否为续算模式
            backup_info: (stat_backup, time_backup) 备份信息元组
        
        Returns:
            脚本内容字符串
        """
        # 获取路径（使用绝对路径，但保留符号链接，不解析到物理路径）
        # 优先从配置文件读取，否则使用默认路径
        abacus_py_default = (Path(__file__).parent.parent / 'abacus.py').absolute()
        abacus_py_conf = get('ENV', 'ABACUS_PY', None)
        abacus_py = Path(abacus_py_conf).absolute() if abacus_py_conf else abacus_py_default
        
        monitor_py = (Path(__file__).parent / 'monitor.py').absolute()
        stru_file_abs = job_dir / 'STRU.vasp'
        
        # 获取配置
        conda_path = get('ENV', 'CONDA_PATH', '/XYFS01/nscc-gz_pinchen_1/sf_install/miniconda3')
        conda_env = get('ENV', 'CONDA_ENV', 'dftflow')
        abacus_dir = get('ABACUS', 'ABACUS_DIR', '/XYFS01/nscc-gz_pinchen_1/sf_box/abacus-develop-LTSv3.10.0/bin')
        abacus_exe = get('ABACUS', 'ABACUS_EXE', 'abacus')
        abacus_bin = f"{abacus_dir}/{abacus_exe}"
        modules_str = get('MODULE', 'MODULES', 'intel/oneapi2023.2_noimpi mpi/mpich/4.1.2-icc-oneapi2023.2-ch4')
        
        # 生成脚本头部
        if is_resume:
            script_content = [
                "#!/bin/bash",
                f"# Resume script for {job_name}",
                f"# Generated at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                f"# Resuming from: {stages[0]}",
                f"# Stages to run: {', '.join(stages)}",
            ]
            if backup_info and backup_info[0]:
                script_content.append(f"# Previous logs backed up: {backup_info[0]}, {backup_info[1]}")
        else:
            script_content = [
                "#!/bin/bash",
                f"# Job script for {job_name}",
                f"# Generated by abacusflow",
                f"# Run this script from: {job_dir}",
            ]
        
        # 环境初始化
        script_content.extend([
            "",
            "# ===== 环境初始化 =====",
            "export OMP_NUM_THREADS=1",
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
            f"cd {job_dir} || exit 1",
            "",
        ])
        
        # 续算模式添加辅助函数
        if is_resume:
            script_content.extend(self._get_resume_helper_functions())
        
        # 初始化日志
        if is_resume:
            resume_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            script_content.extend([
                "echo '[RESUME] Starting workflow resume...'",
                f"echo '[RESUME] Resuming from stage: {stages[0]}'",
                "",
                "# 使用现有日志文件（已备份）",
                "stat_log=stat.log",
                "time_log=time.log",
                f"echo '' >> $stat_log",
                f"echo '# === Resume from {stages[0]} at {resume_time} ===' >> $time_log",
            ])
        else:
            script_content.extend([
                "echo '[...]TASK START!'",
                f"echo '[...]Structure: {stru_file_abs}'",
                'echo "[...]Working directory: $(pwd)"',
                "",
                f"if [ ! -f STRU.vasp ]; then",
                f"  cp {stru_file_abs} ./STRU.vasp",
                "fi",
                "",
                "# 初始化状态日志和时间统计",
                "stat_log=stat.log",
                "> $stat_log",
                "time_log=time.log",
                "> $time_log",
                "echo '# Stage statistics' > $time_log",
                "echo '# Format: Stage | Duration(s) | Nodes | Cores | Core-hours | Status' >> $time_log",
            ])
        
        script_content.extend([
            "WORKFLOW_START=$(date +%s)",
            ""
        ])
        
        # 遍历阶段
        for stage in stages:
            stage_lines = self._generate_stage_script(
                stage=stage,
                abacus_py=abacus_py,
                abacus_bin=abacus_bin,
                monitor_py=monitor_py,
                is_resume=is_resume
            )
            script_content.extend(stage_lines)
        
        # 脚本结束部分
        script_content.extend(self._get_workflow_summary_script(abacus_py, job_dir))
        
        return '\n'.join(script_content)
    
    def _get_resume_helper_functions(self):
        """返回续算模式需要的bash辅助函数"""
        return [
            "",
            "# ===== 续算辅助函数 =====",
            "",
            "detect_resume_number() {",
            "    local stage_dir=$1",
            "    local max_r=-1",
            "    ",
            "    for f in \"$stage_dir\"/*_R*_*; do",
            "        if [ -f \"$f\" ]; then",
            "            if [[ $(basename \"$f\") =~ _R([0-9]+)_ ]]; then",
            "                r_num=\"${BASH_REMATCH[1]}\"",
            "                if [ \"$r_num\" -gt \"$max_r\" ]; then",
            "                    max_r=$r_num",
            "                fi",
            "            fi",
            "        fi",
            "    done",
            "    ",
            "    echo $((max_r + 1))",
            "}",
            "",
            "detect_previous_try_count() {",
            "    local stage_dir=$1",
            "    local count=0",
            "    ",
            "    while [ -f \"$stage_dir/INPUT${count}\" ]; do",
            "        count=$((count + 1))",
            "    done",
            "    ",
            "    echo $count",
            "}",
            "",
            "should_clean_and_restart() {",
            "    local stage_dir=$1",
            "    local stage_name=$2",
            "    ",
            "    local try_count=$(detect_previous_try_count \"$stage_dir\")",
            "    ",
            "    # 首次尝试就失败 → 清空",
            "    if [ \"$try_count\" -eq 0 ]; then",
            "        echo \"clean|First attempt failed, input likely has errors\"",
            "        return",
            "    fi",
            "    ",
            "    # 多次尝试 → 检查是否有 STRU_ION_D",
            "    if [[ \"$stage_name\" =~ [Rr]elax ]]; then",
            "        local stru_ion_d=$(find \"$stage_dir\"/OUT.* -name \"STRU_ION_D\" 2>/dev/null | head -1)",
            "        if [ -z \"$stru_ion_d\" ]; then",
            "            echo \"clean|Multiple attempts but no STRU_ION_D found\"",
            "            return",
            "        fi",
            "    fi",
            "    ",
            "    echo \"continue|Continue from attempt $try_count\"",
            "}",
            "",
        ]
    
    def _generate_stage_script(self, stage, abacus_py, abacus_bin, monitor_py, is_resume=False):
        """生成单个阶段的脚本内容"""
        stage_config = WORKFLOW[stage]
        template = stage_config.get('template', stage)
        nodes = stage_config.get('node', 1)
        cores = stage_config.get('core', 64)
        try_num = stage_config.get('try_num', 2)
        ignore_error = stage_config.get('ignore_error', "False")
        partition = get('ALLOW', 'PARTITION', 'deimos')
        
        lines = [
            f"echo '[...]start {stage} task'",
            f"{stage.upper()}_START=$(date +%s)",
            f"echo \"[INFO] {stage} started at $(date '+%Y-%m-%d %H:%M:%S')\"",
        ]
        
        # 续算模式的目录处理
        if is_resume:
            lines.extend([
                "",
                f"# === 续算模式：智能处理 {stage} 目录 ===",
                f"if [ -d {stage} ]; then",
                f"  echo '[RESUME] {stage} directory exists, analyzing...'",
                "  ",
                "  # 检测续算次数",
                f"  RESUME_NUM=$(detect_resume_number {stage})",
                "  echo \"[RESUME] This is resume attempt #$RESUME_NUM\"",
                "  ",
                "  # 判断是否应该清空重算",
                f"  DECISION=$(should_clean_and_restart {stage} {stage})",
                "  ACTION=$(echo \"$DECISION\" | cut -d'|' -f1)",
                "  REASON=$(echo \"$DECISION\" | cut -d'|' -f2)",
                "  ",
                "  echo \"[RESUME] Decision: $ACTION\"",
                "  echo \"[RESUME] Reason: $REASON\"",
                "  ",
                "  if [ \"$ACTION\" = \"clean\" ]; then",
                f"    echo '[CLEAN] Removing {stage} directory for fresh start'",
                f"    rm -rf {stage}",
                f"    mkdir {stage} && cd {stage} || exit",
                "    RESUME_START_TRY=0",
                "  else",
                "    echo '[CONTINUE] Continuing from previous attempts'",
                f"    cd {stage}",
                "    ",
                "    # 检测之前尝试了多少次",
                "    PREV_TRY_COUNT=$(detect_previous_try_count .)",
                "    echo \"[CONTINUE] Previous attempts: $PREV_TRY_COUNT\"",
                "    ",
                "    # 续算从下一次尝试开始",
                "    RESUME_START_TRY=$PREV_TRY_COUNT",
                "    ",
                "    # 备份当前文件为 *_R${RESUME_NUM}_pre",
                "    for file in INPUT KPT STRU; do",
                "      if [ -f \"$file\" ]; then",
                "        cp \"$file\" \"${file}_R${RESUME_NUM}_pre\"",
                "        echo \"[BACKUP] $file -> ${file}_R${RESUME_NUM}_pre\"",
                "      fi",
                "    done",
                "    ",
                "    # 如果有 STRU_ION_D，使用它作为新的 STRU",
                "    STRU_ION_D=$(find OUT.* -name \"STRU_ION_D\" 2>/dev/null | head -1)",
                "    if [ -n \"$STRU_ION_D\" ]; then",
                "      echo \"[CONTINUE] Using optimized structure from $STRU_ION_D\"",
                "      cp \"$STRU_ION_D\" STRU",
                "    fi",
                "  fi",
                "else",
                f"  mkdir {stage} && cd {stage} || exit",
                "  RESUME_START_TRY=0",
                "fi",
            ])
        else:
            lines.extend([
                f"if [ ! -d {stage} ];then",
                f"  mkdir {stage} && cd {stage} || exit",
                "else",
                f"  cd {stage}",
                "fi",
            ])
        
        lines.extend([
            "",
            f"echo '[...]prepare {stage} inputs.'",
            f"python {abacus_py} generate --work_dir . --stage {template}",
        ])
        
        if ignore_error == "True":
            lines.extend([
                "",
                f"# Create ignore.txt marker for {stage}",
                "touch ignore.txt",
            ])
        
        # 循环迭代部分
        if is_resume:
            lines.extend([
                "",
                "# 续算模式：从 RESUME_START_TRY 开始",
                f"END_TRY=$((RESUME_START_TRY + {try_num}))",
                "for try_num in $(seq $RESUME_START_TRY $((END_TRY - 1)))",
            ])
        else:
            lines.extend([
                "",
                f"for try_num in $(seq 0 {try_num-1})",
            ])
        
        lines.extend([
            "  do",
            f"  echo \"[...]task {stage} round: $try_num on {nodes} node {cores} core\"",
            "  ",
            "  # 运行 ABACUS 并实时监控",
            f"  if command -v stdbuf >/dev/null 2>&1 && [ -f {monitor_py} ]; then",
            f"    stdbuf -oL -eL yhrun -N {nodes} -n {cores} -p {partition} {abacus_bin} 2>&1 | tee running.log | python -u {monitor_py}",
            "  else",
            f"    yhrun -N {nodes} -n {cores} -p {partition} {abacus_bin} > running.log 2>&1",
            "  fi",
            "  ",
            "  # 检查计算结果",
            "  echo \"[...]checking iteration $try_num result...\"",
            f"  python {abacus_py} errors --work_dir .",
            f"  python {abacus_py} converge --work_dir .",
            "  ",
            "  # 如果收敛成功，退出循环",
            "  if [ -f \"converge.txt\" ]; then",
            "    echo \"[...]iteration $try_num: converged successfully!\"",
        ])
        
        if 'spin' in stage.lower() or 'spin' in template.lower():
            lines.append(f"    python {abacus_py} spin --work_dir .")
        
        lines.extend([
            "    break",
            "  fi",
            "  ",
        ])
        
        # 备份逻辑
        if is_resume:
            lines.extend([
                "  # 未收敛且不是最后一次，准备下一次迭代",
                "  if [ $try_num -lt $((END_TRY - 1)) ]; then",
                "    echo '[...]calculation not done, prepare to next loop'",
                "    ",
                "    # 续算模式：使用新的备份命名",
                "    for file in INPUT KPT STRU; do",
                "      if [ -f \"$file\" ]; then",
                "        cp \"$file\" \"${file}_R${RESUME_NUM}_${try_num}\"",
                "        echo \"[BACKUP] $file -> ${file}_R${RESUME_NUM}_${try_num}\"",
                "      fi",
                "    done",
                "    ",
                f"    python {abacus_py} update --work_dir . --try_num $try_num --stage {stage}",
            ])
        else:
            lines.extend([
                f"  if [ $try_num -lt {try_num-1} ]; then",
                "    echo '[...]calculation not done, prepare to next loop'",
                f"    python {abacus_py} update --work_dir . --try_num $try_num --stage {stage}",
            ])
        
        lines.extend([
            "  else",
            "    echo '[...]last iteration completed'",
            "  fi",
            "done",
            "",
            "if [ ! -f 'converge.txt' ]; then",
            "  echo '[...]Job failed or not converged'",
            "  if [ -f 'error.txt' ]; then",
            "    echo '[...]Error detected (error.txt exists), cannot continue even with ignore.txt'",
            f"    echo '{stage} failed' >> ../stat.log",
            "    exit 1",
            "  elif [ ! -f 'ignore.txt' ]; then",
            "    echo '[...]Not converged and no ignore marker, exiting...'",
            f"    echo '{stage} failed' >> ../stat.log",
            "    exit 1",
            "  else",
            "    echo '[...]Not converged but can be ignored (ignore.txt exists), continuing...'",
        ])
        
        if 'spin' in stage.lower() or 'spin' in template.lower():
            lines.append(f"    python {abacus_py} spin --work_dir .")
        
        lines.extend([
            f"    echo '{stage} success (ignored)' >> ../stat.log",
            "  fi",
            "else",
            f"  echo '[...]{stage} job done!'",
        ])
        
        if 'spin' in stage.lower() or 'spin' in template.lower():
            lines.append(f"  python {abacus_py} spin --work_dir .")
        
        lines.extend([
            f"  echo '{stage} success' >> ../stat.log",
            "fi",
            "",
            f"# 计算 {stage} 阶段用时和资源消耗",
            f"{stage.upper()}_END=$(date +%s)",
            f"{stage.upper()}_DURATION=$(( ${stage.upper()}_END - ${stage.upper()}_START ))",
            f"{stage.upper()}_HOURS=$(awk \"BEGIN {{printf \\\"%.4f\\\", ${stage.upper()}_DURATION/3600}}\")",
            f"{stage.upper()}_CORE_HOURS=$(awk \"BEGIN {{printf \\\"%.2f\\\", {nodes}*{cores}*${stage.upper()}_HOURS}}\")",
            f"echo \"[INFO] {stage} completed in ${{{stage.upper()}_DURATION}}s (${{{stage.upper()}_HOURS}}h)\"",
            f"echo \"[INFO] {stage} consumed ${{{stage.upper()}_CORE_HOURS}} Core-hours ({nodes}nodes * {cores}cores * ${{{stage.upper()}_HOURS}}h)\"",
            f"echo \"{stage}|${{{stage.upper()}_DURATION}}|{nodes}|{cores}|${{{stage.upper()}_CORE_HOURS}}|$(tail -1 ../stat.log 2>/dev/null || echo 'unknown')\" >> ../time.log",
            "",
            "cd ..",
            ""
        ])
        
        return lines
    
    def _get_workflow_summary_script(self, abacus_py, job_dir):
        """返回工作流总结部分的脚本"""
        return [
            "",
            "# ===== 计算总用时和资源消耗 =====",
            "WORKFLOW_END=$(date +%s)",
            "TOTAL_DURATION=$((WORKFLOW_END - WORKFLOW_START))",
            "TOTAL_HOURS=$(awk \"BEGIN {printf \\\"%.4f\\\", $TOTAL_DURATION/3600}\")",
            "",
            "echo ''",
            "echo '======================================================================'",
            "echo '                    WORKFLOW SUMMARY'",
            "echo '======================================================================'",
            "echo ''",
            "echo '[INFO] Individual Stage Statistics:'",
            "echo '----------------------------------------------------------------------'",
            "printf '%-15s %-12s %-8s %-8s %-15s %-15s\\n' 'Stage' 'Duration(s)' 'Nodes' 'Cores' 'Core-hours' 'Status'",
            "echo '----------------------------------------------------------------------'",
            "tail -n +3 time.log | while IFS='|' read -r stage duration nodes cores core_hours status; do",
            "  hours=$(awk \"BEGIN {printf \\\"%.2f\\\", $duration/3600}\")",
            "  printf '%-15s %-12s %-8s %-8s %-15s %-15s\\n' \"$stage\" \"${duration}s (${hours}h)\" \"$nodes\" \"$cores\" \"$core_hours\" \"$status\"",
            "done",
            "echo '----------------------------------------------------------------------'",
            "",
            "TOTAL_CORE_HOURS=$(tail -n +3 time.log | awk -F'|' '{sum+=$5} END {printf \"%.2f\", sum}')",
            "",
            "echo ''",
            "echo '[INFO] Total Workflow Statistics:'",
            "echo \"  Total Duration    : ${TOTAL_DURATION}s (${TOTAL_HOURS}h)\"",
            "echo \"  Total Core-hours  : ${TOTAL_CORE_HOURS}\"",
            "echo \"  Workflow Start    : $(date -d @$WORKFLOW_START '+%Y-%m-%d %H:%M:%S')\"",
            "echo \"  Workflow End      : $(date -d @$WORKFLOW_END '+%Y-%m-%d %H:%M:%S')\"",
            "echo ''",
            "echo '======================================================================'",
            "",
            f"python {abacus_py} summary --root {job_dir}",
            "",
            "echo '[...]ALL STAGES COMPLETED!'",
            "echo 'Workflow completed' >> $stat_log",
            ""
        ]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ABACUS Flow Manager')
    parser.add_argument('action', choices=['init', 'submit', 'status'], help='Action to perform')
    args = parser.parse_args()
    
    # 默认工作目录为 work_cal
    manager = AbacusFlowManager('work_cal')
    
    if args.action == 'init':
        manager.generate_scripts()
