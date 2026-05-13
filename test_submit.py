#!/usr/bin/env python3
"""
测试自动提交系统

这个脚本用于测试提交系统的各个组件
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from abacus.submit_manager import (
    ConfigManager,
    SlurmMonitor,
    ScriptParser,
    JobQueue,
    JobInfo,
    AutoSubmitter
)
import logging

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('TestSubmit')


def test_config_manager():
    """测试配置管理器"""
    print("=" * 70)
    print("测试 ConfigManager")
    print("=" * 70)
    
    config = ConfigManager("config/condor.ini")
    
    print(f"分区: {config.partition}")
    print(f"总节点数: {config.total_nodes}")
    print(f"最大任务数: {config.max_jobs}")
    print(f"提交间隔: {config.interval_time}秒")
    print(f"检查间隔: {config.check_time}秒")
    print(f"默认节点数: {config.default_nodes}")
    print(f"每节点核心数: {config.cores_per_node}")
    print()


def test_slurm_monitor():
    """测试Slurm监控器"""
    print("=" * 70)
    print("测试 SlurmMonitor")
    print("=" * 70)
    
    config = ConfigManager("config/condor.ini")
    monitor = SlurmMonitor(config, logger)
    
    total_jobs, running_jobs, used_nodes = monitor.get_queue_status()
    
    print(f"总任务数: {total_jobs}")
    print(f"运行中任务: {running_jobs}")
    print(f"使用节点数: {used_nodes}")
    print()
    
    # 测试是否可以提交
    test_job = JobInfo(
        script_path="test.sh",
        job_dir=".",
        job_name="test_job",
        nodes=1,
        cores=16,
        partition=config.partition
    )
    
    can_submit, reason = monitor.can_submit(test_job)
    print(f"是否可以提交测试作业: {can_submit}")
    print(f"原因: {reason}")
    print()


def test_script_parser():
    """测试脚本解析器"""
    print("=" * 70)
    print("测试 ScriptParser")
    print("=" * 70)
    
    # 测试解析实际脚本
    test_script = "work_cal7/hmat_0/hmat_0.sh"
    
    if os.path.exists(test_script):
        nodes, cores = ScriptParser.parse_script(test_script)
        job_name = ScriptParser.get_job_name(test_script)
        
        print(f"脚本: {test_script}")
        print(f"作业名: {job_name}")
        print(f"节点数: {nodes}")
        print(f"核心数: {cores}")
    else:
        print(f"测试脚本不存在: {test_script}")
    print()


def test_job_queue():
    """测试作业队列"""
    print("=" * 70)
    print("测试 JobQueue")
    print("=" * 70)
    
    queue = JobQueue(logger)
    
    # 添加测试作业
    for i in range(5):
        job = JobInfo(
            script_path=f"test_{i}.sh",
            job_dir=".",
            job_name=f"test_job_{i}",
            nodes=1,
            cores=16,
            partition="deimos"
        )
        queue.add_job(job)
    
    print(f"待提交作业数: {queue.pending_count}")
    
    # 提取作业
    job = queue.get_next_job()
    if job:
        print(f"提取作业: {job.job_name}")
        queue.mark_submitted(job, "12345")
    
    print(f"已提交作业数: {queue.submitted_count}")
    print(f"待提交作业数: {queue.pending_count}")
    
    # 标记失败
    job = queue.get_next_job()
    if job:
        queue.mark_failed(job, "测试失败")
    
    print(f"失败作业数: {queue.failed_count}")
    print()


def test_scan_scripts():
    """测试脚本扫描"""
    print("=" * 70)
    print("测试脚本扫描")
    print("=" * 70)
    
    submitter = AutoSubmitter("config/condor.ini")
    
    work_dir = "work_cal7"
    if os.path.exists(work_dir):
        jobs = submitter.scan_scripts(work_dir)
        print(f"找到 {len(jobs)} 个作业")
        
        if jobs:
            print("\n前5个作业:")
            for job in jobs[:5]:
                print(f"  - {job.job_name}: 节点={job.nodes}, 核心={job.cores}")
    else:
        print(f"工作目录不存在: {work_dir}")
    print()


def main():
    """运行所有测试"""
    print("\n")
    print("=" * 70)
    print(" 自动提交系统测试")
    print("=" * 70)
    print("\n")
    
    try:
        test_config_manager()
        test_slurm_monitor()
        test_script_parser()
        test_job_queue()
        test_scan_scripts()
        
        print("=" * 70)
        print("所有测试完成")
        print("=" * 70)
        
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()


