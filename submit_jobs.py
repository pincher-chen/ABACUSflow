#!/usr/bin/env python3
"""
自动提交作业到Slurm系统

使用方法:
    # 提交新作业（workflow模式生成的脚本）
    python submit_jobs.py work_cal7
    
    # 提交续算作业
    python submit_jobs.py batch_work/ --resume
    
    # 提交单步计算作业（single模式生成的脚本）
    python submit_jobs.py batch_work_single/ --single

该脚本会：
1. 扫描目录下的所有.sh脚本（或 *_resume.sh / *_single.sh）
2. 解析脚本中的节点数配置
3. 根据config/condor.ini的配置自动提交作业
4. 监控队列状态，避免超出资源限制
5. 间隔提交，降低调度系统负载
"""

import sys
import os
from pathlib import Path
from abacus.submit_manager import AutoSubmitter

# 获取 abacusflow 根目录（submit_jobs.py 所在目录）
ABACUSFLOW_ROOT = Path(__file__).parent.absolute()
DEFAULT_CONFIG = str(ABACUSFLOW_ROOT / "config" / "condor.ini")


def print_help():
    """打印帮助信息"""
    print("使用方法: python submit_jobs.py <work_directory> [--resume | --single]")
    print()
    print("选项:")
    print("  --resume    提交续算作业（*_resume.sh）")
    print("  --single    提交单步计算作业（*_single.sh）")
    print("  --help, -h  显示此帮助信息")
    print()
    print("示例:")
    print("  # 提交新作业（workflow模式）")
    print("  python submit_jobs.py work_cal7")
    print()
    print("  # 提交续算作业（*_resume.sh）")
    print("  python submit_jobs.py batch_work/ --resume")
    print()
    print("  # 提交单步计算作业（*_single.sh）")
    print("  python submit_jobs.py batch_work_single/ --single")
    print()
    print("配置文件: config/condor.ini")
    print("  - TOTAL_NODE: 最大使用节点数")
    print("  - MAX_JOBS: 最大队列任务数")
    print("  - INTERVAL_TIME: 提交间隔（秒）")
    print("  - CHECK_TIME: 检查队列间隔（秒）")


def main():
    # 检查帮助参数
    if '--help' in sys.argv or '-h' in sys.argv:
        print_help()
        sys.exit(0)
    
    if len(sys.argv) < 2 or sys.argv[1].startswith('-'):
        print_help()
        sys.exit(1)
    
    work_dir = sys.argv[1]
    
    # 检查模式
    include_resume = '--resume' in sys.argv
    include_single = '--single' in sys.argv
    
    # 创建自动提交器（使用绝对路径）
    submitter = AutoSubmitter(config_file=DEFAULT_CONFIG)
    
    # 运行提交流程
    try:
        submitter.run(work_dir, max_retries=3, include_resume=include_resume, include_single=include_single)
    except KeyboardInterrupt:
        print("\n\n用户中断，正在停止...")
        submitter.stop()


if __name__ == '__main__':
    main()

