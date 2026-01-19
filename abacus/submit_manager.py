#!/usr/bin/env python3
"""
Slurm Auto Submission Manager
自动提交管理器，支持大规模作业批量提交到Slurm系统

功能特性：
1. 自动扫描生成的作业脚本
2. 内存队列管理，支持10w+作业
3. 智能资源监控，避免超出节点/任务限制
4. 间隔提交，降低调度器负载
5. 自动重试和错误处理
"""

import os
import re
import time
import subprocess
import configparser
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from collections import deque
import logging


@dataclass
class JobInfo:
    """作业信息"""
    script_path: str        # 脚本路径
    job_dir: str           # 作业目录
    job_name: str          # 作业名称
    nodes: int             # 节点数
    cores: int             # 核心数
    partition: str         # 分区
    job_id: Optional[str] = None  # 提交后的作业ID
    status: str = 'pending'       # 状态: pending, submitted, failed
    retry_count: int = 0          # 重试次数


class ConfigManager:
    """配置管理器"""
    
    def __init__(self, config_file: str = "config/condor.ini"):
        self.config_file = config_file
        self.config = configparser.ConfigParser()
        self._last_mtime = 0
        self.load_config()
    
    def load_config(self):
        """加载或重新加载配置文件"""
        try:
            self.config.read(self.config_file)
            self._last_mtime = os.path.getmtime(self.config_file)
        except Exception as e:
            logging.warning(f"无法读取配置文件 {self.config_file}: {e}")
    
    def check_and_reload(self) -> bool:
        """
        检查配置文件是否有更新，如果有则重新加载
        
        Returns:
            是否重新加载了配置
        """
        try:
            current_mtime = os.path.getmtime(self.config_file)
            if current_mtime > self._last_mtime:
                old_config = {
                    'total_node': self.total_nodes,
                    'max_jobs': self.max_jobs,
                    'interval_time': self.interval_time,
                    'check_time': self.check_time
                }
                
                self.load_config()
                
                new_config = {
                    'total_node': self.total_nodes,
                    'max_jobs': self.max_jobs,
                    'interval_time': self.interval_time,
                    'check_time': self.check_time
                }
                
                # 检查是否有实质性变化
                if old_config != new_config:
                    logging.info("检测到配置文件更新，重新加载配置:")
                    for key in old_config:
                        if old_config[key] != new_config[key]:
                            logging.info(f"  {key}: {old_config[key]} -> {new_config[key]}")
                    return True
        except Exception as e:
            logging.warning(f"检查配置文件更新时出错: {e}")
        
        return False
        
    def get(self, section: str, key: str, default=None):
        """获取配置项"""
        try:
            return self.config.get(section, key)
        except:
            return default
    
    def getint(self, section: str, key: str, default: int = 0):
        """获取整数配置"""
        try:
            return self.config.getint(section, key)
        except:
            return default
    
    def getfloat(self, section: str, key: str, default: float = 0.0):
        """获取浮点数配置"""
        try:
            return self.config.getfloat(section, key)
        except:
            return default
    
    @property
    def partition(self) -> str:
        """分区名称"""
        return self.get('ALLOW', 'PARTITION', 'default')
    
    @property
    def total_nodes(self) -> int:
        """总节点数限制"""
        return self.getint('ALLOW', 'TOTAL_NODE', 10)
    
    @property
    def max_jobs(self) -> int:
        """最大任务数限制"""
        return self.getint('ALLOW', 'MAX_JOBS', 1000)
    
    @property
    def interval_time(self) -> float:
        """提交间隔时间(秒)"""
        return self.getfloat('ALLOW', 'INTERVAL_TIME', 0.5)
    
    @property
    def check_time(self) -> float:
        """检查队列状态时间(秒)"""
        return self.getfloat('ALLOW', 'CHECK_TIME', 60)
    
    @property
    def default_nodes(self) -> int:
        """默认节点数"""
        return self.getint('ALLOW', 'NODES', 1)
    
    @property
    def cores_per_node(self) -> int:
        """每节点核心数"""
        return self.getint('ALLOW', 'CORES_PER_NODE', 64)
    
    @property
    def config_reload_interval(self) -> float:
        """配置重载检查间隔(秒)"""
        return self.getfloat('ALLOW', 'CONFIG_RELOAD_INTERVAL', 600)


class SlurmMonitor:
    """Slurm队列监控器"""
    
    def __init__(self, config: ConfigManager, logger: logging.Logger):
        self.config = config
        self.logger = logger
        
    def get_queue_status(self) -> Tuple[int, int, int]:
        """
        获取当前队列状态
        
        Returns:
            (total_jobs, running_jobs, used_nodes)
        """
        try:
            # 使用squeue获取当前用户的作业
            # 只使用 partition 的第一个部分（如果包含额外参数则忽略）
            partition_name = self.config.partition.split()[0]
            cmd = f"squeue -u $USER -h -o '%T %D' -p {partition_name}"
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                self.logger.warning(f"squeue command failed: {result.stderr}")
                return 0, 0, 0
            
            total_jobs = 0
            running_jobs = 0
            used_nodes = 0
            
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                    
                parts = line.split()
                if len(parts) >= 2:
                    status = parts[0]
                    nodes = int(parts[1])
                    
                    total_jobs += 1
                    if status in ['RUNNING', 'R']:
                        running_jobs += 1
                        used_nodes += nodes
            
            self.logger.debug(
                f"Queue status: total={total_jobs}, "
                f"running={running_jobs}, nodes={used_nodes}"
            )
            
            return total_jobs, running_jobs, used_nodes
            
        except subprocess.TimeoutExpired:
            self.logger.error("squeue command timeout")
            return 0, 0, 0
        except Exception as e:
            self.logger.error(f"Error getting queue status: {e}")
            return 0, 0, 0
    
    def can_submit(self, job: JobInfo) -> Tuple[bool, str]:
        """
        检查是否可以提交作业
        
        Args:
            job: 作业信息
            
        Returns:
            (can_submit, reason)
        """
        total_jobs, running_jobs, used_nodes = self.get_queue_status()
        
        # 检查任务数限制
        if total_jobs >= self.config.max_jobs:
            return False, f"达到最大任务数限制 {self.config.max_jobs}"
        
        # 检查节点数限制
        if used_nodes + job.nodes > self.config.total_nodes:
            return False, f"节点资源不足 (当前:{used_nodes}, 需要:{job.nodes}, 限制:{self.config.total_nodes})"
        
        return True, "资源充足"


class ScriptParser:
    """脚本解析器，从脚本中提取节点和核心信息"""
    
    @staticmethod
    def parse_script(script_path: str) -> Tuple[int, int]:
        """
        解析脚本获取节点数和核心数
        
        优先从SBATCH指令中提取，如果没有则从yhrun/srun/mpirun命令中推断
        
        Returns:
            (nodes, cores)
        """
        nodes = 1  # 默认值
        cores = 16  # 默认值
        
        try:
            with open(script_path, 'r') as f:
                content = f.read()
                
                # 尝试从SBATCH指令中提取
                # #SBATCH -N 2
                # #SBATCH --nodes=2
                node_match = re.search(r'#SBATCH\s+(?:-N|--nodes=?)(\d+)', content)
                if node_match:
                    nodes = int(node_match.group(1))
                
                # #SBATCH -n 32
                # #SBATCH --ntasks=32
                task_match = re.search(r'#SBATCH\s+(?:-n|--ntasks=?)(\d+)', content)
                if task_match:
                    cores = int(task_match.group(1))
                
                # 尝试从mpirun命令中提取
                # mpirun -np 16
                mpirun_match = re.search(r'mpirun\s+-np\s+(\d+)', content)
                if mpirun_match:
                    cores = int(mpirun_match.group(1))
                
                # 尝试从yhrun/srun命令中提取
                # yhrun -N 1 -n 32 或 srun -N 1 -n 32
                yhrun_match = re.search(r'(?:yhrun|srun)\s+.*?-N\s+(\d+)\s+.*?-n\s+(\d+)', content)
                if yhrun_match:
                    nodes = int(yhrun_match.group(1))
                    cores = int(yhrun_match.group(2))
                
                # 尝试从脚本注释中提取（single脚本格式）
                # # Resources: 1 node(s), 32 processes, OMP=1
                resource_match = re.search(r'#\s*Resources:\s*(\d+)\s*node\(s\),\s*(\d+)\s*processes', content)
                if resource_match:
                    nodes = int(resource_match.group(1))
                    cores = int(resource_match.group(2))
                
        except Exception as e:
            logging.warning(f"解析脚本失败 {script_path}: {e}")
        
        return nodes, cores
    
    @staticmethod
    def get_job_name(script_path: str) -> str:
        """从脚本路径提取作业名称"""
        path = Path(script_path)
        # 使用父目录名称作为作业名称
        return path.parent.name


class JobQueue:
    """作业队列管理器"""
    
    def __init__(self, logger: logging.Logger):
        self.queue = deque()
        self.submitted_jobs = {}  # job_id -> JobInfo
        self.failed_jobs = []
        self.logger = logger
        
    def add_job(self, job: JobInfo):
        """添加作业到队列"""
        self.queue.append(job)
        
    def add_jobs(self, jobs: List[JobInfo]):
        """批量添加作业"""
        self.queue.extend(jobs)
        
    def get_next_job(self) -> Optional[JobInfo]:
        """获取下一个待提交作业"""
        if self.queue:
            return self.queue.popleft()
        return None
    
    def mark_submitted(self, job: JobInfo, job_id: str):
        """标记作业已提交"""
        job.job_id = job_id
        job.status = 'submitted'
        self.submitted_jobs[job_id] = job
        self.logger.info(f"作业已提交: {job.job_name} (ID: {job_id})")
        
    def mark_failed(self, job: JobInfo, reason: str):
        """标记作业提交失败"""
        job.status = 'failed'
        job.retry_count += 1
        self.failed_jobs.append((job, reason))
        self.logger.error(f"作业提交失败: {job.job_name}, 原因: {reason}")
        
    @property
    def pending_count(self) -> int:
        """待提交作业数"""
        return len(self.queue)
    
    @property
    def submitted_count(self) -> int:
        """已提交作业数"""
        return len(self.submitted_jobs)
    
    @property
    def failed_count(self) -> int:
        """失败作业数"""
        return len(self.failed_jobs)


class AutoSubmitter:
    """自动提交器"""
    
    def __init__(self, config_file: str = "config/condor.ini"):
        self.config = ConfigManager(config_file)
        self.logger = self._setup_logger()
        self.monitor = SlurmMonitor(self.config, self.logger)
        self.queue = JobQueue(self.logger)
        self.running = False
        
    def _setup_logger(self) -> logging.Logger:
        """设置日志"""
        logger = logging.getLogger('AutoSubmitter')
        logger.setLevel(logging.INFO)
        
        # 控制台处理器
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # 文件处理器
        log_dir = Path('logs')
        log_dir.mkdir(exist_ok=True)
        file_handler = logging.FileHandler(
            log_dir / f'submit_{time.strftime("%Y%m%d_%H%M%S")}.log'
        )
        file_handler.setLevel(logging.DEBUG)
        
        # 格式化
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        console_handler.setFormatter(formatter)
        file_handler.setFormatter(formatter)
        
        logger.addHandler(console_handler)
        logger.addHandler(file_handler)
        
        return logger
    
    def scan_scripts(self, work_dir: str, pattern: str = "*.sh", include_resume: bool = False, include_single: bool = False) -> List[JobInfo]:
        """
        扫描工作目录中的作业脚本
        
        Args:
            work_dir: 工作目录
            pattern: 脚本文件模式
            include_resume: 是否包含续算脚本 (*_resume.sh)
            include_single: 是否包含单步计算脚本 (*_single.sh)
            
        Returns:
            作业信息列表
        """
        self.logger.info(f"扫描作业脚本: {work_dir}")
        if include_resume:
            self.logger.info("  模式: 续算脚本 (*_resume.sh)")
        elif include_single:
            self.logger.info("  模式: 单步计算脚本 (*_single.sh)")
        else:
            self.logger.info("  模式: 工作流脚本 (排除 *_resume.sh 和 *_single.sh)")
        
        jobs = []
        work_path = Path(work_dir)
        
        if not work_path.exists():
            self.logger.error(f"工作目录不存在: {work_dir}")
            return jobs
        
        # 递归查找所有脚本文件
        script_files = list(work_path.rglob(pattern))
        
        # 过滤脚本
        if include_resume:
            # 只包含续算脚本
            script_files = [
                f for f in script_files 
                if f.name.endswith('_resume.sh')
            ]
        elif include_single:
            # 只包含单步计算脚本
            script_files = [
                f for f in script_files 
                if f.name.endswith('_single.sh')
            ]
        else:
            # 排除续算脚本和单步计算脚本（默认行为，用于workflow）
            script_files = [
                f for f in script_files 
                if not f.name.endswith('_resume.sh') and not f.name.endswith('_single.sh')
            ]
        
        self.logger.info(f"找到 {len(script_files)} 个脚本文件")
        
        for script_file in script_files:
            try:
                # 解析脚本信息
                nodes, cores = ScriptParser.parse_script(str(script_file))
                job_name = ScriptParser.get_job_name(str(script_file))
                
                job = JobInfo(
                    script_path=str(script_file),
                    job_dir=str(script_file.parent),
                    job_name=job_name,
                    nodes=nodes,
                    cores=cores,
                    partition=self.config.partition
                )
                
                jobs.append(job)
                
            except Exception as e:
                self.logger.error(f"解析脚本失败 {script_file}: {e}")
        
        return jobs
    
    def submit_job(self, job: JobInfo) -> Tuple[bool, str]:
        """
        提交单个作业
        
        Args:
            job: 作业信息
            
        Returns:
            (success, job_id_or_error)
        """
        try:
            # 确保使用绝对路径（但不解析符号链接）
            script_path = Path(job.script_path).absolute()
            job_dir = Path(job.job_dir).absolute()
            
            # 构建sbatch命令
            # 使用脚本中解析出的实际核心数，而不是独占整个节点
            # 这样可以更高效地利用集群资源
            cmd = [
                'sbatch',
                '-N', str(job.nodes),
                '-n', str(job.cores),  # 使用脚本中实际需要的核心数
            ]
            
            # 处理 partition 参数（可能包含额外选项如 --reservation）
            # 支持两种格式：
            # 1. PARTITION = deimos
            # 2. PARTITION = allsys --reservation=hpc_test
            partition_parts = job.partition.split()
            if partition_parts:
                cmd.extend(['-p', partition_parts[0]])
                # 添加额外的参数（如 --reservation）
                if len(partition_parts) > 1:
                    cmd.extend(partition_parts[1:])
            
            cmd.extend([
                '--job-name', job.job_name,
                str(script_path)
            ])
            
            cmd_str = ' '.join(cmd)
            self.logger.debug(f"提交命令: {cmd_str}")
            self.logger.debug(f"工作目录: {job_dir}")
            
            # 执行提交
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30,
                cwd=str(job_dir)
            )
            
            if result.returncode != 0:
                error_msg = result.stderr.strip()
                self.logger.error(f"作业 {job.job_name} 提交失败: returncode={result.returncode}")
                self.logger.error(f"  命令: {cmd_str}")
                self.logger.error(f"  标准输出: {result.stdout.strip()}")
                self.logger.error(f"  标准错误: {error_msg}")
                return False, error_msg
            
            # 提取作业ID
            # 输出格式: "Submitted batch job 12345"
            output = result.stdout.strip()
            match = re.search(r'Submitted batch job (\d+)', output)
            
            if match:
                job_id = match.group(1)
                self.logger.debug(f"作业 {job.job_name} 提交成功: ID={job_id}")
                return True, job_id
            else:
                error_msg = f"无法解析作业ID: {output}"
                self.logger.error(f"作业 {job.job_name}: {error_msg}")
                return False, error_msg
                
        except subprocess.TimeoutExpired:
            error_msg = "提交命令超时"
            self.logger.error(f"作业 {job.job_name}: {error_msg}")
            return False, error_msg
        except Exception as e:
            error_msg = str(e)
            self.logger.error(f"作业 {job.job_name} 提交异常: {error_msg}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return False, error_msg
    
    def run(self, work_dir: str, max_retries: int = 3, include_resume: bool = False, include_single: bool = False):
        """
        运行自动提交流程
        
        Args:
            work_dir: 工作目录
            max_retries: 最大重试次数
            include_resume: 是否提交续算脚本
            include_single: 是否提交单步计算脚本
        """
        self.running = True
        
        self.logger.info("=" * 70)
        self.logger.info("Slurm 自动提交器启动")
        self.logger.info("=" * 70)
        self.logger.info(f"工作目录: {work_dir}")
        self.logger.info(f"分区: {self.config.partition}")
        self.logger.info(f"节点限制: {self.config.total_nodes}")
        self.logger.info(f"任务限制: {self.config.max_jobs}")
        self.logger.info(f"提交间隔: {self.config.interval_time}秒")
        self.logger.info(f"检查间隔: {self.config.check_time}秒")
        if include_resume:
            self.logger.info("模式: 续算脚本提交 (*_resume.sh)")
        elif include_single:
            self.logger.info("模式: 单步计算脚本提交 (*_single.sh)")
        else:
            self.logger.info("模式: 工作流脚本提交")
        self.logger.info("=" * 70)
        
        # 扫描作业脚本
        jobs = self.scan_scripts(work_dir, include_resume=include_resume, include_single=include_single)
        
        if not jobs:
            self.logger.warning("没有找到待提交的作业")
            return
        
        # 添加到队列
        self.queue.add_jobs(jobs)
        
        self.logger.info(f"队列中共有 {self.queue.pending_count} 个作业待提交")
        
        # 开始提交循环
        last_check_time = 0
        last_config_reload_time = time.time()
        
        try:
            while self.running and self.queue.pending_count > 0:
                current_time = time.time()
                
                # 定期检查并重新加载配置
                if current_time - last_config_reload_time >= self.config.config_reload_interval:
                    if self.config.check_and_reload():
                        self.logger.info("✅ 配置已更新并应用")
                        self.logger.info(f"   当前节点限制: {self.config.total_nodes}")
                        self.logger.info(f"   当前任务限制: {self.config.max_jobs}")
                        self.logger.info(f"   当前提交间隔: {self.config.interval_time}秒")
                        self.logger.info(f"   当前检查间隔: {self.config.check_time}秒")
                    last_config_reload_time = current_time
                
                # 定期刷新队列状态
                if current_time - last_check_time >= self.config.check_time:
                    total_jobs, running_jobs, used_nodes = self.monitor.get_queue_status()
                    self.logger.info(
                        f"队列状态 - 总任务:{total_jobs}, "
                        f"运行中:{running_jobs}, "
                        f"使用节点:{used_nodes}/{self.config.total_nodes}, "
                        f"待提交:{self.queue.pending_count}"
                    )
                    last_check_time = current_time
                
                # 获取下一个作业
                job = self.queue.get_next_job()
                if not job:
                    break
                
                # 检查是否可以提交
                can_submit, reason = self.monitor.can_submit(job)
                
                if not can_submit:
                    self.logger.info(f"资源不足，等待中... ({reason})")
                    # 放回队列
                    self.queue.add_job(job)
                    # 等待一段时间后重试
                    time.sleep(self.config.check_time)
                    continue
                
                # 提交作业
                success, result = self.submit_job(job)
                
                if success:
                    job_id = result
                    self.queue.mark_submitted(job, job_id)
                else:
                    error_msg = result
                    
                    # 检查是否需要重试
                    if job.retry_count < max_retries:
                        self.logger.warning(
                            f"作业提交失败，将重试 ({job.retry_count + 1}/{max_retries}): "
                            f"{job.job_name}, 原因: {error_msg[:100]}"
                        )
                        self.queue.add_job(job)
                    else:
                        self.logger.error(
                            f"作业 {job.job_name} 重试次数已达上限，标记为失败"
                        )
                        self.queue.mark_failed(job, error_msg)
                
                # 间隔等待
                time.sleep(self.config.interval_time)
                
        except KeyboardInterrupt:
            self.logger.info("\n用户中断，正在停止...")
            self.running = False
        
        # 输出统计信息
        self.logger.info("=" * 70)
        self.logger.info("提交完成统计:")
        self.logger.info(f"  成功提交: {self.queue.submitted_count}")
        self.logger.info(f"  待提交: {self.queue.pending_count}")
        self.logger.info(f"  失败: {self.queue.failed_count}")
        self.logger.info("=" * 70)
        
        # 输出失败作业详情
        if self.queue.failed_jobs:
            self.logger.info("\n失败作业列表:")
            for job, reason in self.queue.failed_jobs:
                self.logger.info(f"  - {job.job_name}: {reason}")
    
    def stop(self):
        """停止自动提交"""
        self.running = False


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Slurm自动提交管理器',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 提交work_cal7目录下的所有作业（workflow模式）
  python -m abacus.submit_manager work_cal7
  
  # 提交续算脚本
  python -m abacus.submit_manager batch_work/ --resume
  
  # 提交单步计算脚本
  python -m abacus.submit_manager batch_work_single/ --single
  
  # 使用自定义配置文件
  python -m abacus.submit_manager work_cal7 --config my_config.ini
  
  # 设置最大重试次数
  python -m abacus.submit_manager work_cal7 --max-retries 5
        """
    )
    
    parser.add_argument(
        'work_dir',
        help='工作目录，包含待提交的作业脚本'
    )
    
    parser.add_argument(
        '--config',
        default='config/condor.ini',
        help='配置文件路径 (默认: config/condor.ini)'
    )
    
    parser.add_argument(
        '--max-retries',
        type=int,
        default=3,
        help='最大重试次数 (默认: 3)'
    )
    
    parser.add_argument(
        '--pattern',
        default='*.sh',
        help='脚本文件匹配模式 (默认: *.sh)'
    )
    
    parser.add_argument(
        '--resume',
        action='store_true',
        help='提交续算脚本 (*_resume.sh)'
    )
    
    parser.add_argument(
        '--single',
        action='store_true',
        help='提交单步计算脚本 (*_single.sh)'
    )
    
    args = parser.parse_args()
    
    # 创建提交器并运行
    submitter = AutoSubmitter(config_file=args.config)
    submitter.run(args.work_dir, max_retries=args.max_retries, include_resume=args.resume, include_single=args.single)


if __name__ == '__main__':
    main()

