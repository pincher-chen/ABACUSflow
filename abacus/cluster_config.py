# -*- coding: utf-8 -*-
"""
集群相关配置的统一读取层。

设计目标：让 ABACUSflow 不依赖任何特定集群的硬编码路径 / 调度器命令。
所有集群相关信息都从 config/condor.ini 读取，未配置时回落到通用默认值
（而不是作者所在集群的路径）。

配置文件查找顺序（优先级从高到低）：
  1. 环境变量 ABACUSFLOW_CONFIG 指定的文件
  2. 当前工作目录下的 condor.ini
  3. 当前工作目录下的 config/condor.ini
  4. 包内默认文件 <包目录>/config/condor.ini
  5. 环境变量 ABACUSFLOW_HOME 下的 config/condor.ini

每个配置项都支持用环境变量覆盖，便于在不改文件的情况下适配新集群，
例如：ABACUSFLOW_ABACUS_DIR=/opt/abacus/bin
"""

import configparser
import os
import shutil
from pathlib import Path
from typing import Optional

ENV_PREFIX = "ABACUSFLOW_"

# 各配置项在“未设置”时可以视为兼容的旧键名（用于平滑迁移）
_LEGACY_KEYS = {
    ("ABACUS", "ABACUS_BIN"): [("ABACUS", "ABACUS_EXE")],
}

# (section, key) -> 可用的环境变量名列表（ABACUSFLOW_<KEY> 始终可用）
_ENV_ALIASES_FOR = {
    ("SCHEDULER", "TYPE"): ["ABACUSFLOW_SCHEDULER"],
    ("ENV", "SETUP"): ["ABACUSFLOW_ENV_SETUP"],
}


def find_config_file(explicit: Optional[str] = None) -> Optional[str]:
    """按优先级查找配置文件，返回第一个存在的路径。

    顺序：
      1. 显式参数
      2. 环境变量 ABACUSFLOW_CONFIG（指向一个 ini 文件）
      3. 环境变量 ABACUSFLOW_CONFIG_DIR（指向一个配置目录，取其中的
         condor.ini；兼容直接指向 config/ 或指向项目根目录两种写法）
      4. 从当前目录向上逐级查找 condor.ini 或 config/condor.ini
         （允许把配置放在作业目录的父目录，例如 <proj>/config/condor.ini）
      5. 环境变量 ABACUSFLOW_HOME 下
      6. 包内默认 config/condor.ini
    """
    candidates = []
    if explicit:
        candidates.append(explicit)
    if os.environ.get("ABACUSFLOW_CONFIG"):
        candidates.append(os.environ["ABACUSFLOW_CONFIG"])

    config_dir = os.environ.get("ABACUSFLOW_CONFIG_DIR")
    if config_dir:
        # 目录可能直接是配置目录，也可能是项目根目录（其下有 config/）
        candidates.append(os.path.join(config_dir, "condor.ini"))
        candidates.append(os.path.join(config_dir, "config", "condor.ini"))

    cur = Path(os.getcwd()).resolve()
    for parent in [cur, *cur.parents]:
        candidates.append(str(parent / "condor.ini"))
        candidates.append(str(parent / "config" / "condor.ini"))

    if os.environ.get("ABACUSFLOW_HOME"):
        candidates.append(os.path.join(os.environ["ABACUSFLOW_HOME"], "config", "condor.ini"))

    pkg_root = Path(__file__).resolve().parent.parent
    candidates.append(str(pkg_root / "config" / "condor.ini"))

    for cand in candidates:
        if cand and os.path.isfile(cand):
            return os.path.abspath(cand)
    return None


def load_config(explicit: Optional[str] = None) -> configparser.ConfigParser:
    """加载配置，缺文件时返回空配置（各项走默认值）。"""
    cfg = configparser.ConfigParser()
    path = find_config_file(explicit)
    if path:
        cfg.read(path, encoding="utf-8")
    return cfg


class ClusterConfig:
    """集群配置访问器：ini -> 环境变量覆盖 -> 默认值。"""

    def __init__(self, parser: Optional[configparser.ConfigParser] = None,
                 path: Optional[str] = None):
        self.parser = parser if parser is not None else load_config(path)
        self.path = path or find_config_file()

    # ---- 底层读取 --------------------------------------------------------
    def get(self, section: str, key: str, fallback: Optional[str] = None) -> Optional[str]:
        # 环境变量覆盖：别名（如 ABACUSFLOW_SCHEDULER）+ 通用 ABACUSFLOW_<KEY>
        names = _ENV_ALIASES_FOR.get((section, key), []) + [ENV_PREFIX + key]
        for name in names:
            val = os.environ.get(name)
            if val not in (None, ""):
                return val

        try:
            val = self.parser.get(section, key)
            if val.strip() != "":
                return val.strip()
        except (configparser.NoSectionError, configparser.NoOptionError):
            pass

        for sec, k in _LEGACY_KEYS.get((section, key), []):
            try:
                val = self.parser.get(sec, k)
                if val.strip() != "":
                    return val.strip()
            except (configparser.NoSectionError, configparser.NoOptionError):
                pass
        return fallback

    def get_required(self, section: str, key: str, purpose: str) -> str:
        val = self.get(section, key)
        if not val:
            raise SystemExit(
                f"[ERROR] 缺少配置项 [{section}] {key}（用于{purpose}）。\n"
                f"        请在 {self.path or 'config/condor.ini'} 中设置，"
                f"或设置环境变量 {ENV_PREFIX}{key}。"
            )
        return val

    # ---- 环境初始化 ------------------------------------------------------
    @property
    def env_setup(self) -> Optional[str]:
        """自定义环境初始化脚本片段，例如 'module load ...' 或 'source ...'。"""
        return self.get("ENV", "SETUP")

    @property
    def conda_path(self) -> Optional[str]:
        return self.get("ENV", "CONDA_PATH")

    @property
    def conda_env(self) -> Optional[str]:
        return self.get("ENV", "CONDA_ENV")

    @property
    def python_cmd(self) -> str:
        """脚本中调用 Python 的命令，默认为 'python'。"""
        return self.get("ENV", "PYTHON", "python")

    @property
    def abacus_py(self) -> Optional[str]:
        return self.get("ENV", "ABACUS_PY")

    # ---- ABACUS ---------------------------------------------------------
    @property
    def abacus_bin(self) -> str:
        """ABACUS 可执行文件。

        优先级：ABACUS_BIN（可直接是 PATH 中的命令名，如容器内的 'abacus'）
                > ABACUS_DIR + ABACUS_EXE。
        两者都未配置时报错，避免静默生成一个不可用的裸命令。
        """
        full = self.get("ABACUS", "ABACUS_BIN")
        if full:
            return full
        abacus_dir = self.get("ABACUS", "ABACUS_DIR")
        if not abacus_dir:
            raise SystemExit(
                "[ERROR] 未配置 ABACUS 可执行文件。请在 "
                f"{self.path or 'config/condor.ini'} 中设置 [ABACUS] ABACUS_BIN "
                "(可执行文件完整路径，容器内可直接写 abacus)，"
                "或设置 [ABACUS] ABACUS_DIR，"
                "或设置环境变量 ABACUSFLOW_ABACUS_BIN / ABACUSFLOW_ABACUS_DIR。"
            )
        exe = self.get("ABACUS", "ABACUS_EXE", "abacus")
        return os.path.join(abacus_dir, exe)

    @property
    def pseudo_dir(self) -> Optional[str]:
        return self.get("ABACUS", "PSEUDO_POTENTIAL_DIR")

    @property
    def orbital_dir(self) -> Optional[str]:
        return self.get("ABACUS", "ORBITAL_DIR")

    # ---- Module ---------------------------------------------------------
    @property
    def modules(self) -> Optional[str]:
        return self.get("MODULE", "MODULES")

    @property
    def module_setup(self) -> Optional[str]:
        """module 命令初始化片段；未配置时只输出注释，不做任何 export。"""
        return self.get("MODULE", "MODULE_SETUP")

    # ---- 调度器 ---------------------------------------------------------
    @property
    def scheduler_type(self) -> str:
        return (self.get("SCHEDULER", "TYPE", "slurm") or "slurm").lower()

    @property
    def submit_cmd(self) -> str:
        """提交作业命令，未配置时按调度器类型推断。"""
        explicit = self.get("SCHEDULER", "SUBMIT_CMD")
        if explicit:
            return explicit
        return "yhbatch" if self.scheduler_type == "yh" else "sbatch"

    @property
    def run_cmd(self) -> str:
        """运行 MPI 任务的命令模板，支持占位符 {nodes} {cores} {partition} {cmd}。"""
        explicit = self.get("SCHEDULER", "RUN_CMD")
        if explicit:
            return explicit
        if self.scheduler_type == "yh":
            return "yhrun -N {nodes} -n {cores} -p {partition} {cmd}"
        return "srun -N {nodes} -n {cores} -p {partition} {cmd}"

    @property
    def queue_cmd(self) -> str:
        explicit = self.get("SCHEDULER", "QUEUE_CMD")
        if explicit:
            return explicit
        return "squeue" if self.scheduler_type == "slurm" else "squeue"

    def build_run_command(self, nodes, cores, partition: Optional[str], cmd: str) -> str:
        template = self.run_cmd
        tmpl = template.replace("{nodes}", str(nodes)).replace("{cores}", str(cores))
        tc = tmpl.replace("{cmd}", cmd)
        if "{partition}" in tc:
            if not partition:
                # 模板需要分区但未配置：去掉 -p {partition} 片段
                tc = tc.replace("-p {partition}", "").replace("--partition={partition}", "")
                tc = tc.replace("--partition {partition}", "")
            else:
                tc = tc.replace("{partition}", partition)
        return " ".join(tc.split())

    @property
    def partition(self) -> Optional[str]:
        return self.get("ALLOW", "PARTITION")

    # ---- 结构文件 -------------------------------------------------------
    @property
    def stru_suffixes(self) -> list:
        raw = self.get("STRU", "SUFFIX", "*.vasp *.poscar *.cif *.mcif") or ""
        return raw.split()

    # ---- 资源限制 -------------------------------------------------------
    def get_int(self, key: str, default: int, section: str = "ALLOW") -> int:
        val = self.get(section, key)
        try:
            return int(val) if val is not None else default
        except (TypeError, ValueError):
            return default

    def get_float(self, key: str, default: float, section: str = "ALLOW") -> float:
        val = self.get(section, key)
        try:
            return float(val) if val is not None else default
        except (TypeError, ValueError):
            return default


def which_or_none(cmd: str) -> Optional[str]:
    return shutil.which(cmd)


def default_python() -> str:
    """当前解释器路径，用于生成脚本时保证与运行时一致。"""
    import sys
    return sys.executable or "python"


# 生成脚本时使用的 Python 命令：优先 condor.ini [ENV] PYTHON，
# 否则使用运行本程序所用的解释器（容器/SIF 内即为镜像自带的 python）。
def runtime_python(ccfg: "ClusterConfig") -> str:
    """生成脚本中使用的 Python 命令。

    优先 condor.ini [ENV] PYTHON。未配置时回落到 'python'，避免把本机
    绝对路径写进脚本，使脚本在 SIF / 其他节点上依然可移植。
    """
    explicit = ccfg.get("ENV", "PYTHON")
    return explicit if explicit else "python"


def render_env_setup(ccfg: "ClusterConfig", omp_threads: Optional[int] = None) -> list:
    """渲染生成脚本中的“环境初始化”片段。

    只输出配置中真实声明过的内容，不再写入任何集群专有的路径或 module
    初始化命令。未配置时脚本依然可用（依赖脚本运行时已有的环境）。
    """
    lines = []
    if omp_threads is not None:
        lines.append(f"export OMP_NUM_THREADS={omp_threads}")

    # 1) 完全自定义：直接插入配置文件中的片段
    setup = ccfg.env_setup
    if setup:
        lines.append("# 自定义环境初始化（condor.ini [ENV] SETUP）")
        lines.extend(setup.replace("\\n", "\n").splitlines())

    # 2) Conda（仅在配置了路径时）
    conda_path = ccfg.conda_path
    if conda_path:
        conda_env = ccfg.conda_env
        lines.append("# conda 环境")
        lines.append(f"source {conda_path.rstrip('/')}/etc/profile.d/conda.sh")
        if conda_env:
            lines.append(f"conda activate {conda_env}")

    # 3) module（仅在配置了 module 列表或初始化片段时）
    modules_str = ccfg.modules
    module_setup = ccfg.module_setup
    if modules_str or module_setup:
        lines.append("")
        lines.append("# 加载环境模块（condor.ini [MODULE]）")
        if module_setup:
            lines.extend(module_setup.replace("\\n", "\n").splitlines())
        if modules_str:
            for mod in modules_str.split():
                lines.append(f"module load {mod} 2>/dev/null || echo \"[WARN] module load {mod} failed\"")

    return lines

