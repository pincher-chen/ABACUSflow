#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path
import configparser

_config_root = Path(os.path.abspath(__file__)).parent
_condor = _config_root / "condor.ini"
_workflow = _config_root / "workflow.json"

# 读取配置文件
CONDOR = configparser.ConfigParser()
if _condor.exists():
    CONDOR.read(_condor)
else:
    # 如果配置文件不存在，创建默认配置
    CONDOR = configparser.ConfigParser()
    CONDOR.add_section('ABACUS')
    CONDOR.add_section('STRU')
    CONDOR.add_section('MODULE')
    CONDOR.add_section('ALLOW')
    CONDOR.add_section('PARAMETERS')

# 读取工作流配置
WORKFLOW = {}
if _workflow.exists():
    import json
    with open(_workflow, 'r') as f:
        WORKFLOW = json.load(f)

# 加载模板
INCAR_TEMPLATE = {}
_temp_dir = _config_root / "template"
if _temp_dir.exists():
    import yaml
    for template_file in _temp_dir.glob("*.yaml"):
        with open(template_file, 'r') as f:
            stage_name = template_file.stem
            try:
                INCAR_TEMPLATE[stage_name] = yaml.safe_load(f)
            except Exception as e:
                print(f"Error loading template {template_file}: {e}")

PACKAGE_ROOT = _config_root.parent

def get(section, key, fallback=None):
    """获取配置值"""
    try:
        return CONDOR.get(section, key, fallback=fallback)
    except (configparser.NoSectionError, configparser.NoOptionError):
        return fallback

if __name__ == '__main__':
    # 测试读取配置
    print("ABACUS_DIR:", get('ABACUS', 'ABACUS_DIR'))
    print("PARTITION:", get('ALLOW', 'PARTITION'))

