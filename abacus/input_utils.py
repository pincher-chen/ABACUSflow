#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INPUT 文件处理工具模块
包含 ABACUS INPUT 文件的解析和更新功能
"""

def parse_input_file(content):
    """解析 INPUT 文件为字典"""
    params = {}
    for line in content.split('\n'):
        line = line.strip()
        if line and not line.startswith('#') and not line.startswith('INPUT_PARAMETERS'):
            parts = line.split()
            if len(parts) >= 2:
                params[parts[0]] = parts[1]
    return params

def write_input_file(params, original_content):
    """根据参数字典更新 INPUT 文件内容"""
    lines = original_content.split('\n')
    new_lines = []
    
    for line in lines:
        stripped = line.strip()
        if stripped and not stripped.startswith('#') and not stripped.startswith('INPUT_PARAMETERS'):
            parts = stripped.split()
            if len(parts) >= 2:
                key = parts[0]
                if key in params:
                    # 保持原有的格式和注释
                    indent = len(line) - len(line.lstrip())
                    comment_parts = line.split('#', 1)
                    comment = f" # {comment_parts[1]}" if len(comment_parts) > 1 else ""
                    new_line = ' ' * indent + f"{key:<18} {params[key]}{comment}"
                    new_lines.append(new_line)
                    continue
        new_lines.append(line)
    
    return '\n'.join(new_lines)


