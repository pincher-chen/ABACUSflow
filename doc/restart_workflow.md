# ABACUS 工作流续算功能设计文档

版本: v1.0  
日期: 2025-12-25  
作者: AbacusFlow Team

---

## 📋 目录

- [1. 概述](#1-概述)
- [2. 设计目标](#2-设计目标)
- [3. 核心功能](#3-核心功能)
- [4. 技术实现](#4-技术实现)
- [5. 使用指南](#5-使用指南)
- [6. 文件命名规范](#6-文件命名规范)
- [7. 故障排查](#7-故障排查)
- [8. 实现细节](#8-实现细节)

---

## 1. 概述

### 1.1 背景

ABACUS 高通量计算工作流包含多个顺序执行的阶段：

```
Test_spin → Coarse_relax → Relax → Scf → Band → Dos
```

在实际计算中，工作流可能因为以下原因中断：
- 某个阶段计算失败
- 参数设置不当导致不收敛
- 系统资源限制
- 作业时间超限

**续算功能**允许用户从失败或中断的阶段继续计算，而不需要从头重新运行整个工作流，节省计算资源和时间。

### 1.2 主要特性

- ✅ **智能状态检测**：自动分析工作流各阶段状态
- ✅ **自动续算起点**：智能确定从哪个阶段开始续算
- ✅ **日志备份**：自动备份历史日志，保留完整记录
- ✅ **错误清理**：自动清理失败标记，准备重新计算
- ✅ **两种续算策略**：
  - 首次尝试就失败 → 清空重算（可能是参数错误）
  - 多次尝试未收敛 → 继续迭代（利用已有结果）
- ✅ **续算编号系统**：通过 `_R{N}_` 区分不同批次的续算
- ✅ **Dry-run 模式**：预览续算计划而不实际执行

---

## 2. 设计目标

### 2.1 功能目标

1. **自动化**：用户只需一条命令即可完成续算
2. **可靠性**：准确判断各阶段状态，避免误判
3. **灵活性**：支持手动指定续算起点
4. **可追溯性**：完整保留所有尝试历史
5. **向后兼容**：不影响现有工作流的正常运行

### 2.2 非功能目标

1. **易用性**：命令简单直观，提示信息清晰
2. **可扩展性**：易于添加新的续算策略
3. **健壮性**：处理各种异常情况
4. **可维护性**：代码模块化，职责清晰

---

## 3. 核心功能

### 3.1 状态检测系统

#### 3.1.1 阶段状态定义

| 状态 | 图标 | 说明 | 判断依据 |
|------|------|------|----------|
| `success` | ✅ | 成功完成 | 存在 `converge.txt` |
| `ignored` | ⚠️  | 未收敛但被忽略 | 存在 `ignore.txt` 但无 `converge.txt` |
| `failed` | ❌ | 计算失败 | 存在 `error.txt` |
| `incomplete` | 🔄 | 运行未完成 | 存在 `running.log` 但无结论文件 |
| `not_started` | ⭕ | 未开始 | 目录不存在或为空 |

#### 3.1.2 检测优先级

```python
1. error.txt 存在 → failed
2. converge.txt 存在 → success  
3. ignore.txt 存在 → ignored
4. running.log 存在 → incomplete
5. 其他 → not_started
```

#### 3.1.3 stat.log 权威性

`stat.log` 是工作流状态的**权威记录**，续算时优先参考其内容：

```
Test_spin success
Coarse_relax success (ignored)
Relax failed
```

### 3.2 续算起点确定

#### 3.2.1 自动检测逻辑

```python
for stage in workflow_stages:
    if stat.log 记录该阶段 failed:
        return stage  # 从这里开始
    
    if 目录不存在:
        return stage  # 从这里开始
    
    if 目录存在但 incomplete:
        return stage  # 从这里开始
    
    if success or ignored:
        continue  # 继续检查下一个阶段

return None  # 所有阶段都完成
```

#### 3.2.2 手动指定

用户可以使用 `--from-stage` 参数强制指定续算起点：

```bash
python abacus.py resume work_cal/hmat_0/ --from-stage Scf
```

### 3.3 续算策略分析

对于失败的阶段，系统会智能判断应该如何处理：

#### 策略 A：清空重算

**适用场景**：首次尝试就失败

**判断依据**：
- 目录中没有 `INPUT0`, `INPUT1` 等备份文件
- 说明只运行了一次就失败，可能是输入参数有严重错误

**处理方式**：
```bash
rm -rf Relax/
mkdir Relax
# 重新生成输入文件并计算
```

#### 策略 B：继续迭代

**适用场景**：多次尝试未收敛

**判断依据**：
- 存在 `INPUT0`, `INPUT1`, `INPUT2` 等多个备份文件
- 对于 Relax 类型计算，存在 `OUT.*/STRU_ION_D` 优化结构

**处理方式**：
```bash
# 保留目录和已有结果
cd Relax/
# 备份当前文件
cp INPUT INPUT_R0_pre
cp KPT KPT_R0_pre  
cp STRU STRU_R0_pre
# 使用优化结构继续
cp OUT.Relax/STRU_ION_D STRU
# 从 try_num=3 继续尝试
```

### 3.4 日志备份系统

#### 3.4.1 备份规则

续算前自动备份 `stat.log` 和 `time.log`：

```bash
stat.log → stat0.log   # 首次续算
stat.log → stat1.log   # 第二次续算
stat.log → stat2.log   # 第三次续算
...
```

自动递增编号，避免覆盖历史记录。

#### 3.4.2 日志追加模式

续算脚本不会覆盖现有日志，而是追加新的记录：

```bash
# 续算脚本中
stat_log=stat.log  # 不覆盖
time_log=time.log
echo '' >> $stat_log  # 追加
echo '# === Resume from Relax at 2025-12-25 15:30:00 ===' >> $time_log
```

---

## 4. 技术实现

### 4.1 模块架构

```
abacusflow/
├── abacus.py                 # CLI 入口，包含 resume 命令
├── abacus/
│   ├── submit.py             # 工作流脚本生成
│   │   ├── generate_scripts()          # 生成完整工作流
│   │   ├── generate_resume_script()    # 生成续算脚本
│   │   ├── _create_workflow_script()   # 统一的脚本生成逻辑
│   │   ├── _generate_stage_script()    # 单个阶段脚本
│   │   └── _get_resume_helper_functions() # Bash 辅助函数
│   └── resume_utils.py       # 续算辅助函数
│       ├── detect_stage_status()       # 检测阶段状态
│       ├── parse_stat_log()            # 解析 stat.log
│       ├── determine_resume_point()    # 确定续算起点
│       ├── backup_logs()               # 备份日志
│       ├── clean_stage_errors()        # 清理错误
│       ├── detect_resume_number()      # 检测续算次数
│       ├── detect_previous_try_count() # 检测尝试次数
│       └── should_clean_and_restart()  # 判断续算策略
```

### 4.2 核心流程

```
用户执行: python abacus.py resume work_cal/hmat_0/
    ↓
读取工作目录，检查 STRU.vasp 是否存在
    ↓
遍历所有阶段，检测状态
    ├─ 读取目录中的状态文件 (converge.txt, error.txt, etc.)
    ├─ 解析 stat.log
    └─ 显示状态表格
    ↓
确定续算起点
    ├─ 自动检测（默认）
    └─ 手动指定（--from-stage）
    ↓
分析续算策略
    ├─ 检测续算次数 (RESUME_NUM)
    ├─ 检测之前尝试次数 (PREV_TRY_COUNT)
    └─ 判断是清空还是继续
    ↓
备份日志 (stat.log → stat{N}.log)
    ↓
清理错误标记 (error.txt)
    ↓
生成续算脚本 ({job_name}_resume.sh)
    ├─ 脚本头部（续算信息）
    ├─ 环境初始化
    ├─ Bash 辅助函数（续算模式）
    ├─ 从续算起点开始的各阶段
    │   ├─ 目录智能处理
    │   ├─ 输入文件生成
    │   ├─ ABACUS 运行
    │   ├─ 结果检查
    │   └─ 续算备份 (*_R{N}_{try_num})
    └─ 工作流总结
    ↓
提示用户提交作业
```

### 4.3 关键数据结构

#### 4.3.1 阶段状态映射

```python
{
    'Test_spin': 'success',
    'Coarse_relax': 'ignored',
    'Relax': 'failed',
    'Scf': 'not_started',
    'Band': 'not_started',
    'Dos': 'not_started'
}
```

#### 4.3.2 续算策略决策

```python
(should_clean, reason) = should_clean_and_restart(stage_dir, stage_name)
# 返回:
# (True, "First attempt failed, input likely has errors")
# (False, "Continue from attempt 3")
```

---

## 5. 使用指南

### 5.1 基本用法

#### 5.1.1 自动续算（推荐）

```bash
# 进入工作目录
cd /XYFS01/nscc-gz_pinchen_1/work/dev/abacusflow/work_cal7/hmat_0

# 执行续算（自动检测起点）
python /path/to/abacus.py resume .

# 输出示例：
# ======================================================================
# 📊 Current Workflow Status
# ======================================================================
#   ✅ Test_spin      - File: success      | Log: success
#   ⚠️  Coarse_relax  - File: ignored      | Log: ignored
#   ❌ Relax          - File: failed       | Log: failed
#   ⭕ Scf            - File: not_started  | Log: not_started
#   ⭕ Band           - File: not_started  | Log: not_started
#   ⭕ Dos            - File: not_started  | Log: not_started
# ======================================================================
# 
# 🎯 Resume Plan: Relax failed in previous run
#    Starting from: Relax
#    Stages to run: Relax, Scf, Band, Dos
#
# 📋 Resume Strategy Analysis:
#   📁 Relax:
#      Resume attempt: #0
#      Previous tries: 2
#      Action: 🔄 Continue
#      Reason: Continue from attempt 2
#
# 📦 Backup Plan:
#    stat.log will be backed up
#    time.log will be backed up
#
# 🧹 Cleanup Plan:
#    Relax/error.txt will be removed
#
# 🔧 Generating resume script...
# [BACKUP] stat.log -> stat0.log
# [BACKUP] time.log -> time0.log
# [CLEAN] Removed Relax/error.txt
# Generated resume script: /path/to/hmat_0/hmat_0_resume.sh
#
# ✅ Resume script generated successfully!
#    Script: /path/to/hmat_0/hmat_0_resume.sh
#
# 📝 Next steps:
#    cd /path/to/hmat_0
#    yhbatch hmat_0_resume.sh

# 提交续算作业
yhbatch hmat_0_resume.sh
```

#### 5.1.2 预览模式（Dry-run）

```bash
# 仅查看续算计划，不生成脚本
python abacus.py resume . --dry-run

# 输出续算计划后退出，不执行任何操作
```

#### 5.1.3 指定续算起点

```bash
# 强制从 Scf 阶段开始（忽略自动检测）
python abacus.py resume . --from-stage Scf
```

#### 5.1.4 不备份日志

```bash
# 跳过日志备份（覆盖原有日志）
python abacus.py resume . --no-backup
```

#### 5.1.5 不清理错误

```bash
# 保留 error.txt（用于调试）
python abacus.py resume . --no-clean
```

### 5.2 典型场景

#### 场景 1：首次尝试就失败

```bash
# 初始状态
Relax/
├── INPUT              # 参数有误，导致立即失败
├── KPT
├── STRU
├── error.txt
└── running.log

# 执行续算
python abacus.py resume .

# 续算策略：清空重算
# [CLEAN] Removing Relax directory for fresh start

# 续算后
Relax/
├── INPUT              # 重新生成（请先手动修复参数）
├── KPT
├── STRU
└── running.log
```

**注意**：如果是输入参数错误，建议先手动修改 template 或 workflow.json，再续算。

#### 场景 2：多次尝试未收敛

```bash
# 初始状态（已尝试 3 次）
Relax/
├── INPUT, KPT, STRU   # 最后一次
├── INPUT0, KPT0, STRU0  # 第 1 次备份
├── INPUT1, KPT1, STRU1  # 第 2 次备份
├── INPUT2, KPT2, STRU2  # 第 3 次备份
├── OUT.Relax/STRU_ION_D
└── error.txt

# 执行续算
python abacus.py resume .

# 续算策略：继续迭代
# [CONTINUE] Previous attempts: 3
# [BACKUP] INPUT -> INPUT_R0_pre
# [CONTINUE] Using optimized structure from OUT.Relax/STRU_ION_D

# 续算后（从 try_num=3 开始，再尝试 2 次）
Relax/
├── INPUT, KPT, STRU
├── INPUT0~2, KPT0~2, STRU0~2     # 原始备份
├── INPUT_R0_pre, KPT_R0_pre, STRU_R0_pre  # 续算前状态
├── INPUT_R0_3, KPT_R0_3, STRU_R0_3  # 续算第 4 次尝试
├── INPUT_R0_4, KPT_R0_4, STRU_R0_4  # 续算第 5 次尝试
└── OUT.Relax/STRU_ION_D
```

#### 场景 3：二次续算

```bash
# 首次续算后仍未收敛，再次续算
python abacus.py resume .

# 续算编号自动递增
# [RESUME] This is resume attempt #1
# [BACKUP] INPUT -> INPUT_R1_pre

# 二次续算后
Relax/
├── INPUT, KPT, STRU
├── INPUT0~2                      # 原始备份
├── INPUT_R0_pre, INPUT_R0_3~4    # 首次续算
├── INPUT_R1_pre, INPUT_R1_5~6    # 二次续算
└── OUT.Relax/STRU_ION_D
```

#### 场景 4：批量续算

```bash
# 找出所有失败的计算
cd work_cal7

for dir in */; do
    if [ -f "$dir/stat.log" ] && grep -q "failed" "$dir/stat.log"; then
        echo "Generating resume script for $dir"
        python /path/to/abacus.py resume "$dir"
    fi
done

# 批量提交
for script in */*_resume.sh; do
    echo "Submitting $script"
    yhbatch "$script"
done
```

### 5.3 续算后的修改

#### 5.3.1 调整 K 点密度

续算前修改 `config/workflow.json`：

```json
{
  "Relax": {
    "kval": [0.04, 0.02, 0.015],  # 添加更密的 K 点
    "try_num": 3
  }
}
```

然后续算，系统会自动使用新的 K 点密度。

#### 5.3.2 调整收敛参数

修改 `config/template/Relax.yaml`：

```yaml
scf_nmax: 150      # 增加 SCF 步数
mixing_beta: 0.3   # 降低混合参数
force_thr_ev: 0.02 # 放宽力收敛标准
```

然后续算。

---

## 6. 文件命名规范

### 6.1 原始工作流备份

在原始工作流中，每次尝试的文件会被备份：

```
INPUT   → INPUT0  (第 1 次尝试的备份)
INPUT   → INPUT1  (第 2 次尝试的备份)
INPUT   → INPUT2  (第 3 次尝试的备份)
...
```

同理适用于 `KPT` 和 `STRU`。

### 6.2 续算模式备份

续算模式使用新的命名格式：

```
格式: {文件名}_R{续算编号}_{状态}

示例:
INPUT_R0_pre    # 首次续算前的状态
INPUT_R0_0      # 首次续算的第 1 次尝试
INPUT_R0_1      # 首次续算的第 2 次尝试

INPUT_R1_pre    # 二次续算前的状态
INPUT_R1_2      # 二次续算的第 1 次尝试
INPUT_R1_3      # 二次续算的第 2 次尝试
```

### 6.3 完整示例

假设 Relax 阶段经历了以下过程：

```
1. 原始运行（try_num=3）
   INPUT0, INPUT1, INPUT2

2. 首次续算（try_num=2，从 3 开始）
   INPUT_R0_pre, INPUT_R0_3, INPUT_R0_4

3. 二次续算（try_num=2，从 5 开始）
   INPUT_R1_pre, INPUT_R1_5, INPUT_R1_6
```

最终目录结构：

```
Relax/
├── INPUT, KPT, STRU           # 最新文件
├── INPUT0~2                   # 原始尝试
├── INPUT_R0_pre               # 首次续算前
├── INPUT_R0_3~4               # 首次续算尝试
├── INPUT_R1_pre               # 二次续算前
├── INPUT_R1_5~6               # 二次续算尝试
├── KPT0~2, STRU0~2
├── KPT_R0_pre~4, STRU_R0_pre~4
├── KPT_R1_pre, KPT_R1_5~6
├── STRU_R1_pre, STRU_R1_5~6
└── OUT.Relax/STRU_ION_D
```

### 6.4 日志备份

```
stat.log   # 当前日志
stat0.log  # 首次续算备份
stat1.log  # 二次续算备份
stat2.log  # 三次续算备份
...

time.log   # 当前日志
time0.log  # 首次续算备份
time1.log  # 二次续算备份
time2.log  # 三次续算备份
...
```

---

## 7. 故障排查

### 7.1 常见问题

#### Q1: 续算脚本生成失败

```
[ERROR] Failed to generate resume script: Invalid stage: Relax
```

**原因**：指定的阶段名称不在 workflow.json 中

**解决**：检查 `--from-stage` 参数是否正确，必须是以下之一：
- Test_spin
- Coarse_relax
- Relax
- Scf
- Band
- Dos

#### Q2: 找不到 STRU.vasp

```
[ERROR] Invalid work directory (no STRU.vasp found)
```

**原因**：工作目录路径不正确

**解决**：确保在正确的工作目录中（应该包含 STRU.vasp 文件）

```bash
# 正确路径示例
cd work_cal/hmat_0/
python /path/to/abacus.py resume .
```

#### Q3: 续算后仍然失败

如果续算后计算仍然失败，可能的原因：

1. **参数设置不当**
   - 检查 K 点密度是否太粗
   - 检查 SCF 参数是否合理
   - 检查结构是否有问题

2. **收敛标准太严格**
   - 适当放宽 `force_thr_ev`
   - 增加 `scf_nmax` 和 `relax_nmax`

3. **系统本身难以收敛**
   - 考虑更换 `mixing_type`
   - 尝试不同的 `smearing_method`
   - 检查赝势文件是否合适

#### Q4: 备份文件太多占用空间

可以定期清理旧的备份文件：

```bash
# 只保留最近的备份
cd work_cal/hmat_0/Relax/
rm -f INPUT_R0_* KPT_R0_* STRU_R0_*  # 删除首次续算备份

# 清理所有续算备份
find . -name "*_R*_*" -delete
```

**注意**：确保不需要这些历史记录后再删除。

### 7.2 调试技巧

#### 使用 --dry-run 预览

```bash
# 先预览续算计划
python abacus.py resume . --dry-run

# 确认无误后再实际执行
python abacus.py resume .
```

#### 检查续算脚本

```bash
# 查看生成的续算脚本
cat hmat_0_resume.sh

# 检查辅助函数
grep -A 20 "detect_resume_number" hmat_0_resume.sh
```

#### 手动测试阶段

```bash
# 进入失败的阶段目录
cd Relax/

# 手动生成输入文件
python /path/to/abacus.py generate --work_dir . --stage Relax

# 检查生成的文件
ls -lh INPUT KPT STRU

# 手动运行一次测试
yhrun -N 1 -n 16 -p deimos /path/to/abacus > test.log 2>&1
```

---

## 8. 实现细节

### 8.1 Bash 辅助函数

续算脚本中包含以下 Bash 辅助函数：

```bash
detect_resume_number() {
    # 检测当前是第几次续算
    # 返回: 0, 1, 2, ...
}

detect_previous_try_count() {
    # 检测之前尝试了多少次
    # 返回: 0, 1, 2, ...
}

should_clean_and_restart() {
    # 判断是否应该清空重算
    # 返回: "clean|reason" 或 "continue|reason"
}
```

### 8.2 续算脚本结构

```bash
#!/bin/bash
# Resume script for hmat_0
# Generated at 2025-12-25 15:30:45
# Resuming from: Relax
# Stages to run: Relax, Scf, Band, Dos
# Previous logs backed up: stat0.log, time0.log

# ===== 环境初始化 =====
# (与原始脚本相同)

# ===== 续算辅助函数 =====
detect_resume_number() { ... }
detect_previous_try_count() { ... }
should_clean_and_restart() { ... }

# ===== 日志初始化 =====
stat_log=stat.log
time_log=time.log
echo '' >> $stat_log
echo '# === Resume from Relax at 2025-12-25 15:30:45 ===' >> $time_log

# ===== Relax 阶段 =====
if [ -d Relax ]; then
    # 智能处理现有目录
    RESUME_NUM=$(detect_resume_number Relax)
    DECISION=$(should_clean_and_restart Relax Relax)
    
    if [ "$ACTION" = "clean" ]; then
        rm -rf Relax
        mkdir Relax && cd Relax
        RESUME_START_TRY=0
    else
        cd Relax
        PREV_TRY_COUNT=$(detect_previous_try_count .)
        RESUME_START_TRY=$PREV_TRY_COUNT
        
        # 备份当前文件
        cp INPUT INPUT_R${RESUME_NUM}_pre
        cp KPT KPT_R${RESUME_NUM}_pre
        cp STRU STRU_R${RESUME_NUM}_pre
        
        # 使用优化结构
        STRU_ION_D=$(find OUT.* -name "STRU_ION_D" | head -1)
        cp "$STRU_ION_D" STRU
    fi
else
    mkdir Relax && cd Relax
    RESUME_START_TRY=0
fi

# 生成输入文件
python /path/to/abacus.py generate --work_dir . --stage Relax

# 循环尝试
END_TRY=$((RESUME_START_TRY + 2))
for try_num in $(seq $RESUME_START_TRY $((END_TRY - 1))); do
    # 运行 ABACUS
    yhrun -N 1 -n 16 -p deimos /path/to/abacus ...
    
    # 检查结果
    python /path/to/abacus.py errors --work_dir .
    python /path/to/abacus.py converge --work_dir .
    
    if [ -f "converge.txt" ]; then
        break
    fi
    
    # 备份（续算模式）
    if [ $try_num -lt $((END_TRY - 1)) ]; then
        cp INPUT INPUT_R${RESUME_NUM}_${try_num}
        cp KPT KPT_R${RESUME_NUM}_${try_num}
        cp STRU STRU_R${RESUME_NUM}_${try_num}
        
        python /path/to/abacus.py update --work_dir . --try_num $try_num --stage Relax
    fi
done

# 错误处理
if [ ! -f 'converge.txt' ]; then
    if [ -f 'error.txt' ]; then
        echo 'Relax failed' >> ../stat.log
        exit 1
    fi
fi

cd ..

# ===== Scf 阶段 =====
# (与 Relax 类似)

# ===== Band 阶段 =====
# (与 Relax 类似)

# ===== Dos 阶段 =====
# (与 Relax 类似)

# ===== 工作流总结 =====
# (与原始脚本相同)
```

### 8.3 Python API

#### 命令行接口

```python
@cli.command()
@click.argument('work_dir', type=click.Path(exists=True))
@click.option('--from-stage', help='强制从指定阶段开始')
@click.option('--dry-run', is_flag=True, help='仅显示续算计划')
@click.option('--no-backup', is_flag=True, help='不备份日志')
@click.option('--no-clean', is_flag=True, help='不清理错误')
def resume(work_dir, from_stage, dry_run, no_backup, no_clean):
    """智能续算工作流"""
    # ... 实现 ...
```

#### 核心函数

```python
# abacus/resume_utils.py

def detect_stage_status(stage_dir, stage_name):
    """检测阶段状态"""
    # 返回: 'success', 'ignored', 'failed', 'incomplete', 'not_started'

def parse_stat_log(work_dir):
    """解析 stat.log"""
    # 返回: {'Test_spin': 'success', 'Relax': 'failed', ...}

def determine_resume_point(work_dir, workflow_stages):
    """确定续算起点"""
    # 返回: (resume_stage_index, reason) 或 (None, reason)

def backup_logs(work_dir):
    """备份日志"""
    # 返回: (stat_backup_name, time_backup_name)

def clean_stage_errors(work_dir, start_stage):
    """清理错误标记"""
    # 删除 error.txt 和 stat.log 中的失败记录

def detect_resume_number(stage_dir):
    """检测续算次数"""
    # 返回: 0, 1, 2, ...

def detect_previous_try_count(stage_dir):
    """检测之前的尝试次数"""
    # 返回: 0, 1, 2, ...

def should_clean_and_restart(stage_dir, stage_name):
    """判断续算策略"""
    # 返回: (should_clean: bool, reason: str)
```

```python
# abacus/submit.py

class AbacusFlowManager:
    def generate_resume_script(self, work_dir, start_stage, clean_errors=True, no_backup=False):
        """生成续算脚本"""
        # 1. 备份日志
        # 2. 清理错误
        # 3. 生成续算脚本
        # 返回: 脚本路径
    
    def _create_workflow_script(self, job_dir, job_name, stages, is_resume=False, backup_info=None):
        """统一的脚本生成逻辑"""
        # 完整工作流和续算工作流共用此方法
        # 返回: 脚本内容字符串
    
    def _generate_stage_script(self, stage, abacus_py, abacus_bin, monitor_py, is_resume=False):
        """生成单个阶段的脚本"""
        # 返回: 脚本行列表
    
    def _get_resume_helper_functions(self):
        """返回 Bash 辅助函数"""
        # 返回: 脚本行列表
```

---

## 9. 未来扩展

### 9.1 计划功能

- [ ] Web UI 可视化续算控制
- [ ] 续算历史记录数据库
- [ ] 自动参数优化建议
- [ ] 续算性能分析报告
- [ ] 远程续算触发（API）

### 9.2 改进方向

- [ ] 更智能的续算策略（机器学习）
- [ ] 分布式续算支持
- [ ] 更详细的续算日志
- [ ] 续算失败原因分析
- [ ] 自动生成续算建议

---

## 10. 版本历史

### v1.0 (2025-12-25)

- ✅ 初始版本
- ✅ 智能状态检测
- ✅ 自动续算起点确定
- ✅ 日志备份系统
- ✅ 两种续算策略（清空/继续）
- ✅ 续算编号系统
- ✅ Dry-run 模式
- ✅ 完整文档

---

## 11. 参考资料

- [ABACUS 官方文档](http://abacus.ustc.edu.cn/)
- [工作流配置说明](../README.md)
- [输入文件说明](../config/template/README.md)
- [错误处理指南](./error_handling.md)

---

## 12. 贡献指南

欢迎提交 Issue 和 Pull Request！

### 报告 Bug

请提供以下信息：
- 工作目录结构
- `stat.log` 内容
- 续算命令和输出
- 期望行为 vs 实际行为

### 功能建议

请描述：
- 使用场景
- 预期功能
- 可能的实现方案

---

## 联系方式

- GitHub: [abacusflow](https://github.com/username/abacusflow)
- Email: support@abacusflow.org

---

**Happy Computing! 🚀**


