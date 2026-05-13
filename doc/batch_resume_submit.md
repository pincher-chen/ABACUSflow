# 批量续算与自动提交

批量续算和自动提交可以组合使用，处理大批量失败或未收敛作业的重跑需求。

## 完整流程

```bash
# 1. 查看哪些作业需要续算（不生成脚本）
python abacus.py batch-resume work_cal/ --dry-run

# 2. 生成所有续算脚本（*_resume.sh）
python abacus.py batch-resume work_cal/

# 3. 自动提交续算脚本
python submit_jobs.py work_cal/ --resume
```

步骤 2+3 可以合并为一条命令：

```bash
python abacus.py batch-resume work_cal/ --auto-submit
```

## 常用选项

```bash
# 所有作业从指定阶段开始续算
python abacus.py batch-resume work_cal/ --from-stage Relax

# 后台运行（大批量时推荐）
nohup python abacus.py batch-resume work_cal/ --auto-submit > resume.log 2>&1 &

# 不备份日志文件
python abacus.py batch-resume work_cal/ --no-backup

# 不清理已有的 error.txt
python abacus.py batch-resume work_cal/ --no-clean
```

## 资源配额调整

提交过程中可以直接修改 `config/condor.ini`，系统会在下次检查时自动应用，无需重启：

```ini
[ALLOW]
TOTAL_NODE = 20      # 从 10 调高到 20
MAX_JOBS   = 2000
```

## 手动批量提交（不用自动提交系统）

```bash
# 生成续算脚本
python abacus.py batch-resume work_cal/

# 手动逐个提交
for script in work_cal/*/*_resume.sh; do
    yhbatch "$script"
done
```

## 监控进度

```bash
# 查看提交日志
tail -f logs/submit_*.log

# 查看队列状态
watch -n 10 'squeue -u $USER'
```

## 常见问题

**Q: 中断后如何继续？**  
直接重新运行相同命令，系统会跳过已在队列中的作业。

**Q: 如何只对部分目录续算？**  
对每个子目录分别运行 `python abacus.py resume <subdir>`，再用 `submit_jobs.py` 提交。

**Q: 续算后仍然失败怎么办？**  
检查参数设置（K 点密度、收敛标准、mixing 参数），修改对应的 `config/template/*.yaml` 或
`config/workflow.json` 后重新续算。
