#python abacus.py workflow InputPoscar/ work_cal7/
python abacus.py workflow  batch_poscar/ batch_work4/
python submit_jobs.py batch_work4/
#python abacus.py batch-resume batch_work3/ --dry-run
#python abacus.py batch-resume batch_work3/ --auto-submit

# 使用Scf模板，kval=0.02，自动设置资源
#python abacus.py single batch_poscar/ batch_work_single/ -t Scf -k 0.02

# 手动指定节点和进程数
#python abacus.py single batch_poscar/ batch_work_single/ -t Relax -k 0.04 -n 1 -p 64

# 不自动设置资源，必须手动指定
#python abacus.py single batch_poscar/ work_out/ -t Scf --no-auto-resource -n 2 -p 128
#python submit_jobs.py batch_work_single/ --single
