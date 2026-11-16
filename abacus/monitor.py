#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ABACUS 计算实时监控脚本
实时提取并显示离子步和电子步的关键信息
"""

import sys
import re
from datetime import datetime

def format_time():
    """格式化当前时间"""
    return datetime.now().strftime("%H:%M:%S")

def monitor_abacus_output():
    """
    实时监控 ABACUS 输出，提取关键信息
    - 电子步迭代（ITER / ELEC）
    - 离子步信息（ION）
    - 总能量（ETOT / E_KohnSham）
    - 收敛信息（Density error, Energy diff）
    支持 PW 和 LCAO 两种基组
    """
    current_scf_step = None
    current_ion_step = None
    last_energy = None
    in_scf_iteration = False
    
    # LCAO 相关变量
    current_lcao_ion = None
    current_lcao_elec = None
    lcao_energy_ev = None
    lcao_prev_energy = None
    lcao_density_error = None
    lcao_magnetism = None
    in_lcao_block = False
    scf_threshold = 1e-7  # 默认 SCF 收敛阈值
    
    # 用于跟踪是否是 Relax 计算
    is_relax = False
    
    try:
        for line in sys.stdin:
            line_stripped = line.strip()
            
            # 检测是否是 Relax 计算
            if 'RELAX' in line or 'Ion relaxation' in line or 'relax_nmax' in line:
                is_relax = True
            
            # ============ LCAO 格式解析 ============
            # 匹配 "LCAO ALGORITHM --------------- ION=   1  ELEC=   7"
            lcao_match = re.search(r'LCAO ALGORITHM.*?ION=\s*(\d+)\s+ELEC=\s*(\d+)', line)
            if lcao_match:
                ion_num = int(lcao_match.group(1))
                elec_num = int(lcao_match.group(2))
                
                # 如果是新的 ION 步
                if current_lcao_ion is None or ion_num != current_lcao_ion:
                    current_lcao_ion = ion_num
                    print(f"\n{'='*70}")
                    print(f"[{format_time()}] 🔄 ION STEP {ion_num}")
                    print(f"{'='*70}")
                
                current_lcao_elec = elec_num
                in_lcao_block = True
                lcao_density_error = None
                lcao_magnetism = None
                lcao_energy_ev = None
                sys.stdout.flush()
                continue
            
            # 在 LCAO block 中解析信息
            if in_lcao_block:
                # 磁矩
                mag_match = re.search(
                    r'total\s+magnetism\s+\(Bohr\s+mag/cell\)\s*=\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)',
                    line, re.IGNORECASE)
                if mag_match:
                    lcao_magnetism = float(mag_match.group(1))
                
                # Density error（这是收敛判据）
                density_match = re.search(r'Density error is\s+([-+]?[\d.eE+-]+)', line)
                if density_match:
                    lcao_density_error = float(density_match.group(1))
                
                # E_KohnSham (eV)
                energy_match = re.search(r'E_KohnSham\s+([-+]?[\d.]+)\s+([-+]?[\d.]+)', line)
                if energy_match:
                    lcao_energy_ev = float(energy_match.group(2))
                
                # LCAO block 结束
                if line.strip().startswith('---') and lcao_energy_ev is not None:
                    in_lcao_block = False
                    
                    # 计算能量变化
                    energy_diff = None
                    if lcao_prev_energy is not None:
                        energy_diff = lcao_energy_ev - lcao_prev_energy
                    
                    # 判断是否需要显示（每5步或密度误差小于1e-4时）
                    should_display = (
                        current_lcao_elec % 5 == 0 or 
                        current_lcao_elec <= 2 or
                        (lcao_density_error is not None and lcao_density_error < 1e-4)
                    )
                    
                    if should_display:
                        # 显示这个电子步的信息
                        print(f"  ELEC{current_lcao_elec:3d}  E={lcao_energy_ev:>14.6f} eV", end="")
                        
                        if energy_diff is not None:
                            print(f"  ΔE={energy_diff:>12.6f} eV", end="")
                        
                        if lcao_density_error is not None:
                            # 显示密度误差（这是收敛判据，与 scf_thr 比较）
                            print(f"  ΔDens={lcao_density_error:.3e}", end="")
                            
                            # 与 scf_thr 比较
                            if lcao_density_error < scf_threshold:
                                print(f"  ✅converged", end="")
                            elif lcao_density_error < 1e-4:
                                print(f"  🔸approaching", end="")
                        
                        if lcao_magnetism is not None and abs(lcao_magnetism) > 1e-4:
                            print(f"  mag={lcao_magnetism:.6f}", end="")
                        
                        print()  # 换行
                        sys.stdout.flush()
                    
                    # 更新上一步能量
                    lcao_prev_energy = lcao_energy_ev
                    continue
            
            # ============ 离子步信息（PW格式） ============
            # 匹配 "STEP OF ION RELAXATION : 3"
            ion_match = re.search(r'STEP OF ION RELAXATION\s*:\s*(\d+)', line)
            if ion_match:
                current_ion_step = int(ion_match.group(1))
                if current_lcao_ion is None:  # 避免重复显示
                    print(f"\n{'='*70}")
                    print(f"[{format_time()}] 🔄 ION STEP {current_ion_step}")
                    print(f"{'='*70}")
                    sys.stdout.flush()
            
            # ============ SCF 迭代信息 ============
            # 匹配 "ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)"
            if 'ITER' in line and 'ETOT' in line and 'EDIFF' in line:
                in_scf_iteration = True
                if current_ion_step is None:
                    print(f"\n{'='*70}")
                    print(f"[{format_time()}] ⚡ SCF Iterations")
                    print(f"{'='*70}")
                print(f"  {line_stripped}")
                sys.stdout.flush()
                continue
            
            # 匹配 SCF 迭代行
            # 格式1: "CG1     -2.169e+03  0.000e+00  3.052e-01  2.1"
            # 格式2: " GE1      6.48e-01   6.48e-01  -1.28908902e+03   0.00000000e+00   2.2072e-01  14.62"
            if in_scf_iteration:
                # 格式2（带 TMAG AMAG）
                scf_line_match2 = re.match(r'^\s*(CG|DA|GE|GC|BFGS|SD)(\d+)\s+([-+]?\d+\.\d+e[+-]\d+)\s+([-+]?\d+\.\d+e[+-]\d+)\s+([-+]?\d+\.\d+e[+-]\d+)\s+([-+]?\d+\.\d+e[+-]\d+)\s+([-+]?\d+\.\d+e[+-]\d+)\s+([\d.]+)', line_stripped)
                if scf_line_match2:
                    method = scf_line_match2.group(1)
                    step = scf_line_match2.group(2)
                    tmag = scf_line_match2.group(3)
                    amag = scf_line_match2.group(4)
                    etot = scf_line_match2.group(5)
                    ediff = scf_line_match2.group(6)
                    drho = scf_line_match2.group(7)
                    time_s = scf_line_match2.group(8)
                    
                    current_scf_step = int(step)
                    last_energy = float(etot)
                    
                    # 格式化输出（每5步显示一次，或者收敛时显示）
                    ediff_val = float(ediff)
                    drho_val = float(drho)
                    if current_scf_step % 5 == 0 or current_scf_step <= 2 or abs(drho_val) < 1e-4:
                        print(f"  {method}{step:>3}  E={etot:>14}eV  ΔE={ediff:>12}  ΔDens={drho:>10}  t={time_s:>6}s", end="")
                        if abs(drho_val) < 1e-7:
                            print(f"  ✅", end="")
                        elif abs(drho_val) < 1e-4:
                            print(f"  🔸", end="")
                        print()
                        sys.stdout.flush()
                    continue
                
                # 格式1（不带 TMAG AMAG）
                scf_line_match = re.match(r'^\s*(CG|DA|GE|GC|BFGS|SD)(\d+)\s+([-+]?\d+\.\d+e[+-]\d+)\s+([-+]?\d+\.\d+e[+-]\d+)\s+([-+]?\d+\.\d+e[+-]\d+)\s+([\d.]+)', line_stripped)
                if scf_line_match:
                    method = scf_line_match.group(1)
                    step = scf_line_match.group(2)
                    etot = scf_line_match.group(3)
                    ediff = scf_line_match.group(4)
                    drho = scf_line_match.group(5)
                    time_s = scf_line_match.group(6)
                    
                    current_scf_step = int(step)
                    last_energy = float(etot)
                    
                    # 格式化输出（每5步显示一次，或者 ediff 很小时显示）
                    ediff_val = float(ediff)
                    if current_scf_step % 5 == 0 or abs(ediff_val) < 1e-5:
                        print(f"  {method}{step:>3}  E={etot:>12}eV  ΔE={ediff:>10}  ρ={drho:>10}  t={time_s}s")
                        sys.stdout.flush()
                    continue
                
                # SCF 结束
                if line_stripped.startswith('---') or 'convergence' in line_stripped.lower():
                    in_scf_iteration = False
            
            # ============ 收敛信息 ============
            if 'charge density convergence is achieved' in line:
                print(f"  ✅ SCF converged!")
                if last_energy is not None:
                    print(f"     Final Energy: {last_energy:.6e} eV")
                sys.stdout.flush()
            
            if 'convergence has NOT been achieved' in line:
                print(f"  ⚠️  SCF NOT converged")
                sys.stdout.flush()
            
            # ============ 力和应力信息（Relax） ============
            if is_relax and 'LARGEST GRAD' in line:
                force_match = re.search(r'LARGEST GRAD\s*\(eV/Angstrom\)\s*=\s*([\d.]+)', line)
                if force_match:
                    max_force = float(force_match.group(1))
                    print(f"  🔧 Max Force: {max_force:.6f} eV/Å")
                    sys.stdout.flush()
            
            # ============ Relax 收敛信息 ============
            if 'Relaxation is converged!' in line:
                print(f"\n{'='*70}")
                print(f"[{format_time()}] ✅ RELAXATION CONVERGED!")
                print(f"{'='*70}")
                sys.stdout.flush()
            
            if 'Relaxation is not converged yet' in line:
                print(f"  ⚠️  Relaxation not converged yet")
                sys.stdout.flush()
            
            if 'Relaxation is converged, but the SCF is unconverged' in line:
                print(f"\n{'='*70}")
                print(f"[{format_time()}] ⚠️  RELAXATION CONVERGED (but SCF unconverged - unreliable!)")
                print(f"{'='*70}")
                sys.stdout.flush()
            
            # ============ 总时间信息 ============
            if 'Total  Time' in line or 'Total time' in line:
                time_match = re.search(r':\s*([\d.]+)', line)
                if time_match:
                    total_time = float(time_match.group(1))
                    print(f"\n{'='*70}")
                    print(f"[{format_time()}] ⏱️  Calculation Finished - Total Time: {total_time:.2f} s")
                    print(f"{'='*70}")
                    sys.stdout.flush()
            
            # ============ 错误信息 ============
            if 'ERROR' in line or 'Error' in line or 'NOTICE' in line:
                if 'NOTICE' not in line or '!!!' in line:
                    print(f"\n{'='*70}")
                    print(f"[{format_time()}] ❌ ERROR DETECTED:")
                    print(f"  {line_stripped}")
                    print(f"{'='*70}")
                    sys.stdout.flush()
    
    except KeyboardInterrupt:
        print(f"\n[{format_time()}] Monitoring stopped by user")
    except Exception as e:
        print(f"\n[{format_time()}] Monitor error: {e}", file=sys.stderr)

if __name__ == '__main__':
    monitor_abacus_output()

