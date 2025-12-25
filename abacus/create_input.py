# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 17:02:37 2018

@author: shenzx
"""
from __future__ import print_function
import warnings
import shutil
from os.path import join
import numpy as np
# 修复导入：ASE 3.26 中没有 ase.io.abacus，使用本地 abacus/stru_utils.py 模块中的函数
import os
import sys
# 从 abacus.stru_utils 导入 STRU 写入函数
from abacus.stru_utils import write_input_stru
from ase.calculators.calculator import FileIOCalculator
# copyright © Key Lab of Quantum Information, CAS, China
"""This module defines an ASE interface to ABACUS

Developed on the basis of modules by Zhen-Xiong Shen.
 The path of the directory containing the
 pseudopotential and basis directories (LDA, PBE, SG15, ORBITAL, ...)
 should be set by the enviromental flag $ABACUS_PP_PATH, $ABACUS_ORBITAL_PATH.

The user should also set the enviroment flag
 $ABACUS_SCRIPT pointing to a python script looking

like::
    import os
    exitcode = os.system('abacus')
http://abacus.ustc.edu.cn/
"""

# Parameters list that can be set in INPUT.  -START-
# 1
general_keys = [
        'suffix',           # the name of main output directory
        'latname',          # the name of lattice name
        'atom_file',        # the filename of file containing atom positions
        'kpoint_file',      # the name of file containing k points
        'pseudo_dir',       # the directory containing pseudo files
        'orbital_dir',      # the directory containing orbital basis files
        'pseudo_type',      # the type pseudo files
        'dft_functional',   # exchange correlation functional
        'calculation',      # test; scf; relax; nscf; ienvelope; istate;
        'ntype',            # atom species number
        'nspin',            # 1: single spin; 2: up and down spin;
        'nbands',           # number of bands
        'nbands_istate',    # number of bands around Fermi level for istate calulation
        'symmetry',         # turn symmetry on or off
        'kspacing',         # spacing between k points
        'nelec'             # input number of electrons
        ]
# 2
pw_keys = [
        'ecutwfc',         # energy cutoff for wave functions
        'ethr',            # threshold for eigenvalues is cg electron iterations
        'dr2',             # charge density error
        'scf_thr',         # SCF convergence threshold (alias for dr2)
        'scf_nmax',        # maximum number of SCF iterations (alias for niter)
        'start_wfc',       # start wave functions are from 'atomic' or 'file'
        'start_charge',    # start charge is from 'atomic' or file
        'charge_extrap',   # atomic; first-order; second-order; dm:coefficients of SIA
        'out_charge',      # >0 output charge density for selected electron steps
        'out_chg',         # output charge density (alias for out_charge)
        'init_chg',        # read charge density from file
        'out_potential',   # output realspace potential
        'out_wf',          # output wave functions
        'out_dos',         # output energy and dos
        'out_band',        # output energy and band structure
        'nx',              # number of points along x axis for FFT grid
        'ny',              # number of points along y axis for FFT grid
        'nz'               # number of points along z axis for FFT grid
        ]
# 3
relaxation_keys = [
        'ks_solver',       # cg; david; lapack; genelpa; hpseps;
        'niter',           # number of electron iterations
        'vna',             # use the vna or not
        'grid_speed',      # 1:normal 2:fast
        'force_set',       # output the force_set or not
        'force',           # calculate the force
        'nstep',           # number of ion iteration steps
        'out_stru',        # output the structure files after each ion step
        'force_thr',       # force threshold,  unit: Ry/Bohr
        'force_thr_ev',    # force threshold,  unit: eV/Angstrom
        'force_thr_ev2',   # force invalid threshold,  unit: eV/Angstrom
        'stress_thr',      # stress threshold
        'press1',          # target pressure, unit: KBar
        'press2',          # target pressure, unit: KBar
        'press3',          # target pressure, unit: KBar
        'bfgs_w1',         # wolfe condition 1 for bfgs
        'bfgs_w2',         # wolfe condition 2 for bfgs
        'trust_radius_max', # maximal trust radius,  unit: Bohr
        'trust_radius_min', # minimal trust radius,  unit: Bohr
        'trust_radius_ini', # initial trust radius,  unit: Bohr
        'stress',           # calculate the stress or not
        'fixed_axes',       # which axes are fixed
        'move_method',      # bfgs; sd; cg; cg_bfgs;
        'relax_method',     # alias for move_method
        'relax_nmax',       # maximum number of ionic steps
        'out_level',        # ie(for electrons); i(for ions);
        'out_dm'            # >0 output density matrix
        ]
# 4
lcao_keys = [
        'basis_type',       # PW; LCAO in pw; LCAO
        'search_radius',    # input search radius (Bohr)
        'search_pbc',       # input periodic boundary condition
        'lcao_ecut',        # energy cutoff for LCAO
        'lcao_dk',          # delta k for 1D integration in LCAO
        'lcao_dr',          # delta r for 1D integration in LCAO
        'lcao_rmax',        # max R for 1D two-center integration table
        'out_hs',           # output H and S matrix
        'out_lowf',         # ouput LCAO wave functions
        'bx',               # division of an element grid in FFT grid along x
        'by',               # division of an element grid in FFT grid along y
        'bz'                # division of an element grid in FFT grid along z
        ]
# 5
smearing_keys = [
        'smearing',         # type of smearing: gauss; fd; fixed; mp; mp2
        'smearing_method',  # type of smearing (alias for smearing): gaussian, gauss, fd, fixed, mp, mp2
        'sigma',            # energy range for smearing
        'smearing_sigma'    # energy range for smearing (alias for sigma)
        ]
# 6
charge_mixing_keys = [
        'mixing_type',      # plain; kerker; pulay; pulay-kerker
        'mixing_beta',      # mixing parameter: 0 means no new charge
        'mixing_ndim',      # mixing dimension in pulay
        'mixing_gg0'        # mixing parameter in kerker
        ]
# 7
dos_keys = [
        'dos_emin_ev',      # minimal range for dos
        'dos_emax_ev',      # maximal range for dos
        'dos_edelta_ev',    # delta energy for dos
        'dos_sigma'         # gauss b coefficeinet(default=0.07)        
        ]
# 8
technique_keys = [
        'gamma_only',       # gamma only
        'diago_proc',       # number of proc used to diago
        'npool',            # number of pools for k points,  pw only
        'sparse_matrix',    # use sparse matrix, in DMM
        'atom_distribution', # distribute atoms, in DMM
        'mem_saver',         # memory saver for many k points used
        'printe'             # print band energy for selectively ionic steps
        ]
# 9
siao_keys = [
        'selinv_npole',     # number of selected poles
        'selinv_temp',      # temperature for Fermi-Dirac distribution
        'selinv_gap',       # supposed gap in the calculation
        'selinv_deltae',    # expected energy range
        'selinv_mu',        # chosen mu as Fermi energy
        'selinv_threshold', # threshold for calculated electron number
        'selinv_niter',     # max number of steps to update mu
        ]
# 10
molecular_dynamics_keys = [
        'md_mdtype',        # choose ensemble
        'md_dt',            # time step
        'md_nresn',         # parameter during integrater
        'md_nyosh',         # parameter during integrater
        'md_qmass',         # mass of thermostat
        'md_tfirst',        # temperature first
        'md_tlast',         # temperature last
        'md_dumpmdfred',    # The period to dump MD information for monitoring and restarting MD
        'md_mdoutpath',     # output path of md
        'md_domsd',         # whether compute <r(t)-r(0)>
        'md_domsdatom',     # whether compute msd for each atom
        'md_rstmd',         # whether restart
        'md_fixtemperature', # period to change temperature
        'md_ediff',          # parameter for constraining total energy change
        'md_ediffg',         # parameter for constraining max force change
        'md_msdstarttime'    # choose which step that msd be calculated
        ]
# 11
efield_keys = [
        'efield',           # add electric field
        'edir',             # add electric field
        'emaxpos',          # maximal position of efield [0, 1
        'eopreg',           # where sawlike potential decrease
        'eamp',             # amplitute of the efield,  unit is a.u.
        'eamp_v'            # amplitute of the efield,  unit is V/A
        ]
# 12
bfield_keys = [
        'bfield',           # add magnetic field
        'bfield_teslax',    # magnetic field strength
        'bfield_teslay',    # magnetic field strength
        'bfield_teslaz',    # magnetic field strength
        'bfield_gauge_x',   # magnetic field gauge origin
        'bfield_gauge_y',   # magnetic field gauge origin
        'bfield_gauge_z'    # magnetic field gauge origin
        ]
# 13
test_keys = [
        'out_alllog',       # output information for each processor,  when parallel
        'nurse',            # for coders  
        'colour',           # for coders,  make their live colourful
        't_in_h',           # calculate the kinetic energy or not
        'vl_in_h',          # calculate the local potential or not
        'vnl_in_h',         # calculate the nonlocal potential or not
        'zeeman_in_h',      # calculate the zeeman term or not
        'test_force',       # test the force
        'test_stress'       # test the force
        ]
# 14
other_methods_keys = [
        'mlwf_flag',        # turn MLWF on or off
        'opt_epsilon2',     # calculate the dielectic function
        'opt_nbands'        # number of bands for optical calculation
        ]
# 15
vdw_d2_keys = [
        'vdwD2',            # calculate vdw-D2 or not
        'vdwD2_scaling',    # scaling of vdw-D2
        'vdwD2_d',          # damping parameter
        'vdwD2_C6_file',    # filename of C6
        'vdwD2_C6_unit',    # unit of C6,  Jnm6/mol or eVA6
        'vdwD2_R0_file',    # filename of R0
        'vdwD2_R0_unit',    # unit of R0,  A or Bohr
        'vdwD2_model',      # expression model of periodic structure,  radius or period
        'vdwD2_radius',     # radius cutoff for periodic structure
        'vdwD2_radius_unit', # unit of radius cutoff for periodic structure
        'vdwD2_period'       # periods of periodic structure
        ]
# 16
spectrum_keys = [
        'spectral_type',    # the type of the calculated spectrum
        'spectral_method',  # 0: tddft(linear response)
        'kernel_type',      # the kernel type: rpa,  tdlda ...
        'eels_method',      # 0: hilbert_transform method; 1: standard method
        'absorption_method', # 0: vasp's method  1: pwscf's method
        'system',           # the calculate system
        'eta',              # eta(Ry)
        'domega',           # domega(Ry)
        'nomega',           # nomega
        'ecut_chi',         # the dimension of chi matrix
        'q_start',          # the position of the first q point in direct coordinate
        'q_direction',      # the q direction
        'nq',               # the total number of qpoints for calculation
        'out_epsilon',      # output epsilon or not
        'out_chi',          # output chi or not
        'out_chi0',         # output chi0 or not
        'fermi_level',      # the change of the fermi_level(Ry)
        'coulomb_cutoff',   #  turn on the coulomb_cutoff or not
        'kmesh_interpolation', # calculting <i, 0|j, R>
        'qcar',             # (unit: 2PI/lat0)
        'lcao_box',         # the scale for searching the existence of the overlap <i, 0|j, R>
        'intrasmear',       # Eta
        'shift',            # shift
        'metalcalc',        # metal or not
        'eps_degauss',      # degauss in calculating epsilon0
        'noncolin',         # using non-collinear-spin
        'lspinorb',         # consider the spin-orbit interaction
        'starting_spin_angle'  # starting_spin_angle

        ]
# 17
tddft_keys = [
        'tddft',            # calculate tddft or not
        'td_dr2',           # threshold for electronic iteration of tddft
        'td_dt',            # time of ion step
        'td_force_dt',      # time of force change
        'val_elec_01',      # val_elec_01
        'val_elec_02',      # val_elec_02
        'val_elec_03',      # val_elec_03
        'vext',             # add extern potential or not
        'vext_dire'         # extern potential direction
] 

# Parameters list that can be set in INPUT.  -END-

class AbacusInput(object):
    
    # Initialize internal dictionary of input parameters to None  -START-
    def __init__(self, restart=None):
        """
        self.directory = './'        # shenzx v20200724
        self.stru_filename = 'STRU'  # shenzx v20200724
        self.pseudo_dir = './'   # shenzx v20200724
        self.potential_name = None   # shenzx v20200724
        self.basis_dir = './'        # shenzx v20200724
        self.basis_name = None       # shenzx v20200724
        self.fix = 1                 # shenzx v20200724
        self.coordinates_type = 'Cartesian'      # shenzx v20200724
        """
        self.general_params = {}
        self.pw_params = {}
        self.relaxation_params = {}
        self.lcao_params = {}
        self.smearing_params = {}
        self.charge_mixing_params = {}
        self.dos_params = {}
        self.technique_params = {}
        self.siao_params = {}
        self.molecular_dynamics_params = {}
        self.efield_params = {}
        self.bfield_params = {}
        self.test_params = {}
        self.other_method_params = {}
        self.vdw_d2_params = {}
        self.spectrum_params = {}
        self.tddft_params = {}
        
        for key in general_keys:
            self.general_params[key] = None
        for key in pw_keys:
            self.pw_params[key] = None
        for key in relaxation_keys:
            self.relaxation_params[key] = None
        for key in lcao_keys:
            self.lcao_params[key] = None
        for key in smearing_keys:
            self.smearing_params[key] = None
        for key in charge_mixing_keys:
            self.charge_mixing_params[key] = None
        for key in dos_keys:
            self.dos_params[key] = None
        for key in technique_keys:
            self.technique_params[key] = None
        for key in siao_keys:
            self.siao_params[key] = None
        for key in molecular_dynamics_keys:
            self.molecular_dynamics_params[key] = None
        for key in efield_keys:
            self.efield_params[key] = None
        for key in bfield_keys:
            self.bfield_params[key] = None
        for key in test_keys:
            self.test_params[key] = None
        for key in other_methods_keys:
            self.other_method_params[key] = None
        for key in vdw_d2_keys:
            self.vdw_d2_params[key] = None
        for key in spectrum_keys:
            self.spectrum_params[key] = None
        for key in tddft_keys:
            self.tddft_params[key] = None
        # Initialize internal dictionary of input parameters to None  -END-

        # Appoint the KPT parameters which are not INPUT parameters  -START-
        self.kpt_params = {
                'knumber': 0,           # The number of K points
                'kmode': 'Gamma',       # Mode of K points, can be Gamma, MP, Line, Direct, Cartesian
                'kpts': [1, 1, 1, 0, 0, 0]  # Give the K points
                }
        # Appoint the KPT parameters which are not INPUT parameters  -END-

    # Set the INPUT and KPT parameters  -START-
    def set(self, **kwargs):
        for key in kwargs:
            if key in self.general_params:
                self.general_params[key] = kwargs[key]
            elif key in self.pw_params:
                self.pw_params[key] = kwargs[key]
            elif key in self.relaxation_params:
                self.relaxation_params[key] = kwargs[key]
            elif key in self.lcao_params:
                self.lcao_params[key] = kwargs[key]
            elif key in self.smearing_params:
                self.smearing_params[key] = kwargs[key]
            elif key in self.charge_mixing_params:
                self.charge_mixing_params[key] = kwargs[key]
            elif key in self.dos_params:
                self.dos_params[key] = kwargs[key]
            elif key in self.technique_params:
                self.technique_params[key] = kwargs[key]
            elif key in self.siao_params:
                self.siao_params[key] = kwargs[key]
            elif key in self.molecular_dynamics_params:
                self.molecular_dynamics_params[key] = kwargs[key]
            elif key in self.efield_params:
                self.efield_params[key] = kwargs[key]
            elif key in self.bfield_params:
                self.bfield_params[key] = kwargs[key]
            elif key in self.test_params:
                self.test_params[key] = kwargs[key]
            elif key in self.other_method_params:
                self.other_method_params[key] = kwargs[key]
            elif key in self.vdw_d2_params:
                self.vdw_d2_params[key] = kwargs[key]
            elif key in self.spectrum_params:
                self.spectrum_params[key] = kwargs[key]
            elif key in self.tddft_params:
                self.spectrum_params[key] = kwargs[key]
            elif key in self.kpt_params:
                self.kpt_params[key] = kwargs[key]
            else:
                raise TypeError('Parameter not defined:  ' + key)
    # Set the INPUT and KPT parameters  -END-

    # Write INPUT file  -START-
    def write_input_input(self, directory='./', **kwargs):
        with open(join(directory, 'INPUT'), 'w') as input_file:
            input_file.write('INPUT_PARAMETERS\n')
            input_file.write('# Created by Atomic Simulation Enviroment\n')
            for key, val in self.general_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
            
            for key, val in self.pw_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.relaxation_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.lcao_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.smearing_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.charge_mixing_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.dos_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.technique_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.siao_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.molecular_dynamics_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.efield_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.bfield_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.test_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.other_method_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.vdw_d2_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
                    
            for key, val in self.spectrum_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
            input_file.write('\n')

            for key, val in self.tddft_params.items():
                if val is not None:
                    params = str(key) + ' ' * (20 - len(key)) + str(val)
                    input_file.write('%s\n' % params)
            input_file.write('\n')
    # Write INPUT file  -END-
    # Read  INPUT  file  --START-

    def read_input_input(self,
                         filename='INPUT',
                         directory='./',
                         **kwargs):
        with open(join(directory, filename), 'r') as file:
            file.readline()
            lines = file.readlines()

        for line in lines: 
            try:
                line = line.replace("# ", "#  ")
                data = line.split()
                if len(data) == 0:
                    continue
                elif data[0][0] == "# ":
                    continue
                
                key = data[0]
                if key in general_keys:
                    self.general_params[key] = data[1]
                elif key in pw_keys:
                    self.pw_params[key] = data[1]
                elif key in relaxation_keys:
                    self.relaxation_params[key] = data[1]
                elif key in lcao_keys:
                    self.lcao_params[key] = data[1]
                elif key in smearing_keys:
                    self.smearing_params[key] = data[1]
                elif key in charge_mixing_keys:
                    self.charge_mixing_params[key] = data[1]
                elif key in dos_keys:
                    self.dos_params[key] = data[1]
                elif key in technique_keys:
                    self.technique_params[key] = data[1]
                elif key in siao_keys:
                    self.siao_params[key] = data[1]
                elif key in molecular_dynamics_keys:
                    self.molecular_dynamics_params[key] = data[1]
                elif key in efield_keys:
                    self.efield_params[key] = data[1]
                elif key in bfield_keys:
                    self.bfield_params[key] = data[1]
                elif key in test_keys:
                    self.test_params[key] = data[1]
                elif key in other_methods_keys:
                    self.other_method_params[key] = data[1]
                elif key in vdw_d2_keys:
                    if key == 'vdwD2_period':
                        self.vdw_d2_params[key] = (data[1] + '  '
                        + data[2] + '  ' + data[3])
                    else:
                        self.vdw_d2_params[key] = data[1]
                elif key in spectrum_keys:
                    if key in ['q_start', 'q_direction', 'qcar', 'lcao_box']:
                        self.spectrum_params[key] = (data[1] + '  '
                        + data[2] + '  ' + data[3])
                    else:
                        self.spectrum_params[key] = data[1]
                elif key in tddft_keys:
                    self.tddft_params[key] = data[1]

                return 'ok'  
            
            except  KeyError:
                raise IOError('keyword "%s" in INPUT is'
                                  'not know by calculator.' % key)
                    
            except IndexError:
                raise IOError('Value missing for keyword "%s" .' % key)
    # Read  INPUT  file  --END- 

    # Write KPT  -START-
    def write_input_kpt(self,
                        directory='./',
                        filename='KPT',
                        **kwargs):
        k = self.kpt_params
        if self.technique_params['gamma_only'] is None:
            return warnings.warn(" 'gamma_only' parameter has not been set, "
                                 "please set it to 0 or 1")

        elif self.technique_params['gamma_only'] == 1:
            with open(join(directory, filename), 'w') as kpoint:
                kpoint.write('K_POINTS\n')
                kpoint.write('0\n')
                kpoint.write('Gamma\n')
                kpoint.write('1 1 1 0 0 0')

        elif self.technique_params['gamma_only'] == 0:
            with open(join(directory, filename), 'w') as kpoint:
                kpoint.write('K_POINTS\n')
                kpoint.write('%s\n' % str(k['knumber']))
                kpoint.write('%s\n' % str(k['kmode']))
                if k['kmode'] in ['Gamma', 'MP']:
                    for n in range(len(k['kpts'])):
                        kpoint.write('%s  ' % str(k['kpts'][n]))
                            
                elif k['kmode'] in ['Direct', 'Cartesian', 'Line']:
                    for n in range(len(k['kpts'])):
                        for i in range(len(k['kpts'][n])):
                            kpoint.write('%s  ' % str(k['kpts'][n][i]))
                        kpoint.write('\n')

                else:
                    raise ValueError("The value of kmode is not right, set to "
                                     "Gamma, MP, Direct, Cartesian, or Line.")
        else:
            return warnings.warn("The value of gamma_only is not right, "
                                 "please set it to 0 or 1")
    # Write KPT  -END-

    # Read KPT file  -START-

    def read_kpt(self,
                 filename='KPT',
                 directory='./',
                 **kwargs):
        with open(filename, 'r') as file:
            lines = file.readlines()

        if lines[2][-1] == '\n':
            kmode = lines[2][:-1]
        else:
            kmode = lines[2]

        if kmode in ['Gamma', 'MP']:
            self.kpt_params['kmode'] = lines[2][:-1]
            self.kpt_params['knumber'] = lines[1].split()[0]
            self.kpt_params['kpts'] = np.array(lines[3].split())

        elif kmode in ['Cartesian', 'Direct', 'Line']:
            self.kpt_params['kmode'] = lines[2][:-1]
            self.kpt_params['knumber'] = lines[1].split()[0]
            self.kpt_params['kpts'] = np.array([list(map(float, line.split())) 
                                     for line in lines[3:]])

        else:
            raise ValueError("The value of kmode is not right, set to "
                                     "Gamma, MP, Direct, Cartesian, or Line.")
    # Read KPT file  -END-

    # Write and read potential  -START-
    def write_potential(self,
                        pseudo_dir='./',
                        potential_name=None,
                        directory='./',
                        **kwargs):
        if pseudo_dir == directory:
            return 'It is ok,  pseudo_dir is in work directory'
        else:
            if self.potential_name == None:
                raise ValueError('The value of "potential_name" is not right, '
                                 'please set it to a list')
            else:
                self.pseudo_dir = pseudo_dir
                self.potential_name = potential_name
                for name in self.potential_name:
                    shutil.copyfile(join(self.pseudo_dir, name),
                                    join(directory, name))

    def read_potential(self,
                       pseudo_dir='./',
                       potential_name=None,
                       **kwargs):
        if self.potential_name is None:
            raise ValueError('The value of "potential_name" is not right, '
                             'please set it to a list')
        else:
            self.pseudo_dir = pseudo_dir
            self.potential_name = potential_name
    # Write and read potential  -END-

    # Write and read orbital basis  -START-
    def write_basis(self,
                    basis_dir='./',
                    basis_name=None,
                    directory='./',
                    **kwargs):
        if basis_dir == directory:
            print('It is ok,  basis_dir is in work directory')
        else:
            if self.basis_name is None:
                raise ValueError('The value of "basis_name" is not right, '
                                 'please set it to a list')
            else:
                self.basis_dir = basis_dir
                self.basis_name = basis_name
                for name in self.basis_name:
                    shutil.copyfile(join(self.basis_dir, name),
                                    join(directory, name))

        print("basis_dir = %s, basis_name = %s"%(basis_dir, basis_name))

    def read_basis(self,
                   basis_dir='./',
                   basis_name=None,
                   directory='./',
                   **kwargs):
        if self.basis_name is None:
            raise ValueError('The value of "basis_name" is not right, '
                             'please set it to a list')
        else:
            self.basis_dir = basis_dir
            self.basis_name = basis_name

    # Write and read orbital basis  -START-
    def write_input(self,
                    atoms,
                    properties=None,
                    system_changes=None):
        """Write input parameters to files-file."""

        FileIOCalculator.write_input(self,
                                     atoms,
                                     properties,
                                     system_changes)

        if(system_changes is None):  # shenzx v20200724
            system_changes = '  '    # shenzx v20200724
        if ('numbers' in system_changes or 
                'initial_magmoms' in system_changes):
            self.initialize(atoms)

        write_input_stru(stru=atoms,
                         filename=self.stru_filename,
                         pseudo_dir=self.pseudo_dir,
                         potential_name=self.potential_name,
                         basis_dir=self.basis_dir,
                         basis_name=self.basis_name,
                         fix=self.fix,
                         directory=self.directory,
                         coordinates_type=self.coordinates_type)
        self.write_input_input(directory=self.directory)
        self.write_input_kpt(directory=self.directory)
    # Write all input file -END-


if __name__ == "__main__":
    # Test a writing...
    import os
    print(os.getcwd())
    from ase import Atoms  # just use when test this module
    test = AbacusInput()
    test.set(gamma_only = 1)
    test.write_input(atoms = Atoms("CO"))
