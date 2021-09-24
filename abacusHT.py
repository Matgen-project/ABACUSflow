#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 导入模块
import os,sys
from math import ceil, pi
import numpy as np
from ase.io import read
from ase.io.abacus import write_input_stru
from ase.calculators.abacus.create_input import AbacusInput
args = sys.argv



# 指定结构文件、赝势库、轨道库路径，以及K点密度相关的KSPACING值
STRUPATH= args[1]
POTPATH='/WORK/nscc-gz_material_1/ICSD_vasp/abacus/abacus/AbacusHighThroughput/PotSG15/'
ORBPATH='/WORK/nscc-gz_material_1/ICSD_vasp/abacus/abacus/AbacusHighThroughput/OrbSG15std/'
KSPACING=0.13
CALC_PATH = args[2]
PRESENT_DIR = os.getcwd()
if not os.path.exists(CALC_PATH):
    os.makedirs(CALC_PATH)

for root, dirs, files in os.walk(top=STRUPATH, topdown=True):
    for file in files:
        # 读取vasp格式的结构文件，获取结构元素种类，结构文件名字
        stru = read(os.path.join(root, file), format='vasp')
        ntype = len(set(stru.get_chemical_symbols()))
        atom_file = file.split('.')[0]
        
        # 根据结构文件名字创建abacus运行文件夹
        cpth = os.path.join(CALC_PATH, atom_file)
        print(cpth) 
        if not os.path.exists(cpth):
            os.makedirs(cpth)
        os.chdir(cpth)
        
        # 在abacus运行文件夹下创建INPUT
        input = AbacusInput()
        input.set(atom_file=atom_file,
                  ntype=ntype,
                  kpoint_file='KPT',
                  pseudo_dir='./',
                  calculation='scf',
                  nspin=2,
                  ecutwfc=60,
                  dr2=1e-06,
                  ks_solver='genelpa',
                  niter=100,
                  basis_type='lcao',
                  smearing='gauss',
                  sigma=0.002,
                  mixing_type='pulay',
                  mixing_beta=0.4,
                  symmetry=1,
                  gamma_only=0)
        input.write_input_input()
        
        # 在abacus运行文件夹下创建KPT
        Kpoints = [int(ceil(2 * pi / KSPACING * np.linalg.norm(stru.get_reciprocal_cell()[i]))) for
                  i in range(3)]
        Kpoints += [0, 0, 0]
        input.set(knumber=0,
                  kmode='Gamma',
                  kpts=Kpoints)
        input.write_input_kpt()

        # 在abacus运行文件夹下创建abacus结构文件(与vasp结构文件同名)，
        # 并复制相应元素赝势和轨道到运行文件夹下
        write_input_stru(stru=stru,
                        pseudo_dir=POTPATH,
                        basis_dir=ORBPATH,
                        potential_name='PotSG15',
                        basis_name='SG15std',
                        coordinates_type='Direct',
                        spin=2,
                        filename=atom_file)
        
        # 复制abacus运行脚本并提交到服务器
        os.system('cp {}/abacus.sh ./'.format(PRESENT_DIR))
        os.system('yhbatch -N 1 abacus.sh')
        os.chdir(PRESENT_DIR)
