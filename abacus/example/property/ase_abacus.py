#!/usr/bin/env python 
# -*- coding: utf-8 -*-

#PBS -N job_name
#PBS -o job.log
#PBS -e job.err
#PBS -q gold5120
#PBS -l nodes=1:ppn=8
#PBS -l walltime=22:00:00

import os
os.system("module load intelmpi/2018.1.163")
os.system("export OMP_NUM_THREADS=1")
os.system("cd  $PBS_O_WORKDIR")
import ase
import ase.io as aio
import ase.build as abd
import ase.optimize as aopt
import ase.visualize as av
import ase.constraints as ac
from ase import Atoms
from ase.units import Ry

from ase.calculators.abacus.abacus_out import *
from ase.calculators.abacus.create_input import *

#os.environ['ASE_ABACUS_COMMAND'] = "mpijob /opt/ABACUS/1.0.1dev_intel-mpi-mkl-2017_genELPA/bin/ABACUS>>PREFIX.log"
os.environ['ASE_ABACUS_COMMAND']="mpirun -machinefile $PBS_NODEFILE -np 8 /home/shenzx/software/abacus/abacus_v2.0/bin/ABACUS.mpi.2.0 >> PREFIX.log"

bulk = aio.read( "/home/shenzx/project/python_20190718/ase_20191211/example/property/SiC_mp-8062_conventional_standard.cif",
 format='cif') 

calc = Abacus(  label = "/home/shenzx/project/python_20190718/ase_20191211/example/property/SiC/test",
                atoms = bulk,
                pseudo_dir = "/home/shenzx/software/abacus/SG15_ONCV_PBE_1.0",
                potential_name = "potential_pbe_sg15" , 
                basis_dir = "/home/shenzx/software/abacus/Orb_DZP_E100_Standard_v1.0" ,
                basis_name = [ "Si_gga_8au_100Ry_2s2p1d.orb", "C_gga_8au_100Ry_2s2p1d.orb" ] ,  
                niter = 1000, 
                dr2 = "5.0e-7", 
                ecutwfc = 100, 
                calculation = "scf",
                nspin = 1,
                force = 1,
                ks_solver = 'genelpa', 
                basis_type =  'lcao',  
                gamma_only = 0, 
                knumber = 0,                 
                kmode = 'Gamma',            
                kpts = [ 3, 3, 3, 0, 0, 0],  
                ) 

bulk.set_calculator(calc)

print("Potential Energy: eV")
print(bulk.get_potential_energy())
print("Force:")
print(bulk.get_forces())
print("Fermi Level:")
print(calc.get_fermi_level())

quit()
