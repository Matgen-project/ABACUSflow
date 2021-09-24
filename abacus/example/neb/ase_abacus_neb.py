#!/usr/bin/env python 
# -*- coding: utf-8 -*-

#PBS -N job_name
#PBS -o job.log
#PBS -e job.err
#PBS -q gold5120
#PBS -l nodes=1:ppn=12
#PBS -l walltime=22:00:00

import os
#os.system("module load intelmpi/2018.1.163")
#os.system("export OMP_NUM_THREADS=1")
#os.system("cd  $PBS_O_WORKDIR")
from ase import Atoms,Atom

from ase.calculators.abacus.abacus_out import *
from ase.calculators.abacus.create_input import *

from ase.io import write
from ase.optimize import MDMin
from ase.neb import NEB 
#import matplotlib.pyplot as plt 
from ase.neb import NEBTools

os.environ['ASE_ABACUS_COMMAND']="mpirun -machinefile $PBS_NODEFILE -np 12 /home/shenzx/software/abacus/abacus_v2.0/bin/ABACUS.mpi.2.0>>PREFIX.log"
#Create a structure
al110_cell = (4.0614,2.8718,2.8718)
al110_positions = [(0.0,0.0,0.0),
                   (2.0307,1.4359,-1.4359)]
initial = Atoms('Al2',
                positions=al110_positions,
                cell=al110_cell,
                pbc=(1,1,0))
#initial *= (2, 2, 1)
initial.append(Atom('C', (2.0307, 1.4359, 4.3077)))
initial.center(vacuum=4.0, axis=2)

final = initial.copy()
final.positions[-1][1] += 2.8718
final.positions[-1][0] += 4.0614

# Construct a list of images:
images = [initial]
for i in range(3):
    images.append(initial.copy())
images.append(final)
print('Create images successfully')

for image in images:
    # Let all images use an abacus calculator:
    image.set_calculator(Abacus(
          label = "/home/shenzx/project/python_20190718/ase_20191211/example/neb/AlC/neb",
          atoms=image,
          pseudo_dir = "/home/shenzx/software/abacus/SG15_ONCV_PBE_1.0",
          potential_name = "potential_pbe_sg15" , 
          basis_dir = "/home/shenzx/software/abacus/Orb_DZP_E100_Standard_v1.0" ,
          basis_name = [ "Al_gga_9au_100Ry_4s4p1d.orb", "C_gga_8au_100Ry_2s2p1d.orb" ] ,
          calculation='scf',
          ntype=2,
          nbands=20,
          ecutwfc=50,
          dr2="1.0e-6",
          niter=100,
          force=1,
          smearing='gaussian',
          sigma=0.02,
          mixing_type='pulay-kerker',
          ks_solver = 'genelpa',
          mixing_beta=0.4,
          basis_type='lcao',
          atom_file='STRU',
          gamma_only = 0,
          knumber = 0,    
          kmode = 'Gamma',    
          kpts = [ 1, 1, 1, 0, 0, 0]
          ))

# Create a Nudged Elastic Band:
neb = NEB(images)
# Make a starting guess for the minimum energy path (a straight line
# from the initial to the final state):
neb.interpolate()
# Relax the NEB path:
minimizer = MDMin(neb)
print('Relax the NEB path, please wait!!!')
minimizer.run(fmax=0.05)
# Write the path to a trajectory:
write('neb.traj', images)
print('Write the path to a trajectory successfully')

#images = read('neb.traj@:')
#view(images1)
#view(images2)
#view(images)

nebtools = NEBTools(images)

# Get the calculated barrier and the energy change of the reaction.
Ef, dE = nebtools.get_barrier()
print(Ef,'   ',dE)
# Get the barrier without any interpolation between highest images.
Ef, dE = nebtools.get_barrier(fit=False)
print(Ef,'   ',dE)
# Get the actual maximum force at this point in the simulation.
max_force = nebtools.get_fmax()
print('max_force:   ',max_force)
# Create a figure like that coming from ASE-GUI.
#print('Create a figure,wait !!!')
#fig = nebtools.plot_band()
#fig.savefig('diffusion-barrier.png')

# Create a figure with custom parameters.
#fig = plt.figure(figsize=(5.5, 4.0))
#ax = fig.add_axes((0.15, 0.15, 0.8, 0.75))
#nebtools.plot_band(ax)
#fig.savefig('diffusion-barrier.png')
#print('Create a figure successfully')
