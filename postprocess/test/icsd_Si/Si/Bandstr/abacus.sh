#!/bin/sh
source /WORK/app/toolshs/cnmodule.sh
module load intel-compilers/15.0.1 MPI/Intel/MPICH/3.1-icc15-dyn ELPA/2016.05.004-intel-15 BLAS/3.5.0-icc15 LAPACK/3.5.0-icc15 fftw/3.3.5 scalapack/2.0.2 

export PATH=/WORK/nscc-gz_material_1/gaosen/apps/ABACUS2.1.0/bin:$PATH
export OMP_NUM_THREADS=1

yhrun -N 1 -n 1 ABACUS.mpi.2.1.0
