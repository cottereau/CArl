#!/bin/bash
#PBS -l walltime=3:59:00
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -P omaha
#PBS -q defaultq

module purge
module load intel-compilers/2018.2
module load intel-mpi/2018.2  
module load intel-mkl/2018.2
module load libgmp/6.1.2
module load libmpfr/3.1.5
module load cmake/3.14.4 

export PETSC_DIR="/home/gattif/srclib/petsc"
export PETSC_ARCH="arch-linux2-cxx-opt"
export CGAL_DIR="/home/gattif/LOCAL/cgal"
export LIBMESH_DIR="/home/gattif/LOCAL/libmesh"
export BOOST_ROOT="/home/gattif/LOCAL/boost"
cd $PBS_O_WORKDIR
