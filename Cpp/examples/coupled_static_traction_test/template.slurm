#!/bin/bash

#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpu_short

module purge
module load git/2.25.0/gcc-9.2.0
module load gmp/6.1.2/gcc-9.2.0
module load mpfr/3.1.6/gcc-9.2.0
module load cmake/3.16.2/gcc-9.2.0
module load openmpi/4.0.2/gcc-9.2.0
module load openblas/0.3.8/gcc-9.2.0
module load netlib-lapack/3.8.0/gcc-9.2.0
module load gcc/9.2.0/gcc-4.8.5

export srclib=/workdir/invsem01/srclib
export locinst=/workdir/invsem01/LOCAL

export BOOST_ROOT=$locinst/boost

export CGAL_DIR=$locinst/cgal

export LIBMESH_DIR=$locinst/libmesh

export PETSC_DIR=$locinst/petsc
export PETSC_ARCH="linux-openmpi-gcc"

export HDF5_DIR=$locinst/hdf5

export CARL_DIR=$locinst/CArl/Cpp/
export CARLBUILD=$locinst/CArl/Cpp/bin

# To compute in the submission directory
cd ${SLURM_SUBMIT_DIR}
