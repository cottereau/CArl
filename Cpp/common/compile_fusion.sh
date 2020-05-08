#!/bin/bash
module purge
module load git/2.22.0
module load openmpi/4.0.0
module load hdf5/1.8.20
module load cmake/3.16.0
module load blas/3.8.0 
module load lapack/3.8.0
module load libgmp/6.1.2
module load libmpfr/3.1.5

export locinst=/home/user/LOCAL

export BOOST_ROOT=$locinst/boost
export CGAL_DIR=$locinst/cgal
export LIBMESH_DIR=$locinst/libmesh
export PETSC_DIR=$locinst/petsc
export PETSC_ARCH="linux-openmpi-gcc"

cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -DBOOST_ROOT=${BOOST_ROOT} ..
make -j24
