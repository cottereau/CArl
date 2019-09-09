#!/bin/bash
module purge
module load intel-compilers/2018.2
module load intel-mpi/2018.2  
module load intel-mkl/2018.2
module load libgmp/6.1.2
module load libmpfr/3.1.5
module load cmake/3.14.4 
#module load python/2.7.13

export PETSC_DIR="/home/gattif/srclib/petsc"
export PETSC_ARCH="arch-linux2-cxx-opt"
export CGAL_DIR="/home/gattif/LOCAL/cgal"
export LIBMESH_DIR="/home/gattif/LOCAL/libmesh"
export BOOST_ROOT="/home/gattif/LOCAL/boost"
export CC=icc CXX=icpc 
export I_MPI_CC=icc I_MPI_CXX=icpc 

cmake -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_C_FLAGS="-std=gnu99" ../

make -j 24 #>& ../log-2018.txt
