# REQUIREMENTS
Some libraries are already available on __FUSION__ and can be uploaded via `module` command:

`module purge`
`module load git/2.22.0`
`module load openmpi/4.0.0`
`module load hdf5/1.8.20`
`module load cmake/3.16.0`
`module load blas/3.8.0`
`module load lapack/3.8.0`
`module load valgrind/3.15.0`
`module load libgmp/6.1.2`
`module load libmpfr/3.1.5`

For non-native libraries, source files can be stored, for instance, in `~/srclib` and respective binaries/libraries in `~\LOCAL`:

`mkdir -p ~/srclib`
`mkdir -p ~/LOCAL`

Before install and compile, please export the following environment variables:

`export srclib=/home/user/srclib`
`export locinst=/home/user/LOCAL`
`export BOOST_ROOT=$locinst/boost`
`export CGAL_DIR=$locinst/cgal`
`export LIBMESH_DIR=$locinst/libmesh`
`export PETSC_DIR=$locinst/petsc`
`export PETSC_ARCH="linux-openmpi-gcc"


# INSTALL ___BOOST___

Download ___BOOST___ from [source (version 1.65.1)](https://sourceforge.net/projects/boost/files/boost/1.65.1/) in $srclib

Configure:

`cd $srclib`
`tar -xvf boost_1_65_1.tar.bz2`
`cd boost_1_65_1` 

`./bootstrap.sh --prefix=${locinst}/boost/ --with-toolset=gcc variant=release`

Before compile, modify the file `project-config.jam` by adding the line "using mpi;" at the end of file.

Install:

`./b2`
`./b2 install` 


# INSTALL ___CGAL___

Download ___CGAL___ (version 4.7) from [github repository](https://github.com/CGAL/cgal.git cgal):

`cd $srclib`
`git clone https://github.com/CGAL/cgal.git cgal`

Configure: 

`mkdir -p $locinst/cgal`
`cd $locinst/cgal`

`cmake -DCMAKE_INSTALL_PREFIX=${locinst}/cgal -DCMAKE_BUILD_TYPE=Release -DGMP_LIBRARIES=/gpfs/opt/libraries/libgmp/6.1.2/lib/libgmp.so -DMPFR_LIBRARIES=/gpfs/opt/libraries/libmpfr/3.1.5/lib/libmpfr.so -DGMP_INCLUDE_DIR=/gpfs/opt/libraries/libgmp/6.1.2/include -DMPFR_INCLUDE_DIR=/gpfs/opt/libraries/libmpfr/3.1.5/include/ $srclib/cgal`

Install:

`make -j24`
`make install -j24`


# INSTALL ___PETSc___

Download ___PETSc___ (version 3.13.1) from [github repository](https://gitlab.com/petsc/petsc.git):

`cd $locinst`
`git clone -b v3.13.1 https://gitlab.com/petsc/petsc.git petsc`

Configure:

`./configure PETSC_DIR=${locinst}/petsc PETSC_ARCH=linux-openmpi-gcc CXXOPTFLAGS=-O3 COPTFLAGS=-O3 --with-mpi-dir=${MPI_ROOT} --with-batch --known-mpi-shared-libraries=0 --with-debugging=0 --with-boost-dir=${BOOST_ROOT}`

Install (with check): 

`make PETSC_DIR=$locinst/petsc PETSC_ARCH=linux-openmpi-gcc all`
`make PETSC_DIR=$locinst/petsc PETSC_ARCH=linux-openmpi-gcc check`

# INSTALL ___LIBMESH___

Download ___LIBMESH___ (version 1.5.1) from the [github repository](git://github.com/libMesh/libmesh.git):

`cd $srclib`
`git clone git://github.com/libMesh/libmesh.git libmesh`

Configure:

`./configure --prefix=${locinst}/libmesh/ CXX=mpicxx CC=mpicc --with-cxx=mpicxx --with-cc=mpicc --with-vtk=yes --enable-unique-ptr --with-mpi=/gpfs/opt/libraries/openmpi/4.0.0/bin --disable-strict-lgpl --with-hdf5=/gpfs/opt/libraries/hdf5/1.8.20 --with-methods=opt --with-boost=${BOOST_ROOT}`

Make (with check): 

`make -j24`
`make -j24 check`
`make install`

# INSTALL __CArl__

Download ___CArl___ from [github repository](git@github.com:cottereau/CArl.git)

Configure:

`cd $locinst/CArl/Cpp/bin`
`cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -DBOOST_ROOT=${BOOST_ROOT} ..`

Install:

`make -j24`
