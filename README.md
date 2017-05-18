# CArl

## PRESENTATION

This project is focused on the development of a software based on the [Arlequin multi-model coupling method](https://www.sciencedirect.com/science/article/pii/S0045782508003630). The main interest of this software is to allow, by its specific structure, the easy interfacing of different third-party softwares (developed and maintained outside of this project), and adapted to each of the models appearing in the coupling. Currently, this project includes two implementations of the CArl sofware:

1. a [MATLAB](http://www.mathworks.fr/products/matlab/) implementation.
2. a parallel C++ / MPI implementation, based on [libMesh](https://libmesh.github.io) and [CGAL](http://www.cgal.org).

This software is mainly developed at laboratoire MSSMat (Ecole Centrale Paris - CNRS).

* contact : [Regis Cottereau](mailto:regis.cottereau@ecp.fr)
* contributors (by order of first commit): R. Cottereau, C. Zaccardi, Y. Le Guennec, D. Neron, T. M. Schlittler

## MATLAB IMPLEMENTATION

The MATLAB implementation of the CArl software can be found at the directory `MATLAB`. Currently, the software to which it is interfaced includes :

1. a 1D/2D FEM acoustic code,
1. a Timoschenko beam code, 
1. an elastic code, and 
1. Comsol (http://www.comsol.com).

### INSTALLATION

Before using the software, you should make sure that you update the matlab path with the appropriate directories. In matlab, run
`>> addpath( genpath('install_dir_CArl/'));`
where you replace `install_dir_CArl` by the full path to the directory `CArl/`

Additionally, you might want to write this line in the startup file (generally located in `~/matlab/startup.m`)

This code has not been extensively tested, but should run on Matlab versions R2013a and newer (mainly because the objects triangulation and delaunayTriangulation are required)

To use the option __FE2D__ (optional), one should install and make available in the path the functions downloadable at http://www.mathworks.in/matlabcentral/fileexchange/27826-fast-assembly-of-stiffness-and-matrices-in-finite-element-method (by Talal Rahman and Jan Valdman)

To use the option __beam__ (optional), one should install and make available in the path the functions downloadable at https://github.com/wme7/aero-matlab/tree/master/FEM/Timoshenko_beam (by Manuel Diaz and A Ferreira)

To use the option __Comsol__ (optional), one should install and make available in the path COMSOL Multiphysics (see http://www.comsol.com/)

### USE

The main calling routine is `CArl.m`

Some examples (in 1D and 2D) can be launched through use of the routine `Test.m` (see the corresponding help). The tests should run without any problem on computers with around 2Go.
 
## C++ / MPI IMPLEMENTATION

The C++ / MPI implementation of the CArl software can found at the directory `Cpp`. Currently, it is capable of interfacing with external solvers based on the [PETSc](http://www.mcs.anl.gov/petsc/) toolkit (including the libMesh solvers, when compiled with PETSc support).

The usage of this implementation will be added in the near future, together with a general documentation and examples.

### PRE-REQUISITES

This code implementation use the following third party libraries: (the numbers indicate the oldest version for which they were tested)

1. [Boost](http://www.boost.org) (version 1.60.0)
1. [CGAL](http://www.cgal.org) (version 4.7)
1. [PETSc](http://www.mcs.anl.gov/petsc/) (version 3.6.2)
1. [libMesh](https://libmesh.github.io) (version 1.1.0, installed with [TetGen](http://wias-berlin.de/software/tetgen/) support for now)

The following compiler combinations were tested:

1. OS X / macOS : Clang (version 7.0.0) and OpenMPI (version 1.10.0)
1. Linux : Intel C++ compilers (version 16.0.3) and Intel MPI (version 5.1.2)

This code was not tested or compiled with other operational systems.

### INSTALLATION

The installation is done using [CMake](https://cmake.org) (version 3.4.2), and the following commands:

`cd [CArl root directory]/Cpp/bin`

`cmake ..`

`make`

This will compile the CArl software using the default system compilers and with the `Release` optimization flags (`-O3 -DNDEBUG`). If you want to change these, use the appropriate flags or an interface such as `ccmake` or `cmake-gui`.

The CMake script will search for the Boost and CGAL libraries at the default include paths. For the libMesh installation, it will search for a `LIBMESH_DIR` environement variable. If the environement variable is not found, it will set it as `/usr/local`. In both cases, the script will search the `libmesh-config` binary at the `$LIBMESH_DIR/bin` directory. 

## REFERENCES

Original references for the theory are (among others)

1. H. Ben Dhia. Multiscale mechanical problems: the Arlequin method, _Comptes Rendus de l'Academie des Sciences - Series IIB 326_ (1998), pp. 899-904.
1. H. Ben Dhia, G. Rateau. The Arlequin method as a flexible engineering design tool, _Int. J. Numer. Meths. Engr._ 62 (2005), pp. 1442-1462.

Many other papers make use of the method or describe specific aspects and implementation details. Some of the scientific papers that made use of the CArl software are included in [references.rtf](references.rtf)
