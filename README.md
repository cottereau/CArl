# CArl

## PRESENTATION

This project is focused on the development of a software based on the [Arlequin multi-model coupling method](https://www.sciencedirect.com/science/article/pii/S0045782508003630). The main interest of this software is to allow, by its specific structure, the easy interfacing of different third-party softwares (developed and maintained outside of this project), and adapted to each of the models appearing in the coupling. Currently, the software to which CArl is interfaced includes :

1. a 1D/2D FEM acoustic code,
1. a Timoschenko beam code, 
1. an elastic code, and 
1. Comsol (http://www.comsol.com).

This software is mainly developed at laboratoire MSSMat (Ecole Centrale Paris - CNRS).

* contact : [Regis Cottereau](mailto:regis.cottereau@ecp.fr)
* contributors (by order of first commit): R. Cottereau, C. Zaccardi, Y. Le Guennec, D. Neron, T. M. Schlittler

It is developed for [MATLAB](http://www.mathworks.fr/products/matlab/). This may however evolve in future versions.
 
## REFERENCES

Original references for the theory are (among others)

1. H. Ben Dhia. Multiscale mechanical problems: the Arlequin method, _Comptes Rendus de l'Academie des Sciences - Series IIB 326_ (1998), pp. 899-904.
1. H. Ben Dhia, G. Rateau. The Arlequin method as a flexible engineering design tool, _Int. J. Numer. Meths. Engr._ 62 (2005), pp. 1442-1462.

Many other papers make use of the method or describe specific aspects and implementation details. Some of the scientific papers that made use of the CArl software are included in [references.rtf](references.rtf)

## INSTALLATION

Before using the software, you should make sure that you update the matlab path with the appropriate directories. In matlab, run
`>> addpath( genpath('install_dir_CArl/'));`
where you replace `install_dir_CArl` by the full path to the directory `CArl/`

Additionally, you might want to write this line in the startup file (generally located in `~/matlab/startup.m`)

This code has not been extensively tested, but should run on Matlab versions R2013a and newer (mainly because the objects triangulation and delaunayTriangulation are required)

To use the option __FE2D__ (optional), one should install and make available in the path the functions downloadable at http://www.mathworks.in/matlabcentral/fileexchange/27826-fast-assembly-of-stiffness-and-matrices-in-finite-element-method (by Talal Rahman and Jan Valdman)

To use the option __beam__ (optional), one should install and make available in the path the functions downloadable at https://github.com/wme7/aero-matlab/tree/master/FEM/Timoshenko_beam (by Manuel Diaz and A Ferreira)

To use the option __Comsol__ (optional), one should install and make available in the path COMSOL Multiphysics (see http://www.comsol.com/)

## USE

The main calling routine is `CArl.m`

Some examples (in 1D and 2D) can be launched through use of the routine `Test.m` (see the corresponding help). The tests should run without any problem on computers with around 2Go.

## C++ (WORK IN PROGRESS)

A C++ implementation of CArl, saved in the folder "Cpp", is under development. Its main goal is to reproduce and expand the capabilities of the current MATLAB implementation. In the meantime, certain functions of the latter might be substituted by the former, using MATLAB’s [MEX-functions](http://fr.mathworks.com/help/matlab/matlab_external/introducing-mex-files.html).

### IMPLEMENTED CODE

* Efficient identification (linear, in most cases) of the intersections between the triangles/tetrahedrons of two meshings using an algorithm developed by [Gander *et al.*, 2009](http://link.springer.com/chapter/10.1007/978-3-642-02677-5_19)). The test executables `intersections_2D` and `intersections_3D` apply this algorithm to test 2D and 3D meshes, respectively.

### PRE-REQUISITES

The C++ version of this software depends on the following libraries

* [Boost](http://www.boost.org), tested with version 1.58: `random` and `variant` libraries.
* [CMake](http://www.cmake.org), tested with version 3.2.3: compilation.
* [CGAL](http://www.cgal.org), tested with version 4.6: data structures and geometry methods.

### INSTALLATION / COMPILATION

1. `cd CArl/Cpp/Mesh_intersection`
2. `cmake -DCGAL_DIR=[directory with the CGAL cmake files] .`
3. `make`

This will build the executables `intersections_2D` and `intersections_3D` inside the `Mesh_intersection`. Usually, the CGAL cmake files are found at the folder `${CGAL_ROOT}/cmake`. **DO NOT** use CGAL’s `cgal_create_CMakeLists` script - the `CMakeLists.txt` included with this code is already configured to compile the executables.