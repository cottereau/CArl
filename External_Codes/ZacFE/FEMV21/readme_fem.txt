Finite Element Toolbox Ver. 2.1
===============================
Date: 2002-04-01
Author: B. Rasmus Anthin.



WHAT IS IT
----------

This is a toolbox for computing ODEs or PDEs in BVPs using FEM in
1D, 2D and 3D.

* FEM1
   this is the main routine for solving ODE BVPs using any combination
   of Dirichlet, Neumann or Robin conditions.
* GENMAT1
   called by FEM1 and is the core of the program. It generates the
   matrices used for solving the linear equation system for the ODE.
* REFINE1
   with this routine you can refine the mesh over certain critical
   grid points. The the gridpoints will become nonuniformly linearly
   spaced.
* TEST1
   this is a test for the 1D case of FEM. Look through this example
   carefully in order to fully understand how FEM1 works.
* FEM2.mat
   contains examples over meshes/triangulations in 2D.
* PLOTGRID2
   plot mesh/triangulation in 2D and put a number in each corresponding
   element (triangle).
* QUADSPACE
   generates quadratically spaced vectors. That is, the spaces are
   linearly decreasing/increasing.
* FEM2
   this is the main routine for solving PDE BVPs using any combination
   of Dirichlet, Neumann or Robin conditions.
* GENMAT2
   called by FEM2 and is the core of the program. It generates the
   matrices used for solving the linear equation system for the PDE.
* TEST2
   test for the 2D case of FEM. Test this for better learning how
   to use FEM2 and other utilities.


INSTALLATION
------------

0. UNINSTALLATION: If you already have one of the earlier distributions,
then the best thing would be if you could remove the earlier archive and
then go proceed to the next step (2) below. To remove the archive simply
enter the "fem" directory and in Unix/Linux write "rm -rf". Please be
sure that you REALLY are in the current directory or else other valuable
files may be lost FOREVER! So use "rm -rf" with great caution!!! You
should not have any of your own files located in this directory if you are
running Matlab under Unix/Linux!
If you are running Matlab under Windows then use either "Windows
Commander" or the Windows built in "File Explorer" in order to delete the
earier installed archive.

1. Enter your matlab directory or the directory where you want to put the
toolbox.

2. Place the fem.#.#.zip archive file in this directory.

3. Linux/Unix: Install by writing "gunzip fem.#.#.zip" (you can
also use "unzip" as well).
   Windows: Copy the entire archive to the directory using
Windows Commander available from "http://www.ghisler.com",
or use WinZip.
You should then end up with a file structure like this:
/fem/
  fem1.m
  genmat1.m
  .
  .

4. Update your matlab path using the command "pathtool" or "editpath" if
you have not done this already.

If you would encounter any problems during the above steps, the please
contact me "e8rasmus@etek.chalmers.se" and I'll do my best to help you
through.



GETTING STARTED
---------------

Look at the included examples, such as TEST1.
And use help for a function for more info on how it works.


UPDATES AND VERSIONS
--------------------

Ver 1.0:
Only 1D case implemented. Works good. However the routines are very
unefficient due to that I'm integrating over each pair of basis functions.

Ver 1.1:
Using sparse matrices in the implentation of GENMAT1 described above.
This yields more efficient calculations and a grid size of 5000 only
takes a couple of seconds on a unix-machine.

Ver 1.2:
Error discovered in GENMAT1: A division with three was not done on the
vector at the right side of the equation.
The refinement algorithm now works properly (REFINE1). Earlier the
refinement took place only between the adjacent points of the refinement-
points and the refinement was exponential. Now the refinement is done
with linearly increasingly/decreasingly spaced grids and is much more smoth.
See REFINE1 and QUADSPACE for more info.

Ver 1.3:
Removed the Neumann special case. This is otherwise easily acheived
by setting GAMMA to zero in the Robin boundary conditions.
Added GENMAT2 and FEM2. GENMAT2 is generally completed but will be improved
later on. FEM2 is not completed and should not be used.

Ver 2.0:
Major changes of syntax.
Completed the 2D FEM solver. However the boundary orientation checking is
not yet implemented and satisfactory. However it should work well for most
cases if you design the Robin boundaries having the desired orientation.

Ver 2.1:
No this is no aprils fool joke. I think I have removed all bugs in the
FEM2 program. TEST2 is updated accordingly (take a glance on it in
order to better understand how the syntax is supposed to be).
The orientation is the same as how the order of nodes are entered in
the corresponding boundary matrix. The boundary-matrices are inputed to
the FEM2 as a cell-array of boundary condition matrices (see FEM2).


EOF.