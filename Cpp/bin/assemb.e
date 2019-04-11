python(3):ERROR:105: Unable to locate a modulefile for 'python/2.7.13'
Currently Loaded Modulefiles:
  1) cmake/3.12.2             6) intel-mpi/2018.2
  2) git/2.19.0               7) paraview/5.6.0
  3) gmsh/4.0.4               8) phdf5/1.10.3-intel-mpi
  4) intel-compilers/2018.2   9) texlive/2017
  5) intel-mkl/2018.2
[1]PETSC ERROR: [2]PETSC ERROR: ------------------------------------------------------------------------
[2]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation, probably memory access out of range
[2]PETSC ERROR: Try option -start_in_debugger or -on_error_attach_debugger
[2]PETSC ERROR: [3]PETSC ERROR: ------------------------------------------------------------------------
[3]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation, probably memory access out of range
[3]PETSC ERROR: Try option -start_in_debugger or -on_error_attach_debugger
[3]PETSC ERROR: or see http://www.mcs.anl.gov/petsc/documentation/faq.html#valgrind
[3]PETSC ERROR: or try http://valgrind.org on GNU/linux and Apple Mac OS X to find memory corruption errors
[3]PETSC ERROR: configure using --with-debugging=yes, recompile, link, and run 
[0]PETSC ERROR: ------------------------------------------------------------------------
[0]PETSC ERROR: or see http://www.mcs.anl.gov/petsc/documentation/faq.html#valgrind
[2]PETSC ERROR: or try http://valgrind.org on GNU/linux and Apple Mac OS X to find memory corruption errors
[2]PETSC ERROR: configure using --with-debugging=yes, recompile, link, and run 
[2]PETSC ERROR: to get more information on the crash.
[2]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[3]PETSC ERROR: to get more information on the crash.
[3]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[3]PETSC ERROR: Signal received
[3]PETSC ERROR: [2]PETSC ERROR: Signal received
[2]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
[2]PETSC ERROR: Petsc Release Version 3.8.4, unknown 
[2]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
[3]PETSC ERROR: Petsc Release Version 3.8.4, unknown 
[3]PETSC ERROR: ./libmesh_assemble_lin_homogeneous__min_x_clamped_dyn on a arch-linux2-cxx-opt named fusion-shm18 by gattif Thu Jan 24 19:23:42 2019
[3]PETSC ERROR: ./libmesh_assemble_lin_homogeneous__min_x_clamped_dyn on a arch-linux2-cxx-opt named fusion-shm18 by gattif Thu Jan 24 19:23:42 2019
[2]PETSC ERROR: Configure options --known-level1-dcache-size=32768 --known-level1-dcache-linesize=64 --known-level1-dcache-assoc=8 --known-sizeof-char=1 --known-sizeof-void-p=8 --known-sizeof-short=2 --known-sizeof-int=4 --known-sizeof-long=8 --known-sizeof-long-long=8 --known-sizeof-float=4 --known-sizeof-double=8 --known-sizeof-size_t=8 --known-bits-per-byte=8 --known-memcmp-ok=1 --known-sizeof-MPI_Comm=4 --known-sizeof-MPI_Fint=4 --known-mpi-long-double=1 --known-mpi-int64_t=1 --known-mpi-c-double-complex=1 --known-sdot-returns-double=0 --known-snrm2-returns-double=0 --known-mklspblas-supports-zero-based=0 --known-has-attribute-aligned=1 --with-mpi-dir=/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/ CXXOPTFLAGS=-O3 COPTFLAGS=-O3 --with-debugging=0 --with-batch --known-mpi-shared-libraries=0 --with-blaslapack-dir=/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mkl/ --with-boost-dir=/home/gattif/LOCAL/boost
[2]PETSC ERROR: Configure options --known-level1-dcache-size=32768 --known-level1-dcache-linesize=64 --known-level1-dcache-assoc=8 --known-sizeof-char=1 --known-sizeof-void-p=8 --known-sizeof-short=2 --known-sizeof-int=4 --known-sizeof-long=8 --known-sizeof-long-long=8 --known-sizeof-float=4 --known-sizeof-double=8 --known-sizeof-size_t=8 --known-bits-per-byte=8 --known-memcmp-ok=1 --known-sizeof-MPI_Comm=4 --known-sizeof-MPI_Fint=4 --known-mpi-long-double=1 --known-mpi-int64_t=1 --known-mpi-c-double-complex=1 --known-sdot-returns-double=0 --known-snrm2-returns-double=0 --known-mklspblas-supports-zero-based=0 --known-has-attribute-aligned=1 --with-mpi-dir=/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/ CXXOPTFLAGS=-O3 COPTFLAGS=-O3 --with-debugging=0 --with-batch --known-mpi-shared-libraries=0 --with-blaslapack-dir=/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mkl/ --with-boost-dir=/home/gattif/LOCAL/boost
[3]PETSC ERROR: #1 User provided function() line 0 in  unknown file
#1 User provided function() line 0 in  unknown file
------------------------------------------------------------------------
application called MPI_Abort(MPI_COMM_WORLD, 59) - process 2
application called MPI_Abort(MPI_COMM_WORLD, 59) - process 3
Caught signal number 11 SEGV: Segmentation Violation, probably memory access out of range
[0]PETSC ERROR: Try option -start_in_debugger or -on_error_attach_debugger
[0]PETSC ERROR: or see http://www.mcs.anl.gov/petsc/documentation/faq.html#valgrind
[0]PETSC ERROR: or try http://valgrind.org on GNU/linux and Apple Mac OS X to find memory corruption errors
[0]PETSC ERROR: configure using --with-debugging=yes, recompile, link, and run 
[0]PETSC ERROR: to get more information on the crash.
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Signal received
[0]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
[0]PETSC ERROR: Petsc Release Version 3.8.4, unknown 
[0]PETSC ERROR: ./libmesh_assemble_lin_homogeneous__min_x_clamped_dyn on a arch-linux2-cxx-opt named fusion-shm18 by gattif Thu Jan 24 19:23:42 2019
[0]PETSC ERROR: Configure options --known-level1-dcache-size=32768 --known-level1-dcache-linesize=64 --known-level1-dcache-assoc=8 --known-sizeof-char=1 --known-sizeof-void-p=8 --known-sizeof-short=2 --known-sizeof-int=4 --known-sizeof-long=8 --known-sizeof-long-long=8 --known-sizeof-float=4 --known-sizeof-double=8 --known-sizeof-size_t=8 --known-bits-per-byte=8 --known-memcmp-ok=1 --known-sizeof-MPI_Comm=4 --known-sizeof-MPI_Fint=4 --known-mpi-long-double=1 --known-mpi-int64_t=1 --known-mpi-c-double-complex=1 --known-sdot-returns-double=0 --known-snrm2-returns-double=0 --known-mklspblas-supports-zero-based=0 --known-has-attribute-aligned=1 --with-mpi-dir=/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/ CXXOPTFLAGS=-O3 COPTFLAGS=-O3 --with-debugging=0 --with-batch --known-mpi-shared-libraries=0 --with-blaslapack-dir=/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mkl/ --with-boost-dir=/home/gattif/LOCAL/boost
[0]PETSC ERROR: #1 User provided function() line 0 in  unknown file
application called MPI_Abort(MPI_COMM_WORLD, 59) - process 0
[1]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation, probably memory access out of range
[1]PETSC ERROR: Try option -start_in_debugger or -on_error_attach_debugger
[1]PETSC ERROR: or see http://www.mcs.anl.gov/petsc/documentation/faq.html#valgrind
[1]PETSC ERROR: or try http://valgrind.org on GNU/linux and Apple Mac OS X to find memory corruption errors
[1]PETSC ERROR: configure using --with-debugging=yes, recompile, link, and run 
[1]PETSC ERROR: to get more information on the crash.
[1]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[1]PETSC ERROR: Signal received
[1]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
[1]PETSC ERROR: Petsc Release Version 3.8.4, unknown 
[1]PETSC ERROR: ./libmesh_assemble_lin_homogeneous__min_x_clamped_dyn on a arch-linux2-cxx-opt named fusion-shm18 by gattif Thu Jan 24 19:23:42 2019
[1]PETSC ERROR: Configure options --known-level1-dcache-size=32768 --known-level1-dcache-linesize=64 --known-level1-dcache-assoc=8 --known-sizeof-char=1 --known-sizeof-void-p=8 --known-sizeof-short=2 --known-sizeof-int=4 --known-sizeof-long=8 --known-sizeof-long-long=8 --known-sizeof-float=4 --known-sizeof-double=8 --known-sizeof-size_t=8 --known-bits-per-byte=8 --known-memcmp-ok=1 --known-sizeof-MPI_Comm=4 --known-sizeof-MPI_Fint=4 --known-mpi-long-double=1 --known-mpi-int64_t=1 --known-mpi-c-double-complex=1 --known-sdot-returns-double=0 --known-snrm2-returns-double=0 --known-mklspblas-supports-zero-based=0 --known-has-attribute-aligned=1 --with-mpi-dir=/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/ CXXOPTFLAGS=-O3 COPTFLAGS=-O3 --with-debugging=0 --with-batch --known-mpi-shared-libraries=0 --with-blaslapack-dir=/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mkl/ --with-boost-dir=/home/gattif/LOCAL/boost
[1]PETSC ERROR: #1 User provided function() line 0 in  unknown file
application called MPI_Abort(MPI_COMM_WORLD, 59) - process 1
