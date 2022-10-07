/*
 *
 *  Created on: May 29，2022
 *      Author: Chensheng Luo
 */
#include "libmesh_solve_linear_system.h"

/** \brief Solve \f$ AX=B \f$ by applying directly \f$ X=A^{-1}B \f$, with \f$ A^{-1} \f$ already calculated
 *  
 * This file use the same parser as the libmesh solver parser, but several parameters are useless, only following parameters are used.
 *  Required parameters:
 *    - `SysMatrix` : path to the inversion of system matrix file, *i.e.* \f$ A^{-1}\f$.
 *    - `SysRHSVector` : path to the system RHS vector file, *i.e.* \f$ B\f$.
 *    - `OutputBase` : output filename base.
 *
 */

int main(int argc, char** argv) {

  // --- Initialize libMesh
  libMesh::LibMeshInit init(argc, argv);

  // Do performance log?
  libMesh::PerfLog perf_log("Main program");

  // libMesh's C++ / MPI communicator wrapper
  libMesh::Parallel::Communicator& WorldComm = init.comm();

  // Number of processors and processor rank.
  int rank = WorldComm.rank();
  int nodes = WorldComm.size();

  // --- Set up inputs

  // Command line parser
  GetPot command_line(argc, argv);

  // File parser
  GetPot field_parser;

  // If there is an input file, parse it to get the parameters. Else, parse the command line
  std::string input_filename;
  if (command_line.search(2, "--inputfile", "-i")) {
    input_filename = command_line.next(input_filename);
    field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
  } else {
    field_parser = command_line;
  }

  carl::libmesh_solve_linear_system_input_params input_params;
  carl::get_input_params(field_parser, input_params);
  unsigned int n_timesteps = input_params.n_timesteps;
  // Check libMesh installation dimension
  const unsigned int dim = 3;

  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // --- Set the matrix and vectors

  // Set up the PETSC versions
  Mat sys_mat_PETSC;
  Vec sys_rhs_vec_PETSC;
  Vec sys_sol_vec_PETSC;
  PetscInt M,N,Nvec;

  MatCreate(WorldComm.get(),&sys_mat_PETSC);
  VecCreate(WorldComm.get(),&sys_rhs_vec_PETSC);
  VecCreate(WorldComm.get(),&sys_sol_vec_PETSC);
  
  // Read
  carl::read_PETSC_matrix(sys_mat_PETSC, input_params.sys_matrix_file, WorldComm.get());
  carl::read_PETSC_vector(sys_rhs_vec_PETSC, input_params.sys_rhs_vec_file, WorldComm.get());

    MatGetSize(sys_mat_PETSC,&M,&N);
    VecGetSize(sys_rhs_vec_PETSC,&Nvec);
    if(Nvec!=N){
        homemade_error_msg("Matrix dimension doesn't match!");
    }

    VecSetSizes(sys_sol_vec_PETSC,PETSC_DECIDE,M);
    VecSetFromOptions(sys_sol_vec_PETSC);

    MatMult(sys_mat_PETSC,sys_rhs_vec_PETSC,sys_sol_vec_PETSC);

    libMesh::PetscVector<libMesh::Number> sys_sol_vec(sys_sol_vec_PETSC,WorldComm);

// Print MatLab debugging output? Variable defined at "carl_headers.h"
    #ifdef PRINT_MATLAB_DEBUG
        sys_sol_vec.print_matlab(input_params.output_base + "_sys_sol_vec.m");
    #endif

  // Export the solution vector
  carl::write_PETSC_vector(sys_sol_vec, input_params.output_base + "_sys_sol_vec.petscvec");

  // --- Cleanup!
  MatDestroy(&sys_mat_PETSC);
  VecDestroy(&sys_rhs_vec_PETSC);

  return 0;
}