/*
 *
 *  Created on: April 26, 2022
 *      Author: Chensheng Luo
 */
#include "CArl_loop_dyn.h"

/** \file CArl_loop_dyn_coupling_setup_finish.cpp

\brief  **DYN-CG** Program responsible

This program's input file description can be found at the documentation of the function carl::get_input_params(GetPot& field_parser, feti_loop_dyn_params& input_params). 

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

  // --- Set up input_params

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

  carl::feti_loop_dyn_params input_params;
  get_input_params(field_parser, input_params);

  // Object containing the FETI operations
  carl::FETI_Operations CG_op(WorldComm,input_params.scratch_folder_path+"/CG_solver",input_params.coupling_folder_path);

  // --- Define if the rb modes will be used or not
  CG_op.using_rb_modes(input_params.bUseRigidBodyModes);

  // --- Read the files!

  // Read up the coupling matricesconst std::string& filename)
  CG_op.set_coupling_matrix_R_micro();
  CG_op.set_coupling_matrix_R_BIG();

    // Read the decoupled solutions, K_i * u_0,i  = F_i 
  CG_op.read_decoupled_solutions();

  // Read the solutions of K_i * x_i(0)  = C_i^t * phi(0)
  CG_op.read_ext_solver_output();
  // Read operations needed if we are using the rigid body modes
  if(input_params.bUseRigidBodyModes)
  {
    // Read the RB-related vectors and matrices
    CG_op.read_null_space_vecs(input_params.RB_vectors_base,input_params.nb_of_rb_vectors);
    CG_op.read_null_space_inv_RITRI_mat();
  }

  // --- Set up any matrices or vectors needed before calculating the outputs
  // Set up the preconditioner
  CG_op.set_preconditioner(input_params.CG_precond_type /*, initial_set = true */ );

  // --- Calculate the output vectors! All are saved internaly inside the object
  // Calculate r(0)
  CG_op.calculate_initial_r();

  // Calculate p(0)
  CG_op.calculate_initial_p();

  // Calculations needed if we are using the rigid body modes
  if(input_params.bUseRigidBodyModes)
  {
      // Calculate the rigid body modes correction RB_corr
      CG_op.calculate_rb_correction();
  }

  // --- Export output vectors!
  // Export 'r(0)' and 'p(0)'
  CG_op.export_inital_vecs();

  // Export the Ct_i * p(0) vectors
  CG_op.export_ext_solver_rhs_initial();

  // Export the scalar data, rho(0) and, if pertinent, |RB_corr|
  CG_op.export_initial_scalar_data();

  if(WorldComm.rank() == 0)
    {
        std::string iter_script_command = ". " + input_params.scratch_folder_path + "/CG_solver/coupling_iterate.sh";
        if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
        {
            std::cout << " !!! LOCAL test: MPI commands cannot be launched recursivelly !!! " << std::endl;
            std::cout << "     Run the following program by hand: " << std::endl;
            std::cout << iter_script_command << std::endl;
        } else {
            carl::exec_command(iter_script_command);
        }
    }

    return 0;
}