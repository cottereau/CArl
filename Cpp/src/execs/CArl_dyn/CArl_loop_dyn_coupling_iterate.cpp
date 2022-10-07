/*
 *
 *  Created on: April 26, 2022
 *      Author: Chensheng Luo
 */
#include "CArl_loop_dyn.h"
#include "FETI_operations.h"

/** \file CArl_loop_dyn_coupling_iterate.cpp

\brief **DYN-CG** Program responsible 

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

  // --- Read the common files: coupling matrices, null space vectors ...
  // Read up the coupling matricesconst std::string& filename)
  CG_op.set_coupling_matrix_R_micro();
  CG_op.set_coupling_matrix_R_BIG();

  // --- Set up any matrices or vectors needed before calculating the outputs
  // Set up the preconditioner
  CG_op.set_preconditioner(input_params.CG_precond_type, /* initial_set = */ false);

  // Read operations needed if we are using the rigid body modes
  if(input_params.bUseRigidBodyModes)
  {
    // Read the RB-related vectors and matrices
    CG_op.read_null_space_vecs(input_params.RB_vectors_base,input_params.nb_of_rb_vectors);
    CG_op.read_null_space_inv_RITRI_mat();
  }

  // --- Read the previous iteration files: scalar data, iteration vectors ...
  /* Read the scalar data from the previous iterations '0 ... kkk'
   * We now have: 'kkk'
   *        'rho(0)'
   *        'rho(kkk)'
   *        '| RB_corr(kkk) |'
   *        'p(0 ... kkk - 1).q(0 ... kkk - 1)'
   */
  CG_op.read_scalar_data();

  // Read the vector data from the previous iterations '0 ... kkk'
  /* We now have: 'r(kkk)'
   *        'phi(kkk)'
   *        'p(0 ... kkk)'
   *        'q(0 ... kkk - 1)'
   */
  CG_op.read_vector_data();

  // Read the previous iteration 'kkk' external solver output, 'x_i(kkk)'
  CG_op.read_ext_solver_output();

  // --- Iterate!
  // Calculate 'q(kkk) = C_1 * x_1(kkk) + C_2 * x_2(kkk)' and 'p(kkk).q(kkk)'
  CG_op.calculate_q();

  // Calculate 'phi(kkk + 1) = phi(kkk) + gamma * p(kkk)'
  // gamma = rho(kkk) / ( p(kkk).q(kkk) )
  CG_op.calculate_phi();

  // Calculate 'r(kkk + 1) = r(kkk) - gamma * q(kkk)'
  CG_op.calculate_r();

  // Calculate 'z(kkk + 1)' (formula depends on preconditioner and projection settings)
  CG_op.calculate_z();

  // Calculate 'p(kkk + 1)'
  CG_op.calculate_p();

  // Calculate 'RB_corr(kkk+1)'
  if(input_params.bUseRigidBodyModes)
  {
    CG_op.calculate_rb_correction();
  }

  // Calculate the scalar data, 'rho(kkk+1)' and '| RB_corr(kkk+1) |'
  CG_op.calculate_scalar_data();

  /* Export scalar data
   * Data to export: 'kkk+1'
   *           'rho(0)'
   *           'rho(kkk+1)'
   *           '| RB_corr(kkk+1) |'
   *           'p(kkk).q(kkk)'
   */
  CG_op.export_scalar_data();

  /* Export the iteration vectors
   * Vectors to export: 'r(kkk+1)'
   *                    'phi(kkk+1)'
   *                    'p(kkk+1)'
   *                    'q(kkk)'
   */
  CG_op.export_iter_vecs();

  // // --- Check the convergence
  carl::IterationStatus current_iteration_status = carl::IterationStatus::ITERATING;
  
  if(input_params.bUseRigidBodyModes)
  {
    current_iteration_status = CG_op.check_convergence(input_params.CG_coupled_conv_rel, input_params.CG_coupled_conv_abs, input_params.CG_coupled_conv_max, input_params.CG_coupled_div, input_params.CG_coupled_conv_corr);
  } else {
    current_iteration_status = CG_op.check_convergence(input_params.CG_coupled_conv_rel, input_params.CG_coupled_conv_abs, input_params.CG_coupled_conv_max, input_params.CG_coupled_div);
  }

  // Print the current values of the convergence parameters
  CG_op.print_previous_iters_conv( /* nb. of iterations = 5 */);

  switch (current_iteration_status)
  {
    case carl::IterationStatus::ITERATING :
        // --- Continue the iteration

        // Export the Ct_i * p(kkk+1) vectors
        CG_op.export_ext_solver_rhs_Ct_p();

        // --- Launch the "iter_script.sh" script --- ONLY ON THE FIRST PROC!
        if(WorldComm.rank() == 0)
        {
          std::string iter_script_command = ". " + input_params.scratch_folder_path + "/CG_solver/coupling_iterate.sh";
          if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
          {
            std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
            std::cout << iter_script_command << std::endl << std::endl;
          } else {
            carl::exec_command(iter_script_command);
          }
        }
        break;
    case carl::IterationStatus::CONVERGED :
        // --- Well ... converged!

        // Export the Ct_i * phi(kkk+1) vectors
        CG_op.export_ext_solver_rhs_Ct_phi();

        if(input_params.bUseRigidBodyModes){
          // Export the rigid body modes correction vector
          CG_op.export_rb_correction_vector();
        }
        // --- Launch the "sol_script.sh" script --- ONLY ON THE FIRST PROC!
        if(WorldComm.rank() == 0)
        {
          std::string sol_script_command = ". " + input_params.scratch_folder_path + "/CG_solver/coupling_solution.sh";
          if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
          {
            std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
            std::cout << sol_script_command << std::endl << std::endl;
          } else {
            carl::exec_command(sol_script_command);
          }
        }
        break;
    case carl::IterationStatus::DIVERGED :
        std::cout << " !!! Divergence detected!!! Stopped automatically! " << std::endl;
        break;
  }

  return 0;
}