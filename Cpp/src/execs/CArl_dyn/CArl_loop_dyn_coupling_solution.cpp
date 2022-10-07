/*
 *
 *  Created on: April 26, 2022
 *      Author: Chensheng Luo
 */
#include "CArl_loop_dyn.h"


/** \file CArl_loop_dyn_coupling_solution.cpp

\brief  **DYN-CG** Program responsible to 
 
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
  carl::FETI_Operations CG_op(WorldComm,input_params.scratch_folder_path+"/CG_solver");

  // --- Define if the rb modes will be used or not
  CG_op.using_rb_modes(input_params.bUseRigidBodyModes);

  // Read the decoupled solutions, K_i * u_0,i  = F_i 
  CG_op.read_decoupled_solutions();

  // Read the solutions of K_i * x_i(FINAL)  = C_i^t * phi(FINAL)
  CG_op.read_ext_solver_output();

  // Read the rb modes correction, 'RB_corr(FINAL)'
  if(input_params.bUseRigidBodyModes)
  {
    CG_op.read_rb_corr();
   }

  // Calculate the solution
  CG_op.calculate_dynamic_coupling_solution(1/(input_params.newmark_A.beta*input_params.newmark_A.deltat*input_params.newmark_A.deltat),
    1/(input_params.newmark_B.beta*input_params.newmark_B.deltat*input_params.newmark_B.deltat));

  // Export it (finallly!)
  CG_op.export_dynamic_coupled_solution(input_params.scratch_folder_path);

  if(WorldComm.rank() == 0){
     std::string init_script_command = "sbatch " + input_params.scratch_folder_path + "/inner_ope_Blink.slurm";
     if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
     {
        std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
        std::cout << init_script_command << std::endl ;
     } else {
        carl::exec_command(init_script_command);
     }
  }

  return 0;
  
}