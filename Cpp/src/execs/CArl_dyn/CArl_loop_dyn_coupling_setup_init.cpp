/*
 *
 *  Created on: April 26, 2022
 *      Author: Chensheng Luo
 */
#include "CArl_loop_dyn.h"

/** \file CArl_loop_dyn_coupling_setup_init.cpp

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

  carl::FETI_Dyn_Operations Dyn_op(WorldComm,input_params.scratch_folder_path,input_params.result_folder_path);

  // get jjj and m
  std::string progression_filename;
  progression_filename = input_params.scratch_folder_path+"/iteration_progression.txt";
  field_parser.parse_input_file(progression_filename, "#", "\n", " \t\n");
  carl::feti_loop_dyn_iteration_progression_params progression_params;
  get_input_params(field_parser, progression_params);
  int jjj=progression_params.inner_loop_progression+1;
  int m=input_params.inner_loop_times;
  int index=progression_params.inner_loop_progression+progression_params.outer_loop_progression*input_params.inner_loop_times+1;

  Dyn_op.interpolate_A_disp(jjj,m);

  Dyn_op.prepare_CG_free_result(input_params.scratch_folder_path,
    input_params.scratch_folder_path+"/CG_solver");
  Dyn_op.prepare_CG_scaled_matrix(input_params.matrix_A.mass_tilde,input_params.matrix_B.mass_tilde,
    input_params.scratch_folder_path+"/CG_solver",
    &input_params.newmark_A,
    &input_params.newmark_B);
  
  if(input_params.bUseRigidBodyModes)
       {
  Dyn_op.prepare_CG_scaled_vector(input_params.RB_vectors_base,input_params.nb_of_rb_vectors,
    input_params.scratch_folder_path+"/force_B/force_"+std::to_string(index)+".petscvec",
    input_params.scratch_folder_path+"/CG_solver",
    &input_params.newmark_B);
    }
  
  
   carl::FETI_Operations CG_op(WorldComm,input_params.scratch_folder_path+"/CG_solver",input_params.coupling_folder_path);

   // --- Define if the rb modes will be used or not
   CG_op.using_rb_modes(input_params.bUseRigidBodyModes);

   // Read the matrices
   CG_op.set_coupling_matrix_R_micro();
   CG_op.set_coupling_matrix_R_BIG();
   if(input_params.bUseRigidBodyModes)
   {// --- If rigid body: calculate phi(0), if needed
    // Set the null space
    CG_op.set_null_space(input_params.scratch_folder_path+"/CG_solver/prepared_rb_vector",input_params.nb_of_rb_vectors);
    
   // Calculate the inital solution
   CG_op.calculate_null_space_phi_0(input_params.scratch_folder_path+"/CG_solver/current_force.petscvec");
   }else{
    CG_op.calculate_no_null_space_phi_0();
   }
   // Export the vectors
   CG_op.export_phi();
   CG_op.export_ext_solver_rhs_Ct_phi();
   
 

  if(WorldComm.rank() == 0)
  {
       std::string init_script_command = ". " + input_params.scratch_folder_path + "/CG_solver/coupling_init.sh";
       if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
       {
        std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
        std::cout << init_script_command << std::endl << std::endl;
       } else {
        carl::exec_command(init_script_command);
       }
  }
  return 0;
}