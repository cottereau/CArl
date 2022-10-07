/*
 *
 *  Created on: Nov 17，2021
 *      Author: Chensheng Luo
 */
#include "CArl_loop_dyn.h"
#include <unistd.h>
/** \file CArl_loop_dyn_setup.cpp

\brief  **DYN-DI/DYN-CG** Program responsible to set up all environments, which will create all necesarry input files.

This program's input file description can be found at the documentation of the function carl::get_input_params(GetPot& field_parser, feti_loop_dyn_params& input_params). 

In this step, it will generate all setup files, as well as prepare all forces for futur calculation.
Initial conditions will equally be considered to generate the initial RHS vectors.
*/


int main(int argc, char* argv[]) {

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



  std::cout << " !!! Finish parser " << std::endl;

  // --- Crete the files / folders needed
  if(input_params.dyn_solver == carl::DynamicSolver::DI){
        carl::Dyn_DI_Solver_Files_Setup DI_files_setup(WorldComm,input_params);
        DI_files_setup.set_scratch_folder();
        DI_files_setup.generate_libmesh_external_solver_inputs();
        DI_files_setup.generate_libmesh_external_solver_script();
        DI_files_setup.generate_inner_operation_script();
        DI_files_setup.generate_combined_scripts();
        DI_files_setup.generate_progression_inputs();
        std::cout << " !!! Script Files Generated " << std::endl;
    }else if(input_params.dyn_solver == carl::DynamicSolver::CG ){
        carl::Dyn_CG_Solver_Files_Setup CG_files_setup(WorldComm,input_params);
        //Prepare work for FETI_files_setup
        CG_files_setup.set_scratch_folder();
        CG_files_setup.generate_libmesh_external_solver_inputs();
        CG_files_setup.generate_libmesh_external_solver_script();
        CG_files_setup.generate_inner_operation_script();
        CG_files_setup.generate_combined_scripts();
        CG_files_setup.generate_progression_inputs();
        std::cout << " !!! Script Files Generated " << std::endl;
    }else{
        homemade_error_msg("Invalid dynamic solver, please configure to DI or CG!");
    }

  //Generate script file
  carl::FETI_Dyn_Operations Dyn_op(WorldComm,input_params.scratch_folder_path,input_params.result_folder_path);
  
  #ifdef TEST_COUPLING_RESULT
    Dyn_op.init_test_coupling();
  #endif
  
  #ifdef PRINT_CALCULATION_TIME
    Dyn_op.export_calculation_time(0,"Begin");
  #endif

  Dyn_op.init_prepare_force(input_params.force_prepare_method,
    input_params.force_prepare_params,
    input_params.inner_loop_times,
    input_params.newmark_B.deltat,
    &Dyn_op.vector_A,
    &Dyn_op.vector_B);

  Dyn_op.set_initial_condition(&input_params.initial_A,
    &input_params.initial_B,
    &input_params.matrix_A,
    &input_params.matrix_B,
    &input_params.newmark_A,
    &input_params.newmark_B,
    &Dyn_op.vector_A,
    &Dyn_op.vector_B,
    input_params.result_folder_path);

  // Calculate A_free_acc

  if(WorldComm.rank() == 0){
     std::string init_script_command = ". " + input_params.scratch_folder_path + "/Afree.sh";
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