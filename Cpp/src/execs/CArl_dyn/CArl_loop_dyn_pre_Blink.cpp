/*
 *  Created on: Nov 17，2021
 *      Author: Chensheng Luo
 */
#include "CArl_loop_dyn.h"

/** \file CArl_loop_dyn_pre_Blink.cpp

\brief  **DYN-DI** Program responsible to calculate \f$ C^{B} \Lambda \f$.
 
This program's input file description can be found at the documentation of the function carl::get_input_params(GetPot& field_parser, feti_loop_dyn_params& input_params). 


In this step, it will use the following files in scratch foler:
     

It will create:

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
  
  // Calculate the rhs vector for solving B_link_acceleration
  Dyn_op.rhs_link(&Dyn_op.vector_B,Dyn_op.m_coupling_vector,&input_params.matrix_B);
      

  //execute next file

  if(WorldComm.rank() == 0){
     std::string init_script_command = ". " + input_params.scratch_folder_path + "/Blink.sh";
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