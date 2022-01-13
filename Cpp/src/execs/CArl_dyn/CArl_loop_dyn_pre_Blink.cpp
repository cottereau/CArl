#include "CArl_loop_dyn.h"

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

  carl::FETI_Dyn_Operations feti_op(WorldComm,input_params.scratch_folder_path);
  
  // Calculate the rhs vector for solving B_link_acceleration
  bool first_calculate = true; 
  feti_op.rhs_link(input_params.coupling_matrix_B,
    input_params.scratch_folder_path+"/coupling_interpolation_vec.petscvec",
    input_params.scratch_folder_path+"/rhs_vec_A_link.petscvec");
      

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