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

  //total_iterations = (int)(input_params.simulation_duration/input_params.deltatA);
  //inner_iterations = (int)(input_params.deltatA/input_params.deltatB);

  //Generate script file
  carl::Dyn_Solver_Files_Setup FETI_files_setup(WorldComm,input_params);
  //Prepare work for FETI_files_setup
  FETI_files_setup.set_scratch_folder();

  FETI_files_setup.generate_libmesh_external_solver_inputs();
  FETI_files_setup.generate_libmesh_external_solver_script();
  FETI_files_setup.generate_inner_operation_script();
  FETI_files_setup.generate_combined_scripts();

  carl::FETI_Dyn_Operations feti_op(WorldComm,input_params.scratch_folder_path);
  feti_op.init_prepare_vector();

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