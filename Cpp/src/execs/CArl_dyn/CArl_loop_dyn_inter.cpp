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

  // get jjj and m
  std::string progression_filename;
  progression_filename = input_params.scratch_folder_path+"/iteration_progression.txt";
  field_parser.parse_input_file(progression_filename, "#", "\n", " \t\n");
  carl::feti_loop_dyn_iteration_progression_params progression_params;
  get_input_params(field_parser, progression_params);
  int jjj=progression_params.inner_loop_progression+1;
  int m=input_params.inner_loop_times;

  feti_op.interpolate_A_acceleration(jjj,m);
  feti_op.rhs_interpolation(input_params.coupling_matrix_A,input_params.coupling_matrix_B);



}